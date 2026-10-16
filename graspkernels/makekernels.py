"""Build GRASP text kernels for homogeneous spheres with Mishchenko's spher.f.

Entry point: fun(ifn, dest), with ifn a JSON file in input_files/kernelgeneration/
(run via runkernelgeneration.py). Keys:

  name          name of the kernel set; outputs go to <dest>/<name>/
  fnpre         file name preface of the kernels (e.g. "kernels")
  ratioName     shape id used in the file names (e.g. "100")
  ratioValue    value written in the "element, ratio" lines (e.g. 1.0)
  wavelength    wavelength [um] of the kernels (must match the grid file)
  gridFile      GRASP grid1.dat providing the radius nodes and scattering angles;
                it is copied verbatim into the output directory
  riFromKernel  GRASP kernel file whose refractive-index nodes are used, or
  mr, mi        explicit lists of real and imaginary (>= 0) refractive-index nodes
  elems         scattering-matrix elements to write, e.g. ["11","12","22","33","34","44"]
  edgeNodes     treatment of the first and last node:
                  "grasp": inner half of the hat only, with weight h (as in the standard
                           GRASP v1.1.3 sphere and spheroid kernels)
                  "half":  inner half of the hat only, with weight h/2 (exact for a
                           piecewise-linear size distribution that ends at the grid edge)
                  "full":  hat extended beyond the grid by one step
  quadNK        Gaussian points per quadrature subinterval in spher.f
  quadPtsPerHalf  quadrature points per half hat, as a schedule [[x1, P1], [x2, P2], ...]:
                P1 points if the half ends at size parameter <= x1, else P2 if <= x2, etc.
                (rounded up to a multiple of quadNK). The points needed are set by the
                Mie resonance widths, roughly a fixed number per half hat; see README.md.
  quadMiRef, quadMiMinFactor  for mi > quadMiRef the schedule is scaled by
                max(quadMiMinFactor, sqrt(quadMiRef/mi)), since resonances broaden with
                absorption (default: no scaling)
  ddelt         spher.f accuracy parameter DDELT (expansion truncation)
  mode          7 (NDISTR=7 triangle, default) or 6 (NDISTR=6 check construction)
  nproc         number of spher_tri processes run at the same time
  nice          nice increment for those processes
  keepWork      keep the per-refractive-index job and output files
  writeConversionJson  also write <name>_conversion.json for runkernelconversion.py
  quadRefine    optional {"miMax": m, "xMin": x, "factor": f}: nodes with apex size parameter
                > x, for imaginary indices <= m, use f times the points of the schedule
                (more points where weak absorption leaves narrow resonances)
  reuseFrom     optional raw_kernels.npz of an earlier build with the same grid, refractive
                indices and angles. Its values are copied, and only part of the kernel is
                computed, chosen by reuseRecompute: "edges" (first and last node, e.g. to
                change edgeNodes; default), "refine" (the quadRefine region, optionally
                limited to apex size parameters within reuseXRange = [xlo, xhi]) or "none"
                (nothing; only rewrite the text files, e.g. with other digits). Without
                reuseFrom the whole kernel is computed, with the same result.
  digits        significant digits of the text kernel values (default 7, as KERNELS_BASE; max 8)

Each refractive index is one spher_tri job; finished jobs are skipped on a rerun,
so an interrupted build can be resumed.
"""
import json
import os
import shutil
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import graspkerneltext as gkt

HERE = os.path.dirname(os.path.abspath(__file__))
EXE = os.path.join(HERE, 'spher_tri')
ELEM_COLS = {'11': 0, '12': 1, '22': 2, '33': 3, '34': 4, '44': 5}  # spher_tri output order
ELEM_SIGN = {'11': 1, '12': 1, '22': 1, '33': 1, '34': -1, '44': 1}  # GRASP sign vs Mishchenko
DEFAULTS = {'fnpre': 'kernels', 'ratioName': '100', 'ratioValue': 1.0,
            'elems': ['11', '12', '22', '33', '34', '44'], 'edgeNodes': 'grasp',
            'quadNK': 10, 'quadPtsPerHalf': [[1e9, 400]], 'quadMiRef': 0.0, 'quadMiMinFactor': 1.0,
            'ddelt': 1e-10, 'mode': 7,
            'nproc': 1, 'nice': 0, 'keepWork': False, 'writeConversionJson': True, 'digits': 7}


def load(ifn):
    with open(ifn) as fp:
        data = json.load(fp)
    cfg = dict(DEFAULTS)
    cfg.update(data)
    return cfg


def nodes(r, edge):
    """Hat edges for every node and the weight factor applied to edge nodes.

    Returns rlo, rhi (rlo == r or rhi == r means that half is omitted) and fac."""
    lr = np.log(r)
    rlo = np.exp(np.r_[2*lr[0]-lr[1], lr[:-1]])
    rhi = np.exp(np.r_[lr[1:], 2*lr[-1]-lr[-2]])
    fac = np.ones(len(r))
    if edge in ('grasp', 'half'):
        rlo[0], rhi[-1] = r[0], r[-1]
        if edge == 'grasp':
            fac[0] = fac[-1] = 2.0
    elif edge != 'full':
        raise ValueError('edgeNodes must be "grasp", "half" or "full"')
    return rlo, rhi, fac


def refined(m, xk, cfg):
    """True if node with apex size parameter xk uses the quadRefine schedule for index m."""
    qr = cfg.get('quadRefine')
    return bool(qr) and m.imag <= qr['miMax'] and xk > qr['xMin']


def nsub(x_hi, mi, cfg, refine=False):
    """Quadrature subintervals for a half hat ending at size parameter x_hi."""
    fac = 1.0
    if cfg['quadMiRef'] > 0 and mi > cfg['quadMiRef']:
        fac = max(cfg['quadMiMinFactor'], np.sqrt(cfg['quadMiRef'] / mi))
    if refine:
        fac *= cfg['quadRefine']['factor']
    for xmax, npts in cfg['quadPtsPerHalf']:
        if x_hi <= xmax:
            return max(1, int(np.ceil(fac * npts / cfg['quadNK'])))
    raise ValueError('quadPtsPerHalf does not cover x = %g' % x_hi)


def write_job(fn, cfg, angles, m, r, rlo, rhi):
    wl = cfg['wavelength']
    with open(fn, 'w') as f:
        f.write('%.8f %d %.3e %d\n%d\n' % (wl, cfg['quadNK'], cfg['ddelt'], cfg['mode'], len(angles)))
        f.write(' '.join('%.6f' % a for a in angles) + '\n')
        f.write('1\n%.9f %.9e\n%d\n' % (m.real, m.imag, len(r)))
        for rk, lo, hi in zip(r, rlo, rhi):
            ref = refined(m, 2*np.pi*rk/wl, cfg)
            nlo = 0 if lo >= rk else nsub(2*np.pi*rk/wl, m.imag, cfg, ref)
            nhi = 0 if hi <= rk else nsub(2*np.pi*hi/wl, m.imag, cfg, ref)
            f.write('%.10e %.10e %.10e %d %d\n' % (rk, lo, hi, nlo, nhi))


def read_out(fn, nr, nang):
    hdr = np.fromfile(fn, dtype='<i4', count=3)
    assert tuple(hdr) == (1, nr, nang), 'unexpected spher_tri output header %s in %s' % (hdr, fn)
    v = np.fromfile(fn, dtype='<f8', offset=12)
    assert v.size == nr*(2+6*nang), 'truncated spher_tri output %s' % fn
    return v.reshape(nr, 2+6*nang)


def run_jobs(cfg, work, angles, jobs, r, rlo, rhi):
    """jobs: list of (n, m, ks) = output index, refractive index, node indices to compute."""
    os.makedirs(work, exist_ok=True)
    exe = os.path.join(work, 'spher_tri')  # private copy: rebuilding the repo cannot touch running jobs
    shutil.copy2(EXE, exe)
    todo = []
    for n, m, ks in sorted(jobs, key=lambda j: j[1].imag):   # most expensive (least absorbing) first
        job, out = os.path.join(work, 'job%04d.txt' % n), os.path.join(work, 'out%04d.bin' % n)
        if os.path.exists(out) and os.path.getsize(out) == 12 + 8*len(ks)*(2+6*len(angles)):
            continue
        write_job(job, cfg, angles, m, r[ks], rlo[ks], rhi[ks])
        todo.append((n, job, out))
    print('%d of %d jobs to compute with %d process(es)' % (len(todo), len(jobs), cfg['nproc']), flush=True)
    t0 = time.time()

    def one(item):
        n, job, out = item
        tmp = out + '.part'
        with open(job + '.log', 'w') as log:
            subprocess.run(['nice', '-n', str(cfg['nice']), exe, job, tmp], check=True, stderr=log)
        os.replace(tmp, out)
        return n

    with ThreadPoolExecutor(max_workers=cfg['nproc']) as ex:
        futs = [ex.submit(one, it) for it in todo]
        for i, f in enumerate(as_completed(futs)):
            n = f.result()
            el = time.time() - t0
            print('  done %d/%d (m index %d), elapsed %.1f min, est. remaining %.1f min'
                  % (i+1, len(todo), n, el/60, el/60*(len(todo)-i-1)/(i+1)), flush=True)


def fun(ifn, dest):
    cfg = load(ifn)
    if not os.path.exists(EXE):
        raise FileNotFoundError('%s not found; build it with compile_all.sh' % EXE)
    r, wl_grid, angles = gkt.read_grid(cfg['gridFile'])
    assert abs(wl_grid - cfg['wavelength']) < 1e-6, 'wavelength differs from the grid file'
    if 'riFromKernel' in cfg:
        mr, mi = gkt.read_ri_nodes(cfg['riFromKernel'])
    else:
        mr, mi = np.sort(np.array(cfg['mr'], float)), np.sort(np.array(cfg['mi'], float))
    mlist = [complex(a, b) for a in mr for b in mi]
    rlo, rhi, fac = nodes(r, cfg['edgeNodes'])
    nr, nang = len(r), len(angles)
    x = 2*np.pi*r/cfg['wavelength']

    outdir = os.path.join(dest, cfg['name'])
    work = os.path.join(outdir, 'work')
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, 'generation_settings.json'), 'w') as fp:
        json.dump(cfg, fp, indent=4)

    allk = np.arange(nr)
    if cfg.get('reuseFrom'):
        old = np.load(cfg['reuseFrom'])
        for key, val in (('r', r), ('angles', angles), ('mr', mr), ('mi', mi)):
            assert np.array_equal(old[key], val), 'reuseFrom: %s differs from this configuration' % key
        what = cfg.get('reuseRecompute', 'edges')
        if what == 'edges':
            jobs = [(n, m, np.array([0, nr-1])) for n, m in enumerate(mlist)]
        elif what == 'refine':
            xlo, xhi = cfg.get('reuseXRange', [-np.inf, np.inf])
            jobs = [(n, m, allk[[refined(m, xk, cfg) and xlo < xk <= xhi for xk in x]]) for n, m in enumerate(mlist)]
            jobs = [j for j in jobs if len(j[2])]
        elif what == 'none':
            jobs = []
        else:
            raise ValueError('reuseRecompute must be "edges", "refine" or "none"')
        work = work + '_' + what
        print('reusing %s; recomputing %s (%d jobs)' % (cfg['reuseFrom'], what, len(jobs)), flush=True)
        ext, sca = old['ext'].copy(), old['sca'].copy()
        pm = {el: old['P'+el].copy() for el in cfg['elems']}
    else:
        jobs = [(n, m, allk) for n, m in enumerate(mlist)]
        ext, sca = np.zeros((len(mr), len(mi), nr)), np.zeros((len(mr), len(mi), nr))
        pm = {el: np.zeros((len(mr), len(mi), nr, nang)) for el in cfg['elems']}
    if jobs:
        run_jobs(cfg, work, angles, jobs, r, rlo, rhi)

    for n, m, ks in jobs:
        i, j = divmod(n, len(mi))
        v = read_out(os.path.join(work, 'out%04d.bin' % n), len(ks), nang) * fac[ks, None]
        ext[i, j, ks], sca[i, j, ks] = v[:, 0], v[:, 1]
        for el in cfg['elems']:
            c = ELEM_COLS[el]
            pm[el][i, j, ks, :] = ELEM_SIGN[el] * v[:, 2+c*nang:2+(c+1)*nang]
    np.savez_compressed(os.path.join(outdir, 'raw_kernels.npz'), r=r, angles=angles, mr=mr, mi=mi,
                        ext=ext, sca=sca, **{'P'+el: P for el, P in pm.items()})

    print('writing GRASP text kernels to %s' % outdir)
    shutil.copyfile(cfg['gridFile'], os.path.join(outdir, 'grid1.dat'))
    gkt.write_kernels(outdir, cfg['fnpre'], cfg['ratioName'], cfg['ratioValue'], cfg['wavelength'],
                      r, mr, mi, ext, ext - sca, angles, pm, digits=cfg['digits'])
    if cfg['writeConversionJson']:
        cj = gkt.conversion_json(cfg['name'], os.path.abspath(outdir), cfg['fnpre'], cfg['ratioName'],
                                 nr, len(mr), len(mi), nang, cfg['elems'])
        with open(os.path.join(outdir, cfg['name'] + '_conversion.json'), 'w') as fp:
            json.dump(cj, fp, indent=4)
    if not cfg['keepWork'] and os.path.isdir(work):
        shutil.rmtree(work)
