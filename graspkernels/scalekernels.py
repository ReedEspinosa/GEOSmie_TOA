"""Multiply a GRASP text-kernel set by a constant factor and rewrite it.

Used to correct the TAMUdust2020 hexahedra kernels, which were built as Q_A*pi*r_v^2/V
instead of C/V = Q_A*A/V and are therefore low by the sphericity Psi at every size
(see the V5 notes in README.md).

Entry point: fun(ifn, dest), with ifn a JSON file in input_files/kernelgeneration/
(run via runkernelscaling.py). Keys:

  name          name of the output folder; outputs go to <dest>/<name>/
  source        folder holding the kernels to scale (and grid1.dat)
  fnpre         file name preface of the kernels (e.g. "kernels")
  ratioName     shape id in the file names (e.g. "299")
  elems         scattering-matrix elements, e.g. ["11","12","22","33","34","44"]
  factor        multiplier for EXTINCTION, ABSORPTION and every absolute P_ij
  factorNote    free text recorded with the factor (where it comes from)
  copyGrid      copy <source>/grid1.dat to the output folder unless one is already there
                (default true)
  digits        significant digits written (default 7, as KERNELS_BASE; use 8 for sources with
                8 digits, e.g. the TAMUdust2020 kernels, so that ext - abs keeps its precision)

The output uses the KERNELS_BASE layout of graspkerneltext.write_kernels. A file <fnpre>_<ratioName>_scaling.json records the source files (with md5
checksums), the factor and the GEOSmie_TOA commit.
"""
import datetime
import hashlib
import json
import os
import shutil
import subprocess
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import graspkerneltext as gkt
import grasppositivity

HERE = os.path.dirname(os.path.abspath(__file__))


def md5(fn, block=1 << 24):
    h = hashlib.md5()
    with open(fn, 'rb') as f:
        for chunk in iter(lambda: f.read(block), b''):
            h.update(chunk)
    return h.hexdigest()


def git_commit(path=HERE):
    try:
        rev = subprocess.run(['git', '-C', path, 'rev-parse', 'HEAD'], capture_output=True, text=True).stdout.strip()
        dirty = subprocess.run(['git', '-C', path, 'status', '--porcelain', '--untracked-files=no'],
                               capture_output=True, text=True).stdout.strip()
        return rev + ('-dirty' if dirty else '')
    except OSError:
        return 'unknown'


def fun(ifn, dest):
    with open(ifn) as fp:
        cfg = json.load(fp)
    cfg.setdefault('fnpre', 'kernels')
    cfg.setdefault('elems', ['11', '12', '22', '33', '34', '44'])
    cfg.setdefault('copyGrid', True)
    cfg.setdefault('digits', 7)
    src, fnpre, ratio, factor = cfg['source'], cfg['fnpre'], cfg['ratioName'], float(cfg['factor'])
    outdir = os.path.join(dest, cfg['name'])
    os.makedirs(outdir, exist_ok=True)

    print('reading %s/%s_%s_*' % (src, fnpre, ratio))
    k = gkt.read_kernels(src, fnpre, ratio, cfg['elems'])
    r, wl, angles = gkt.read_grid(os.path.join(src, 'grid1.dat'))
    if len(r) != k['ext'].shape[-1] or abs(wl - k['wavelength']) > 1e-6:
        raise ValueError('grid1.dat does not match the kernels (radii or wavelength)')
    if not np.allclose(angles, k['angles']):
        raise ValueError('grid1.dat angles do not match the kernels')

    print('writing %s/%s_%s_* (factor %.6f)' % (outdir, fnpre, ratio, factor))
    pm = {el: factor * k['pmat'][el] for el in cfg['elems']}
    gkt.write_kernels(outdir, fnpre, ratio, k['ratio_value'], k['wavelength'], r, k['mr'], k['mi'],
                      factor * k['ext'], factor * k['abs'], k['angles'], pm, digits=cfg['digits'])
    if cfg['copyGrid'] and not os.path.exists(os.path.join(outdir, 'grid1.dat')):
        shutil.copyfile(os.path.join(src, 'grid1.dat'), os.path.join(outdir, 'grid1.dat'))

    names = ['%s_%s_%s.txt' % (fnpre, ratio, el) for el in ['00'] + list(cfg['elems'])]
    prov = {'source': os.path.abspath(src),
            'source_md5': {n: md5(os.path.join(src, n)) for n in names + ['grid1.dat']},
            'factor': factor, 'factorNote': cfg.get('factorNote', ''),
            'output_md5': {n: md5(os.path.join(outdir, n)) for n in names},
            'grasp_positivity_violations': grasppositivity.violations(pm, k['angles']),
            'geosmie_commit': git_commit(), 'created': datetime.datetime.now().isoformat(timespec='seconds'),
            'settings': cfg}
    with open(os.path.join(outdir, '%s_%s_scaling.json' % (fnpre, ratio)), 'w') as fp:
        json.dump(prov, fp, indent=4)
    return outdir
