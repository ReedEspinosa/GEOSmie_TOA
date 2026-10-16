"""Read and write GRASP text kernels (kernels_<ratio>_<elem>.txt and grid1.dat).

The writer reproduces the layout of the standard GRASP v1.1.3 kernels in
grasp_open/src/retrieval/internal_files/KERNELS_BASE (Fortran-style E16.7 values,
10 per line). The header and block strings follow those files and Greema Regmi's
GRASPhexKernel.py (GSFC-GRASP-Python-Interface/TAMU_GRASP_Kernel_Generator),
generalized to any ratio, wavelength, grid and angle set.

Conventions of the values (all in 1/um for radii and wavelength in um):
  ext[k], abs[k]  = INT T_k(ln r) * C_ext/V (resp. C_abs/V) d ln r
  P_ij[k, angle]  = INT T_k(ln r) * C_sca/V * F_ij d ln r
with T_k the hat function of node k and F_ij the scattering matrix normalized so that
0.5 * INT F11 sin(theta) dtheta = 1. Element 34 uses GRASP's sign, which is opposite to
Mishchenko's b2.
"""
import math
import os
import numpy as np

PER_LINE = 10  # values per line, as in KERNELS_BASE
FOOTER = ['', 'Shape distribution:', '   1   nr', '   #              r             rd']


def fe(v, w=16, d=7):
    """Fortran-style Ew.d string (0.1 <= |mantissa| < 1), right-justified to width w."""
    v = float(v)
    if v == 0.0 or not math.isfinite(v):
        if not math.isfinite(v):
            raise ValueError('non-finite kernel value %r' % v)
        return ('0.' + '0'*d + 'E+00').rjust(w)
    s = '%.*e' % (d-1, abs(v))           # e.g. 7.654016e-07
    mant, ex = s.split('e')
    digits = mant.replace('.', '')
    out = '%s0.%sE%+03d' % ('-' if v < 0 else '', digits, int(ex)+1)
    return out.rjust(w)


def _lines(vals, per=PER_LINE, d=7):
    vals = list(vals)
    return [''.join(fe(v, 16, d) for v in vals[i:i+per]) for i in range(0, len(vals), per)]


def _ri_header(mr, mi):
    return ['%16.7f%16.7f range of real refr. index' % (min(mr), max(mr)),
            '%s%s range of imag refr. index' % (fe(-min(mi)), fe(-max(mi))),
            '%12d%12d  number of nodes for real and imag parts of refr. index' % (len(mr), -len(mi))]


def _footer(ratio_value):
    return FOOTER + ['   1%15.6f%s' % (ratio_value, fe(1.0, 15, 6))]


def write_grid(path, r, wavelength, angles):
    """grid1.dat in the KERNELS_BASE layout."""
    with open(os.path.join(path, 'grid1.dat'), 'w') as f:
        f.write('%4d%8.3f\n' % (len(r), wavelength))
        for ri in r:
            f.write('%s\n' % fe(ri, 15, 7))
        f.write('%12d\n' % len(angles))
        for a in angles:
            f.write('%6.2f\n' % a)


def write_kernels(path, fnpre, ratio_name, ratio_value, wavelength, r, mr, mi, ext, absb,
                  angles=None, pmat=None, digits=7):
    """Write <fnpre>_<ratio_name>_00.txt and, if pmat is given, one file per element.

    digits: significant digits of the kernel values (7 as in KERNELS_BASE; up to 8 fits the
    16-character field). Scattering is ext - abs in GRASP, so where the single-scattering albedo
    is small its precision is limited by this.

    r: radius nodes [um]; mr, mi: refractive-index nodes (mi >= 0);
    ext, absb: arrays (len(mr), len(mi), len(r)); pmat: dict elem -> array
    (len(mr), len(mi), len(r), len(angles)), elem in '11','12','22','33','34','44'.
    """
    if not 1 <= digits <= 8:
        raise ValueError('digits must be 1-8 to fit the 16-character field')
    os.makedirs(path, exist_ok=True)
    head = ['%s%s%16.5f  rmin, rmax, ratio' % (fe(min(r), 16, 6), fe(max(r), 16, 6), ratio_value),
            '%12d   number of grid radii' % -len(r)]
    blk = lambda el, i, j: ['%3d%14.5f  element, ratio' % (int(el), ratio_value),
                            '%14.5f%16.6f%s  wavel, rreal, rimag' % (wavelength, mr[i], fe(-mi[j]))]
    out = head + _ri_header(mr, mi)
    for i in range(len(mr)):
        for j in range(len(mi)):
            out += blk('0', i, j) + ['  EXTINCTION (1/um):'] + _lines(ext[i, j], d=digits) \
                   + ['  ABSORPTION (1/um):'] + _lines(absb[i, j], d=digits)
    _write(os.path.join(path, '%s_%s_00.txt' % (fnpre, ratio_name)), out + _footer(ratio_value))
    if pmat is None:
        return
    angle_lines = [''.join('%12.3f' % a for a in angles[i:i+PER_LINE]) for i in range(0, len(angles), PER_LINE)]
    for el, P in pmat.items():
        out = head + ['%12d   number of scattering angles' % len(angles)] + angle_lines + _ri_header(mr, mi)
        for i in range(len(mr)):
            for j in range(len(mi)):
                out += blk(el, i, j)
                for k in range(len(r)):
                    out += _lines(P[i, j, k], d=digits)
        _write(os.path.join(path, '%s_%s_%s.txt' % (fnpre, ratio_name, el)), out + _footer(ratio_value))


def _write(fn, lines):
    with open(fn, 'w') as f:
        f.write('\n'.join(lines) + '\n')


def conversion_json(name, path, fnpre, ratio_name, nr, nmr, nmi, nang, elems,
                    renormUpperX=0.3, angFwdPeak=2):
    """Settings for GEOSmie's convertkernels.py (input_files/kernelconversion/*.json) that match
    files written by write_kernels."""
    nl = math.ceil(nr/PER_LINE)
    na = math.ceil(nang/PER_LINE)
    return {"name": name, "path": path, "fnpre": fnpre, "ratios": [ratio_name],
            "contlen": 2 + 2*(1 + nl), "mrlen": nmr, "milen": nmi,
            "scathdrlen": 3 + na + 3, "scatcontlen": 2 + nr*na,
            "elems": elems, "scatelemlen": na,
            "renormUpperX": renormUpperX, "angFwdPeak": angFwdPeak, "kernelScaleFact": 1}


# ---------------------------------------------------------------- readers
def read_grid(fn):
    """Radii, wavelength and scattering angles from a GRASP grid1.dat."""
    tok = open(fn).read().split()
    n, wl = int(tok[0]), float(tok[1])
    r = np.array(tok[2:2+n], float)
    nang = int(tok[2+n])
    ang = np.array(tok[3+n:3+n+nang], float)
    return r, wl, ang


def read_ri_nodes(fn):
    """Real and imaginary (>= 0) refractive-index nodes, in file order, from any GRASP kernel file."""
    mr, mi = [], []
    with open(fn) as f:
        for line in f:
            if 'wavel, rreal, rimag' in line:
                t = line.split()
                mr.append(float(t[1]))
                mi.append(abs(float(t[2])))
    nmi = len(set(mi))
    return np.array(mr[::nmi]), np.array(mi[:nmi])


def read_kernels(path, fnpre, ratio_name, elems=('11', '12', '22', '33', '34', '44')):
    """Read a complete GRASP text-kernel set for one ratio, in any of the layouts in use
    (KERNELS_BASE, 10 values per line; or one row per radius as in the TAMUdust2020 kernels).

    Returns a dict with ratio_value, wavelength, mr, mi (>= 0), ext and abs (len(mr), len(mi), nr),
    angles, and pmat (dict elem -> array (len(mr), len(mi), nr, len(angles))). Values are as stored
    (i.e. multiplied by the node weight dlnr).
    """
    out = {}
    fn = os.path.join(path, '%s_%s_00.txt' % (fnpre, ratio_name))
    with open(fn) as f:
        lines = f.read().split('\n')
    out['ratio_value'] = float(lines[0].split()[2])
    nr = abs(int(lines[1].split()[0]))
    wl, mr, mi, ext, absb = [], [], [], [], []
    cur = None
    for ln in lines[2:]:
        if 'wavel, rreal, rimag' in ln:
            t = ln.split()
            wl.append(float(t[0])); mr.append(float(t[1])); mi.append(abs(float(t[2])))
        elif 'EXTINCTION' in ln:
            cur = []; ext.append(cur)
        elif 'ABSORPTION' in ln:
            cur = []; absb.append(cur)
        elif 'element, ratio' in ln or 'Shape distribution' in ln:
            cur = None
        elif cur is not None:
            cur.append(ln)
    nmi = len(set(mi))
    nmr = len(mr) // nmi
    out['wavelength'] = wl[0]
    out['mr'] = np.array(mr[::nmi])
    out['mi'] = np.array(mi[:nmi])
    for key, blocks in (('ext', ext), ('abs', absb)):
        arr = np.array(' '.join(' '.join(b) for b in blocks).split(), float)
        if arr.size != nmr * nmi * nr:
            raise ValueError('%s: expected %d %s values, found %d' % (fn, nmr*nmi*nr, key, arr.size))
        out[key] = arr.reshape(nmr, nmi, nr)
    out['angles'], out['pmat'] = None, {}
    for el in elems:
        fn = os.path.join(path, '%s_%s_%s.txt' % (fnpre, ratio_name, el))
        with open(fn) as f:
            lines = f.read().split('\n')
        i = next(k for k, ln in enumerate(lines) if 'number of scattering angles' in ln)
        nang = int(lines[i].split()[0])
        j = next(k for k, ln in enumerate(lines) if 'range of real refr. index' in ln)
        angles = np.array(' '.join(lines[i+1:j]).split()[:nang], float)
        if out['angles'] is None:
            out['angles'] = angles
        elif not np.array_equal(angles, out['angles']):
            raise ValueError('%s: scattering angles differ from the other elements' % fn)
        data, take = [], False
        for ln in lines[j:]:
            if 'wavel, rreal, rimag' in ln:
                take = True
            elif 'element, ratio' in ln or 'Shape distribution' in ln:
                take = False
            elif take:
                data.append(ln)
        arr = np.array(' '.join(data).split(), float)
        if arr.size != nmr * nmi * nr * nang:
            raise ValueError('%s: expected %d values, found %d' % (fn, nmr*nmi*nr*nang, arr.size))
        out['pmat'][el] = arr.reshape(nmr, nmi, nr, nang)
    return out
