"""Compare two sets of GRASP text sphere kernels on the same grid and refractive indices.

  python graspkernels/tests/compare_kernels.py <new dir> <reference dir> [ratio] [--skip-angles 13,...]

Reports, for every imaginary refractive index, the median and maximum over real index,
radius node and angle of
  ext, abs : |new/ref - 1|
  P11      : |new/ref - 1|
  P_ij     : |new - ref| / ref_P11   (other elements cross zero)
separately for interior and edge nodes. Angles listed in --skip-angles are left out
(the "13.000" column of the GRASP v1.1.3 kernels was computed at ~12.94 deg).
"""
import os
import sys
import numpy as np

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
sys.path.insert(0, ROOT)
import convertkernels as ck
from graspkernels import graspkerneltext as gkt

args = [a for a in sys.argv[1:] if not a.startswith('--')]
skip = []
if '--skip-angles' in sys.argv:
    skip = [float(v) for v in sys.argv[sys.argv.index('--skip-angles')+1].split(',')]
    args = [a for a in args if a != sys.argv[sys.argv.index('--skip-angles')+1]]
new, ref = args[0], args[1]
ratio = args[2] if len(args) > 2 else '100'
ELEMS = ['11', '12', '22', '33', '34', '44']


def load(d):
    r, wl, ang = gkt.read_grid(os.path.join(d, 'grid1.dat'))
    mr, mi = gkt.read_ri_nodes(os.path.join(d, 'kernels_%s_00.txt' % ratio))
    cj = gkt.conversion_json('x', d, 'kernels', ratio, len(r), len(mr), len(mi), len(ang), [])
    sh = (len(mr), len(mi), len(r))
    e, a, _ = ck.read00(d, ratio, 'kernels', cj['contlen'])
    out = {'r': r, 'ang': ang, 'mr': mr, 'mi': mi, 'ext': np.array(e).reshape(sh), 'abs': np.array(a).reshape(sh)}
    for el in ELEMS:
        s = ck.readScatEle(d, ratio, el, 'kernels', cj['scathdrlen'], cj['scatcontlen'], cj['scatelemlen'])
        out[el] = np.array(s).reshape(sh + (len(ang),))
    return out


A, B = load(new), load(ref)
assert np.allclose(A['r'], B['r']) and np.allclose(A['mr'], B['mr']) and np.allclose(A['mi'], B['mi'])
keep = ~np.isin(B['ang'], skip)
x = 2*np.pi*B['r']/0.34
nodesets = {'interior': np.arange(1, len(x)-1), 'edges': np.array([0, len(x)-1])}
print('new: %s\nref: %s\nangles skipped: %s' % (new, ref, skip))
for nm, ks in nodesets.items():
    print('\n== %s nodes (x = %.3g .. %.3g)' % (nm, x[ks].min(), x[ks].max()))
    print('     mi        ext med/max        abs med/max        P11 med/max        other P_ij max/P11')
    for j, mij in enumerate(B['mi']):
        de = np.abs(A['ext'][:, j, ks]/B['ext'][:, j, ks] - 1)
        da = np.abs(A['abs'][:, j, ks]/B['abs'][:, j, ks] - 1)
        p11 = B['11'][:, j, ks][..., keep]
        d11 = np.abs(A['11'][:, j, ks][..., keep]/p11 - 1)
        dother = max(np.max(np.abs(A[el][:, j, ks][..., keep] - B[el][:, j, ks][..., keep])/p11) for el in ELEMS[1:])
        print('  %9.3e   %.1e / %.1e    %.1e / %.1e    %.1e / %.1e    %.1e'
              % (mij, np.median(de), de.max(), np.median(da), da.max(), np.median(d11), d11.max(), dother))
if skip:
    i = np.isin(B['ang'], skip)
    d = np.abs(A['11'][..., i]/B['11'][..., i] - 1)
    print('\nskipped angle(s) %s, P11 |new/ref-1|: median %.1e, max %.1e' % (skip, np.median(d), d.max()))
