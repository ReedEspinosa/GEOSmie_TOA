"""Quadrature convergence spot check of a kernel set built by makekernels.py.

Reruns selected refractive indices with the quadrature schedule multiplied by FACTOR and
compares with the build's raw_kernels.npz.

  python graspkernels/tests/convergence_spotcheck.py <generation json> <build dir> <scratch dir> \
         [FACTOR] [m indices, e.g. 0,15,30] [node stride]

With a node stride N, only every N-th node above x = 30 and the last node are recomputed.

m indices count over (mr, mi) in file order (mi inner). Reported per x range: ext and sca
|rel. diff| (median/max), P11 |rel. diff| (median / 95th percentile / max over angles
and nodes), and the sin-weighted P11 difference |dP11| sin / P11 sin.
"""
import json
import os
import subprocess
import sys
import numpy as np

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
sys.path.insert(0, ROOT)
from graspkernels import makekernels as mk

ifn, build, tmp = sys.argv[1:4]
factor = float(sys.argv[4]) if len(sys.argv) > 4 else 2.0
cfg = mk.load(ifn)
d = np.load(os.path.join(build, 'raw_kernels.npz'))
r, ang, mr, mi = d['r'], d['angles'], d['mr'], d['mi']
mlist = [complex(a, b) for a in mr for b in mi]
sel = [int(v) for v in sys.argv[5].split(',')] if len(sys.argv) > 5 else [0, len(mi)//2, len(mi)-1]
cfg['quadPtsPerHalf'] = [[x, p*factor] for x, p in cfg['quadPtsPerHalf']]
rlo, rhi, fac = mk.nodes(r, cfg['edgeNodes'])
xall = 2*np.pi*r/cfg['wavelength']
ks = np.arange(len(r))
if len(sys.argv) > 6:
    stride = int(sys.argv[6])
    ks = np.union1d(np.nonzero(xall > 30)[0][::stride], [len(r)-1])
os.makedirs(tmp, exist_ok=True)
nang = len(ang)
x = xall[ks]
ranges = [(0, 30), (30, 300), (300, 1000), (1000, 1e9)]
sw = np.sin(np.radians(ang))
for n in sel:
    m = mlist[n]
    job, out = os.path.join(tmp, 'sc%04d.txt' % n), os.path.join(tmp, 'sc%04d.bin' % n)
    mk.write_job(job, cfg, ang, m, r[ks], rlo[ks], rhi[ks])
    subprocess.run(['nice', '-n', str(cfg['nice']), mk.EXE, job, out], check=True, stderr=subprocess.DEVNULL)
    v = mk.read_out(out, len(ks), nang) * fac[ks, None]
    i, j = divmod(n, len(mi))
    old_ext, old_sca, old11 = d['ext'][i, j, ks], d['sca'][i, j, ks], d['P11'][i, j, ks]
    print('m = %.6f + %.4e i (index %d), quadrature x%g, %d nodes' % (m.real, m.imag, n, factor, len(ks)))
    for lo, hi in ranges:
        k = (x > lo) & (x <= hi)
        if not k.any():
            continue
        de = np.abs(old_ext[k]/v[k, 0] - 1)
        ds = np.abs(old_sca[k]/v[k, 1] - 1)
        p = v[k, 2:2+nang]
        dp = np.abs(old11[k]/p - 1)
        dw = np.sum(np.abs(old11[k]-p)*sw, 1)/np.sum(p*sw, 1)
        print('   %6g < x <= %-6g ext %.1e/%.1e  sca %.1e/%.1e  P11 %.1e / %.1e / %.1e  sin-weighted max %.1e'
              % (lo, hi, np.median(de), de.max(), np.median(ds), ds.max(), np.median(dp),
                 np.percentile(dp, 95), dp.max(), dw.max()))
