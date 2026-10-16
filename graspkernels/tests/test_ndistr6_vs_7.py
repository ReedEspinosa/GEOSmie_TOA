"""Check of the NDISTR=7 addition to spher.f against unmodified upstream code.

spher_tri MODE=6 computes each half hat with the upstream bimodal volume log-normal
distribution (NDISTR=6), GAM=-1 and two very broad, offset modes whose difference is
linear in ln r (see spher_tri.f, subroutine HALF). MODE=7 uses the new NDISTR=7
triangle. Both must give the same kernels to ~1e-8.

  python graspkernels/tests/test_ndistr6_vs_7.py <scratch dir>
"""
import os
import subprocess
import sys
import numpy as np

EXE = os.path.join(os.path.dirname(__file__), '..', 'spher_tri')
tmp = sys.argv[1]
os.makedirs(tmp, exist_ok=True)
ang = np.arange(181.)
rk = np.array([0.003, 0.05, 0.4, 2.0, 9.0])     # x = 0.06 ... 166 at 0.34 um
h = 0.0463215
ms = [complex(1.37, 1e-4), complex(1.5, 1e-2), complex(1.67, 0.1)]
res = {}
for mode in (7, 6):
    job, out = os.path.join(tmp, 'm%d.job' % mode), os.path.join(tmp, 'm%d.bin' % mode)
    with open(job, 'w') as f:
        f.write('0.34 10 1e-10 %d\n181\n' % mode + ' '.join('%g' % a for a in ang) + '\n%d\n' % len(ms))
        f.writelines('%.6f %.6e\n' % (m.real, m.imag) for m in ms)
        f.write('%d\n' % len(rk))
        for r in rk:
            n = max(4, int(np.ceil(4 * 2*np.pi*r/0.34 * h)))
            f.write('%.9e %.9e %.9e %d %d\n' % (r, r*np.exp(-h), r*np.exp(h), n, n))
    subprocess.run([EXE, job, out], check=True, stderr=subprocess.DEVNULL)
    res[mode] = np.fromfile(out, '<f8', offset=12).reshape(len(ms), len(rk), -1)
a, b = res[7], res[6]
scale = np.abs(a).max(axis=-1, keepdims=True)          # per (m, node) scale of each record
rel_ext = np.abs(b[..., :2] / a[..., :2] - 1).max()
rel_pm = (np.abs(b[..., 2:] - a[..., 2:]) / np.abs(a[..., 2:2+181]).max(axis=-1, keepdims=True)).max()
print('max rel. diff ext/sca: %.2e ; max diff of P_ij relative to max P11: %.2e' % (rel_ext, rel_pm))
ok = rel_ext < 1e-6 and rel_pm < 1e-6
print('PASS' if ok else 'FAIL')
sys.exit(0 if ok else 1)
