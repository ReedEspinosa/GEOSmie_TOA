"""Independent check of spher_tri kernels with other Mie codes and an independent quadrature.

The hat-function integrals are recomputed here with Gauss-Legendre quadrature in ln r
over each half hat, using
  - W. Wiscombe's MIEV0 (compiled in double precision from the user's copy), and
  - pymiecoated (GEOSmie's Python Mie code), for x <= 50 only, where it is accurate to
    better than 1e-6 (its error grows to ~1e-4 at 100 < x < 3000, so it cannot referee
    1e-4 differences there).
Cases use mi >= 1e-3, for which both quadratures converge, so the differences measure
the Mie codes, the triangle definition, the normalization and the sign conventions.

Sign of element 34: Mishchenko's b2 = -Im(S2 S1*) with MIEV0's amplitudes; pymiecoated's
amplitudes are the complex conjugates of MIEV0's (opposite time convention), so with them
b2 = +Im(S2 S1*). GRASP's element 34 is -b2 (see makekernels.ELEM_SIGN).

  python graspkernels/tests/test_independent_mie.py <MIEV0 dir> <scratch dir>

<MIEV0 dir> must contain MIEV0.f and ErrPack.f.
"""
import os
import subprocess
import sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
EXE = os.path.join(HERE, '..', 'spher_tri')
mdir, tmp = sys.argv[1], sys.argv[2]
os.makedirs(tmp, exist_ok=True)
LAM, H = 0.34, 0.0463215
ANG = np.arange(181.)
K = 2*np.pi/LAM
NQ = 400                                   # Gauss-Legendre points per half (both codes)

# ---- build the MIEV0 driver
env = dict(os.environ)
if sys.platform == 'darwin' and 'SDKROOT' not in env:
    env['SDKROOT'] = subprocess.run(['xcrun', '--show-sdk-path'], capture_output=True, text=True).stdout.strip()
drv = os.path.join(tmp, 'miev0_driver')
subprocess.run(['gfortran', '-O2', '-fdefault-real-8', '-o', drv, os.path.join(HERE, 'miev0_driver.f'),
                os.path.join(mdir, 'MIEV0.f'), os.path.join(mdir, 'ErrPack.f')], check=True, env=env,
               stderr=subprocess.DEVNULL)


def miev0(xs, m):
    inp = '%d\n%s\n%d\n' % (len(ANG), ' '.join('%g' % a for a in ANG), len(xs))
    inp += ''.join('%.17g %.17g %.17g\n' % (x, m.real, m.imag) for x in xs)
    out = np.array(subprocess.run([drv], input=inp, capture_output=True, text=True, check=True).stdout.split(), float)
    out = out.reshape(len(xs), 3 + 4*len(ANG))
    return out[:, 1], out[:, 2], out[:, 3:].reshape(len(xs), len(ANG), 4)


def pymie(xs, m):
    sys.path.insert(0, os.path.join(ROOT, 'pymiecoated'))
    from pymiecoated import Mie
    qe, qs, S = [], [], []
    for x in xs:
        mie = Mie(x=x, m=m)
        qe.append(mie.qext()); qs.append(mie.qsca())
        row = []
        for u in np.cos(np.radians(ANG)):
            s1, s2 = mie.S12(u)
            a1, a2 = abs(s1)**2, abs(s2)**2
            row.append([0.5*(a1+a2), 0.5*(a2-a1), (s2*np.conj(s1)).real, (s2*np.conj(s1)).imag])
        S.append(row)
    return np.array(qe), np.array(qs), np.array(S)


def hat(rk, code, m):
    """K_ext, K_sca and K_11, K_12, K_33, K_34 (Bohren-Huffman S_ij) for node rk."""
    z, w = np.polynomial.legendre.leggauss(NQ)
    uk = np.log(rk)
    tot = None
    for sgn in (-1, 1):                   # left and right half
        u = uk + sgn * H * (z + 1) / 2     # maps [-1,1] onto the half
        T = 1 - np.abs(u - uk) / H
        wt = w * H / 2 * T
        r = np.exp(u)
        qe, qs, S = code(K * r, m)
        cv = 3 / (4 * r)                   # (C/V) = Q * 3/(4r)
        ext, sca = np.sum(wt * qe * cv), np.sum(wt * qs * cv)
        pm = np.einsum('i,ijk->jk', wt * 3 / (K**2 * r**3), S)   # 4 pi S_ij /(k^2 V)
        part = np.r_[ext, sca, pm.T.ravel()]
        tot = part if tot is None else tot + part
    return tot


def spher(rk_list, m):
    job, out = os.path.join(tmp, 'ind.job'), os.path.join(tmp, 'ind.bin')
    n = int(np.ceil(NQ / 10))
    with open(job, 'w') as f:
        f.write('%g 10 1e-12 7\n%d\n%s\n1\n%.17g %.17g\n%d\n' % (LAM, len(ANG), ' '.join('%g' % a for a in ANG),
                                                               m.real, m.imag, len(rk_list)))
        for rk in rk_list:
            f.write('%.17e %.17e %.17e %d %d\n' % (rk, rk*np.exp(-H), rk*np.exp(H), n, n))
    subprocess.run([EXE, job, out], check=True, stderr=subprocess.DEVNULL)
    return np.fromfile(out, '<f8', offset=12).reshape(len(rk_list), -1)


na = len(ANG)
ok = True
worst = 0.0
for m in (complex(1.5, 0.01), complex(1.37, 0.1), complex(1.6, 1e-3)):
    xn = np.array([0.5, 5., 30., 150., 600.])
    rk = xn / K
    s = spher(rk, m)
    for i, x in enumerate(xn):
        sv = s[i]
        mish = {'ext': sv[0], 'sca': sv[1], '11': sv[2:2+na], '12': sv[2+na:2+2*na], '22': sv[2+2*na:2+3*na],
                '33': sv[2+3*na:2+4*na], '34': sv[2+4*na:2+5*na], '44': sv[2+5*na:2+6*na]}
        codes = [('MIEV0', miev0, -1)] + ([('pymiecoated', pymie, 1)] if x <= 50 else [])
        for name, code, s34 in codes:
            ref = hat(rk[i], code, m)
            bh = {'ext': ref[0], 'sca': ref[1], '11': ref[2:2+na], '12': ref[2+na:2+2*na],
                  '33': ref[2+2*na:2+3*na], '34': ref[2+3*na:2+4*na]}
            scale = np.abs(bh['11']).max()
            d = {'ext': abs(mish['ext']/bh['ext']-1), 'sca': abs(mish['sca']/bh['sca']-1)}
            for el, bel, sg in [('11', '11', 1), ('22', '11', 1), ('12', '12', 1), ('33', '33', 1), ('44', '33', 1),
                                ('34', '34', s34)]:
                d[el] = np.abs(mish[el] - sg*bh[bel]).max() / scale
            worst = max(worst, max(d.values()))
            print('m=%-12s x=%6.1f %-11s ext %.1e sca %.1e | max|dP|/maxP11: 11 %.1e 22 %.1e 12 %.1e 33 %.1e 44 %.1e 34 %.1e'
                  % (m, x, name, d['ext'], d['sca'], d['11'], d['22'], d['12'], d['33'], d['44'], d['34']))
            ok &= max(d.values()) < 1e-5
print('worst difference %.1e -> %s' % (worst, 'PASS' if ok else 'FAIL'))
sys.exit(0 if ok else 1)
