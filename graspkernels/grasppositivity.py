"""Make scattering-matrix kernels usable by the standard GRASP distribution.

GRASP stores element 33 as log(K33) for scattering angles <= 40 deg and element 44 as log(K44) for
angles <= 50 deg (dls_fixget_bin.f90, "LOG(U...)"), and interpolates those logs in refractive index
(dls_intrpl_orgn_bin.f90). Every node value there must therefore be positive. Size-resolved sphere
kernels violate this near Mie resonances (F33/F11 and F44/F11 dip to about -0.8 at 15-50 deg for
x ~ 5-19); kernels averaged over wide size bins, orientations or shapes usually do not.
(The standard GRASP v1.1.3 sphere kernel holds the value 1000.0 at exactly the 23 places where its true
value is negative.)

fix() applies a local, sum-conserving correction: for each refractive index, element and angle in the
log range, around each run of size nodes with K <= eps*K11 it takes the smallest window (grown one node
on each side at a time) in which sum(K) >= eps*sum(K11), and sets F/F11 in that window to its
K11-weighted mean, K = K11 * sum(K)/sum(K11). sum(K) over the window is unchanged, so the result is
exact for a size distribution that is flat across the window. Only the GRASP text files should be
written this way; the exact values belong in the netCDF files.
"""
import numpy as np

LIMITS = {'33': 40.0, '44': 50.0}   # GRASP takes log(K) at angles <= these


def violations(pmat, angles, limits=LIMITS, eps=0.0):
    """Number of node values K <= eps*K11 where GRASP takes log(K), per element."""
    out = {}
    for el, amax in limits.items():
        if el in pmat:
            sel = np.asarray(angles) <= amax
            out[el] = int(np.sum(pmat[el][..., sel] <= eps * pmat['11'][..., sel]))
    return out


def fix(pmat, angles, limits=LIMITS, eps=1e-3):
    """Return (corrected copy of pmat, statistics). pmat: dict elem -> array (nmr, nmi, nr, nang),
    must contain '11'. Elements not in limits are returned unchanged (same objects)."""
    out = dict(pmat)
    stats = {'eps': eps, 'limits_deg': {k: v for k, v in limits.items() if k in pmat}}
    for el, amax in limits.items():
        if el not in pmat:
            continue
        K = pmat[el].copy()
        K11 = pmat['11']
        widths = []
        for i, j, ia in zip(*np.nonzero(np.any(K <= eps * K11, axis=2) & (np.asarray(angles) <= amax))):
            k, k11 = K[i, j, :, ia], K11[i, j, :, ia]
            bad = np.flatnonzero(k <= eps * k11)
            while bad.size:
                lo = hi = bad[0]
                while hi + 1 < k.size and k[hi + 1] <= eps * k11[hi + 1]:
                    hi += 1
                while k[lo:hi+1].sum() < eps * k11[lo:hi+1].sum():
                    if lo == 0 and hi == k.size - 1:
                        raise ValueError('element %s: no positive window at mr #%d, mi #%d, angle %g' % (el, i, j, angles[ia]))
                    lo, hi = max(lo - 1, 0), min(hi + 1, k.size - 1)
                k[lo:hi+1] = k11[lo:hi+1] * (k[lo:hi+1].sum() / k11[lo:hi+1].sum())
                widths.append(hi - lo + 1)
                bad = np.flatnonzero(k <= eps * k11)
        out[el] = K
        nchanged = int(np.sum(K != pmat[el]))
        w = np.array(widths) if widths else np.zeros(1, int)
        stats[el] = {'violations_before': int(np.sum(pmat[el][..., np.asarray(angles) <= amax] <= 0)),
                     'windows': len(widths), 'window_nodes_median': float(np.median(w)), 'window_nodes_max': int(w.max()),
                     'values_changed': nchanged, 'fraction_of_element_changed': nchanged / K.size}
    return out, stats
