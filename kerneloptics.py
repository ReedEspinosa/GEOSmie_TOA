"""
Evaluate bulk optical efficiencies for an arbitrary PSD directly from the
GRASP-kernel netCDF files produced by convertkernels.py (e.g.
kernel-grasp-v1.1.3-integrated_V4.nc and the single-aspect-ratio spheroid /
hexahedra files in the same format).

Usage conventions established by ground-truth validation against exact Mie
(ASCENT-ACP scripts/mie_ground_truth.py, 2026-08-31; sphere slice agrees to
<=0.2% median, <=0.4% worst case on measured and synthetic PSDs):

1. The stored ext/sca variables are volume-normalized kernels (1/um) at the
   reference wavelength lambda0 = 0.340 um; exact rescale to wavelength lam
   is K_lam = K_stored * (0.340/lam) [K = 3*pi*Q/(2*x*lam)].
2. Each stored value is the triangular-in-ln(x) basis AVERAGE of the volume
   kernel, not a point sample. Before interpolating in x we therefore
   de-smooth: solve the tridiagonal Galerkin mass-matrix system so the
   piecewise-linear reconstruction is the L2 projection of the true kernel.
   Skipping this step biases coefficients low by ~1.7% and lets Mie ripples
   alias against the node spacing for narrow PSDs (up to ~4%).
3. The PSD is resampled to a fine ln(r) grid (piecewise-linear dN/dlnr) and
   integrated against the de-smoothed, log-x-interpolated kernel. Never
   sample the PSD at the kernel x-nodes (under-resolves the PSD).
4. mi is stored negative in the files; pass a non-negative mi here.
"""
import numpy as np
import netCDF4 as nc
from scipy.linalg import solve_banded

LAMBDA0 = 0.340  # um; reference wavelength of the stored kernels

_cache = {}


def _loadKernels(kernelFile, ratioInd):
  key = (kernelFile, ratioInd)
  if key not in _cache:
    with nc.Dataset(kernelFile) as ds:
      _cache[key] = {
        'x': ds['x'][:].filled(np.nan),
        'mr': ds['mr'][:].filled(np.nan),
        'mi': np.abs(ds['mi'][:].filled(np.nan)),  # stored negative
        'ext': ds['ext'][ratioInd].filled(np.nan),  # (mr, mi, x), 1/um at lambda0
        'sca': ds['sca'][ratioInd].filled(np.nan),
      }
  return _cache[key]


def _interpMrMi(kern3d, allmr, allmi, mr, mi):
  """Bilinear in (mr, log mi) -> kernel vs x. mi=0 maps to the first node."""
  i = int(np.clip(np.searchsorted(allmr, mr) - 1, 0, len(allmr) - 2))
  fmr = (mr - allmr[i]) / (allmr[i + 1] - allmr[i])
  logmi = np.log(np.maximum(allmi, 1e-10))
  lmi = np.log(max(mi, 1e-10))
  j = int(np.clip(np.searchsorted(logmi, lmi) - 1, 0, len(allmi) - 2))
  fmi = (lmi - logmi[j]) / (logmi[j + 1] - logmi[j])
  return ((1 - fmr) * (1 - fmi) * kern3d[i, j]
          + fmr * (1 - fmi) * kern3d[i + 1, j]
          + (1 - fmr) * fmi * kern3d[i, j + 1]
          + fmr * fmi * kern3d[i + 1, j + 1])


def _desmooth(kbar, lnx):
  """
  Invert the triangular basis-averaging: solve M g = kbar where
  M_ii-1 = h_{i-1}/6, M_ii = (h_{i-1}+h_i)/3, M_ii+1 = h_i/6, each row
  normalized by the basis integral (h_{i-1}+h_i)/2 (half-triangles at the
  boundaries). The result g is the node vector of the L2 projection of the
  true kernel onto the piecewise-linear space.
  """
  n = len(kbar)
  h = np.diff(lnx)
  ab = np.zeros((3, n))
  norm = np.empty(n)
  norm[0] = h[0] / 2
  norm[-1] = h[-1] / 2
  norm[1:-1] = (h[:-1] + h[1:]) / 2
  ab[1, 0] = (h[0] / 3) / norm[0]
  ab[0, 1] = (h[0] / 6) / norm[0]
  ab[1, -1] = (h[-1] / 3) / norm[-1]
  ab[2, -2] = (h[-1] / 6) / norm[-1]
  ab[1, 1:-1] = ((h[:-1] + h[1:]) / 3) / norm[1:-1]
  ab[2, :-2] = (h[:-1] / 6) / norm[1:-1]
  ab[0, 2:] = (h[1:] / 6) / norm[1:-1]
  return solve_banded((1, 1), ab, kbar)


def calculateEfficiencies(kernelFile, lam, ratioInd, r, psd, mr, mi,
                          psdConvention='dndr', nFine=8000):
  """
  PSD-averaged qext, qsca, qabs from a GRASP-kernel netCDF file.

  kernelFile - path to a kernel-*.nc file written by convertkernels.py
  lam        - wavelength [um]
  ratioInd   - index on the file's shape/ratio dimension (e.g. 0 = spheres
               in kernel-grasp-v1.1.3-integrated_V4.nc)
  r          - 1-D radii [um], ascending, arbitrary spacing
  psd        - size distribution at r; normalization is arbitrary since
               efficiencies are ratios
  mr, mi     - refractive index (mi >= 0)
  psdConvention - 'dndr' (default) or 'dndlnr'
  nFine      - points of the fine ln(r) quadrature grid

  Returns dict with 'qext', 'qsca', 'qabs' (dimensionless; total cross
  section per total projected area) and 'area' (projected area per particle
  in PSD units, i.e. sum(pi r^2 dN) / sum(dN)).
  """
  r = np.asarray(r, float)
  psd = np.asarray(psd, float)
  if r.ndim != 1 or r.shape != psd.shape:
    raise ValueError("r and psd must be 1-D arrays of equal length")
  if np.any(np.diff(r) <= 0):
    raise ValueError("r must be strictly ascending")
  if psdConvention == 'dndr':
    dndlnr = psd * r
  elif psdConvention == 'dndlnr':
    dndlnr = psd.copy()
  else:
    raise ValueError("psdConvention must be 'dndr' or 'dndlnr'")

  k = _loadKernels(kernelFile, ratioInd)
  lnxNodes = np.log(k['x'])

  # fine quadrature grid, piecewise-linear dN/dlnr between the input radii
  lnrF = np.linspace(np.log(r[0]), np.log(r[-1]), nFine)
  rF = np.exp(lnrF)
  nF = np.interp(lnrF, np.log(r), dndlnr)
  xF = 2 * np.pi * rF / lam
  if xF[0] < k['x'][0] * 0.999 or xF[-1] > k['x'][-1] * 1.001:
    raise ValueError("size parameter range [%.4g, %.4g] exceeds kernel "
                     "coverage [%.4g, %.4g]" % (xF[0], xF[-1],
                                                k['x'][0], k['x'][-1]))

  out = {}
  volF = 4. / 3. * np.pi * rF ** 3 * nF     # dV/dlnr (PSD units * um3)
  areaTot = np.trapz(np.pi * rF ** 2 * nF, lnrF)
  numTot = np.trapz(nF, lnrF)
  for key in ('ext', 'sca'):
    kbar = _interpMrMi(k[key], k['mr'], k['mi'], mr, mi)
    g = _desmooth(kbar, lnxNodes)
    kFine = np.interp(np.log(xF), lnxNodes, g) * (LAMBDA0 / lam)  # 1/um at lam
    out[key] = np.trapz(kFine * volF, lnrF)  # total cross section [um2 * PSD units]
  out['abs'] = out['ext'] - out['sca']

  return {'qext': out['ext'] / areaTot,
          'qsca': out['sca'] / areaTot,
          'qabs': out['abs'] / areaTot,
          'area': areaTot / numTot}
