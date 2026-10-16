# graspkernels: GRASP text kernels for spheres from Mishchenko's `spher.f`

Builds GRASP-format sphere kernels (`kernels_<ratio>_00.txt` and `_11`, `_12`, `_22`,
`_33`, `_34`, `_44`) on any radius grid and set of refractive indices. Its first use is
the sphere kernel for the TAMUdust2020 hexahedra grid (273 radii), replacing the
interpolated sphere kernel shipped with the hexahedra, which is 6% too high and invalid
above x ≈ 600.

## Contents

| File | Role |
|---|---|
| `spher.f` | M. I. Mishchenko's Lorenz–Mie code for polydisperse spheres (version of 08/06/2005) as a subroutine library. Changes from upstream are listed in its header: the demo main program was moved out, the reference output was moved to `tests/`, and `DISTRB` gained `NDISTR=7`, a size distribution that is triangular in ln r in volume. |
| `spher_main.f` | Upstream demo main program, verbatim. With `spher.f` it rebuilds the original program. |
| `spher_tri.f` | Driver: loops over refractive indices and radius nodes, calls `SPHER` once per half of each hat function, and evaluates the scattering matrix at the requested angles (`MATR_ANG`, copied from `MATR` in `gsf/spher_expan.f`). |
| `makekernels.py` | `fun(ifn, dest)`: JSON config → parallel `spher_tri` jobs → GRASP text files, `raw_kernels.npz` and a matching `kernelconversion` JSON. |
| `graspkerneltext.py` | GRASP text-kernel writer (layout of the GRASP v1.1.3 `KERNELS_BASE` files; header and block strings follow Greema Regmi's `GRASPhexKernel.py`) and small readers. |
| `tests/` | Validation scripts (below). |

Build with `compile_all.sh` in the repository root. Run with
```
python runkernelgeneration.py --filename input_files/kernelgeneration/<config>.json --dest <output dir>
```
Configs: `grasp_sphere_base41.json` rebuilds the standard GRASP sphere kernel (validation);
`grasp_sphere_tamu273.json` builds the sphere kernel on the TAMUdust2020 hexahedra grid.
The JSON keys are documented at the top of `makekernels.py`. Each refractive index is one
job, and finished jobs are skipped on a rerun, so an interrupted build resumes where it
stopped.

## Rescaling kernels (`scalekernels.py`)

`scalekernels.fun(ifn, dest)` (run with `python runkernelscaling.py --filename <config>.json --dest <dir>`) multiplies EXTINCTION, ABSORPTION and every scattering-matrix element of a kernel set by one factor. It rewrites the files with `graspkerneltext.write_kernels`, keeps the ratio label, grid and angles, and writes a provenance file `<fnpre>_<ratio>_scaling.json` (source checksums, factor, commit). `graspkerneltext.read_kernels` reads any of the layouts in use, KERNELS_BASE or one row per radius.

Its first use (October 2026) corrects the TAMUdust2020 hexahedra kernels. They had been built as Q_A·π·r_v²/V instead of C/V = Q_A·A/V and were low by the sphericity:
- configs `input_files/kernelgeneration/scale_TAMUdust2020_hex_V5.json` and `scale_TAMUdust2020_hex_varAngRes_V5.json`;
- α = A/(π·r_v²) = 1.42726 (Ψ = 0.70064, the ensemble TAMUdust2020 delivers at nominal S = 0.7);
- written with 8 significant digits like the source files (`digits`), so that ext − abs keeps its precision.

The sphere kernel was rewritten with 8 digits from its `raw_kernels.npz` (`grasp_sphere_tamu273_8digits.json`, `reuseRecompute = "none"`).

## GRASP positivity (`grasppositivity.py`)

The standard GRASP distribution stores element 33 as log(K33) for angles ≤ 40° and element 44 as log(K44) for angles ≤ 50° (`dls_fixget_bin.f90`), and interpolates those logs in refractive index. Size-resolved sphere kernels have negative node values there near Mie resonances: on the TAMUdust2020 grid, 1,479 K33 and 1,950 K44 values, at x ≈ 5–19 and 15–50°. GRASP then stops with a NaN.

`graspPositivity: true` in a generation config applies `grasppositivity.fix` to the text files only; `raw_kernels.npz` keeps the exact values. Around each negative run it takes the smallest window of size nodes with a positive sum (median 3, at most 7 nodes) and sets F/F11 there to its K11-weighted mean. Σ K over the window is conserved, so a size distribution that is flat across the window is unaffected. Maximum change in size-integrated F33/F11 and F44/F11 for lognormal volume modes, over all refractive indices and angles: 0.12 (σ = 0.1), 0.019 (0.2), 0.008 (0.3), 0.004 (0.5), 0.003 (0.7).

The standard GRASP v1.1.3 sphere kernel holds 1000.0 at exactly the 23 places (all at m_i = 0.5) where its true K33 (≤ 40°) or K44 (≤ 50°) is negative.

Config: `grasp_sphere_tamu273_8digits_grasppositive.json`. netCDF files should be built from the exact values (`convertkernels` `ratioPaths`).

## What is computed

For node k with hat function T_k(u), u = ln r (0 at the neighbouring nodes, 1 at r_k):

    ext_k  = ∫ T_k (C_ext/V) du         abs_k = ext_k − ∫ T_k (C_sca/V) du
    P_ij,k = ∫ T_k (C_sca/V) F_ij du,   with ½ ∫ F11 sin θ dθ = 1

in 1/µm. This is the GRASP kernel definition. It already contains the node spacing (∫T_k du = dlnr),
so no dlnr rescaling is involved. Each half hat is one `SPHER` call: with the triangle as the
volume distribution, ⟨C⟩/⟨V⟩ from `SPHER` is the T-weighted mean of C/V over that half,
and the node value is Σ_halves (h/2)·⟨C⟩/⟨V⟩. Splitting at the apex keeps the kink out of
the Gaussian quadrature.

Conventions, all checked against the GRASP v1.1.3 kernels:
- F22 = F11 and F44 = F33 (spheres); they are computed independently and written as such.
- Element 34 in GRASP files is −b2 in Mishchenko's notation, i.e. Im(S2 S1*) with
  Wiscombe's (MIEV0) amplitudes. Element 12 is b1 (negative in the Rayleigh limit).
- `edgeNodes = "grasp"`: the first and last node use only the inner half hat, with weight h
  instead of h/2. This is what the standard GRASP kernels do. It explains the last node of
  the GRASP sphere kernel being ~9% above a symmetric triangle average at x = 626.
  `"half"` (weight h/2, exact for a piecewise-linear size distribution ending at the grid
  edge) and `"full"` (hat extended one step beyond the grid) are also available.
  GRASP sums node kernels with equal weights (`OPTCHAR_bin`) and log-interpolates them onto
  the retrieval grid (`MATRIX_FIX_bin`), so an edge node needs the full-hat weight. The
  `"grasp"` choice approximates the missing outer half by the inner half. This
  was presumably needed because the spheroid database ends at the grid edge. It is
  biased by about (2h/3)·d ln(C/V)/d ln r: +9% at x = 626 on the 41-node grid, about +1.5%
  at the last node of the 273-node grid. For spheres the outer half can simply be computed,
  so the TAMUdust2020-grid kernel uses `"full"`, which is also closer to the hexahedra
  kernels (point values × dlnr, symmetric about each node). With `reuseFrom`, an existing
  build can be switched to another edge treatment by recomputing only the two end nodes.

## Quadrature

Accuracy is limited by Mie resonances, whose width grows with absorption and roughly in
proportion to x. The points needed are therefore roughly a fixed number per half hat
(`quadPtsPerHalf`), fewer for more absorbing particles (`quadMiRef`, `quadMiMinFactor`),
and more where needed (`quadRefine`).
Convergence tests (half width 0.0463 in ln r, relative to much denser references):

- m_i ≥ 1e-3: 100 points per half converge ext, sca and P11 to ≤ 1e-6 for x ≤ 200.
- m_i = 1e-4: 400 points: ≤ 1e-5; 200 points: P11 median 3e-4, 95th percentile ≤ 2e-3.
- m_i = 1e-10: 1600 points still leave P11 errors of ~1e-3 (95th percentile).
- Large x (370–2800): ext and sca converge to ≤ 2e-4 even at 0.5–1 point per unit x. P11
  needs ~4 points per unit x for a median error ~1e-4 at m_i ≤ 1e-4.

Cost is about 1.15e-8 · x² s per quadrature point (the 2x Gauss angles `SPHER` uses to build
the expansion), so the largest nodes dominate the run time.

## Tests

| Script | Checks |
|---|---|
| `tests/test_spher_upstream.py` | `spher_main` reproduces the 366 numbers of Mishchenko's reference output exactly. |
| `tests/test_ndistr6_vs_7.py` | `NDISTR=7` agrees with a construction that uses only upstream `NDISTR=6` (difference of two broad log-normals, linear in ln r). |
| `tests/test_writer_roundtrip.py` | Reading and re-writing the GRASP v1.1.3 kernels reproduces all eight files line by line. |
| `tests/test_independent_mie.py` | Kernels recomputed with Wiscombe's MIEV0 and with `pymiecoated`, using an independent Gauss–Legendre quadrature in Python. |
| `tests/compare_kernels.py` | Rebuilt kernels vs. a reference set (e.g. the GRASP v1.1.3 kernels). |
| `tests/convergence_spotcheck.py` | Reruns selected refractive indices of a build with more quadrature points. |

## Validation results (2026-10-15)

- `test_spher_upstream.py`: 366 of 366 numbers identical.
- `test_ndistr6_vs_7.py`: max relative difference 6.5e-12.
- `test_writer_roundtrip.py`: all eight `KERNELS_BASE` files (grid1.dat and kernels_100_*)
  reproduced line by line.
- `test_independent_mie.py`: m = 1.5+0.01i, 1.37+0.1i, 1.6+0.001i, x = 0.5–600: ext, sca
  and all six elements agree with MIEV0 to ≤ 6.4e-9 (and with `pymiecoated`, x ≤ 30, to
  ≤ 3e-9).
- Rebuild of the GRASP v1.1.3 sphere kernel (`grasp_sphere_base41.json`; 41 nodes, 22 × 16
  refractive indices, `edgeNodes = "grasp"`, 40 min on 4 cores), compared with
  `compare_kernels.py --skip-angles 13`: median |difference| ~1e-6 (the 7-digit print
  precision) for ext, abs and all elements at every refractive index, interior and edge
  nodes, including the "9% high" last node. Larger maxima occur only for weak absorption,
  where narrow resonances limit both quadratures: m_i = 5e-4–1.3e-3, ≤ 4.6e-5 in ext and
  ≤ 1.6e-2 in P11 at single angles; m_i = 1e-10, ≤ 3.7e-4 in ext and ≤ 7e-2 in P11. The
  m_i = 1e-10 absorption is resonance-dominated and not converged in either kernel. The GEOSmie netCDF made from the
  rebuild with `runkernelconversion.py` agrees with `kernel-grasp-v1.1.3-integrated_V4.nc`
  to a median of 4e-7 in ext and qext, and to ≤ 4e-4 in g for x ≥ 1.

- Production kernel on the TAMUdust2020 hexahedra grid (`grasp_sphere_tamu273.json`; 273
  nodes, 8 × 31 refractive indices, `edgeNodes = "full"`). All values are finite and
  positive. P22 = P11 and P44 = P33 to 2e-5. ½∫K11 sin θ dθ = K_sca to 6e-6 for x < 10. The
  x grid and refractive indices are identical to the hexahedra kernels.
  - How it was built: the base schedule on 4 cores (~5 h), then `quadRefine` (8× the points
    for m_i ≤ 2e-3, all nodes) recomputed via `reuseFrom` on 8 cores (~8.5 h). A from-scratch
    run of the JSON gives the same values.
  - Convergence, from `convergence_spotcheck.py` with 2× the final points (i.e. 16× the base
    schedule in the refined region; every 6th node above x = 30 plus the last, and all
    nodes x ≤ 30). Changes:
    - m_i = 1e-4 (m_r = 1.37, 1.67): ext and sca ≤ 5e-8 everywhere. P11 ≤ 6e-8 for x ≤ 1000;
      for x > 1000 its 95th percentile is ≤ 1.6e-6 and its maximum 2.8e-5.
    - m_i = 4e-4 and 2e-3: ext and sca ≤ 3e-12; P11 ≤ 2e-11 for x ≤ 1000 and ≤ 1.7e-7 above.
    - m_i ≥ 3e-3 (not refined, 2× check of the base schedule): ext and sca ≤ 1e-12; P11
      ≤ 1e-8 for x < 1000 and ≤ 8e-6 above.
  - What the refinement corrected (base vs final, m_i ≤ 2e-3): ext and sca ≤ 1.4e-3 and P11
    sin-weighted ≤ 1.1% at x > 1000; at 300 < x ≤ 1000, P11 95th percentile ≤ 1.2e-3; below
    x = 300, P11 ≤ 2.5e-3 at single angles.
- Against the old interpolated sphere kernel shipped with the hexahedra (2024 folder):
  old/new = 1.0604 in absorption and 1.0605 in extinction at x < 0.1 (the dlnr error).
  Old P11 is ~12% high at small x, i.e. 1.0604 × 1.057, the bias of the coarse (0.27 in
  ln r) triangle smoothing for Rayleigh scattering. Above x = 583 the old values are 20
  to 1500 times too large.

## Notes found along the way

- The 13° column of the GRASP v1.1.3 sphere kernel is not computed for 15 of its 16
  imaginary indices (all except m_i = 1e-10). For every real index and radius node it is
  exactly a two-point interpolation of the stored 12° and 14° values: linear in ln F for
  F11, F22, F33 and F44 (the geometric mean), and linear for F12 and F34. P11 there differs from the
  true 13° value by a median of 5e-4–9e-4 and up to ~17%. The m_i = 1e-10 column is
  computed and agrees with this code to 1.7e-6, as at all other angles. Among the
  single-axis-ratio DLS spheroid kernels only the ratio-1 (sphere) member shows the
  interpolation, so it reaches the spheroid mixture (`kernels_299_*`) only through that
  member. The ln/linear split matches the angle-interpolation convention in GRASP's
  `dls_fixget_bin.f90`.
- The GRASP v1.1.3 sphere kernel contains 23 fill values of exactly 1000.0 (2 in
  `kernels_100_33.txt`, 21 in `kernels_100_44.txt`), all at m_i = 0.5, x = 3.6–4.7 and
  scattering angles 40–50°, where F33 and F44 cross zero. The correct values are ~1e-3.
- `pymiecoated` (GEOSmie's Python Mie code) agrees with MIEV0 to ~1e-9 for x ≤ 50, but its
  Q_ext and Q_sca are off by up to ~5e-4 for 100 < x < 3000 at low absorption. Its
  amplitudes S1, S2 are the complex conjugates of MIEV0's, so its Im(S2 S1*) has the
  opposite sign.
