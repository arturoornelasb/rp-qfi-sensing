# Q3 / DR-C -- Sierpinski gasket QFI (RESULTS)

Geometry-only fractality, V-ensemble (8 draws), fits over g = [6, 7, 8] (largest 3 sizes).

| g | N | usable/window | degen frac | ln chi_d | ln chi_g | ln IPR |
|---|---|---|---|---|---|---|
| 5 | 366 | 11/91 | 0.88 | +8.590 | +7.919 | -4.732 |
| 6 | 1095 | 23/274 | 0.92 | +12.218 | +10.861 | -5.556 |
| 7 | 3282 | 43/821 | 0.95 | +13.946 | +12.895 | -6.494 |
| 8 | 9843 | 512/2460 | 0.00 | +4.082 | +4.152 | -6.649 |

- alpha_diag = -3.706; alpha_goe = -3.056; enh = -0.650; D2_meas = 0.498; 1 - D2_meas = 0.502; dev = 1.152
- **R1**: DEGENERACY-LIMITED (usable states below floor at a fit size) (margin 0.15)
- **R2**: DOES NOT EXTEND (report the measured pair) (tol 0.15)

**Reading: DEGENERACY-LIMITED (usable states below floor at a fit size); DOES NOT EXTEND (report the measured pair).**

Wall time: 114.8 s.

## Post-run instrument diagnosis (2026-07-07) -- NO physics conclusions from this run

Two instrument defects, both diagnosed by direct test (scratchpad
diagnostics, reproduced independently):

1. **The g = 8 cell is INVALID.** `np.linalg.eigh` (with vectors) at
   N = 9843 returns eigenvalues off by up to **31** on a width-6 spectrum
   (vs `eigvalsh`, which matches the g = 7 cross-check to 4e-14). The
   true g = 8 window has 89.3% exact-zero gaps (clusters of 683 at
   E = 1.0), fully consistent with g = 7 -- the runner's "0.00 degeneracy"
   and the chi collapse were computed from garbage eigenpairs. Validated
   clean on gRP matrices at N = 2048/4096/6144 (agreement ~ 8e-14,
   residuals ~ 1e-14) -> **Q1/Q1b are unaffected**.
2. **The exclusion design starves the statistics.** 84-95% of window
   states sit in exact multiplets at every g, so excluding them leaves
   11/23/43 states at g = 5/6/7 -- below the preregistered floor by
   construction, not by bad luck.

Fixes are frozen in `../q03b_certified_multiplets/` (certified eigenpairs
+ degenerate-multiplet subspace treatment). Lesson recorded: eigenpair
certification is mandatory in this repo too, exactly as in plates-rp-fem.
