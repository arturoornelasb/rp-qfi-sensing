# Q3b -- certified multiplet instrument (RESULTS)

Valid sizes: [5, 6, 7, 8] (fit over [6, 7, 8]); V-ensemble 8 draws; full-window multiplet treatment.

| g | N | states/multiplets (max) | ln chi_d | ln chi_g | ln IPR | cert agree |
|---|---|---|---|---|---|---|
| 5 | 366 | 105/31 (40) | +5.178 | +5.306 | -3.841 | 1.2e-14 |
| 6 | 1095 | 319/66 (121) | +6.442 | +6.604 | -4.287 | 2.4e-14 |
| 7 | 3282 | 541/26 (364) | +6.883 | +7.403 | -3.795 | 3.9e-14 |
| 8 | 9843 | 738/51 (363) | +3.017 | +3.450 | -4.962 | 6.2e-14 |

- alpha_diag = -1.560; alpha_goe = -1.437; enh = -0.123; D2_meas = 0.308; 1 - D2_meas = 0.692; dev = 0.816
- **R1**: NO EFFECT (margin 0.15)
- **R2**: DOES NOT EXTEND (report the measured pair) (tol 0.15)

**Reading: NO EFFECT; DOES NOT EXTEND (report the measured pair).**

Wall time: 524.3 s.

## Interpretation note (post-hoc, 2026-07-07)

- The instrument worked: all four sizes CERTIFIED (agreement 1.2e-14 to
  6.2e-14) -- scipy driver='evr' handles N = 9843 where numpy's dsyevd
  returned garbage, and the multiplet subspace treatment recovers
  full-window statistics (105/319/541/738 states vs 11/23/43/-- under
  exclusion).
- Caveat on the exponents: ln chi and ln IPR are NON-MONOTONE across
  sizes because the gasket window mixes giant exact multiplets (up to
  ~364 sampled states from one ladder) with regular states in
  size-varying proportions -- the fits carry composition scatter, so the
  robust statement is the ABSENCE of a consistent enhancement, not a
  precise exponent pair.
- **Combined Q1b + Q3b statement (the finding):** the chi-enhancement
  law with exponent 1 - D2 is a property of the RP non-ergodic-extended
  phase (Q1b: SUPPORTS, 4/4 within 0.06), and does NOT transfer naively
  to geometry-only fractality (Q3b: enh = -0.12, NO EFFECT, on one
  degeneracy-heavy lattice family). This addresses report #6's open
  question 3 (the enhancement criterion): what enhances is the RP
  phase's structure, not fractal dimension per se -- and the plate
  program's RP-like regime is therefore the metrologically relevant one.
- Optional refinement if the paper needs the gasket leg sharpened
  (registered design-level as Q3c): stratified analysis -- singleton vs
  multiplet states separately, energy-resolved windows pinned away from
  the giant ladders.
