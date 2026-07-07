# Q1b -- the enhancement-exponent gate, fresh seed (preregistered)

Preregistered 2026-07-07, before its run; frozen AFTER Q1's data taught us
the correct normalization but on an INDEPENDENT seed, so the gate is
clean. Q1 (seed 7) found: alpha_diag = 2 - D2_meas, alpha_goe ~ 1.00 flat
(D2-blind dense reference), so the report-#6 / Skvortsov-Amini-Kravtsov
"sensitivity exponent ~ 1 - D2" is the ENHANCEMENT exponent
alpha_diag - alpha_goe.

## Design

Identical ensemble, sizes, realizations, window, and fits as Q1
(`../q01_chi_vs_d2/README.md`), with seed = 8 (independent draws).

## Preregistered readings

- **R1' (enhancement gate):** at the interior gammas {1.2, 1.4, 1.6, 1.8},
  |(alpha_diag - alpha_goe) - (1 - D2_meas)| <= 0.15 for >= 3 of 4
  points -> SUPPORTS the D2 figure of merit; monotone tracking with
  offset > tol -> PARTIAL; else CHALLENGES.
- **R2 (unchanged):** alpha_diag monotone across the interior grid, total
  change >= 0.3.
- SUPPORTS overall requires R1' AND R2.
- Controls (reported, ungated): ergodic 0.5 / localized 2.5 endpoints;
  P(chi) upper tail at largest N.

## Budget

~8-10 min thread-capped (Q1 took 499 s).
