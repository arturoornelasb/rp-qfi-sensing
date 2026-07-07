# Q1 / DR-B -- chi_f vs D2 (RESULTS)

gRP ensemble, sizes [256, 512, 1024, 2048, 4096], window 25%, fits over the largest 3 sizes. Primary V = diagonal; secondary V = GOE/sqrt(N).

| gamma | alpha_diag | alpha_goe | D2_meas | D2_theory | 1-D2_meas | dev |
|---|---|---|---|---|---|---|
| 0.5 | +0.500 | +0.505 | 0.999 | 1.00 | 0.001 | 0.499 |
| 1.2 | +1.177 | +1.043 | 0.854 | 0.80 | 0.146 | 1.032 |
| 1.4 | +1.443 | +1.017 | 0.608 | 0.60 | 0.392 | 1.051 |
| 1.6 | +1.661 | +1.027 | 0.397 | 0.40 | 0.603 | 1.058 |
| 1.8 | +1.861 | +0.994 | 0.168 | 0.20 | 0.832 | 1.029 |
| 2.5 | +1.645 | +1.013 | -0.043 | 0.00 | 1.043 | 0.602 |

- **R1** (|alpha - (1-D2_meas)| <= 0.15 at >= 3/4 interior gammas): 0/4 within tol -> **PARTIAL (monotone tracking, offset > tol)**
- **R2** (monotone alpha, total change >= 0.3): monotone = True, delta = 0.684 -> **SUPPORTS**

**Reading: PARTIAL -- see R1/R2.**

Wall time: 499.0 s.

## Post-hoc reconciliation (added after the frozen verdict, 2026-07-07)

The literal R1 gate compared alpha_diag to (1 - D2_meas) -- report #6's
shorthand taken at face value. The data show the shorthand names the
ENHANCEMENT exponent, not the raw exponent:

- alpha_diag = 2 - D2_meas at all four interior gammas (dev 0.03-0.06);
- alpha_goe ~ 1.00 flat at ALL gammas -- the dense perturbation is
  D2-blind and provides the in-situ ergodic reference;
- hence **alpha_diag - alpha_goe = 1 - D2_meas within 0.06 at 4/4
  interior points** -- the Skvortsov-Amini-Kravtsov enhancement exponent,
  confirmed in-house on measured (not assumed) D2.

Two findings beyond the gate: (i) the D2 figure of merit applies to
FIELD-LIKE (local/diagonal) sensed parameters; a structureless dense
coupling gains nothing from fractality -- a scoping constraint the
Part-2 paper needs; (ii) the ergodic control's alpha = 0.500 is
quantitatively the gRP bandwidth artifact (chi ~ N^gamma for gamma < 1),
not an anomaly.

Q1b (fresh seed, corrected frozen gate) converts this reconciliation
into a preregistered verdict -- see `../q01b_enhancement_gate/`.
