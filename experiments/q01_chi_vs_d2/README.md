# Q1 / DR-B -- chi_f scaling vs D2 at fixed criticality distance (preregistered)

Preregistered 2026-07-07, before any run. Source claim (deep-research
report #6, from Skvortsov-Amini-Kravtsov PRB 106, 054208 (2022)): in the
RP non-ergodic-extended phase the fidelity susceptibility is enhanced with
**sensitivity exponent ~ 1 - D2**. The referee crux: show the enhancement
is a function of fractality (D2) at fixed distance from criticality, not a
critical-point artifact.

## Ensemble and observables

- Generalized RP: H = diag(eps) + N^{-gamma/2} W, eps_i ~ N(0,1), W GOE
  (offdiag var 1, diag var 2). Phases: gamma < 1 ergodic; 1 < gamma < 2
  fractal NEE with D2 -> 2 - gamma; gamma > 2 localized.
- gamma grid: interior {1.2, 1.4, 1.6, 1.8} (each >= 0.2 from both
  edges = "fixed distance from the transition"); controls {0.5, 2.5}.
- Sizes N = {256, 512, 1024, 2048, 4096}; realizations
  {96, 48, 24, 12, 6}; eigenstate window = central 25% of the spectrum.
- chi_f per window eigenstate n: sum_{m != n} |V_mn|^2 / (E_m - E_n)^2,
  with (a) PRIMARY: V = diag(v), v_i ~ N(0,1) (field-like sensed
  parameter); (b) SECONDARY: V = GOE / sqrt(N) (norm-O(1) dense
  perturbation).
- Typical value chi_typ = exp(mean ln chi) (heavy tails make the mean
  unstable); alpha(gamma) = slope of ln chi_typ vs ln N over the three
  largest sizes.
- D2 MEASURED per gamma: slope of -mean ln IPR vs ln N over the same
  sizes (not assumed 2 - gamma).

## Preregistered readings

- **R1 (figure of merit):** at the interior gammas,
  |alpha_diag - (1 - D2_meas)| <= 0.15 for >= 3 of 4 points -> SUPPORTS
  the D2 figure of merit. Monotone tracking with a systematic offset
  > 0.15 -> PARTIAL (report the calibration). No tracking -> CHALLENGES.
- **R2 (fractality, not criticality):** alpha_diag varies monotonically
  across the interior grid with total change >= 0.3 (theory: 1 - D2 goes
  0.2 -> 0.8). SUPPORTS overall requires R1 AND R2.
- Controls (reported, ungated): ergodic (0.5) and localized (2.5)
  endpoints; the GOE-V secondary exponents; the P(chi) tail exponent at
  the largest N per gamma (context for the platform paper's x^-2 GOE-tail
  must-prove).

## Budget

Pure NumPy eigh + windowed matrix elements; ~30-60 min thread-capped
(OMP/MKL <= 4) so it coexists with the plates-rp-fem queue.
