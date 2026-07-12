# rp-qfi-sensing (private)

In-house de-risking experiments for the **RP ↔ QFI sensitivity-map program**
— the two-part paper defined by deep-research reports #5 and #6 (dossier
repo `the program's private research record, deep_research_interferometria_atomica/`):

- **Part 1 (platform):** interacting tilted+disordered optical lattice;
  observable = fidelity susceptibility chi_f(h) = QFI/4; RP structure
  sourced from interactions/disorder (attaches to Bloch/lattice gravimetry,
  not free-fall readout).
- **Part 2 (figure of merit):** the D2 sensitivity exponent — in the RP
  non-ergodic-extended phase, chi is enhanced with exponent ~ (1 - D2)
  (Skvortsov-Amini-Kravtsov, PRB 106, 054208 (2022)); the plate's
  eigenvector diagnostic becomes metrological (D2 ~ 0.76 -> exponent
  ~ 0.24).

This repo runs the **directly simulable tests** those reports queued,
before any atomic experiment — the same de-risking pattern that settled
Gap A (E14) and the mechanical GOE->GUE (E15) in `plates-rp-fem`.

## Experiments

- `experiments/q01_chi_vs_d2/` — DR-B: chi_f scaling vs D2 at fixed
  distance from the transition, generalized-RP ensemble (PREREGISTERED).
- DR-A (tilted+disordered interacting lattice: Poisson->GOE window vs
  chi_f) and DR-C (Sierpinski NEM: QFI under fractality-without-disorder)
  are registered at design level in `docs/RESEARCH_PLAN.md`.

## Discipline

Same as `plates-rp-fem`: each experiment folder carries a preregistered
README (readings frozen BEFORE execution), runners save raw JSON, executed
experiments are FROZEN, one commit per unit of work.

## Environment

Pure NumPy (MKL-backed); run with the `plates-fem` conda env's python
(a conda environment with numpy/scipy, e.g. `plates-fem`) or any
numpy >= 1.26.
