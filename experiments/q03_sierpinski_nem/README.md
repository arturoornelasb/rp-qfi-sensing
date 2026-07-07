# Q3 / DR-C -- QFI on the Sierpinski gasket: fractality without disorder (preregistered)

Preregistered 2026-07-07, before any run. Deep-research report #6, open
question 2: in the geometry-only Sierpinski non-ergodic-extended system,
**no metrological quantity has ever been computed** -- is chi_f enhanced,
and does it track D2? This is the decisive isolation leg: fractal
eigenstates with ZERO disorder, so any enhancement cannot be a
disorder-criticality effect.

## System and observables

- Tight-binding H on the Sierpinski gasket graph, hopping t = 1, no
  onsite terms, no disorder; generations g = 5..8
  (N = 366, 1095, 3282, 9843), single deterministic lattice per size.
- Statistics come from a V-ENSEMBLE (the lattice is fixed): 8 independent
  draws per size of (a) PRIMARY V = diag(v), v_i ~ N(0,1) and
  (b) SECONDARY V = GOE/sqrt(N), pooling ln chi over window states x
  draws. chi_f as in Q1 (`sum |V_mn|^2 / (E_m - E_n)^2`).
- Window: central 25% of the spectrum. The gasket spectrum carries EXACT
  macroscopic degeneracies; window states with a neighbor gap
  < 1e-10 x spectral width are EXCLUDED from both chi and D2 (IPR in a
  degenerate subspace is basis-arbitrary); if that leaves fewer than 64
  usable states at any size the run reports DEGENERACY-LIMITED. Up to 512
  usable window states per size enter chi (declared cap).
- D2_meas: slope of -mean ln IPR (non-degenerate window states) vs ln N
  over the largest 3 sizes; alpha_diag / alpha_goe: same-size fits of
  ln chi_typ.

## Preregistered readings

- **R1 (does geometric fractality enhance chi?):**
  enh = alpha_diag - alpha_goe over the largest 3 sizes:
  enh >= +0.15 -> ENHANCED; |enh| < 0.15 -> NO EFFECT;
  enh <= -0.15 -> SUPPRESSED (the 3D-Anderson-like sign -- directly
  relevant to report #6's open question 3 about the sign criterion).
- **R2 (exploratory -- does the RP formula extend?):**
  |enh - (1 - D2_meas)| <= 0.15 -> the RP figure of merit EXTENDS to
  geometry-only fractality; otherwise DOES NOT EXTEND (report the
  measured pair; either outcome is a finding -- nobody has computed
  this).
- Reported ungated: D2_meas, degeneracy fractions per size, usable-state
  counts, spectral width.

## Budget

~20-35 min thread-capped (OMP/MKL <= 4); peak RAM ~3 GB (dense eigh at
N = 9843). Runs alongside E17 (which is I/O-light and 8-thread-capped).
