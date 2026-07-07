# Q2 / DR-A -- the platform coexistence test (preregistered)

Preregistered 2026-07-07, before any run. Deep-research report #5's
decisive open question (Q2 there): does an interacting disordered lattice
show its Poisson->GOE crossover **in the same parameter window where the
fidelity susceptibility to a tilt is large** -- and does the chi_f
distribution carry the predicted GOE x^-2 heavy tail? This confirms or
kills the minimal model of the Part-1 (platform) paper: RP structure from
interactions+disorder, tilt = the sensed force only, chi_f(h) = QFI/4.

## Model and observables

- Spinless fermions, t-V chain with onsite disorder and a tilt:
  H = -J sum (c+_i c_{i+1} + h.c.) + V sum n_i n_{i+1} + sum w_i n_i
      + h sum (i - c) n_i,   J = V = 1, w_i ~ U[-W, W], half filling.
  The canonical interacting Poisson<->GOE (ergodic<->MBL) crossover;
  disorder breaks all lattice symmetries (no sector bookkeeping); the
  tilt enters ONLY as the sensed parameter, evaluated at the working
  point h = 0. (The literal tilted Bose-Hubbard / interacting-AAH
  variants of the source preprints are registered alternates, Q2b, if
  the paper needs them; the coexistence question is model-generic.)
- Disorder grid W in {0.5, 1, 2, 3, 4, 5, 6, 8}; sizes L = 12 (dim 924,
  200 realizations) and L = 14 (dim 3432, 64 realizations).
- Per realization: full spectrum + eigenvectors (dims well below the
  validated 6144 eigh ceiling; the FIRST realization of every (L, W)
  cell is additionally certified against eigvalsh at 1e-12 x width --
  a cell that fails is INSTRUMENT-FAIL and excluded).
- Mid-spectrum window: central 25%. Per (L, W): pooled r-tilde over
  window spacings; chi_f per window eigenstate for the centered dipole
  V_op = sum (i - (L-1)/2) n_i; chi_typ = exp(mean ln chi); Hill tail
  index alpha_H of the pooled chi distribution (top 5%).

## Preregistered readings

- **R1 (coexistence -- the decisive gate):** at L = 14, let
  W* = grid-argmax of chi_typ. SUPPORTS if the level statistics at W*
  are mid-crossover: pooled r-tilde(W*) in [0.42, 0.50].
  CHALLENGES-ERGODIC if r-tilde(W*) >= 0.52 (sensitivity peaks deep in
  the chaotic phase); CHALLENGES-LOCALIZED if <= 0.40. L = 12 reported
  for size consistency (ungated).
- **R2 (the GOE x^-2 tail):** at W*, Hill alpha_H in [1.6, 2.4] ->
  SUPPORTS the predicted x^-2 tail (density exponent ~2); outside ->
  report the measured index. Ergodic control (W = 0.5) reported
  alongside, ungated.
- Ungated: the full chi_typ(W) and r-tilde(W) curves (the paper's
  sensitivity-map figure), mean-vs-typical chi (tail dominance), size
  drift 12 -> 14.
- Registered extension (not run today): L = 16 (dim 12,870) with a
  windowed-chi instrument (energy-truncated sums validated against the
  full sum at L = 14) and certified evr eigenpairs.

## Budget

~1-1.5 h thread-capped (OMP/MKL <= 4), ~1 GB RAM; coexists with E17.
