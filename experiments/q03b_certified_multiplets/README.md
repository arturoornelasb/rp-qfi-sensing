# Q3b -- Sierpinski gasket QFI with a certified instrument (preregistered)

Preregistered 2026-07-07, before its run. Q3-as-run was DEGENERACY-LIMITED
and its largest size was invalidated by an eigensolver failure (see
`../q03_sierpinski_nem/RESULTS.md`). Q3b asks the same physics questions
with a fixed instrument. Readings R1/R2 are UNCHANGED from Q3.

## Instrument fixes (both mandatory)

1. **Certified eigenpairs.** Reference eigenvalues from the validated
   no-vector path (`np.linalg.eigvalsh`); vectors from
   `scipy.linalg.eigh(driver='evr')` (a different LAPACK algorithm than
   numpy's dsyevd). Per-size gates: max |E_evr - E_ref| <= 1e-12 x width,
   and spot residual max ||H u - lambda u||_inf <= 1e-10 x width over 64
   sampled columns. A gate failure marks the size INSTRUMENT-FAIL and it
   is dropped (reported; fits then use the remaining sizes and the run is
   COVERAGE-limited if fewer than 3 remain).
2. **Degenerate-multiplet subspace treatment** (replaces exclusion).
   Exact multiplets = clusters with internal gaps <= 1e-10 x width. For
   each V-draw and each in-window multiplet C: diagonalize the projected
   perturbation P_C V P_C (k x k), rotate the multiplet basis by its
   eigenvectors (the correct zeroth-order basis of degenerate perturbation
   theory), then chi for each rotated state sums |<m|V|n'>|^2/(E_m-E_n')^2
   over ALL states m outside C. Source multiplets need no rotation (the
   summed coupling into a subspace is rotation-invariant). This uses the
   FULL window -- no starvation.
   - D2 (primary): IPR of the V-rotated window states, pooled over draws
     (the rotation selects a physical basis in each multiplet);
     (secondary, reported): non-degenerate states only.

## Design (otherwise as Q3)

Generations g = 5..8, window central 25%, V-ensemble 8 draws (diagonal
primary, GOE/sqrt(N) secondary), cluster sampling capped at 512 window
states per size (whole multiplets, declared), fits over the largest 3
valid sizes.

## Preregistered readings (unchanged from Q3)

- **R1:** enh = alpha_diag - alpha_goe: >= +0.15 ENHANCED; |enh| < 0.15
  NO EFFECT; <= -0.15 SUPPRESSED (Anderson-like sign).
- **R2 (exploratory):** |enh - (1 - D2_meas)| <= 0.15 -> the RP formula
  EXTENDS to geometry-only fractality; else DOES NOT EXTEND.
- Ungated: certification numbers per size, multiplet-size distribution,
  D2 secondary.

## Budget

~25-40 min thread-capped (<= 4); peak RAM ~5 GB at g = 8.
