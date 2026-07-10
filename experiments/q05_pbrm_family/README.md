# Q5 -- does the 1 - D2 law extend to PBRM? (preregistered)

Frozen 2026-07-09, before any run. The enhancement-criterion
classification so far: RP non-ergodic-extended YES (Q1b, 4/4); geometric
gasket NO (Q3b). PBRM (critical power-law banded random matrices, a = 1)
is the third family: multifractal eigenstates with D2 tunable by the
bandwidth b, but a DIFFERENT miniband structure than RP. If the law
extends, "1 - D2" is about multifractality per se; if not, it is
RP-miniband-specific -- either answer completes the criterion the
Part-2 paper must state.

## Design

Critical PBRM: H_ij ~ N(0, [1 + (|i-j|/b)^2]^{-1}) symmetrized, at
b in {0.5, 1.5, 4.0} (D2_meas expected ~0.3 / ~0.55 / ~0.8);
N in {512..4096} (48/24/12/6 realizations); mid-window 25%; chi_f for
the diagonal (field-like) primary and dense (GOE/sqrt(N)) reference;
D2 measured from IPR scaling; enh = alpha_diag - alpha_dense.

## Frozen reading

- EXTENDS: |enh - (1 - D2_meas)| <= 0.15 at >= 2 of 3 bandwidths.
- DOES NOT EXTEND: otherwise (report the measured pairs) -- the
  enhancement is then RP-miniband-specific, sharpening the criterion.

## Budget

~20-30 min thread-capped.
