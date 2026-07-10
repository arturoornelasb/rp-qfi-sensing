# Q4 -- which sensed operators inherit the D2 enhancement? (preregistered)

Frozen 2026-07-09, before any run. Q1b established the chi-enhancement
exponent 1 - D2 for a FIELD-LIKE (diagonal) sensed operator, with a
dense (GOE) operator D2-BLIND (alpha ~ 1 flat) -- the scoping finding.
Q4 maps the middle: banded (short-range) and low-rank operators, the
structures real sensed parameters take (local fields, contact
interactions, collective couplings). This is the "which observables
inherit the advantage" table the platform paper needs.

## Design

gRP ensemble at interior gamma in {1.4, 1.8}; N in {512..4096}
(realizations 48/24/12/6); mid-window 25%; for each realization chi_f
per window state for FOUR operator structures:
- diag: V = diag(v), v ~ N(0,1)  [the confirmed field-like case]
- banded: V_ij ~ N(0,1) for |i-j| <= 8, symmetrized, else 0
- lowrank: V = sum_{k=1..4} v_k v_k^T / sqrt(4N)  [collective]
- dense: V = GOE/sqrt(N)  [the blind reference]
Exponents alpha_V over the largest 3 sizes; enhancement
enh_V = alpha_V - alpha_dense; D2 measured from IPR as in Q1.

## Frozen reading (per structure, per gamma)

- INHERITS: |enh_V - (1 - D2_meas)| <= 0.15 (like diag).
- BLIND: |enh_V| < 0.15 (like dense).
- PARTIAL: between -- report the fraction of the full enhancement.
The deliverable is the classification table (any outcome is the
result); consistency across the two gammas required to label a row.

## Budget

~20-30 min thread-capped.
