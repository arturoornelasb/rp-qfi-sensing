# Q3c — stratified/energy-resolved gasket null (preregistered)

Frozen 2026-07-12, before any run. The Q3b follow-up its close-out
registered: Q3b's NO EFFECT verdict (enh = −0.123) carried the caveat
that giant multiplets and regular states mix size-dependently
(composition scatter). Q3c decomposes the statistic so the null (or a
hidden effect) is attributed to a population, not an average.

## Design (Q3b instrument, unchanged)

Same certified eigensolver (eigvalsh reference + `evr` vectors +
agreement/residual gates; INSTRUMENT-FAIL drops the size), same
degenerate-multiplet subspace rotation, same generations g = 5..8 and
8 V-draws. New: every (chi, IPR) sample is tagged by

- **class stratum**: `reg` (non-degenerate states) vs `mult`
  (multiplet members), inside the mid window; and
- **energy stratum**: four bulk windows centered at spectrum fractions
  {0.30, 0.43, 0.57, 0.70}, width 0.10 each (pooled over classes),
  capped at 256 chi-states per window per size.

Per stratum: enhancement exponent alpha_d − alpha_g over the largest
3 sizes, and measured D2 from the IPR ladder of the same states.

## Frozen readings

- **R1 — STRATIFIED-NULL-CONFIRMED** if every stratum has
  |enh| < 0.15: the Q3b null is not a composition artifact; geometric
  fractality gives no enhancement in any population. (Sharpens the
  paper's criterion sentence.)
- **R2 — HIDDEN-STRATUM** if any stratum has |enh| ≥ 0.15: report the
  sign (negative = the Anderson-like SUPPRESSED branch of report #6's
  open-Q3) and whether it is monotone across sizes.
- **R3 (exploratory, ungated):** does the `reg` stratum's deviation
  from 1 − D2 improve on Q3b's pooled 0.816?

Either outcome is a finding; R1 feeds the paper's scope section as a
strengthened null, R2 would be a new result on its own.

## Budget

One certified solve per size (the g = 8, N = 9843 solve is the pole),
then stratified chi/IPR accounting: ~2–5 h total, CPU-only, runs
alongside Q2c-L16 overnight.
