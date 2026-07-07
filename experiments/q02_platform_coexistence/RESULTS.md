# Q2 / DR-A -- platform coexistence test (RESULTS)

t-V chain, J = V = 1, half filling; L = 12 (924, 200 real) and L = 14 (3432, 64 real); window central 25%; chi_f w.r.t. the tilt at h = 0.

| L | W | <r> | ln chi_typ | mean chi | Hill (top 5%) | cert |
|---|---|---|---|---|---|---|
| 12 | 0.5 | 0.5275(0.0012) | +10.616 | 1.93e+05 | 2.04 | 8e-14 |
| 12 | 1.0 | 0.5272(0.0012) | +11.120 | 4.73e+05 | 2.01 | 1e-13 |
| 12 | 2.0 | 0.5108(0.0012) | +11.582 | 1.8e+06 | 1.86 | 1e-13 |
| 12 | 3.0 | 0.4496(0.0013) | +10.704 | 8.89e+06 | 1.67 | 1e-13 |
| 12 | 4.0 | 0.4197(0.0013) | +9.153 | 2.21e+07 | 1.55 | 2e-13 |
| 12 | 5.0 | 0.4015(0.0013) | +7.367 | 1.05e+08 | 1.48 | 2e-13 |
| 12 | 6.0 | 0.3960(0.0013) | +5.833 | 1.2e+09 | 1.41 | 2e-13 |
| 12 | 8.0 | 0.3886(0.0013) | +3.423 | 4.84e+06 | 1.41 | 1e-13 |
| 14 | 0.5 | 0.5318(0.0011) | +12.673 | 1.46e+06 | 2.04 | 2e-13 |
| 14 | 1.0 | 0.5327(0.0011) | +13.221 | 2.34e+06 | 2.01 | 2e-13 |
| 14 | 2.0 | 0.5204(0.0011) | +14.220 | 1.04e+07 | 1.96 | 2e-13 |
| 14 | 3.0 | 0.4717(0.0011) | +13.845 | 1.07e+08 | 1.73 | 2e-13 |
| 14 | 4.0 | 0.4190(0.0012) | +11.901 | 2.48e+09 | 1.54 | 3e-13 |
| 14 | 5.0 | 0.4032(0.0012) | +9.750 | 1.12e+09 | 1.45 | 2e-13 |
| 14 | 6.0 | 0.3951(0.0012) | +7.669 | 1.79e+10 | 1.40 | 8e-13 |
| 14 | 8.0 | 0.3864(0.0012) | +4.437 | 6.4e+09 | 1.37 | 5e-13 |

- **R1 (coexistence)**: CHALLENGES-ERGODIC (W* = 2, <r> = 0.5204)
- **R2 (GOE x^-2 tail)**: SUPPORTS (Hill = 1.96 at W*)
- Ergodic control (W = 0.5, L = 14): <r> = 0.5318, Hill = 2.04

**Reading: CHALLENGES the platform minimal model.**

Wall time: 2075.9 s.

## Interpretation note (post-hoc, 2026-07-07)

The frozen R1 gate reads CHALLENGES-ERGODIC by the letter: chi_typ's
argmax (W* = 2, <r> = 0.520) sits one gridpoint on the ergodic side of
the mid-crossover band. The structure underneath is sharper than the
gate:

- chi_typ RISES through the GOE regime toward the crossover (12.67 ->
  14.22 in ln units from W = 0.5 to 2) and the crossover-edge value at
  W = 3 (<r> = 0.472, inside the band) is within 0.38 ln-units of the
  peak -- the peak hugs the crossover's ergodic EDGE, it is not deep in
  the chaotic phase.
- Inside the crossover the tail index falls BELOW 2 (Hill 1.73 at W = 3,
  1.54 at W = 4): infinite-mean territory. The MEAN chi keeps exploding
  into the transition (rare resonances) while the TYPICAL collapses --
  sensitivity claims there are estimator-dependent, and a sensing device
  banks on the typical.
- R2 is a clean positive: Hill = 1.96-2.04 across the GOE regime at both
  sizes -- the predicted GOE x^-2 tail, confirmed.

**Re-scoped platform statement for the Part-1 paper:** maximum TYPICAL
tilt-sensitivity sits at the chaotic edge of the crossover; entering the
intermediate regime trades typical sensitivity for heavy-tailed
rare-event gains (alpha < 2). The minimal model's "large chi_f in the
crossover window" holds for the mean, not the typical.

Registered follow-up (Q2c, design level): does W*(L) TRACK the drifting
crossover position W_c(L) as L grows (a weaker but real coexistence), or
pin to fixed W? Finer W grid around 2-3.5, sizes 12/14/16 (windowed-chi
instrument at 16).
