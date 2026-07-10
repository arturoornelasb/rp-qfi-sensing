# Research plan and experiment registry

Target: the two-part RP<->QFI paper (platform + D2 figure of merit) defined
by deep-research reports #5 (classical->quantum bridge: CONSTRUCT, on an
interacting tilted+disordered lattice, chi_f(h) = QFI/4) and #6
(multifractality-enhanced sensing: Kravtsov's exact RP result -> sensitivity
exponent ~ 1 - D2; must isolate fractality from criticality). Reports live
in the dossier repo
(`lopez-gonzalez-research-analysis/deep_research_interferometria_atomica/`).

Policy: preregister readings before execution; freeze executed experiments;
one commit per unit; runs are serial with the plates-rp-fem queue when the
machine is CPU-bound.

## Q1 -- DR-B: chi_f scaling vs D2 at fixed criticality distance [PREREGISTERED 2026-07-07]

The load-bearing claim of Part 2, and the referee crux report #6 flagged:
show the chi enhancement is a function of the FRACTAL DIMENSION along the
RP non-ergodic-extended phase, not a critical-point effect. Generalized-RP
ensemble H = diag(eps) + N^{-gamma/2} W; interior points gamma in
{1.2, 1.4, 1.6, 1.8} (>= 0.2 from both phase edges), controls at 0.5
(ergodic) and 2.5 (localized); N = 256..4096; chi_f per mid-spectrum
eigenstate for a diagonal (field-like) perturbation (primary) and a
normalized GOE perturbation (secondary); D2 MEASURED from IPR scaling (not
assumed 2 - gamma). Frozen readings in
`experiments/q01_chi_vs_d2/README.md` (R1: alpha tracks 1 - D2_meas; R2:
alpha varies monotonically inside the phase). Budget: ~30-60 min,
thread-capped.

**[DONE 2026-07-07 -- literal PARTIAL; physics CONFIRMS the corrected
normalization.]** D2_meas tracks theory to 0.03-0.05 across the phase
(0.854/0.608/0.397/0.168 vs 0.8/0.6/0.4/0.2); alpha_diag = 2 - D2_meas
(dev 0.03-0.06) with alpha_goe ~ 1.00 flat (a dense perturbation is
D2-BLIND -- the in-situ ergodic reference), so the enhancement
alpha_diag - alpha_goe = 1 - D2_meas within 0.06 at 4/4 interior points:
the SAK exponent, in-house, on measured D2. Literal R1 (raw alpha vs
1 - D2) reads 0/4 -- report #6's shorthand names the ENHANCEMENT
exponent. Scoping finding for the paper: the figure of merit requires
FIELD-LIKE (local) coupling to the sensed parameter. Ergodic control
alpha = 0.500 is the gRP bandwidth artifact (chi ~ N^gamma, gamma < 1),
quantitatively understood.

## Q1b -- Enhancement-exponent gate, fresh seed [PREREGISTERED 2026-07-07]

Converts Q1's post-hoc reconciliation into a frozen-criteria verdict:
identical design, independent seed (8), gate R1' =
|(alpha_diag - alpha_goe) - (1 - D2_meas)| <= 0.15 at >= 3/4 interior
gammas (R2 unchanged). Readings frozen in
`experiments/q01b_enhancement_gate/README.md` BEFORE its run (gate chosen
from Q1's seed-7 data; Q1b's draws are independent).

**[DONE 2026-07-07 -- SUPPORTS, frozen criteria.]** R1' 4/4 within 0.15
(deviations 0.015/0.018/0.009/0.078; corr(enh, 1-D2_meas) = 0.996);
R2 monotone, delta 0.709. The D2 sensitivity figure of merit --
chi-enhancement exponent = 1 - D2 along the RP fractal phase at fixed
criticality distance, for field-like couplings -- is CONFIRMED in-house
under preregistered criteria. This is the Part-2 paper's core numerical
result; next legs: Q2 (platform coexistence) and Q3 (fractality without
disorder).

## Q2 -- DR-A: the platform coexistence test [PREREGISTERED 2026-07-07; runner ready]

Report #5 open-Q2, the decisive test of the Part-1 minimal model: does an
interacting disordered lattice show its Poisson->GOE crossover in the SAME
parameter window where the fidelity susceptibility to the tilt is large --
and does the chi_f distribution carry the predicted GOE x^-2 heavy tail?
Instantiation: t-V chain (spinless fermions, J = V = 1) + onsite disorder
W + tilt as the sensed parameter at h = 0 (the canonical interacting
ergodic<->MBL crossover; Bose-Hubbard/AAH variants registered as Q2b
alternates). L = 12 (dim 924, 200 realizations) and L = 14 (dim 3432, 64),
W in {0.5..8}; first realization of every cell certified vs eigvalsh;
L = 16 registered as a windowed-chi extension. Frozen readings in
`experiments/q02_platform_coexistence/README.md`: R1 = chi_typ's argmax
W* has mid-crossover statistics (<r> in [0.42, 0.50]) vs
CHALLENGES-ERGODIC / -LOCALIZED; R2 = Hill tail index in [1.6, 2.4] at W*.
Confirms or kills the platform paper's minimal model.

**[DONE 2026-07-07 -- CHALLENGES (by the frozen gate); re-scoped finding
recorded.]** Textbook crossover (<r> 0.532 -> 0.386); chi_typ peaks at
the crossover's ergodic EDGE (W* = 2, <r> = 0.520; the in-band W = 3
value is within 0.38 ln-units), NOT mid-crossover -> R1
CHALLENGES-ERGODIC. R2 SUPPORTS: Hill = 1.96-2.04 across the GOE regime
(the predicted x^-2 tail). Key structure: inside the crossover the tail
index drops below 2 (1.73/1.54) -- mean chi explodes on rare resonances
while typical collapses; sensitivity there is estimator-dependent.
Platform statement re-scoped: max TYPICAL sensitivity at the chaotic
edge; the crossover window offers heavy-tailed rare-event gains only.
Q2c registered (design level): does W*(L) track W_c(L) with size, or pin
to fixed W (finer grid 2-3.5, sizes 12/14/16 windowed).

## Q3 -- DR-C: fractality without disorder (Sierpinski NEM) [PREREGISTERED 2026-07-07; runner ready]

Report #6 open-Q2: no metrological quantity has ever been computed in the
geometry-only Sierpinski non-ergodic-extended system. Tight-binding gasket
g = 5..8 (N up to 9843), zero disorder, V-ensemble statistics (8 draws),
exact-degeneracy exclusion (the gasket spectrum is degeneracy-heavy;
DEGENERACY-LIMITED reported if usable states fall below 64). Frozen
readings in `experiments/q03_sierpinski_nem/README.md`: R1 = sign of the
enhancement (ENHANCED / NO EFFECT / SUPPRESSED at margin 0.15 -- the
suppressed branch is the Anderson-like sign of report #6 open-Q3); R2 =
exploratory, does the RP formula enh = 1 - D2 extend to geometric
fractality. Either outcome is a finding.

**[DONE-AS-RUN 2026-07-07 -- DEGENERACY-LIMITED; no physics conclusions.]**
Two instrument defects diagnosed by direct test: (1) np.linalg.eigh WITH
vectors returns eigenvalues off by up to 31 (width-6 spectrum) at
N = 9843 -- the g = 8 cell was garbage; validated clean at N <= 6144 on
gRP matrices, so Q1/Q1b are unaffected; (2) the degeneracy-EXCLUSION
design starves statistics (84-95% of window states are in exact
multiplets). Both fixed in Q3b. PITFALL recorded: eigenpair certification
is mandatory in this repo (same discipline as plates-rp-fem).

## Q3b -- Sierpinski QFI, certified instrument [PREREGISTERED 2026-07-07; runner ready]

Same physics questions and R1/R2 gates as Q3, fixed instrument:
certified eigenpairs (eigvalsh reference + scipy driver='evr' vectors;
per-size agreement <= 1e-12 x width and residual <= 1e-10 x width gates;
failures drop the size as INSTRUMENT-FAIL) and degenerate-multiplet
SUBSPACE treatment (rotate each in-window multiplet by the projected V's
eigenbasis -- correct degenerate perturbation theory -- then chi over
states outside the multiplet; full-window statistics, no exclusion).
Readings frozen in `experiments/q03b_certified_multiplets/README.md`.

**[DONE 2026-07-07 -- NO EFFECT; DOES NOT EXTEND.]** All four sizes
CERTIFIED (evr driver passes at N = 9843, agreement 6.2e-14 -- the numpy
dsyevd failure is bypassed); full-window statistics recovered
(105/319/541/738 states). Frozen verdict: enh = -0.123 (|enh| < 0.15 ->
NO EFFECT); dev vs 1 - D2 = 0.816 (DOES NOT EXTEND). Caveat: exponents
carry composition scatter (giant multiplets vs regular states mix
size-dependently; ln chi non-monotone) -- the robust statement is the
absence of consistent enhancement. **Combined Q1b + Q3b finding: the
1 - D2 enhancement is RP-PHASE-SPECIFIC, not generic fractality** --
answers report #6 open-Q3's criterion question; the plate program's
RP-like regime is the metrologically relevant one. Q3c (optional,
design level): stratified/energy-resolved gasket analysis if the paper
needs this leg sharpened.

## Q4 -- Operator-structure scope law [DONE 2026-07-09]

Frozen classification delivered: diag INHERITS (both gammas); banded
(range 8) CONDITIONAL -- 58% at gamma = 1.4, full at 1.8, i.e. a
short-range operator inherits once its range <= the fractal support
scale; lowrank (collective) BLIND; dense BLIND. Platform design rule:
SENSE WITH LOCAL COUPLINGS. (Tilt/position = field-like: the Q2 choice
was right by construction.)

## Q5 -- PBRM third family [DONE 2026-07-09 -- EXTENDS, 3/3]

Critical PBRM at b = 0.5/1.5/4.0: enh = 1 - D2_meas within 0.016-0.064
at all three bandwidths. **CRITERION REVISED:** with Q1b (RP yes) and
Q3b (gasket no), the enhancement law is generic to MULTIFRACTAL STATES
OF RANDOM ENSEMBLES -- what the gasket lacked was randomness (exact
deterministic multiplets), not RP miniband structure specifically. This
supersedes the Q3b-era phrasing "RP-phase-specific"; the Part-2 paper's
criterion: multifractal support + random local statistics + field-like
coupling.

## Deferred / watch

- Report #6 open-Q3 (sign criterion: RP enhances, 3D-Anderson suppresses --
  where does the plate sit?) -- theory question, revisit after Q1 data.
- Report #6 open-Q4 (many-body D2 inheritance in the Q2 platform) -- attach
  to Q2 if its crossover window is found.
- P3/P4/P5/P6 of the deep-research queue stay in the dossier repo
  (deep-research API runs, currently paused).
