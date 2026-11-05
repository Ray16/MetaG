# Validation against the evaluated release

The [fixed-pose solvation rotation control](../benchmark/solvation_rotation_validation/README.md)
retains 426 successful pairs and six gas-SCC refusals across all 432 requested
pairs. Dense SASA quadrature reduces orientation dependence in the complete
subcases, while the charged-oxidized energetic endpoint remains withheld.
The original complete-panel gate is not weakened. This numerical improvement
does not establish biochemical Gibbs accuracy or change the production profile.

The [released-ensemble reliability comparison](../benchmark/native_uncertainty_validation/README.md)
now records the native 100-member reaction outputs on the fixed 69-target
development cohort. Its spreads preserve shared-member molecular covariance;
native means are evaluated separately from archived OOF means. Sixty-eight exact
normalized state equations also occur in the shipped training table. The
descriptive risk curves do not validate independent risk reduction or a deployed
reliability rule, and retain the preceding reference-condition qualifications.

The subsequent [minimum-curvature recovery and global-reference study](../benchmark/global_amino_reference_validation/README.md)
adds a numerical recovery path after the final force-converged FIRE/mode-following
attempt still has a strong imaginary mode. The new path requires tighter forces
and a fresh Hessian; it does not soften the vibrational criterion. A complete
alanine parent calculation validates this path, while the original refusal remains
retained. Raw-species source/cache identity changes, and no new full biochemical
sweep or accuracy result is inherited from the earlier release.

That study's molecule-local amino-acid reference prototype is a separate,
unadopted candidate. Its common-target error worsens despite numerical cycle
closure and increased availability from the recorded bank. Reproducing those
results without a quantum runtime does not establish an accurate production
state/reference policy or independent generalization.

The later calculation-condition guard refuses explicit requests outside the
implemented thermodynamic target before molecular work in both APIs. The packaged
364 inputs omit those requests. A bounded source comparison removes only the
new admission guards and docstrings and recovers the preceding runtime ASTs;
see [the compatibility receipt](../benchmark/gpdh_exchange_reference_review/compatibility_audit.json).
This changes the full source identity, without a new physical model or molecular
sweep. The guard does not normalize experimental reference conditions.

The final 31 October submission package additionally rejects atom-map annotations in input SMILES
before routing. This prevents annotated water from bypassing the liquid-water reference. All 364
sweep inputs lack those annotations and pass the new guard; removing exactly that guard and the
module-documentation correction leaves identical runtime abstract syntax trees. The final package
and the actual molecular sweep have distinct full source fingerprints, with unchanged raw-species
physics. See `benchmark/final_review_20261031/` for the source-change proof and original files.
The molecular sweep has not been relabelled as a new run under the later full fingerprint.

The post-audit corrected source has now completed a 364-input evaluation with 270 estimates,
92 workflow refusals and two time-budget outcomes. On the primary reference, 198 estimates
give MAE 8.38 kJ/mol; the preceding source gives 8.37 on the same targets. Coverage fell by
45 previously available primary labels, so the smaller aggregate error is not improvement.
The default benchmark scorer now uses this explicitly fingerprinted corrected-source panel.

Shared-record controls reproduce two equivalent-drawing predictions and one invariant refusal.
A coherent-source audit evaluates 30 of 49 independent nullspace vectors; four do not close,
with residuals up to 8.70 kJ/mol. The default estimator is therefore not established as a
global network potential. Physical admission controls and the new cache identity remain active;
no historical accepted artifact is restored to improve coverage. The earlier validations below
remain separately identified historical records.

## Source identity and conformer-search repeats

The current checkout also provides a separate experimental
[full-molecule state-function profile](state-function-profile.md). Its consistency controls and
accuracy records do not replace the frozen reduction-based benchmark above. Later manuscript
reference/scaffold candidates have their own source and outcome identities. A package version or
unchanged raw molecular engine alone does not establish equivalence of complete estimators.

`seed_offset` shifts the seed supplied to each conformer-generation batch. Batch `j` uses
`j + seed_offset`, starting at `j = 1`. Offsets 0 and 1 therefore reuse seed batches and must not be
described as independent sampling runs. For nonoverlapping batches, choose offsets separated by at
least `conv_max` and record the actual batch seeds. For example, with `conv_max = 8`, offsets 0 and
10000 use disjoint possible ranges 1–8 and 10001–10008. Disjoint seeds alone do not prove complete
gas/aqueous sampling or calibrate total prediction error.

Retain every requested repeat and its failures. Exceeding `thermal_validation_max` is a validation
budget refusal, not evidence that a chemical minimum cannot exist. A higher ceiling defines a
different resource/settings profile and cache identity; do not relabel accepted lower-budget cache
records, discard unresolved representatives, or report a partial ensemble as complete. Report any
additional effort separately from a change to the energy model or minimum-acceptance thresholds.

## Archived refactor and guarded-release validations

MetaG 0.3.0 is a production-only refactor of the development release 0.2.0.dev24 plus four corrections.
The refactor is checked for equivalence first; the corrections are then measured reaction by reaction. Switches for options that were tested and rejected were removed;
the default physics was kept. Equivalence of the refactor (before the corrections) was checked on 2026-10-27 at three levels.

**Reaction level (CPU, exact).** Every reaction of the packaged 364-reaction panel plus the examples
(368 inputs) was scored by both implementations with the same deterministic replay backend. The backend
returns a fixed pseudo-energy for each (charge, canonical SMILES). This exercises the reaction-level
paths invoked by these inputs, including routing, source-backed protonation, hydration microstates, sugar pools, assembly and fail-closed guards. Three variants
were run: the reactions as written, every reaction reversed, and a variant in which some charged species
fail gas-phase identity checks, which exercises the neutral-reference re-route. All 3 × 368 outcomes are
identical, including the unrounded assembled energies, which match bit for bit. The replay detects a
0.01 kJ/mol change in one constant.

**Species level (GPU).** Free energies of methanol, acetic acid, phosphoric acid, tetramethylammonium (+1),
triplet O₂, glycerol and liquid water were computed from empty caches on a V100. They agree with
0.2.0.dev24 to at most 1.2 × 10⁻⁴ kJ/mol. Two dev24 runs of the same species differ by up to
1.7 × 10⁻⁴ kJ/mol through GPU floating-point nondeterminism. Water and O₂ are bit-identical.

**Whole reactions (GPU).** Four panel reactions were scored end to end from cold caches by both
implementations: pyrophosphate hydrolysis (rxn00001), a hydrated carbonyl (rxn00541), hydration plus
truncation (rxn01117), and truncation with a net proton (rxn45677). They give the same ΔrG'°
(−31.5, 13.0, −10.1, 13.1 kJ/mol), identical stage decompositions to 0.01 kJ/mol, the same scored species
and the same wall time.

**Output schema.** The following fields were removed: the historical class-σ interval block (`sigma_pred`,
`ci95`, `ci_info`, `sigma_breakdown`, `sign_status`, `unresolved`; never coverage-calibrated in
production), the always-null `anchor`/`smd`/`water_ref` records, the `explicit` flag, and the route
flags of removed options. `dG`, `dG_raw`, `stages`, `U_samp`, `routes`, `species_scored` and `suspect`
are unchanged. Species-cache keys differ from dev24 by construction (new source identity), so earlier
caches are never reused.

## Corrections in 0.3.0

The corrections are listed in the [CHANGELOG](../CHANGELOG.md). Their footprint was determined on CPU with the
same replay as above:

- direction-invariant routing changes the route of 13 reactions as written; before it, 16 of 368
  reactions violated ΔG(reverse) = −ΔG(forward), afterwards none do;
- the redox-proton guard changes 1 reaction;
- the deduplication fallback changes 14 species, which occur in 10 reactions;
- failing degraded species changes no reaction of the panel.

These 24 reactions were recomputed with 0.3.0 from empty caches on V100 GPUs. Three unchanged control
reactions recomputed the same way reproduce their dev24 values to within 0.1 kJ/mol (rxn00001 -31.5 → -31.5, rxn01329 +8.8 → +8.8, rxn44278 -6.5 → -6.4).
rxn00824 exceeded the 5400 s per-reaction budget in dev24 and remains a time-budget failure. Values are
in kJ/mol:

| Reaction | Correction | dev24 | 0.3.0 |
|---|---|---|---|
| rxn00191 | orientation | +0.8 | +1.6 |
| rxn00260 | orientation | +10.0 | +10.9 |
| rxn00493 | orientation | +9.9 | +10.8 |
| rxn00527 | orientation | +9.9 | +10.7 |
| rxn00753 | dedup | -23.0 | -21.8 |
| rxn00754 | dedup | -22.0 | -23.2 |
| rxn00806 | orientation | +1.2 | +0.9 |
| rxn00824 | dedup | — | — |
| rxn00849 | orientation | -4.6 | -4.1 |
| rxn00906 | orientation | +7.4 | +9.5 |
| rxn00907 | orientation | +7.4 | +9.5 |
| rxn01201 | orientation | +22.0 | +24.1 |
| rxn01336 | dedup | +3.7 | +3.9 |
| rxn01380 | dedup | +3.6 | +3.6 |
| rxn01409 | dedup | -24.1 | -21.4 |
| rxn01532 | dedup | +21.8 | +21.1 |
| rxn01702 | dedup | -40.5 | -40.5 |
| rxn01726 | dedup | -23.8 | -20.4 |
| rxn01999 | orientation | +8.0 | +8.8 |
| rxn02320 | orientation | +3.8 | +4.7 |
| rxn02902 | redox guard | +26.0 | +40.7 |
| rxn02914 | orientation | +8.4 | +9.3 |
| rxn03281 | dedup | -5.8 | -7.1 |
| rxn16141 | orientation | +22.0 | +24.1 |

On the Mg-bounded benchmark (primary policy) the MAE changes from 8.44 to 8.42 kJ/mol on the same 239
labels. `benchmark/metag_0.3.0_predictions.json` holds the combined panel predictions.

### Guarded-release source comparison (29 October 2026)

The 364-reaction rerun used full estimator fingerprint
`6662bffb1b34a2dba2f2a58792ba54dae56bdf246bc1523f81ecbfb7f288eb54`.
A later instrumented repeat of initial default routing found cancelled MCS
searches in two of the 245 scored reactions; the historical executions did not
record cancellation flags. The working package now rejects partial MCS atom
maps, so it is a different estimator from the frozen benchmark source. The
guard was checked on both affected routes and with the full CPU suite. A
label-free cache-only replay across all 364 panel reactions yielded 325
estimates, 35 workflow refusals and four unresolved cache dependencies; all
325 estimates match the earlier frozen results at reported precision. On the
primary reference, 243 estimates give MAE 8.71 kJ/mol; 29 are refusals and two
remain unresolved. This retrospective replay does not assess fresh molecular
sampling or independent accuracy; a guarded-source selection audit remains
outstanding.

The earlier comparison record [current_source_equivalence_20261029.json](../benchmark/history/current_source_equivalence_20261029.json)
applies to the pre-guard working source only. The [guarded source comparison](../benchmark/history/mcs_guard_source_comparison_20261029.json)
records the executable difference in `routing/truncate.py`. Use
`python benchmark/compare_source_trees.py /path/to/evaluated/metag src/metag --output /tmp/source-check.json`
to compare the actual source trees. The frozen evaluated tree must accompany
the study archive; AST comparison is a source-level diagnostic, not a proof
of runtime equivalence.

The guarded replay and M2 study use a frozen source whose package marker is
`0.3.0`; at this checkpoint the unreleased working source used `0.3.1.dev0` to avoid
conflating it with the already released 0.3.0. The [version-marker comparison](../benchmark/history/guarded_version_marker_comparison_20261029.json)
records that all 55 package files matched the guarded evaluated tree byte-for-byte
except the literal version assignment in `__init__.py`, and that the raw-species
fingerprint is unchanged. This source-level finding does not prove runtime or
metadata equivalence. The frozen source and runtime records remain necessary
for reproduction.
