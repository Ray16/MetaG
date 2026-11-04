# Benchmark: Mg-bounded openTECR reference

`reference_mg_bounded.json` holds experimental ΔrG'° labels transformed to pH 7,
I = 0 and pMg 14 for the reactions of the packaged panel
(`src/metag/data/reactions_opentecr_std.json`). They were derived from the openTECR re-curation of TECRDB
by a Legendre transform of every measurement at its reported temperature with the eQuilibrator dissociation
model. A measurement with unknown free Mg²⁺ is admitted only when its label changes by at most 2 kJ/mol
over an assumed pMg interval. Those intervals are assumptions about unreported conditions, not measured
or universally established physical bounds. The rule was recorded before residuals on this reconstructed reference were computed
([PREREGISTRATION.md](PREREGISTRATION.md)). Thresholds 1 and 4 kJ/mol are reported as sensitivities.
These labels are not normalized to 298.15 K, and pMg 14 is a finite approximation to no Mg binding.

```bash
python benchmark/score_checked.py                              # completed corrected-source panel
python benchmark/score_checked.py metag_out/metag_all.json     # a fresh `metag score` run
```

`metag_0.3.0_predictions.json` holds the 0.3.0 panel predictions (320 estimates for 364 inputs). With the
primary policy they give MAE 8.42 kJ/mol on 239 of 274 labels; the zero predictor gives 11.68 on the
same labels. They combine the 0.2.0.dev24 outcomes for reactions the 0.3.0 corrections do not touch
with fresh 0.3.0 calculations for the 24 that they change; see [docs/validation.md](../docs/validation.md).
`metag_dev24_predictions.json` keeps the dev24 outcomes (MAE 8.44) for reference.

Caveats, as reported in the manuscript:

- Discrete method choices were made while looking at TECRDB residuals, so these scores are
  development-exposed rather than independent validation. No parameter fitted to reaction free energies
  enters an estimate. The [development-exposure record](development_exposure.md) describes the
  historical decision audit and the remaining current-version selection gap.
- Labels are medians over measurements at their native temperatures (273–348 K). Many reactions have no
  measurement at 298.15 K.
- Ionic strength is imputed as 0.25 M where it was not reported.
- Fitted comparators (component contribution, dGbyG, dGPredictor) have to be compared by regime: in-sample,
  out-of-fold, or compound held out. Their predictions are not redistributed here.

Data provenance: openTECR (Wagner et al.) and the original TECRDB measurements. Cite them when you use
these labels.

## Reconstructing the reference

`rebuild_mg_bounded.py` recovers the historical label construction. It requires the original frozen
openTECR CSVs, ModelSEED snapshot, prediction-structure file, and the eQuilibrator compound database;
these source datasets are not bundled. `reference_snapshot_provenance.json` records their original hashes
and software versions. With those inputs and the `reference` extra installed, run:

```bash
python benchmark/rebuild_mg_bounded.py --reactions /path/to/prediction_structures.json \
  --source /path/to/snapshot --modelseed /path/to/modelseed \
  --output /tmp/reference.json --audit /tmp/reference_audit.json
```

The audit records the selected measurements, source-row numbers, assumptions, source hashes and current
transformer provenance. It also records each excluded measurement that mapped to the panel, its reason,
mapping candidates and orientation. Rows that could not be mapped to the panel are outside these exclusion
counts and are not listed as panel exclusions.
Compare input and database hashes with the frozen provenance before claiming exact reproduction. This
script reproduces the historical policy; it does not make its assumptions experimentally verified.
The released tables were reproduced exactly using the archived source snapshot and compound database;
all recorded source, reaction-input and transformer hashes matched. The Mg sensitivity is evaluated on
an 11-point grid, so its threshold is a grid-based test, not a proof of a continuous-interval maximum.

An archived out-of-fold dGbyG comparison needs reaction-scale and ionic-strength corrections and an
exact-family identity check. `score_dgbyg_oof.py` defaults to the retained identity metadata and
completed corrected-source predictions; with the complete reconstructed 164-row pairing table it reproduces the manuscript's
69-pair intersection and its 7.79 versus 4.59 kJ/mol MAEs. For the preceding 84-pair result, pass
`--predictions benchmark/metag_guarded_cache_predictions.json` explicitly. Pass `--predictions
benchmark/metag_0.3.0_predictions.json` explicitly for the historical combined-release comparison. See the
[comparison audit](comparison_audit.md) for the calculation, the broader exploratory set and limitations.
The audit also documents `select_dgbyg_cv.py` and `build_dgbyg_pairs.py`, which reconstruct the
available pairs from the authors' archived CV outputs and the retained mapping, with pinned source
and prediction fingerprints. The archive is external; the reconstruction does not independently
verify the original fold assignments or propagate condition-bridge uncertainty.
An exploratory analysis of training-compound and reaction-change similarity is in
[domain audit](domain_audit.md); it does not establish a validated MetaG-favored domain.
The [cycle audit](cycle_audit.md) records development-panel reactions that expose nonadditivity.

## Checking staged manuscript evidence

If the manuscript and historical study sources are staged together, run
`python benchmark/check_manuscript_source_manifest.py /path/to/study/manuscript`.
This strict check requires every original source listed in the manuscript's eight historical
evidence manifests to be present and match its SHA-256. It complements the manuscript's
`tools/check_evidence.py`, which checks macros and outcomes but permits historical sources to be
absent. The strict check covers only those eight source manifests; it does not verify all paper claims,
new analyses, dataset redistribution permissions or the completeness of the final release.

## Frozen-source all-input rerun

`metag_frozen_source_predictions.json` contains the label-free outputs of a separate
364-reaction rerun of the frozen estimator. The original timeouts were retained in the
study record and each had one longer retry; under that rule, 327 reactions have estimates
and 37 have workflow refusals. This is a separate run from the combined 0.3.0 panel above.

```bash
python benchmark/score_checked.py benchmark/metag_frozen_source_predictions.json
# Optional validated entry point: rejects duplicate rows and invalid energies.
python benchmark/score_checked.py benchmark/metag_frozen_source_predictions.json
# Requires matplotlib; writes PDF and PNG beside the output prefix.
python benchmark/plot_frozen_benchmark.py --output-prefix output/fresh_benchmark
```

The original scorer is retained byte-for-byte for the frozen study protocol; the
checked entry point uses the same metric function and returns the same results
for valid inputs.

Against the 2 kJ/mol Mg-bounded reference, the rerun scores 245 of 274 labels with MAE
8.68 kJ/mol; the zero predictor is 11.72 kJ/mol on those same 245. The plotting
command uses the same reference and predictions and plots these 245 pairs.
Rendering can vary with Matplotlib and installed fonts. On the 239 targets
shared with the earlier release, the rerun and earlier release have MAEs of 8.37 and
8.42 kJ/mol respectively. These are descriptive, development-exposed results, not an
independent test or evidence of a meaningful performance improvement. The exact evaluated
source predates the current MCS timeout guard and has a different estimator
fingerprint; see [validation](../docs/validation.md) for the source comparison
and the limits of the earlier benchmark. The original outcomes, retries and source must be deposited in
the final study archive for independent verification.

## Guarded cache-only replay

The evaluated guarded source predates the working package version marker
`0.3.1.dev0`. The version-only source comparison is documented in
[validation](../docs/validation.md); it does not establish runtime equivalence.

`metag_guarded_cache_predictions.json` contains all 364 guarded replay outcomes:
325 estimates, 35 workflow refusals and four unresolved cache dependencies.
No experimental labels were used in the replay. On the primary reference,
243 estimates have MAE 8.71 kJ/mol versus 11.70 for zero; 29 labelled reactions
are refused and two remain unresolved. The 84 matched dGbyG out-of-fold
comparisons have MAEs of 9.07 for MetaG and 4.27 for dGbyG, without fully
propagated bridge uncertainty. All 325 available replay estimates equal the
prior frozen-source results at reported precision. This is a retrospective
cache replay, not new conformer sampling or independent validation. The full
M2 selection audit is complete but retrospective; independently held-out biochemical validation
and the public release of the underlying evidence remain outstanding.

```bash
python benchmark/score_checked.py benchmark/metag_guarded_cache_predictions.json
```

## Completed corrected-source evaluation

`metag_corrected_source_predictions.json` contains all 364 outcomes and identifies the frozen
post-audit source: 270 estimates, 92 workflow refusals and two 5400-s budget outcomes.
The validated scorer and OOF comparison default to this file. The original `score.py` remains
unchanged for explicit historical replay. `corrected_source_metrics.json` reports all Mg policies,
common-target comparisons and the preceding error distribution of newly unavailable reactions.
On the 198 admitted primary labels, current/preceding common-target MAEs are 8.38/8.37 kJ/mol.
No accuracy improvement follows from the smaller aggregate error after coverage falls.

`corrected_cycle_summary.json` records the current-source shared-record cycle audit: 30 of 49
independent input-state vectors are evaluable, of which four fail numerical closure up to
8.70 kJ/mol. The additionally named sucrose and amino-acid controls are not independent basis
vectors. Full raw molecular records, source, unrounded assembly and control-arm evidence remain
in the separate study archive and manuscript evidence, pending public deposition. Keep all
unavailable outcomes; missing cycle values are never zero or proof of thermodynamic closure.
