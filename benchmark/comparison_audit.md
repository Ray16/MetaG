# Audit of an archived dGbyG comparison

This is an exploratory comparison on the development-exposed Mg-bounded benchmark. It is not an
independent test of MetaG and it does not establish out-of-distribution performance. dGbyG's
out-of-fold status applies to the original dGbyG reaction split, not necessarily to its constituent
compounds, related reactions, or the current reference construction. The condition bridge and
reference use eQuilibrator dissociation information; errors may therefore be correlated.

The manuscript's primary matched comparison additionally requires an audited exact compound-family
equation match. After completing the existing physical bridge metadata uniformly across the 164
archived matches, the guarded replay has **84 exact-family reactions**, with MAEs of **9.07 kJ/mol
for MetaG**, **4.27 kJ/mol for dGbyG**, and **8.21 kJ/mol for zero**. Thirty other available pairs have
unresolved identity or stereochemical transfer. The previous 148-row reconstruction omitted 16
records whose bridge fields had not been copied into the compact identity table; four of those
are now eligible and increase the earlier 80-reaction cohort. Metadata completion used existing
identity, signed-extent and condition-bridge records, without reference errors, but was undertaken
after outcome availability had been observed. The 109-pair results below use the earlier
combined-release predictions and are retained as a broader historical exploratory analysis.
If a newly available prediction would require a pair whose audited extent or bridge is absent, the
command stops rather than scoring it with unverified metadata; that pair needs a fresh identity and
condition audit.

The archived `metag_vs_dgbyg_pairs.tsv` was produced by
`tools/thermo_benchmark/paper_reproduction/compare_metag.py` in the parent project. Its `dgbyg` value
is in the dGbyG paper's reaction orientation and extent, at I = 0.25 M. Its `scale` equals paper
extent divided by MetaG extent, and `bridge_delta_I0_minus_I025` is the eQuilibrator correction for
the MetaG reaction, in MetaG's extent. The comparable I = 0 prediction is therefore

```
dGbyG_MetaG_orientation_I0 = dgbyg / scale + bridge_delta_I0_minus_I025
```

The later archived `analysis/jcim_revision_20261028/score.py` instead compares the unmodified
`dgbyg` field to the MetaG-oriented I = 0 reference. This omits the ionic-strength correction and
the orientation/extent conversion. Its reported dGbyG OOF MAE of 6.28 kJ/mol on 109 reactions should
not be used. Ten of those 109 reactions have negative scale. For the earlier combined-release MetaG
predictions, the corrected broader calculation gives:

| Method | n | MAE (kJ/mol) | Median absolute error (kJ/mol) | RMSE (kJ/mol) |
|---|---:|---:|---:|---:|
| MetaG 0.3.0 | 109 | 8.16 | 5.80 | 12.29 |
| dGbyG OOF, harmonized | 109 | 4.64 | 2.76 | 7.68 |
| Zero | 109 | 7.57 | 5.84 | 9.99 |

These numbers use the same 109 labeled reactions with estimates from both methods. They are based on
an archived pairing and prediction ensemble that have not yet received an independent audit of
reaction identity, cross-validation assignments, leakage, and condition-transform applicability. The
correction does not establish that those remaining assumptions are valid. Prediction and reference
uncertainties are not included.

## Exploratory complementarity

MetaG is closer to the reference on 35 of 109 reactions and dGbyG on 74. MetaG's advantage exceeds
5 kJ/mol on 8 reactions, while dGbyG's advantage exceeds 5 kJ/mol on 42. These are descriptive
within-sample counts; they do not provide a rule for selecting a method for a new reaction.

| Enzyme class | n | MetaG MAE | dGbyG MAE | MetaG closer |
|---|---:|---:|---:|---:|
| Transferases (EC 2) | 38 | 7.72 | 4.19 | 11 |
| Hydrolases (EC 3) | 17 | 12.84 | 7.20 | 5 |
| Lyases (EC 4) | 25 | 9.24 | 5.34 | 8 |
| Isomerases (EC 5) | 26 | 4.71 | 2.99 | 10 |

EC 1 (n=1) and EC 6 (n=2) are too small to characterize. Class membership comes from the packaged
reaction annotation; these results are exploratory and not adjusted for chemical dependence. No
major EC class in this paired set demonstrates lower average error for MetaG.

The class summaries and paired wins were examined after the benchmark was known. They do not supply
an outcome-independent rule for choosing a method, and the broader set retains identity uncertainty.

Reproduce with the original pairing table (SHA-256
`ce8905fea7c8d14d2a1a5525d8ee096c214588a2af534e54200c9f6a79f7e932`):

```
python benchmark/score_dgbyg_oof.py --pairs /path/to/metag_vs_dgbyg_pairs.tsv --output /tmp/oof.json
```

The command uses the guarded replay by default. Pass `--predictions benchmark/metag_0.3.0_predictions.json`
to reproduce the earlier combined-release numbers and exploratory breakdown above. The generated JSON
records the exact-family per-reaction pairs, the count of excluded unresolved
pairs, descriptive within-class correlations and hashes of all inputs. The completed identity metadata
in `dgbyg_family_identity_complete.json` is derived from the archived audit and does not include model
predictions; `dgbyg_family_identity.json` preserves the earlier 148-complete record.
The original pairing table is not included in this checkout. The underlying dGbyG cross-validation archive is available at
[Zenodo 13885478](https://doi.org/10.5281/zenodo.13885478) (`results_data.zip`, 35.5 GB); the
authors' original code archive is [Zenodo 15663012](https://doi.org/10.5281/zenodo.15663012).
The comparison requires the historical 675-row CV training table from
[the authors' pinned source revision](https://gitee.com/f-wc/dGbyG/raw/d539d2d4ace8905ae13dc93a18ce6f02493722a7/dGbyG/data/TrainingData.csv)
(SHA-256 `cd8ede797d2f1c088cccf6b6561d258a34e753ff1ea8c43ca8cb36ca28553f3c`).
The retrieval helper requests only the pinned 200 ZIP members (about 6.3 GB uncompressed), verifies
each member's name, length, CRC32 and SHA-256, and resumes from existing verified files. Alternatively,
use locally extracted NPY and CSV files from the same archive. The following builds the selection and
pairing locally:

```
python benchmark/fetch_dgbyg_cv.py --output-dir /tmp/dgbyg-raw
python benchmark/select_dgbyg_cv.py --raw-dir /tmp/dgbyg-raw --output-dir /tmp/dgbyg-selected
python benchmark/build_dgbyg_pairs.py --selected-dir /tmp/dgbyg-selected --training-csv /path/to/TrainingData.csv \
  --output /tmp/dgbyg-pairs.tsv --audit /tmp/dgbyg-pairs-audit.json
python benchmark/score_dgbyg_oof.py --pairs /tmp/dgbyg-pairs.tsv --output /tmp/dgbyg-oof.json
```

The selector checks the pinned archive-member hashes and follows the authors' minimum mean training-loss
epoch for each replicate. The pair builder checks the 100 selected prediction-vector fingerprints in
`dgbyg_cv_manifest.json` and uses the completed mapping and bridge. It reconstructs all 164 archived
rows, including 114 available benchmark pairs and 84 exact-family primary pairs. It does not independently verify the author's
fold assignments or the physical uncertainty of the bridge. The code's license must not be assumed to
resolve redistribution terms for every separate dataset; obtain and verify the author-provided inputs.
