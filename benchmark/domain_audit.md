# Exploratory dGbyG training-domain analysis

This analysis investigates where MetaG might complement dGbyG. Both the MetaG method and the
Mg-bounded reference have been exposed during development. The dGbyG predictions here are from the
published fitted ensemble, **not out-of-fold predictions**. Training data and reference measurements
may overlap. These are hypothesis-generating comparisons, not an independent test or proof of
superiority on novel chemistry.

`analyze_dgbyg_domain.py` parses the released dGbyG training table and the original structures for
each benchmark reaction. It neutralizes molecules with RDKit, compares exact and stereo-free
direction-invariant reaction signatures (water and protons omitted), and measures each compound's
nearest training-compound similarity using radius-2 Morgan fingerprints and Tanimoto similarity.
The minimum over reaction compounds determines the reaction's compound similarity. Unseen means
no exact neutralized isomeric SMILES match in the parsed table; it does not guarantee no related
compounds or related measurements were used in training. Reaction matching is extent-sensitive.

The independent input tables used in this local run are identified by hashes in the generated JSON.
The reaction and compound matching omits 32 training rows from reaction signatures because their
normalized reaction cancels or a structure contains a wildcard; successfully parsed compounds from
cancelled rows still enter the compound set. This can understate training exposure. If removing
stereochemistry collapses a reaction, an exact stereochemical match can still identify a seen reaction;
otherwise stereo-free overlap is indeterminate and is reported separately. Reactions whose exact
signature also collapses remain unclassifiable. These categories are not certified training splits.

On the 239 reactions with earlier combined-release predictions and a primary-policy label:

| Descriptive training-domain category | n | MetaG MAE | dGbyG MAE | MetaG closer |
|---|---:|---:|---:|---:|
| Matched training reaction, exact or stereo-free | 155 | 8.16 | 3.52 | 31 |
| New reaction, all compounds matched | 1 | 14.16 | 29.44 | 1 |
| At least one unmatched compound, minimum similarity ≥ 0.5 | 77 | 9.21 | 6.42 | 29 |
| At least one unmatched compound, minimum similarity < 0.5 | 2 | 4.16 | 5.71 | 1 |
| Stereo-free signature collapses, exact reaction unmatched | 2 | 1.93 | 0.33 | 0 |
| Exact signature also collapses | 2 | 5.71 | 5.71 | 1 |

Errors are in kJ/mol. The two most distant-compound cases are NAD- and NADP-dependent versions of the
same methylenetetrahydrofolate dehydrogenase transformation, not two independent chemistry families.
Their different references also require condition and identity adjudication. They cannot establish
generalization to distant molecules.

## Redox cases with larger individual differences

The published-model comparison suggests individual examples larger than the sugar stereochemistry
effect. D-octopine dehydrogenase (`rxn00402`) has a packaged reference of +36.62 kJ/mol, MetaG +36.3,
and dGbyG +11.28. The reference is one 298.15 K measurement, with ionic strength imputed; it comes from
one source, so it does not establish a repeatable 25 kJ/mol advantage. Glyoxylate reductase
(`rxn00512`) has reference +48.23, MetaG +48.7, and dGbyG +38.61 kJ/mol, all in the same nominal
orientation and conditions.

The glyoxylate reference contains a material source conflict that its overall median obscures. In
the historical Mg-bounded reconstruction, five measurements from Cartwright and Hullin (1966) have
a transformed median of +62.22 kJ/mol, while seven measurements from Zelitch (1955) have a median of
+48.17 kJ/mol. All these values use the historical I=0.25 M imputation and the shared transform model;
the roughly 14 kJ/mol difference is not resolved by these data alone. The complete reaction mapping,
experimental conditions and primary studies need adjudication before using the example as a decisive
accuracy result. The packaged +48.23 value is the median over all 16 admitted measurements, not a
measurement-independent consensus.

As a deliberately narrow exploratory check, species labels exactly equal to NAD, NADH, NADP, NADPH,
NAD+ or NADP+ select 69 paired reactions: MetaG MAE is 8.78 kJ/mol and dGbyG MAE 6.90 kJ/mol; MetaG
is closer on 28. The flag can miss other spellings and structures and is not a validated redox
classifier. It does not establish a MetaG advantage for NAD-dependent reactions. This is the published
dGbyG ensemble, with possible training overlap; the labels and MetaG are development-exposed.

The script also computes a direction-invariant reaction-change similarity: stoichiometric differences
of sparse Morgan radius-2 count fingerprints, compared with generalized Jaccard similarity to each
parsed training reaction. This is a proxy for transformation novelty, not an atom-mapped reaction
center or a guarantee that spectators cancel at every fingerprint radius. For novel-classified
reactions, previously used exploratory bins yield:

| Nearest reaction-change similarity | n | MetaG MAE | dGbyG MAE | MetaG closer |
|---|---:|---:|---:|---:|
| [0, 0.6) | 7 | 7.13 | 7.63 | 3 |
| [0.6, 0.8) | 6 | 17.43 | 5.82 | 1 |
| [0.8, 1] | 66 | 8.69 | 6.77 | 27 |

These small, exploratory groups do not demonstrate a reliable advantageous regime. Among the seven
lowest-similarity cases, three favor MetaG and four favor dGbyG; the 0.5 kJ/mol MAE difference is
especially fragile. Their similarity values agree exactly with the archived reaction-change audit
for all 219 comparable nonempty fingerprint reactions in this run. The fingerprints do not encode
chirality; epimerizations can have an empty change fingerprint and cannot be assessed by this proxy.
The script also reports a separate chirality-aware Morgan change similarity, but its ability to
predict extrapolation errors is not validated. `screen_dgbyg_domain.py` reports nearest matching
training rows for structure-only inputs. See [the rare-sugar candidate screen](candidates/itoh_1994_screen.md).

Reproduce with the archived training and comparator files:

```
python benchmark/analyze_dgbyg_domain.py \
  --training /path/to/TrainingData.csv \
  --comparator /path/to/all_predictions.tsv \
  --output /tmp/metag_dgbyg_domain.json
```

The input files are not included in this release. The run used RDKit 2026.03.5. Before a publication
claim, independently verify the released-weight training-data snapshot, prediction provenance,
reaction mappings, measurement overlap and reference conditions. A suitable next test would freeze a
structural or reaction-change criterion on this development set, then evaluate it on separately
curated public measurements not used in MetaG development, with both methods scored on identical
reactions and the coverage and errors of both reported.
