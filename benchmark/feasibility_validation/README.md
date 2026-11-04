# Independent-benchmark feasibility, 4 November 2026

This study identifies structural holdout and source-curation leads. It admits
**zero experimental references**, performs no prediction or model training, and
does not freeze a test set or validate the production estimator. The enrollment
and success rules in [the validation design](../state_function_validation/design.json)
remain unchanged.

The inspected history contains 383 distinct drawn reaction inputs from a bounded
JSON inventory. A historical input witness does not establish experimental-label
use, and absence from this inventory does not establish independence. Molecular
records, notebooks, other formats, and uninspected files remain unresolved.
Backups and the MetaG-JCIM project were excluded from the inventory.

## Results

| Structural or metadata screen | 671-row corpus | 675-row corpus |
|---|---:|---:|
| Reaction targets | 449 | 453 |
| Equations with concrete structures | 440 | 446 |
| Equations matching inspected development inputs after operational aliases | 328 | 327 |
| Concrete equations without such a match | 112 | 119 |
| Unmatched equations with a resolved broad family absent from inspected history | 30 | 32 |
| Fully resolved equations receiving molecular holdout previews | 22 | 24 |
| Previews with at least one held-out family below similarity 0.5 in retained training support | 4 | 5 |
| Unmatched equations with strict named-equation source leads | 87 | 90 |

These counts are overlapping attributes, not a sequence of independent cohorts.
In particular, an absent resolved family can coexist with an unresolved partner.
The 22/24 previews have no unresolved candidate aliases, but their experimental
identities, ancestry, conditions and development exposure are still unadjudicated.

Protonation, reducing-sugar ring/open drawings, and simple gem-diol/carbonyl
drawings are reconciled only for conservative exposure screening. They are not
evidence that an experimental equilibrium pool has been identified or that its
free energy is modeled correctly. Broad exposure families additionally erase
stereochemistry and canonicalize tautomers; they are not molecular identities or
scaffolds. Enumeration limits remain explicit unknowns.

The previews remove held-out molecular families from **both formation and
reaction labels**, remove exact/reversed/scaled and operationally equivalent test
equations, and exclude every unresolved training row, including all wildcard
rows. The unresolved-row policy excludes 92/93 rows before most family-specific
exclusions, leaving 577–578 or 580–581 training rows. These are limitations of
this conservative audit, not declarations that those experimental labels are
invalid. Comparative retraining would require nested controls using the same
unresolved-row policy and training-only normalization.

Training-support loss matters: the fifth distant case in the 675-row corpus has
nearest similarity 0.4444 after exclusion but 0.5758 in the original concrete
other-family subset. Four cases in each corpus retain a below-0.5 participant
relative to that larger concrete subset. Unresolved chemistry still prevents a
complete original-domain novelty certificate. The two corpora substantially
overlap and cannot be counted as independent replications.

Neither corpus supplies ten cases meeting the retained-domain distance screen,
let alone the ten independently adjudicated source/family clusters required by the existing design. This route
alone therefore cannot support the planned primary unfamiliar-chemistry claim.
It can support focused source review and later exploratory controls.

## Source review

The refreshed OpenTECR file is byte-identical to the 29 October snapshot:
5,735 rows, SHA256 recorded in `nightly_refresh.json`. The exported source table
contains only identifiers, named equations and condition metadata; it omits
equilibrium constants and enthalpies.

Source matching retains phase and stereochemical names and normalizes only
letter case, whitespace, primitive stoichiometric extent and direction. A match
is a lead, not proof that a training target derives from that measurement. Only
one unmatched equation per corpus has any named-source row with all four
temperature/ionic-strength/pH/pMg fields filled. Filled fields still require
validation; blank fields never mean zero ionic strength or absent magnesium.
Primary texts may recover conditions missing from the transcription.

One concrete identity discrepancy is recorded in
`valine_identity_adjudication.json`. The 671-row table encodes a straight-chain
2-oxovalerate product. The original study's
[NIST abstract](https://www.nist.gov/publications/thermodynamics-reactions-catalyzed-branched-chain-amino-acid-transaminase)
instead identifies branched 2-oxoisovalerate. On printed page 1363, the
[2007 evaluated table](https://srd.nist.gov/jpcrdreprint/1.2789450.pdf)
uses the former name in its heading but the latter in its chemical-reference
annotation. This supports a source-name propagation problem; it does not prove
the numerical label's derivation. The original training file and target remain
unchanged. This entry is not eligible for independent enrollment pending review.

`exposure_events.json` records numerical values inadvertently visible during
schema and source inspection. No outcome magnitude, prediction, disagreement or
residual selected these candidates. This is not a claim of complete historical
outcome blindness.

## Replay

From the MetaG checkout, with its dependencies installed:

```bash
PYTHONPATH=src python benchmark/benchmark_feasibility.py \
  --training benchmark/feasibility_validation/training671_structural.csv \
  --panel benchmark/feasibility_validation/historical_inputs.json \
  --source-metadata benchmark/feasibility_validation/opentecr_source_metadata.csv \
  --output /tmp/public671_audit.json \
  --source-output /tmp/public671_source_candidates.json \
  --holdout-output /tmp/public671_holdout_previews.json
```

Use `training675_structural.csv` for the alternate corpus. Outputs contain
zero-based source-row indices. Exact replay is bound to the recorded RDKit
version, audit source and normalization dependencies. Other versions may change
alias construction and must be treated as new audits.

The protocol records original and projected input hashes. `manifest.json`
authenticates all retained files except itself. The structural projections
preserve row order and omit Gibbs-energy targets and errors. Source reports
and holdout previews are reproducible using only these portable inputs and the
MetaG source; the original private working directories are not required.

The dGbyG inputs derive from Fan et al.,
[Cell Systems (2025)](https://doi.org/10.1016/j.cels.2025.101393)
and the authors' repositories. The upstream MIT notice is retained in
`DGBYG_LICENSE`. Exact released-checkpoint correspondence remains unresolved.
No claim of native pretrained-model performance or uncertainty follows from
these projected tables or unexecuted retraining plans.
