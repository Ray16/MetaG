# Released-ensemble uncertainty: development-only comparison

This packet contains all 100 member values for 120 normalized molecular inputs
and 69 reaction-member sums from the released dGbyG model, plus the fixed
condition/extent normalization and complete equal-coverage curves. It is a
separate experiment from the earlier cross-validation-repeat reliability study.
It inherits that cohort's reference-temperature and condition qualifications;
no new experimental temperature normalization is performed.
The corresponding native means are used throughout these curves; native spreads
are not paired with OOF means. Molecular member covariances are retained by
summing each shared checkpoint column before calculating population SD.

All 69 requests returned finite predictions. Sixty-eight exact normalized
reaction equations also occur in the shipped 671-row training table. The other
case does not establish graph, family, measurement or checkpoint independence.
The full cohort was already used during MetaG development. No independent risk
reduction, significance, calibrated total uncertainty or OOD result is claimed.

At equal retained counts, descriptive comparator AURC is 2.0962 for native
spread, 1.7753 for disagreement, 1.7249 for the fixed equal-rank combination and
2.6598 for magnitude. The full-cohort native comparator MAE is 2.9793 kJ/mol;
it must not replace the separately reported OOF MAE. Only one native comparator
error exceeds 20 kJ/mol. No deployment policy is selected from these results.

The first extraction stopped because a NumPy QC-bound scalar was not JSON
serializable. The complete initial producer/log/partial-export evidence is in
infrastructure_attempt1.zip. The revision only converts that exported scalar to
builtin float, retains the scientific protocol and original 900-second budget,
and reruns the fixed cohort. It does not change predictions or select cases.

Portable verification, without Torch or model weights:

    python benchmark/native_uncertainty.py --identity-replay

The arithmetic verifier rebuilds all 6900 reaction-member sums, independently
recomputes means/SDs with decimal arithmetic, and checks all 2484 risk entries.
The optional identity replay uses RDKit to independently normalize the archived
source equations. The original reference source is MIT-licensed and its license
is included in source_implementation.zip. Checkpoint hashes and archived inference
code document model provenance; the scalar audit does not rerun model inference.
