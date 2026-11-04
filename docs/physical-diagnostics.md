# Portable physical diagnostics

The companion `benchmark/thermochemistry_diagnostics.py` operates on user-supplied retained states or a complete set of reaction energies. It does not load UMA, run quantum chemistry, change caches, fit experimental reaction offsets or repair deployed predictions. It is delivered in the source distribution and code archive, rather than the minimal scoring wheel.

Install the declared CPU dependencies (`pip install '.[audit]'` from this source), then run:

```bash
python benchmark/thermochemistry_diagnostics.py transfer examples/diagnostics/transfer.json -o transfer_report.json
python benchmark/thermochemistry_diagnostics.py cycles examples/diagnostics/triangle.json -o cycle_report.json
```

The supplied examples are synthetic analytic controls, not experimental validation. The triangle requires a minimum uniform change of 1/3 kJ/mol; this is a feasibility bound, not an accuracy estimate. Existing output files are refused to preserve prior reports.

For neutral transfer, provide a common gas/solution state basis, explicit RT in kJ/mol, and each state's electronic, thermal and solvent terms. The script returns the partition difference, solvent envelope and minimum gap from a supplied reference center or interval. `reference_interval_kj`, if supplied, must come from the user's own declared reference assessment. No interval is inferred from the observed discrepancy. The tool assumes the state basis and molecular identity/minimum validation are supplied correctly; unique IDs alone do not certify distinct physical basins. Ions are refused because this interface provides no single-ion convention.

For cycles, provide one common conditions record, energies in kJ/mol, and each reaction's stoichiometry keyed by explicit chemical-state identities. State IDs must distinguish protonation, stereochemistry, declared solution pools and reservoirs as appropriate. Missing estimates are refused, never converted to zero. Include only the explicit cohort being assessed and report its availability separately.

The complete exact stoichiometric nullspace defines all closure constraints within that supplied cohort. The report identifies participating reactions; acyclic branches remain untested by closure. The consistency radius minimizes the largest absolute change needed to satisfy those constraints. It is invariant to the choice of cycle basis, at fixed reaction extent/units; an exact rational row-space reduction prevents ill-conditioned basis changes from silently losing constraints. Reaction reversal preserves the result; changing reaction extent requires explicit normalization if a different per-reaction error convention is desired.

The report supplies a primal perturbation and a unit-L1 dual reaction witness. These allow independent verification of the upper and lower bounds without trusting the optimizer alone. The primal vector is a mathematical certificate and is generally nonunique: it is not a physical correction, an attribution of each reaction's error, or a set of replacement predictions. The dual witness localizes an inconsistent relation. A positive radius proves incompatibility on the declared state space. A zero radius, or the absence of testable cycles, establishes no experimental accuracy or prospective uncertainty calibration.

The manuscript's retained-data example shows why this matters: disabling generic truncation removes some exhibited cycle failures while leaving a reference-state inconsistency and the same minimum consistency radius on identical targets. Correcting individual cycles is not sufficient to establish network thermodynamics. Standard cycle constraints, partition inequalities and data reconciliation have prior literature; the companion operationalizes their conditional use alongside molecular component records.

Prior foundations include [GoodVibes](https://doi.org/10.12688/f1000research.22758.1), [CENSO](https://doi.org/10.1021/acs.jpca.1c00971), [component contribution](https://doi.org/10.1371/journal.pcbi.1003098) and [biochemical data reconciliation](https://doi.org/10.1093/bioinformatics/btz741). The diagnostic does not claim these underlying principles as new.
