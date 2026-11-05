# Molecule-local nicotinamide anchor: development study

The fixed approximation is G(full oxidized sibling) + G(methyl reduced core) - G(methyl oxidized core). Its atom/charge-balanced local transfer, source code and single reference choice were frozen before evaluating its new errors; previous development outcomes were already known. NAD and NADP retain separate exact full oxidized references. No source-amino repair or reaction-dependent routing is added.

On 112 common development-exposed labels, MAE changes from 16.6782 to 7.9363 kJ/mol (zero 12.1545). Of 29 affected labels, MAE changes 40.8959 to 7.1341; 83 labels are unchanged. Coverage is 160/364 attempted requests and 123/274 admitted development labels, with 185 missing archived molecular prerequisites and 19 workflow refusals. These are cached-reassembly counts, not a new physical sweep, independent performance estimate or genome-scale coverage.

All 52 successful affected requests cancel their full oxidized reference. Thus this test supports neither the absolute anchor in scaffold synthesis/exchange nor the omission of scaffold-specific redox effects. In particular, NAD/NADP transhydrogenation is forced to zero with zero conditional sampling uncertainty; its development label is -3.92kJ/mol. Locality/model error is excluded from this diagnostic. This policy is not adopted for production.

The immutable dependency is the existing sibling `global_amino_reference_validation` packet. `freeze.json` records its required file hashes; its complete manifest and 596 raw records are reverified. `outcomes.zip` contains all364requests for each arm, with states and primitive sampling derivatives. A virtual mapped intermediate in the inherited ledger is identified by `anchored_transfer` diagnostics; it is not a newly computed raw molecular record. `applicability_diagnosis.json` is explicitly post-outcome and did not select the candidate.

From the repository root, with the audit/test environment:

```sh
python benchmark/audit_anchored_nicotinamide.py --cold-replay
```

The audit reconstructs mapped energies, reaction sums and shared-primitive sampling propagation with 50-digit decimal arithmetic, checks the full baseline against its prior archive and reconstructs every available cycle. Cold replay separately runs both 364-request assemblies and all 640 invariance checks using the archived assembly source, while forbidding quantum imports and network access. This reuses retained records; it does not rerun quantum calculations. No independent primary biochemical labels are added.
