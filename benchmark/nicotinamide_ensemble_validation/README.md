# Complete retained nicotinamide ensembles

All 103 originally retained geometries are included: 41/55 full oxidized/reduced
cofactor minima and 3/4 small-core minima. Nine exact earlier gas/SMD pairs are
reused; 94 new pairs extend the same numerical recipe. Every arm retains the
original gas geometry and RRHO correction. Full reference charges are +1/0 with
protonated phosphate groups. These are not biochemical protonation ensembles.

The exchange is full oxidized + core reduced -> full reduced + core oxidized.
It conserves atoms, charge and molecular count. All four electronic/solvent arms
are reported, including their complete conformer weights and components. No
failed or missing conformer is dropped. The frozen analysis predates the new
DFT outputs; the old nine-point results were already known. Its primary contrast
is the full minus restricted DFT/SMD exchange. The 5 kJ/mol descriptive threshold
measures intervention magnitude, not experimental accuracy or a promotion gate.

Replay the complete numerical evidence with standard-library Python:

```bash
python nicotinamide_ensemble.py .
```

The archive includes immutable scientific inputs, original source records, point
results/SCF logs, numerical fixtures, architecture bridges, producer/ownership
records, frozen analysis and resource-only amendments. Source paths and launcher
settings preserve provenance; they are not portable launch commands. The auditor
reads the archive without extracting or accessing those absolute paths. Independent
50-digit partitions use DFT/SMD water totals directly. This checks saved values,
not recomputation of SCF or an electronic ground-state/stability proof.

Recomputation requires the declared UMA producer and PySCF/GPU4PySCF environment.
The original two-worker plan was enlarged to seven H100 workers under a 30 GPU-hour
ceiling; scientific inputs, tolerances, one-attempt rule and analysis stayed fixed.
Missing-state, aqueous sampling, protonation, Mg binding, ionic-strength and model
errors remain unvalidated. No independent experimental label is added, no empirical
offset is fitted, and no production reference policy is adopted by this control.
