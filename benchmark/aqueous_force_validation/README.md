# ALPB solvent-force validation

This packet checks the gradient of the matched total-energy increment
`E(GFN2/ALPB water) - E(GFN2 gas)` before an experimental solution-geometry
profile is considered. It changes no production scoring path and establishes
no molecular Gibbs accuracy.

The refined-grid run contains 800 native executions on ethanol, acetate,
glycine zwitterion and dihydrogen phosphate. Every Cartesian force component
is compared against central energy differences at two displacement sizes.
Net force, net torque and step-size criteria were declared before evaluation.
All four cases pass with the documented `extreme` SASA grid. The maximum
component difference is below 1.3e-5 eV/angstrom on these fixtures.

The preceding default-grid ethanol trial is retained separately. Its 218
executions pass energy/force finite differences but fail the same torque
criterion. These are four distinct molecules and five tested configurations;
the repeated ethanol trial is not an additional independent chemical case.
Refining the grid changes the numerical potential and must not relabel the
existing production calculations.

Run both numerical replays without xTB, UMA or access to the original cluster:

```bash
python benchmark/audit_alpb_forces.py benchmark/aqueous_force_validation --output /tmp/alpb-force-audit.json
python benchmark/audit_alpb_forces.py benchmark/aqueous_force_validation/default_grid_trial --output /tmp/alpb-default-audit.json
```

The auditor reconstructs every derivative from retained native energies and
gradients, checks phase, charge, spin, displaced coordinates and recorded grid,
and verifies the reported pass or failure. Raw records are compressed without
altering their decompressed bytes. The manifest authenticates files and records
the original JSONL hashes. It is not evidence of physical accuracy.

No Hessian, aqueous entropy, conformer-completeness or large-cofactor accuracy
bound follows from these four fixtures. Calculation source and original paths
remain for provenance; the numerical replay is portable, while the producer
requires the declared xTB binary and runtime. Grid options are documented in
the [official xTB solvation guide](https://xtb-docs.readthedocs.io/en/latest/gbsa.html).
