# Parent-model check of the aqueous geometry intervention

This control holds four previously generated NAD geometries fixed and compares
their electronic/solvent scores. It admits no experimental reference, computes
no reaction Gibbs energy and makes no production-policy change.

For each of the oxidized/reduced ions, charges −1/−2, one geometry was relaxed
against the gas UMA potential and one against
`E_UMA + E_GFN2,ALPB − E_GFN2,gas`. The same ordered atoms, connectivity and
specified stereochemistry are retained within each pair. All geometries were
selected before the new quantum calculations. No new nuclear optimization is
performed in this control.

The parent calculations use restricted-singlet wB97M-V/def2-TZVPD, density
fitting with def2-universal-jkfit, grids 4/3, and SMD water with 590 SWIG surface
points. The overlap gate is 1e-8 and no admitted AO direction is discarded.
Complete orbital SCF convergence is distinct from nuclear-force convergence;
these structures are not parent-DFT minima.

Each aqueous scoring arm evaluates both sets of geometries. The geometry
intervention is `(R_aq − O_aq) − (R_gas − O_gas)`, where the subscripts denote
the geometry-generation potential, not the phase of the final score.

| Aqueous scoring arm | Oxidized geometry change | Reduced geometry change | Intervention contrast |
|---|---:|---:|---:|
| UMA/ALPB | −6.83 | −51.54 | −44.71 |
| DFT/ALPB | −7.52 | −51.02 | −43.50 |
| UMA/SMD | −12.16 | −33.46 | −21.30 |
| DFT/SMD | −12.85 | −32.95 | −20.09 |

Values are kJ/mol. Electronic replacement changes the contrast by +1.21 kJ/mol
at fixed ALPB; replacing ALPB by the parent-density SMD increment changes it
by +23.40 kJ/mol at fixed UMA. The combined change is +24.61 kJ/mol. The response
retains its sign, while its magnitude depends substantially on the solvent
model. This does not establish which model better predicts experimental Gibbs
energies or repair the developmental redox errors.

The SMD increment is shared between the two electronic arms. Their linear
factorial interaction is therefore zero by construction, not a finding of
independent physical errors. Different constrained microgeometries are different
states; this contrast is not a closed cycle of identical equilibrium metabolites.
Common elemental, linear charge and uniform molecular offsets cancel.

`audit_summary.json` reconstructs all four arms from the native logs, per-point
SCF terms and retained UMA/xTB records. It checks nuclear repulsion independently
from atomic numbers and ordered coordinates, full AO counts, dependencies,
solvent conventions and the CPU/GPU water fixture. It never adds the duplicate
`e2` term or the orbital-gap diagnostic to the SCF energy. Missing endpoints,
mixed versions and nonfinite records fail rather than becoming zero.

The inherited RRHO sensitivity is recorded separately and enters none of the
table values. No parent Hessian, complete aqueous conformer ensemble, pH/Mg
transformation, calibrated error interval or biochemical accuracy is established.
Parent electronic agreement alone cannot supply those missing validations.

From the MetaG checkout, with NumPy installed:

```bash
python benchmark/audit_aqueous_geometry_parent.py \
  benchmark/aqueous_geometry_validation \
  --output /tmp/aqueous_geometry_audit.json
```

The retained producer scripts describe the original cluster run and its bounded
resource policy. The audit requires no CUDA, model weights, PySCF, xTB, original
workspace paths or new molecular calculation. Upstream producer metadata is
retained without redistributing model weights. `manifest.json` authenticates all
retained files except itself; the launch manifest separately records the frozen
inputs and numerical producer used before the calculations.
