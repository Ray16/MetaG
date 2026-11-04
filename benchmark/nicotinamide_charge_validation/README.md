# Vertical nicotinamide phosphate-charge control

Two SMD calculations remove two verified P-O-H protons from each previously
selected full-cofactor geometry. Every other coordinate remains fixed. The
oxidized/reduced charges change from+1/0 to-1/-2; electron counts remain344/346.
The endpoint is (charged reduced - charged oxidized) - (reference reduced -
reference oxidized). Atoms, charge and total molecular stoichiometry balance,
so uniform elemental, charge and solution-standard-state offsets cancel.
This exchange is chemically nonempty: a nonzero energy is allowed and does
not constitute a Hess-law cycle failure.

Run the numerical audit without CUDA, MetaG or external workspaces:

```bash
python audit_nicotinamide_charge_control.py . --output /tmp/charge-replay.json
```

The numerical protocol and original point producers are authenticated with
this packet. The original four-arm control is also available separately in
benchmark/nicotinamide_validation. Numerical library versions match those of
the parent points. The CPU/GPU water fixture and both final charged points
pass SCF, overlap-rank, component-sum and geometry-binding checks.

The result is a vertical aqueous electronic-plus-solvent coupling, not a
Gibbs energy, macroscopic pKa or empirical correction. It includes no new
nuclear relaxation, Hessian, RRHO term or aqueous conformational ensemble.
The separate raw H-removal energy differences omit the proton reservoir and
must not be interpreted as deprotonation Gibbs energies. Component sums refer
to independently polarized aqueous densities, rather than independent gas
and solvent interventions. Negative HOMO diagnostics are recorded but do not
prove orbital stability. No experimental label, fitted offset or production
policy is selected here. Startup records retain original cluster paths and
process identities; they are provenance, not a portable cluster launcher.
