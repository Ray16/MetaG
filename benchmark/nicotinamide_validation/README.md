# Fixed-geometry nicotinamide component control

The closed exchange is full NADox + MNAH -> full NADH + MNA+. Full reference
states have charge+1/0 with protonated phosphate groups. This is not a
validated biochemical protonation ensemble. Atom, charge and molecular-count
balance remove elemental, single-ion and uniform solution-standard-state
offsets from this exchange.

All three oxidized and four reduced core minima and one lowest-original-total-G full minimum
per cofactor were fixed before new DFT energies. The original full pools of 41 and 55 minima
remain in raw_records and their original exchange is reported separately.
Every arm retains the original RRHO term. Alternative DFT energies and SMD
shifts are evaluated at unchanged geometries, without new nuclear optimization,
Hessians or aqueous ensemble sampling. SMD includes its CDS term exactly once.

Reconstruct and verify all four arms without MetaG, CUDA or external files:

```bash
python audit_nicotinamide_control.py . --raw-directory raw_records --output /tmp/nicotinamide-replay.json
```

This audits retained numerical results; it does not rerun quantum calculations.
Recomputing those requires the declared PySCF/GPU4PySCF runtime and the original
UMA checkpoint, references and pinned FairChem SDK. Launch/point records and
code retain original paths and process identities for provenance. They are not
a portable cluster launcher. dependency_sources.json lists the copied engine
and SDK hashes; those third-party trees and checkpoint are not redistributed.

The run used one lambda2 V100 with CPUs 2-5 and native thread caps of one. An API
adapter failure and a default overlap-rank truncation were caught and archived
before this final implementation. The same predeclared 1e-8 overlap gate is used
explicitly in the corrected calculation; nMO=nAO is required. The numerical
water fixture compares CPU and GPU density-fitted gas/SMD energies. It does
not bound fitting, functional, grid, missing-state or experimental error.
No orbital-stability proof, total-error interval, fitted offset or production
replacement is claimed. All four planned arms are retained.

Primary redox references were adjudicated in reference_adjudication.json.
They are not condition-matched independent experimental exchange labels.
No experimental accuracy or unfamiliar-chemistry advantage follows from this
component control. The open primary MNA source is
https://doi.org/10.1038/srep37596; the transformed NAD table is
https://doi.org/10.1016/j.bpc.2004.05.003.
