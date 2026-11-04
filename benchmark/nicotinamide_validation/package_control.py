"""Build a self-contained arithmetic replay only after all nine controls pass."""
from pathlib import Path
import hashlib
import importlib.util
import json
import shutil

HERE=Path(__file__).resolve().parent
SOURCE=HERE.parents[1]/'benchmark/audit_nicotinamide_control.py'
spec=importlib.util.spec_from_file_location('control_audit',SOURCE)
auditor=importlib.util.module_from_spec(spec)
spec.loader.exec_module(auditor)
result=auditor.audit(HERE)
assert result['passed'] and result['n_slots']==9
assert json.loads((HERE/'control_status.json').read_text())['status']=='complete'
destination=HERE/'public_bundle'
assert not destination.exists(), 'Preserve existing public bundle; inspect before changing it'
destination.mkdir()
manifest=json.loads((HERE/'launch_manifest.json').read_text())
files=set(manifest['owned_file_sha256'])|{'launch_manifest.json','uma_replay.json',
    'reference_adjudication.json','infrastructure_revision.json','control_status.json',
    'launcher.log','controller_owner.json','launcher_owner.json'}
files.update(str(p.relative_to(HERE)) for p in (HERE/'points').iterdir() if p.is_file())
files.update(p.name for p in HERE.glob('*.console.log'))
files.update(p.name for p in HERE.glob('*.predictor_owner.json'))
for revision in ['infrastructure_attempt1','infrastructure_attempt2']:
    for name in ['reference_point.py','run_control.py','control_status.json','gpu_smoke.console.log','launcher.log']:
        if (HERE/revision/name).exists():files.add(revision+'/'+name)
for name in sorted(files):
    target=destination/name
    target.parent.mkdir(parents=True,exist_ok=True)
    shutil.copy2(HERE/name,target)
raw=destination/'raw_records';raw.mkdir()
for value in json.loads((HERE/'geometry_inputs.json').read_text())['source_records'].values():
    shutil.copy2(value['path'],raw/Path(value['path']).name)
shutil.copy2(SOURCE,destination/SOURCE.name)
shutil.copy2(__file__,destination/Path(__file__).name)
second=auditor.audit(destination,raw)
assert second==result
(destination/'summary.json').write_text(json.dumps(result,indent=2)+'\n')
(HERE/'summary.json').write_text(json.dumps(result,indent=2)+'\n')
macros={
    'nicPointN':str(result['n_slots']),
    'nicRestrictedOriginal':f"{result['arms']['UMA_ALPB']['exchange_kj']:.2f}",
    'nicDFTALPB':f"{result['arms']['DFT_ALPB']['exchange_kj']:.2f}",
    'nicUMASMD':f"{result['arms']['UMA_SMD']['exchange_kj']:.2f}",
    'nicDFTSMD':f"{result['arms']['DFT_SMD']['exchange_kj']:.2f}",
    'nicAllMinimaOriginal':f"{result['original_all_minima_exchange_kj']:.2f}",
    'nicRestrictedEffect':f"{result['restricted_geometry_effect_kj']:.2f}",
    'nicElectronicEffect':f"{result['factorial_contrasts']['DFT_minus_UMA_at_ALPB_kj']:.2f}",
    'nicSolventEffect':f"{result['factorial_contrasts']['SMD_minus_ALPB_at_UMA_kj']:.2f}",
    'nicInteraction':f"{result['factorial_contrasts']['interaction_kj']:.2f}",
}
(destination/'nicotinamide_numbers.tex').write_text('% Generated from audited fixed-geometry control.\n'+
    ''.join('\\newcommand{\\'+name+'}{'+value+'}\n' for name,value in macros.items()))
(destination/'README.md').write_text('''# Fixed-geometry nicotinamide component control

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
''')
owned={str(p.relative_to(destination)):hashlib.sha256(p.read_bytes()).hexdigest()
       for p in destination.rglob('*') if p.is_file()}
(destination/'manifest.json').write_text(json.dumps(dict(files=owned,
    audit_source_sha256=auditor.file_sha(SOURCE),scope=result['scope']),indent=2)+'\n')
print(json.dumps(dict(arms={k:v['exchange_kj'] for k,v in result['arms'].items()},
                     factorial=result['factorial_contrasts'],bundle=str(destination)),indent=2))
