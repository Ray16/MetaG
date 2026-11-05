"""Complete-pool replay and refusal of omitted or inconsistent SCF evidence."""
from pathlib import Path
import hashlib
import importlib.util
import json
import shutil
import zipfile

import pytest

BENCHMARK=Path(__file__).resolve().parents[1]/'benchmark'
spec=importlib.util.spec_from_file_location('cofactor_complete_pool_audit',BENCHMARK/'nicotinamide_ensemble.py')
auditor=importlib.util.module_from_spec(spec);spec.loader.exec_module(auditor)
PACKET=BENCHMARK/'nicotinamide_ensemble_validation'


def test_complete_retained_pool_replay():
    report=auditor.audit(PACKET)
    assert report['verified_pairs']==103 and report['reused_pairs']==9
    assert len(report['independent_decimal_exchange_kj'])==8
    assert not report['production_reference_policy_adopted']


@pytest.mark.parametrize('mutation',['omit_required_pair','inconsistent_water_total'])
def test_incomplete_or_inconsistent_scientific_evidence_refused(tmp_path,mutation):
    root=tmp_path/'packet';shutil.copytree(PACKET,root)
    archive=root/'calculation_records.zip'
    with zipfile.ZipFile(archive) as source:contents={n:source.read(n) for n in source.namelist()}
    key='points/full_red__minimum_012.json'
    assert key in contents
    if mutation=='omit_required_pair':del contents[key]
    else:
        point=json.loads(contents[key]);point['water']['energy_Eh']+=0.001
        contents[key]=json.dumps(point).encode()
    with zipfile.ZipFile(archive,'w',compression=zipfile.ZIP_DEFLATED) as destination:
        for name,value in contents.items():destination.writestr(name,value)
    manifest=json.loads((root/'manifest.json').read_text())
    manifest['files']['calculation_records.zip']=hashlib.sha256(archive.read_bytes()).hexdigest()
    (root/'manifest.json').write_text(json.dumps(manifest))
    # Updating the outer checksum must not make incomplete or internally
    # inconsistent calculations admissible as a complete scientific endpoint.
    with pytest.raises((KeyError,AssertionError)):auditor.audit(root)
