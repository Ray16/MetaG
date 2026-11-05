"""Complete chemical-control replay and refusal of damaged scientific outputs."""
from pathlib import Path
import hashlib,importlib.util,json,shutil,zipfile
import pytest

BENCHMARK=Path(__file__).resolve().parents[2]/'benchmark'
spec=importlib.util.spec_from_file_location('chemical_reference_transfer_audit',BENCHMARK/'reference_transfer.py')
auditor=importlib.util.module_from_spec(spec);spec.loader.exec_module(auditor)
PACKET=BENCHMARK/'reference_transfer_validation'

def test_complete_chemical_control_replays():
    report=auditor.audit(PACKET)
    assert report['n_optimized_starts']==48 and report['n_retained_geometries']==30
    assert report['n_raw_chemical_references']==2 and len(report['predictions_kj'])==8
    assert report['primary_independent_metabolic_references_added']==0

@pytest.mark.parametrize('mutation',['missing_point','inconsistent_SMD_total'])
def test_invalid_SCF_evidence_refused_after_outer_checksum_update(tmp_path,mutation):
    root=tmp_path/'packet';shutil.copytree(PACKET,root);path=root/'calculation_records.zip'
    with zipfile.ZipFile(path) as archive:files={n:archive.read(n) for n in archive.namelist()}
    name='points/acetate_ion__gas__03.json';assert name in files
    if mutation=='missing_point':del files[name]
    else:
        point=json.loads(files[name]);point['water']['energy_Eh']+=.001;files[name]=json.dumps(point).encode()
    with zipfile.ZipFile(path,'w',compression=zipfile.ZIP_DEFLATED) as archive:
        for name,value in files.items():archive.writestr(name,value)
    manifest=json.loads((root/'manifest.json').read_text());manifest['files']['calculation_records.zip']=hashlib.sha256(path.read_bytes()).hexdigest();(root/'manifest.json').write_text(json.dumps(manifest))
    with pytest.raises((FileNotFoundError,AssertionError,ValueError)):auditor.audit(root)
