"""Portable data audit for the global-amino candidate and curvature-recovery case."""
from pathlib import Path
import argparse
import hashlib
import json
import math
import os
import shutil
import subprocess
import sys
import tempfile
import zipfile


def sha(data):return hashlib.sha256(data).hexdigest()

def verify_manifest(packet):
    packet=Path(packet);manifest=json.loads((packet/'manifest.json').read_text())
    actual={str(p.relative_to(packet)) for p in packet.rglob('*') if p.is_file() and '__pycache__' not in p.parts}
    if actual!=set(manifest['files'])|{'manifest.json'}:raise ValueError('Evidence inventory changed')
    for name,entry in manifest['files'].items():
        path=packet/name
        if Path(name).is_absolute() or '..' in Path(name).parts or path.is_symlink():raise ValueError('Invalid evidence path')
        data=path.read_bytes()
        if len(data)!=entry['bytes'] or sha(data)!=entry['sha256']:raise ValueError('Changed evidence bytes: '+name)
    return len(manifest['files'])

def extract(archive_path,destination):
    with zipfile.ZipFile(archive_path) as archive:
        entries=archive.infolist()
        if len({e.filename for e in entries})!=len(entries):raise ValueError('Duplicate archive path')
        if sum(e.file_size for e in entries)>150_000_000:raise ValueError('Archive exceeds study size limit')
        for entry in entries:
            name=Path(entry.filename)
            if name.is_absolute() or '..' in name.parts or (entry.external_attr>>16)&0o170000==0o120000:
                raise ValueError('Invalid archive entry')
            target=Path(destination)/name
            target.parent.mkdir(parents=True,exist_ok=True)
            target.write_bytes(archive.read(entry))

def audit(packet):
    packet=Path(packet);count=verify_manifest(packet)
    catalog=json.loads((packet/'catalog596.json').read_text())
    assert len(catalog['records'])==596
    producer=catalog['source_profiles']['original']['producer']
    recovery=catalog['source_profiles']['curvature_recovery']['producer']
    assert {k for k in set(producer)|set(recovery) if producer.get(k)!=recovery.get(k)}=={'source_sha256'}
    snapshots={}
    for name in ['assembly_source','recovery_source']:
        with zipfile.ZipFile(packet/(name+'.zip')) as archive:
            prefix=name+'/metag/'
            entries={path[len(prefix):]:archive.read(path) for path in archive.namelist() if path.startswith(prefix)}
        snapshots[name]=entries
        selected={path:data for path,data in entries.items() if
            path.startswith('energetics/') and path.endswith('.py') or path in
            ['chem.py','symmetry.py','settings.py','constants.py','timing.py']}
        digest=hashlib.sha256(b'metag-raw-species-source-v1\0')
        for path,data in sorted(selected.items()):
            digest.update(path.encode()+b'\0'+str(len(data)).encode()+b'\0'+data)
        expected=producer if name=='assembly_source' else recovery
        assert digest.hexdigest()==expected['source_sha256']
    assert set(snapshots['assembly_source'])==set(snapshots['recovery_source'])
    changed=sorted(path for path,data in snapshots['assembly_source'].items()
        if snapshots['recovery_source'][path]!=data)
    assert changed==['energetics/minimum_optimizer.py','energetics/species.py']
    maximum=0.;minimum_count=0
    with zipfile.ZipFile(packet/'raw_records.zip') as archive:
        assert set(archive.namelist())=={r['name'] for r in catalog['records']}
        for row in catalog['records']:
            data=archive.read(row['name']);assert sha(data)==row['sha256']
            raw=json.loads(data)
            assert raw['q']==row['charge'] and raw['smi']==row['smiles']
            assert raw['settings']['producer']==catalog['source_profiles'][row['raw_profile']]['producer']
            ensemble=raw['meta']['thermal']['ensemble'];minima=ensemble['validated_minima']
            assert minima and ensemble['n_reference_thermal_approximation']==0
            values=[]
            for point in minima:
                value=point['electronic_kj']+point['solvation_kj']+point['thermal_kj']
                assert abs(value-point['total_kj'])<1e-7
                geometry=dict(point['geometry']);digest=geometry.pop('sha256')
                assert sha(json.dumps(geometry,sort_keys=True,allow_nan=False).encode())==digest
                assert geometry['charge']==row['charge']
                values.append(value)
            low=min(values)
            # Keep the evaluated species-partition convention; pKa/reservoir
            # constants have their own explicitly retained deployed rounding.
            total=low-2.4789*math.log(math.fsum(math.exp(-(v-low)/2.4789) for v in values))
            difference=abs(total-raw['G']);maximum=max(maximum,difference)
            assert difference<1e-7 and raw['G']==row['G_kj'] and raw['sigma']==row['sigma_kj']
            minimum_count+=len(minima)
    report=json.loads((packet/'comparison_report.json').read_text())
    with zipfile.ZipFile(packet/'outcomes.zip') as archive:
        base=json.loads(archive.read('reassembly596/full_extended596.json'))['results']
        new=json.loads(archive.read('reassembly596/global_amino_extended596.json'))['results']
    labels=json.loads((packet/'reference.json').read_text())['2.0']['labels']
    common=sorted(rid for rid in labels if base[rid]['status']==new[rid]['status']=='estimated')
    assert common==report['paired_ids'] and len(common)==112
    independent={}
    for name,rows in [('baseline',base),('candidate',new)]:
        value=math.fsum(abs(rows[rid]['result']['dG_unrounded']-labels[rid]['exp']) for rid in common)/len(common)
        assert abs(value-report['equal_catalog_common_labels'][name]['MAE'])<1e-10
        independent[name]=value
    assert not report['production_reference_policy_adopted'] and report['new_primary_references_admitted']==0
    invariants=json.loads((packet/'invariances596.json').read_text())
    assert (invariants['n_reaction_checks'],invariants['n_protonation_drawings_checked'])==(592,44)
    return dict(status='passed',authenticated_files=count,authenticated_raw_records=596,
        validated_minima_reconstructed=minimum_count,maximum_raw_partition_difference_kj=maximum,
        common_labels=112,independent_common_MAE_kj=independent,
        production_reference_policy_adopted=False,experimental_validation_regime='development_exposed')

def replay(packet,*,physical=False):
    packet=Path(packet)
    with tempfile.TemporaryDirectory(prefix='metag-amino-replay-') as name:
        root=Path(name)
        for path in packet.rglob('*'):
            if path.is_file() and path.suffix in {'.json','.py'} and path.name!='manifest.json':
                target=root/path.relative_to(packet);target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(path,target)
        for archive in ['assembly_source.zip','outcomes.zip','numerical_controls.zip']:extract(packet/archive,root)
        environment=dict(os.environ,PYTHONDONTWRITEBYTECODE='1',OMP_NUM_THREADS='1',OMP_THREAD_LIMIT='1',
            MKL_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',NUMEXPR_NUM_THREADS='1',SPECIES_CACHE='0',METAG_CACHE=str(root/'no_external_cache'))
        if physical:
            (root/'independent_verification.json').unlink()
            result=subprocess.run([sys.executable,str(root/'verify_physical.py')],cwd=root,env=environment,
                capture_output=True,text=True,timeout=180)
            if result.returncode:raise RuntimeError(result.stdout+result.stderr)
            value=json.loads((root/'independent_verification.json').read_text())
            return {k:v for k,v in value.items() if k!='source_states'}
        result=subprocess.run([sys.executable,str(root/'cold_replay.py')],cwd=root,env=environment,
            capture_output=True,text=True,timeout=180)
        if result.returncode:raise RuntimeError(result.stdout+result.stderr)
        return json.loads(result.stdout.strip().splitlines()[-1])

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--packet',type=Path,default=Path(__file__).with_name('global_amino_reference_validation'))
    p.add_argument('--cold-replay',action='store_true');p.add_argument('--physical-replay',action='store_true')
    args=p.parse_args();result=audit(args.packet)
    if args.cold_replay:result['cold_replay']=replay(args.packet)
    if args.physical_replay:result['physical_replay']=replay(args.packet,physical=True)
    print(json.dumps(result,indent=2,allow_nan=False))

if __name__=='__main__':main()
