"""Reconstruct a failed complete-panel invariance control without hiding refusals."""
from collections import defaultdict
from pathlib import Path
import argparse
import hashlib
import json
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile
import zipfile

NUMBER=r'[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][+-]?\d+)?'
def sha(data):return hashlib.sha256(data).hexdigest()


def audit(packet):
    packet=Path(packet);read=lambda name:json.loads((packet/name).read_text())
    manifest=read('manifest.json')
    actual={str(p.relative_to(packet)) for p in packet.rglob('*') if p.is_file() and '__pycache__' not in p.parts}
    assert actual==set(manifest['files'])|{'manifest.json'}
    for name,entry in manifest['files'].items():
        path=packet/name;assert not Path(name).is_absolute() and '..' not in Path(name).parts and not path.is_symlink()
        data=path.read_bytes();assert len(data)==entry['bytes'] and sha(data)==entry['sha256']
    p=read('protocol.json');tasks=read('tasks.json');poses=read('poses.json')
    assert len(tasks)==432 and len({t['key'] for t in tasks})==432 and len(poses)==6
    native_count=0;success=[];failures=[];groups=defaultdict(dict);replay={}
    with zipfile.ZipFile(packet/'raw_native.zip') as native,zipfile.ZipFile(packet/'results.zip') as results:
        names=set(native.namelist());result_names=set(results.namelist())
        assert len(names)==len(native.namelist()) and len(result_names)==len(results.namelist())
        assert all(not Path(name).is_absolute() and '..' not in Path(name).parts for name in names|result_names)
        for task in tasks:
            key=task['key'];prefix='attempts/'+key+'/'
            assert prefix+'owner.json' in names,'Unattempted task must not disappear'
            good=key+'.json' in result_names
            energies=[]
            for phase in ['gas','water']:
                path=prefix+phase+'.json'
                if path not in names:
                    assert not good and phase=='water';continue
                raw=native.read(path);row=json.loads(raw);native_count+=1
                assert row['phase']==phase
                argv=row['argv']
                for flag,value in [('--gfn','2'),('--chrg',str(task['charge'])),('--uhf','0')]:assert argv[argv.index(flag)+1]==value
                assert '--sp' in argv and '--norestart' not in argv
                assert ('--acc' in argv)==(task['accuracy']!='default')
                if '--acc' in argv:assert argv[argv.index('--acc')+1]=='0.001'
                assert ('--alpb' in argv)==(phase=='water')
                lines=row['input_xyz'].splitlines();assert int(lines[0])==len(task['symbols'])
                coordinates=[line.split() for line in lines[2:]]
                assert [r[0] for r in coordinates]==task['symbols']
                assert all(abs(float(v)-target)<=.50000001e-6 for record,xyz in zip(coordinates,task['positions_angstrom']) for v,target in zip(record[1:],xyz))
                match=re.search(r'TOTAL ENERGY\s+('+NUMBER+r')\s+Eh',row['stdout'])
                if good:
                    assert row['returncode']==0 and match and 'convergence criteria satisfied' in row['stdout']
                    energy=float(match[1].replace('D','E'));assert math.isfinite(energy) and energy==row['energy_Eh']
                    energies.append(energy)
                    if phase=='water':
                        assert argv[argv.index('--alpb')+1]=='water'
                        assert re.search(r'Grid points\s+'+('5810' if task['grid']=='extreme' else '230')+r'\s+per atom',row['stdout'])
                else:
                    assert phase=='gas' and row['returncode']!=0 and row['energy_Eh'] is None and match is None
                    assert 'Self consistent charge iterator did not converge' in row['stdout']
                    assert prefix+'failure.json' in names
            if not good:
                failures.append(dict(key=key,pose=task['pose'],transform=task['transform'],grid=task['grid'],accuracy=task['accuracy'],
                    first_causal_failure='Gas-phase self-consistent-charge nonconvergence',unattempted_downstream='water'));continue
            result=json.loads(results.read(key+'.json'))
            assert result['task']==task and result['passed']
            assert result['worker_sha256']==sha((packet/'worker.py').read_bytes())
            assert result['protocol_sha256']==sha((packet/'protocol.json').read_bytes())
            for phase in ['gas','water']:assert sha(native.read(prefix+phase+'.json'))==result['native_record_sha256'][phase]
            delta=(energies[1]-energies[0])*p['Hartree_to_kj_per_mol']
            assert abs(delta-result['solvation_kj'])<1e-10
            groups[(task['grid'],task['accuracy'],task['pose'])][task['transform']]=delta;success.append(key)
    assert len(success)==426 and len(failures)==6 and native_count==858
    assert {r['pose'] for r in failures}=={'charged_ox'}
    diagnostics={}
    for grid in ['default','extreme']:
        for accuracy in ['default','0.001']:
            for pose in poses:
                values=groups[(grid,accuracy,pose)];proper={k:v for k,v in values.items() if k!='inversion_control'}
                complete=len(proper)==17
                span=max(proper.values())-min(proper.values()) if complete else None
                diagnostics[f'{grid}|{accuracy}|{pose}']=dict(n_finite_proper_orientations=len(proper),n_required_proper_orientations=17,
                    orientation_range_kj=span,energetic_endpoint_status='complete_subcase' if complete else 'withheld',
                    rotational_gate_passed=(span<=p['rotational_invariance_tolerance_kj']) if complete else False,
                    inversion_minus_identity_kj=values['inversion_control']-values['identity'] if {'inversion_control','identity'}<=set(values) else None)
                ref=poses[pose]['original_solvation_kj']
                if grid==accuracy=='default' and ref is not None:
                    replay[pose]=values['identity']-ref;assert abs(replay[pose])<=p['native_replay_tolerance_kj']
    cycles={}
    for name,weights in [('full_core',p['fixed_pose_exchange_weights']),('charge',p['vertical_charge_cycle_weights'])]:
        for grid in ['default','extreme']:
            for accuracy in ['default','0.001']:
                if any(len([k for k in groups[(grid,accuracy,pose)] if k!='inversion_control'])!=17 for pose in weights):
                    cycles[f'{name}|{grid}|{accuracy}']=dict(status='withheld_due_failed_native_pairs');continue
                orientations=[k for k in groups[(grid,accuracy,next(iter(weights)))] if k!='inversion_control']
                values={k:math.fsum(coefficient*groups[(grid,accuracy,pose)][k] for pose,coefficient in weights.items()) for k in orientations}
                low=math.fsum(min(coefficient*v for k,v in groups[(grid,accuracy,pose)].items() if k!='inversion_control') for pose,coefficient in weights.items())
                high=math.fsum(max(coefficient*v for k,v in groups[(grid,accuracy,pose)].items() if k!='inversion_control') for pose,coefficient in weights.items())
                cycles[f'{name}|{grid}|{accuracy}']=dict(status='complete_predeclared_subcase',same_rotation_range_kj=max(values.values())-min(values.values()),
                    independently_selected_sampled_orientation_interval_kj=[low,high],sampled_interval_width_kj=high-low,
                    scope='Extrema over sampled orientations only; not a bound over allSO3or total Gibbs error.')
    return dict(status='negative_control_reconstructed',requested_pairs=432,successful_pairs=426,failed_pairs=6,unattempted_pairs=0,
        actual_native_executions=native_count,unattempted_downstream_water_calls=6,full_frozen_complete_panel_endpoint_available=False,
        primary_analysis_gate_weakened=False,failures=failures,complete_subcase_diagnostics=diagnostics,cycle_diagnostics=cycles,
        original_production_replay_difference_kj=replay,new_experimental_references_admitted=0,production_numerical_policy_adopted=False,
        scope='Numerical orientation/control evidence only. Missing charged-oxidized values remain withheld; no full-panel energetic endpoint or experimental accuracy claim.')


def strict_endpoint_replay(packet):
    packet=Path(packet)
    with tempfile.TemporaryDirectory(prefix='metag-rotation-original-gate-') as directory:
        root=Path(directory)
        frozen=json.loads((packet/'freeze.json').read_text())
        for name in list(frozen['sha256'])+['freeze.json']:
            target=root/name;target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(packet/name,target)
        for archive_name,prefix in [('raw_native.zip',''),('results.zip','results/')]:
            with zipfile.ZipFile(packet/archive_name) as archive:
                entries=archive.infolist()
                assert sum(e.file_size for e in entries)<50_000_000
                for entry in entries:
                    path=Path(entry.filename)
                    assert not path.is_absolute() and '..' not in path.parts and not entry.is_dir()
                    assert (entry.external_attr>>16)&0o170000!=0o120000
                    target=root/prefix/path;target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes(archive.read(entry))
        env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1',OMP_NUM_THREADS='1',OMP_THREAD_LIMIT='1',MKL_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',NUMEXPR_NUM_THREADS='1')
        result=subprocess.run([sys.executable,str(root/'analyze.py')],cwd=root,env=env,capture_output=True,text=True,timeout=60)
        assert result.returncode!=0 and 'Complete orientation control withheld; missing' in result.stderr,result.stderr
        assert not (root/'result.json').exists()
        return dict(status='original_complete_panel_endpoint_withheld',original_gate_unchanged=True,quantum_execution=False)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--packet',type=Path,default=Path(__file__).with_name('solvation_rotation_validation'))
    parser.add_argument('--strict-endpoint-replay',action='store_true')
    args=parser.parse_args();result=audit(args.packet)
    if args.strict_endpoint_replay:result['strict_endpoint_replay']=strict_endpoint_replay(args.packet)
    print(json.dumps(result,indent=2))
