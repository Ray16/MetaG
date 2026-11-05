"""Reconstruct orientation diagnostics from all432unchanged native pairs."""
from pathlib import Path
from collections import defaultdict
import hashlib
import json
import math
import re

HERE=Path(__file__).resolve().parent
sha=lambda path:hashlib.sha256(Path(path).read_bytes()).hexdigest()
read=lambda name:json.loads((HERE/name).read_text())
number=r'[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][+-]?\d+)?'
def analyze():
    for name,digest in read('freeze.json')['sha256'].items():assert sha(HERE/name)==digest,name
    p=read('protocol.json');poses=read('poses.json');tasks=read('tasks.json')
    assert len(tasks)==432
    groups=defaultdict(dict);n_calls=0;maximum_energy_reconstruction_error=0.
    for task in tasks:
        key=task['key'];path=HERE/'results'/(key+'.json')
        assert path.exists(),'Complete orientation control withheld; missing '+key
        result=json.loads(path.read_text())
        assert result['passed'] and result['task']==task
        assert result['protocol_sha256']==sha(HERE/'protocol.json')
        assert result['worker_sha256']==sha(HERE/'worker.py')
        energies=[]
        for phase in ['gas','water']:
            path=HERE/'attempts'/key/(phase+'.json');assert sha(path)==result['native_record_sha256'][phase]
            native=json.loads(path.read_text());assert native['returncode']==0
            assert 'convergence criteria satisfied' in native['stdout']
            match=re.search(r'TOTAL ENERGY\s+('+number+r')\s+Eh',native['stdout']);assert match
            value=float(match[1].replace('D','E'));assert math.isfinite(value) and value==native['energy_Eh']
            energies.append(value);n_calls+=1
            argv=native['argv']
            for flag,expected in [('--gfn','2'),('--chrg',str(task['charge'])),('--uhf','0')]:assert argv[argv.index(flag)+1]==expected
            assert '--sp' in argv and '--norestart' not in argv
            assert ('--acc' in argv)==(task['accuracy']!='default')
            if '--acc' in argv:assert argv[argv.index('--acc')+1]=='0.001'
            assert ('--alpb' in argv)==(phase=='water')
            if phase=='water':
                assert argv[argv.index('--alpb')+1]=='water'
                grid=5810 if task['grid']=='extreme' else 230
                assert re.search(r'Grid points\s+'+str(grid)+r'\s+per atom',native['stdout'])
            lines=native['input_xyz'].splitlines();assert int(lines[0])==len(task['symbols'])
            rows=[line.split() for line in lines[2:]]
            assert [row[0] for row in rows]==task['symbols']
            assert all(abs(float(v)-target)<=.50000001e-6 for row,point in zip(rows,task['positions_angstrom']) for v,target in zip(row[1:],point))
        delta=(energies[1]-energies[0])*p['Hartree_to_kj_per_mol']
        maximum_energy_reconstruction_error=max(maximum_energy_reconstruction_error,abs(delta-result['solvation_kj']))
        assert abs(delta-result['solvation_kj'])<1e-10
        groups[(task['grid'],task['accuracy'],task['pose'])][task['transform']]=delta
    assert n_calls==864 and len(groups)==24
    reported={};native_replay={}
    for (grid,accuracy,pose),values in groups.items():
        assert len(values)==18
        proper=[v for name,v in values.items() if name!='inversion_control']
        mean=math.fsum(proper)/len(proper)
        key=f'{grid}|{accuracy}|{pose}'
        reported[key]=dict(n_proper_orientations=17,range_kj=max(proper)-min(proper),
            maximum_absolute_identity_difference_kj=max(abs(v-values['identity']) for v in proper),
            orientation_mean_kj=mean,orientation_SD_kj=math.sqrt(math.fsum((v-mean)**2 for v in proper)/len(proper)),
            inversion_minus_identity_kj=values['inversion_control']-values['identity'],
            rotational_gate_passed=max(proper)-min(proper)<=p['rotational_invariance_tolerance_kj'])
        if grid==accuracy=='default' and poses[pose]['original_solvation_kj'] is not None:
            native_replay[pose]=values['identity']-poses[pose]['original_solvation_kj']
            assert abs(native_replay[pose])<=p['native_replay_tolerance_kj']
    cycles={}
    for label,weights in [('full_core_exchange',p['fixed_pose_exchange_weights']),('vertical_charge_cycle',p['vertical_charge_cycle_weights'])]:
        for grid in ['default','extreme']:
            for accuracy in ['default','0.001']:
                transforms=[name for name in groups[(grid,accuracy,next(iter(weights)))] if name!='inversion_control']
                values={name:math.fsum(coefficient*groups[(grid,accuracy,pose)][name] for pose,coefficient in weights.items()) for name in transforms}
                low=math.fsum(min(coefficient*v for name,v in groups[(grid,accuracy,pose)].items() if name!='inversion_control') for pose,coefficient in weights.items())
                high=math.fsum(max(coefficient*v for name,v in groups[(grid,accuracy,pose)].items() if name!='inversion_control') for pose,coefficient in weights.items())
                cycles[f'{label}|{grid}|{accuracy}']=dict(same_rotation_values_kj=values,
                    same_rotation_range_kj=max(values.values())-min(values.values()),
                    independently_selected_sampled_orientation_interval_kj=[low,high],
                    interval_scope='Extrema over the sampled17orientations; not a bound over allSO3or a physical error/uncertainty interval.')
    return dict(scope=p['scope'],completed_pairs=432,n_native_executions=864,
        maximum_energy_reconstruction_error_kj=maximum_energy_reconstruction_error,
        original_production_replay_difference_kj=native_replay,per_pose=reported,balanced_cycles=cycles,
        no_production_profile_adopted=True,new_experimental_labels_admitted=0,
        interpretation=p['interpretation'])

if __name__=='__main__':
    result=analyze();(HERE/'result.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps(result,indent=2))
