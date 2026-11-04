"""Reconstruct solvent derivative checks from every retained xTB execution."""
from pathlib import Path
import argparse
import gzip
import hashlib
import json
import math
import re
import numpy as np

NUMBER=r'[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][+-]?\d+)?'
def need(condition,message):
    if not condition:raise ValueError(message)
def close(a,b,tolerance,message):need(math.isfinite(a) and math.isfinite(b) and abs(a-b)<=tolerance,message)
def read(path):return json.loads(Path(path).read_text())

def check_record(record,symbols,positions,charge,aqueous,gradient,protocol):
    need(record['returncode']==0,'Native calculation failed')
    need(record['charge']==charge and record['solvent']==aqueous,'Wrong phase or charge')
    need(record['gradient_requested']==gradient,'Wrong run type')
    argv=record['argv']
    for key,value in [('--gfn','2'),('--chrg',str(charge)),('--uhf','0'),('--acc',str(protocol['accuracy']))]:
        need(key in argv and argv[argv.index(key)+1]==value,'Wrong native option:'+key)
    need('--norestart' in argv,'Restart state not isolated')
    need(('--grad' if gradient else '--sp') in argv,'Wrong derivative option')
    if aqueous:
        need('--alpb' in argv and argv[argv.index('--alpb')+1:argv.index('--alpb')+3]==['water','gsolv'],'Wrong solvent convention')
        grid=5810 if protocol.get('sasa_grid')=='extreme' else 230
        need(re.search(r'Grid points\s+'+str(grid)+r'\s+per atom',record['stdout']) is not None,'Wrong recorded solvent grid')
    else:need('--alpb' not in argv,'Gas calculation unexpectedly solvated')
    lines=record['input_xyz'].splitlines();need(int(lines[0])==len(symbols),'Atom count')
    rows=[line.split() for line in lines[2:]]
    need([row[0] for row in rows]==symbols,'Input atom ordering')
    xyz=np.array([[float(v) for v in row[1:]] for row in rows])
    need(xyz.shape==positions.shape and np.max(np.abs(xyz-positions))<1e-14,'Wrong displaced coordinates')
    match=re.search(r'TOTAL ENERGY\s+('+NUMBER+r')\s+Eh',record['stdout'])
    need(match is not None and 'convergence criteria satisfied' in record['stdout'],'Missing converged energy')
    energy=float(match[1].replace('D','E'));close(energy,record['energy_Eh'],1e-12,'Recorded energy mismatch')
    return energy

def extract_gradient(text,symbols,positions,protocol):
    rows=[line.split() for line in text.strip().splitlines()]
    need(rows[0]==['$grad'] and rows[-1]==['$end'],'Gradient data-group boundaries')
    n=len(symbols);coords=rows[2:2+n];grad=rows[2+n:2+2*n]
    need(len(rows)==3+2*n and [r[3].capitalize() for r in coords]==symbols,'Gradient atom ordering')
    xyz=np.array([[float(x.replace('D','E')) for x in row[:3]] for row in coords])*protocol['Bohr_to_A']
    need(np.max(np.abs(xyz-positions))<1e-8,'Gradient coordinates mismatch')
    result=np.array([[float(x.replace('D','E')) for x in row] for row in grad])
    need(result.shape==positions.shape and np.isfinite(result).all(),'Invalid gradient array')
    return result

def audit(root):
    root=Path(root);protocol=read(root/'protocol.json')
    saved=read(root/'result.json') if (root/'result.json').exists() else dict(cases=read(root/'partial_results.json'))
    cases=saved['cases'];count=0;reports=[]
    for case in cases:
        symbols=case['symbols'];positions=np.array(case['positions_A']);charge=case['case']['charge']
        path=root/(case['case']['name']+'.raw.jsonl')
        if path.exists():lines=path.read_text().splitlines()
        else:
            with gzip.open(str(path)+'.gz','rt') as stream:lines=stream.read().splitlines()
        records=[json.loads(line) for line in lines];cursor=0
        def pair(xyz,gradient=False):
            nonlocal cursor,count
            gas,water=records[cursor:cursor+2];cursor+=2;count+=2
            eg=check_record(gas,symbols,xyz,charge,False,gradient,protocol)
            ew=check_record(water,symbols,xyz,charge,True,gradient,protocol)
            return (ew-eg)*protocol['Hartree_to_eV'],gas,water
        energy,gas,water=pair(positions,True)
        force=-(extract_gradient(water['gradient_file'],symbols,positions,protocol)-extract_gradient(gas['gradient_file'],symbols,positions,protocol))*protocol['Hartree_to_eV']/protocol['Bohr_to_A']
        close(energy,case['solvent_delta_eV'],1e-12,'Base solvent energy')
        need(np.max(np.abs(force-np.array(case['solvent_force_eV_A'])))<1e-12,'Base solvent force')
        numeric=[];errors=[]
        for step in protocol['finite_difference_steps_A']:
            array=np.zeros_like(positions)
            for index in range(positions.size):
                displacement=np.zeros_like(positions);displacement.flat[index]=step
                plus,_,_=pair(positions+displacement);minus,_,_=pair(positions-displacement)
                array.flat[index]=(minus-plus)/(2*step)
            numeric.append(array);delta=array-force
            errors.append(dict(step_A=step,max_error_eV_A=float(np.max(np.abs(delta))),RMS_error_eV_A=float(np.sqrt(np.mean(delta**2)))))
        need(cursor==len(records),'Unused or missing native results')
        for calculated,reported in zip(errors,case['finite_differences']):
            for key in calculated:close(calculated[key],reported[key],1e-12,'Derivative statistic mismatch')
        difference=float(np.max(np.abs(numeric[0]-numeric[1])))
        net=float(np.linalg.norm(np.sum(force,axis=0)))
        torque=float(np.linalg.norm(np.sum(np.cross(positions-positions.mean(axis=0),force),axis=0)))
        for value,key in [(difference,'step_difference_eV_A'),(net,'net_force_eV_A'),(torque,'net_torque_eV')]:close(value,case[key],1e-12,key)
        passed=all(r['max_error_eV_A']<=protocol['maximum_force_component_error_eV_A'] and r['RMS_error_eV_A']<=protocol['RMS_force_error_eV_A'] for r in errors) and difference<=protocol['maximum_step_difference_eV_A'] and net<=protocol['maximum_net_force_eV_A'] and torque<=protocol['maximum_net_torque_eV']
        need(passed==case['passed'],'Admission gate mismatch')
        reports.append(dict(name=case['case']['name'],native_executions=len(records),finite_differences=errors,net_torque_eV=torque,passed=passed))
    if 'n_xtb_executions' in saved:need(count==saved['n_xtb_executions'],'Total execution count mismatch')
    return dict(scope='Raw-result numerical replay, including declared failed gates. No backend execution, molecular Gibbs accuracy or production acceptance.',cases=reports,n_native_executions=count,all_cases_passed=all(r['passed'] for r in reports),n_cases=len(reports))

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('root',type=Path);p.add_argument('--output',type=Path,required=True);args=p.parse_args()
    result=audit(args.root);args.output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
