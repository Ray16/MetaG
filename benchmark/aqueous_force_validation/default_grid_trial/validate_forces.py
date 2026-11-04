"""Finite-difference validation of an ALPB total-energy solvent gradient."""
from pathlib import Path
import hashlib
import json
import math
import os
import re
import subprocess
import tempfile
import time
import traceback

HERE=Path(__file__).resolve().parent
protocol=json.loads((HERE/'protocol.json').read_text())
assert not (HERE/'owner.json').exists(),'Inspect the existing run before restarting'
os.sched_setaffinity(0,set(protocol['resource_ceiling']['cpu_affinity']))
assert max(os.getloadavg())+1<=os.cpu_count()
env=dict(os.environ,OMP_NUM_THREADS='1',OMP_THREAD_LIMIT='1',MKL_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',NUMEXPR_NUM_THREADS='1',OMP_STACKSIZE='1G')
os.environ.update({k:v for k,v in env.items() if k.endswith('NUM_THREADS') or k=='OMP_THREAD_LIMIT'})
import numpy as np
from rdkit import Chem,rdBase
from rdkit.Chem import AllChem

def write(path,value):
    p=Path(path);t=p.with_name(p.name+'.tmp');t.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n');t.replace(p)
owner=dict(pid=os.getpid(),host=os.uname().nodename,start_ticks=Path(f'/proc/{os.getpid()}/stat').read_text().rsplit(')',1)[1].split()[19],cpu_affinity=sorted(os.sched_getaffinity(0)),rdkit_version=rdBase.rdkitVersion,source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
assert hashlib.sha256(Path(protocol['xtb']).read_bytes()).hexdigest()==protocol['xtb_sha256']
write(HERE/'owner.json',owner)
start=time.monotonic();deadline=start+protocol['resource_ceiling']['wall_seconds'];count=0
number=r'[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][+-]?\d+)?'

def single(symbols,xyz,charge,solvent,gradient,stream):
    global count
    remaining=deadline-time.monotonic()
    if remaining<=0:raise TimeoutError('Declared CPU validation budget exhausted')
    with tempfile.TemporaryDirectory(prefix='metag-alpb-force-') as directory:
        d=Path(directory);text=str(len(symbols))+'\nmatched solvent derivative\n'+''.join(s+' '+' '.join(format(float(v),'.17g') for v in row)+'\n' for s,row in zip(symbols,xyz))
        (d/'input.xyz').write_text(text)
        argv=[protocol['xtb'],str(d/'input.xyz'),'--gfn','2','--chrg',str(charge),'--uhf','0','--acc',str(protocol['accuracy']),'--norestart','--grad' if gradient else '--sp']
        if solvent:argv+=['--alpb','water','gsolv']
        process=subprocess.Popen(argv,cwd=d,env=env,stdin=subprocess.DEVNULL,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True)
        try:stdout,stderr=process.communicate(timeout=min(20,remaining))
        except subprocess.TimeoutExpired:
            process.kill();stdout,stderr=process.communicate();raise TimeoutError('xTB point exceeded bounded execution time')
        count+=1
        record=dict(index=count,charge=charge,solvent=solvent,gradient_requested=gradient,input_xyz=text,argv=argv,pid=process.pid,returncode=process.returncode,stdout=stdout,stderr=stderr)
        match=re.search(r'TOTAL ENERGY\s+('+number+r')\s+Eh',stdout)
        record['energy_Eh']=float(match.group(1).replace('D','E')) if match else None
        if gradient and (d/'gradient').exists():record['gradient_file']=(d/'gradient').read_text()
        stream.write(json.dumps(record,allow_nan=False)+'\n');stream.flush()
        assert process.returncode==0 and record['energy_Eh'] is not None,(process.returncode,stderr[-500:])
        grad=None
        if gradient:
            lines=record['gradient_file'].splitlines();head=next(i for i,line in enumerate(lines) if 'cycle' in line.lower())
            coordinates=[line.split() for line in lines[head+1:head+1+len(symbols)]]
            assert [row[3].capitalize() for row in coordinates]==symbols
            reported=np.array([[float(v.replace('D','E')) for v in row[:3]] for row in coordinates])*protocol['Bohr_to_A']
            assert np.max(np.abs(reported-xyz))<1e-8,'Gradient coordinate mismatch'
            grad=np.array([[float(v.replace('D','E')) for v in line.split()] for line in lines[head+1+len(symbols):head+1+2*len(symbols)]])
            assert grad.shape==(len(symbols),3) and np.isfinite(grad).all()
        return record['energy_Eh'],grad

def correction(symbols,xyz,charge,gradient,stream):
    eg,gg=single(symbols,xyz,charge,False,gradient,stream)
    ea,ga=single(symbols,xyz,charge,True,gradient,stream)
    energy=(ea-eg)*protocol['Hartree_to_eV']
    force=-(ga-gg)*protocol['Hartree_to_eV']/protocol['Bohr_to_A'] if gradient else None
    return energy,force

results=[]
try:
    for case in protocol['fixtures']:
        molecule=Chem.AddHs(Chem.MolFromSmiles(case['smiles']))
        assert Chem.GetFormalCharge(molecule)==case['charge']
        parameters=AllChem.ETKDGv3();parameters.randomSeed=protocol['conformer_seed']
        assert AllChem.EmbedMolecule(molecule,parameters)==0
        assert AllChem.MMFFOptimizeMolecule(molecule,mmffVariant='MMFF94s',maxIters=1000)==0
        symbols=[a.GetSymbol() for a in molecule.GetAtoms()];xyz=np.array(molecule.GetConformer().GetPositions())
        print('START '+case['name'],flush=True)
        with (HERE/(case['name']+'.raw.jsonl')).open('x') as stream:
            energy,force=correction(symbols,xyz,case['charge'],True,stream)
            estimates=[];reports=[]
            for step in protocol['finite_difference_steps_A']:
                numeric=np.zeros_like(xyz)
                for atom in range(len(symbols)):
                    for axis in range(3):
                        plus=xyz.copy();minus=xyz.copy();plus[atom,axis]+=step;minus[atom,axis]-=step
                        ep,_=correction(symbols,plus,case['charge'],False,stream)
                        em,_=correction(symbols,minus,case['charge'],False,stream)
                        numeric[atom,axis]=-(ep-em)/(2*step)
                difference=numeric-force
                reports.append(dict(step_A=step,max_error_eV_A=float(np.max(np.abs(difference))),RMS_error_eV_A=float(np.sqrt(np.mean(difference**2)))))
                estimates.append(numeric)
        step_difference=float(np.max(np.abs(estimates[0]-estimates[1])))
        net_force=float(np.linalg.norm(force.sum(axis=0)))
        torque=float(np.linalg.norm(np.cross(xyz-xyz.mean(axis=0),force).sum(axis=0)))
        passed=all(r['max_error_eV_A']<=protocol['maximum_force_component_error_eV_A'] and r['RMS_error_eV_A']<=protocol['RMS_force_error_eV_A'] for r in reports) and step_difference<=protocol['maximum_step_difference_eV_A'] and net_force<=protocol['maximum_net_force_eV_A'] and torque<=protocol['maximum_net_torque_eV']
        result=dict(case=case,symbols=symbols,positions_A=xyz.tolist(),solvent_delta_eV=energy,solvent_force_eV_A=force.tolist(),finite_differences=reports,step_difference_eV_A=step_difference,net_force_eV_A=net_force,net_torque_eV=torque,passed=passed)
        results.append(result);write(HERE/(case['name']+'.json'),result);write(HERE/'partial_results.json',results)
        print(case['name'],json.dumps(reports),'passed='+str(passed),flush=True)
        assert passed,'Derivative validation gate failed; preserve all raw outputs'
    write(HERE/'result.json',dict(passed=True,cases=results,n_xtb_executions=count,elapsed_seconds=time.monotonic()-start,scope=protocol['scope'],producer=owner))
except BaseException as error:
    write(HERE/'failure.json',dict(error=str(error),traceback=traceback.format_exc(),completed_cases=results,n_xtb_executions=count,elapsed_seconds=time.monotonic()-start));raise
