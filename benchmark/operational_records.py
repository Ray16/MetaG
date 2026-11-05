"""Authenticate operational molecular records without treating unfinished work as zero energy."""
from collections import Counter
from decimal import Decimal,localcontext
import hashlib
import json
import math
from pathlib import Path

read=lambda p:json.loads(Path(p).read_text())
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()


def finite(value):
    return isinstance(value,(int,float)) and not isinstance(value,bool) and math.isfinite(value)


def execution_state(result):
    """Complete means a declared opportunity ended, not that molecular QC passed."""
    if result is None:return dict(terminal=False,opportunity_complete=False,reason='not_started_or_pending')
    status=result['status']
    execution=result.get('controller_execution',result.get('execution',{}))
    if status=='external_resource_interruption' or execution.get('external_GPU_contention'):
        return dict(terminal=True,opportunity_complete=False,reason='external_resource_interruption')
    if execution.get('common_deadline_censored'):
        return dict(terminal=True,opportunity_complete=False,reason='common_deadline_censored')
    if 'controller_execution' not in result:
        return dict(terminal=False,opportunity_complete=False,reason='controller_has_not_authenticated_completion')
    if status not in {'molecular_record_complete','workflow_no_estimate','execution_failure','execution_timeout','predictor_execution_failure'}:
        raise ValueError('Unknown terminal outcome: '+str(status))
    if status=='molecular_record_complete' and (execution.get('timed_out') or execution.get('exit_code')!=0):
        raise ValueError('Successful molecular record has failed controller execution')
    return dict(terminal=True,opportunity_complete=True,reason=status)


def authenticate(task,result,phase,expected_species_settings,expected_producer_model):
    if not execution_state(result)['opportunity_complete'] or result['status']!='molecular_record_complete':
        raise ValueError('No authenticated complete molecular record')
    phase=Path(phase).resolve();name=Path(result['cache_path']).name
    path=phase/'cache'/task['slot']/name
    if path.is_symlink() or not path.is_file() or sha(path)!=result['cache_sha256']:
        raise ValueError('Raw record path/hash mismatch')
    if not finite(result['energy_kj']) or not finite(result['sampling_sigma_kj']) or result['sampling_sigma_kj']<0:
        raise ValueError('Invalid result energy/sampling value')
    raw=read(path)
    if raw['G']!=result['energy_kj'] or raw['sigma']!=result['sampling_sigma_kj']:
        raise ValueError('Raw/result numerical mismatch')
    producer=raw.get('producer',raw.get('settings',{}).get('producer'))
    if producer!=result['producer'] or producer.get('model')!=expected_producer_model:
        raise ValueError('Raw/result physical producer mismatch')
    if task['kind']=='liquid_water':
        if raw.get('kind')!='liquid_water_reference':raise ValueError('Wrong water reference kind')
        parts=[raw[k] for k in ['electronic_kj','experimental_water_hydration_kj','thermal_kj','concentration_kj']]
        if any(not finite(v) for v in parts) or abs(math.fsum(parts)-raw['G'])>1e-7:
            raise ValueError('Water reference component arithmetic')
        thermo=raw['thermodynamics']
        if thermo['temperature_K']!=298.15 or thermo['pressure_Pa']!=101325. or raw['n_imag']!=0:
            raise ValueError('Water thermochemistry target/curvature mismatch')
        if abs(thermo['enthalpy_correction_kj']+thermo['minus_T_entropy_kj']-raw['thermal_kj'])>1e-7:
            raise ValueError('Water thermal arithmetic')
        if Counter(raw['geometry']['symbols'])!=Counter({'H':2,'O':1}):raise ValueError('Water geometry composition mismatch')
        points=[raw['geometry']]
        maximum=0.;n_minima=1
    else:
        from rdkit import Chem
        if raw['q']!=task['charge'] or raw['smi']!=task['smiles'] or raw['method']!='implicit':
            raise ValueError('Raw molecular identity mismatch')
        if raw['settings'].get('spin',1)!=task['spin']:
            raise ValueError('Raw spin mismatch')
        if any(raw['settings'].get(k)!=v for k,v in expected_species_settings.items()):
            raise ValueError('Raw numerical settings mismatch')
        ensemble=raw['meta']['thermal']['ensemble'];minima=ensemble['validated_minima']
        if not minima or ensemble['n_reference_thermal_approximation']!=0:
            raise ValueError('Missing full validated thermal partition')
        expected_atoms=Counter(a.GetSymbol() for a in Chem.AddHs(Chem.MolFromSmiles(task['smiles'])).GetAtoms())
        if any(Counter(m['geometry']['symbols'])!=expected_atoms for m in minima):
            raise ValueError('Stored geometry element composition mismatch')
        with localcontext() as context:
            context.prec=50;D=lambda x:Decimal(str(x));rt=D('2.4789');totals=[]
            for minimum in minima:
                if any(not finite(minimum[k]) for k in ['electronic_kj','solvation_kj','thermal_kj','total_kj']):
                    raise ValueError('Nonfinite minimum')
                total=sum(D(minimum[k]) for k in ['electronic_kj','solvation_kj','thermal_kj'])
                if abs(total-D(minimum['total_kj']))>D('1e-7'):raise ValueError('Minimum component mismatch')
                totals.append(total)
            low=min(totals);partition=low-rt*sum((-(g-low)/rt).exp() for g in totals).ln()
            maximum=float(abs(partition-D(raw['G'])))
            if maximum>1e-7:raise ValueError('Raw partition mismatch')
        points=[m['geometry'] for m in minima];n_minima=len(minima)
    for point in points:
        body=dict(point);claimed=body.pop('sha256')
        if hashlib.sha256(json.dumps(body,sort_keys=True,allow_nan=False).encode()).hexdigest()!=claimed:
            raise ValueError('Geometry digest mismatch')
        if body['charge']!=task['charge'] or body['spin_multiplicity']!=task['spin']:
            raise ValueError('Geometry charge/spin mismatch')
        xyz=body['positions_angstrom']
        if len(xyz)!=task['n_atoms'] or any(len(v)!=3 or not all(finite(x) for x in v) for v in xyz):
            raise ValueError('Geometry dimensions/coordinates mismatch')
    return raw,dict(slot=task['slot'],kind=task['kind'],source_input_key=task['source_input_key'],
        cache_relative=str(path.relative_to(phase)),cache_sha256=sha(path),n_validated_minima=n_minima,
        maximum_decimal_partition_error_kj=maximum,producer=producer)


def collect(phase,tasks,expected_species_settings,expected_producer_model):
    phase=Path(phase);outcomes={};states={};accepted={};catalog=[]
    for task in tasks:
        path=phase/'results'/(task['slot']+'.json');result=read(path) if path.exists() else None
        state=execution_state(result);states[task['slot']]=state
        if result is not None:outcomes[task['slot']]=result
        if result is not None and state['opportunity_complete'] and result['status']=='molecular_record_complete':
            raw,row=authenticate(task,result,phase,expected_species_settings,expected_producer_model)
            key=task['source_input_key']
            if key in accepted:raise ValueError('Duplicate primitive key')
            accepted[key]=raw;catalog.append(row)
    return dict(outcomes=outcomes,opportunities=states,accepted=accepted,catalog=catalog,
        all_declared_opportunities_complete=all(s['opportunity_complete'] for s in states.values()),
        status_counts=dict(Counter(s['reason'] for s in states.values())))
