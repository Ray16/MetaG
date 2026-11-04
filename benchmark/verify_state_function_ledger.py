"""Independently reassemble retained state potentials and reaction energies.

This verifier checks supplied-state arithmetic and chemical bookkeeping. It
does not establish aqueous-state completeness, parameter accuracy, or prediction
uncertainty. Raw-record checks require the caller's original index and files.
"""
import argparse
from collections import Counter
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path

from rdkit import Chem


def digest(value):
    return hashlib.sha256(json.dumps(value,sort_keys=True,allow_nan=False).encode()).hexdigest()


def canonical(smiles):
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        raise ValueError('Invalid retained molecular structure')
    return Chem.MolToSmiles(molecule,isomericSmiles=True)


def composition(smiles):
    molecule = Chem.AddHs(Chem.MolFromSmiles(smiles))
    return Counter(atom.GetSymbol() for atom in molecule.GetAtoms())


def close(actual, expected, context, tolerance=1e-7):
    if any(isinstance(v,bool) or not isinstance(v,(int,float)) or not math.isfinite(v)
           for v in (actual,expected)) or abs(actual-expected)>tolerance:
        raise ValueError(context+': '+str(actual)+' differs from '+str(expected))


def partition(values, kt):
    low = min(values)
    return low-kt*math.log(math.fsum(math.exp(-(v-low)/kt) for v in values))


def independent_phi(sites, ph, rt_ln10):
    values = []
    for site in sites:
        side,pka = site[:2]
        kind = site[2] if len(site)>2 else 'acid'
        mult = site[3] if len(site)>3 else 1.
        if side!='prod' or kind not in {'acid','base'} or mult<=0:
            raise ValueError('State transform must have owned, product-side terms')
        exponent = ph-pka if kind=='acid' else pka-ph
        values.append(-mult*rt_ln10*(max(exponent,0.)+math.log1p(10.**-abs(exponent))/math.log(10.)))
    return math.fsum(values)


def macro_phi(model, ph, rt_ln10):
    weights = [0.]
    for pka in model['pKas']:
        weights.append(weights[-1]+ph-pka)
    pivot = max(weights)
    log_partition = pivot+math.log10(math.fsum(10.**(v-pivot) for v in weights))
    return -rt_ln10*(log_partition-weights[model['reference_deprotonations']])


def verify(payload, *, raw_root=None, raw_index=None, protocol=None, cycles=None):
    if payload['schema']!='metag-owned-state-function-ledger-v1':
        raise ValueError('Unknown retained-state ledger schema')
    if protocol is not None and payload['arm'] not in protocol['arms']:
        raise ValueError('Intervention absent from declared protocol')
    ledger = payload['ledger']
    if ledger['config'].get('state_function_policy')!='full-molecule-species-potential-v2':
        raise ValueError('Unsupported state-function physics version')
    conditions = payload['conditions']
    if any(conditions.get(k)!=v for k,v in {'T_K':298.15,'pH':7.,'ionic_strength_M':0.,
                                         'pMg':None,'water_activity':1.}.items()):
        raise ValueError('Unsupported conditions for the frozen v2 verifier')
    close(payload['proton_energy_kj'],-26.3-1104.5-2.303*8.314e-3*298.15*7.,'v2 proton convention')
    close(payload['std_state_kj'],8.314e-3*298.15*math.log(.082057*298.15),'v2 solute standard state')
    close(payload['KT_kj'],2.4789,'v2 pool temperature convention',1e-12)
    checked = dict(ledger);claimed = checked.pop('sha256')
    if digest(checked)!=claimed:
        raise ValueError('Ledger content identity mismatch')
    raw = {}
    authenticated = 0
    for row in ledger['raw_molecules']:
        charge,smi = row['charge'],canonical(row['smiles'])
        key = digest({'charge':charge,'smiles':smi})
        if key in raw:
            raise ValueError('Repeated raw molecular identity')
        value,sigma,details = row['record']
        close(value,value,'raw molecular energy')
        close(sigma,sigma,'raw sampling diagnostic')
        if sigma<0:
            raise ValueError('Negative raw sampling diagnostic')
        raw[key] = {'charge':charge,'smiles':smi,'G':value,'sigma':sigma}
        if raw_root is not None:
            if raw_index is None:
                raise ValueError('Original raw index required for raw authentication')
            name = details['source_record']
            if Path(name).name!=name:
                raise ValueError('Raw record path escapes supplied root')
            path = Path(raw_root)/name
            expected = raw_index[name]['sha256']
            if (hashlib.sha256(path.read_bytes()).hexdigest()!=expected
                    or payload['source_records_sha256'][name]!=expected
                    or details['source_sha256']!=expected):
                raise ValueError('Original molecular record fingerprint mismatch')
            source = json.loads(path.read_text())
            if (source['settings']!=payload['raw_settings'] or source['method']!='implicit'
                    or source['q']!=charge or canonical(source['smi'])!=smi):
                raise ValueError('Original molecular record identity/settings mismatch')
            minima = source['meta']['thermal']['ensemble']['validated_minima']
            for minimum in minima:
                close(minimum['total_kj'],minimum['electronic_kj']+minimum['solvation_kj']
                      +minimum['thermal_kj'],'validated minimum components')
            expected_g = source['G']
            if payload['arm']=='lowest_conformer':
                expected_g = min(m['total_kj'] for m in minima)
            elif payload['arm']=='omit_alpb':
                expected_g = partition([m['electronic_kj']+m['thermal_kj'] for m in minima],payload['KT_kj'])
            close(value,expected_g,'raw molecular component intervention')
            close(sigma,source['sigma'],'raw sampling provenance')
            authenticated += 1
    states = ledger['states']
    ph = payload['conditions']['pH']
    rt_ln10 = 2.303*8.314e-3*payload['conditions']['T_K']
    hydration_kt = 8.314462618e-3*payload['conditions']['T_K']
    q_h = payload['proton_energy_kj']
    std = payload['std_state_kj']
    n_pools = 0
    for key,row in states.items():
        record = dict(row);claimed = record.pop('record_sha256')
        if row['state_id']!=key or digest(record)!=claimed:
            raise ValueError('Retained state content identity mismatch')
        plan = row['plan']
        expected_sampling = Counter()
        if plan.get('kind')=='reducing_sugar_equilibrium':
            n_pools += 1
            members = row['states']
            if (len({m['state_id'] for m in members})!=len(members)
                    or {m['smiles'] for m in members}!=set(plan['source_states'])):
                raise ValueError('Incomplete declared pool state basis')
            energies = []
            for member in members:
                child = states[member['state_id']]
                if child['record_sha256']!=member['record_sha256']:
                    raise ValueError('Pool constituent identity mismatch')
                close(member['potential_kj'],child['potential_kj'],'pool constituent energy')
                if child['plan']['input_charge']!=plan['input_charge']:
                    raise ValueError('Pool constituent charge convention mismatch')
                if composition(member['smiles'])!=composition(members[0]['smiles']):
                    raise ValueError('Pool constituents have different elemental formulas')
                energies.append(child['potential_kj'])
            if row.get('pool_intervention')=='lowest_sugar_state':
                if payload['arm']!='lowest_sugar_state':
                    raise ValueError('Lowest-state pool outside its declared intervention')
                expected = min(energies)
                winner = min(members,key=lambda r:(r['potential_kj'],r['smiles']))
                if row['selected_constituent']!=winner['state_id']:
                    raise ValueError('Wrong lowest-state pool constituent')
                weights = [float(m['state_id']==winner['state_id']) for m in members]
            else:
                expected = partition(energies,payload['KT_kj'])
                low = min(energies)
                weights = [math.exp(-(v-low)/payload['KT_kj']) for v in energies]
                norm = math.fsum(weights)
                for member,w in zip(members,weights):
                    close(member['conditional_population'],w/norm,'pool population',1e-12)
                weights = [w/norm for w in weights]
            for member,w in zip(members,weights):
                if w==0. and row.get('pool_intervention')=='lowest_sugar_state':
                    continue
                for sample_key,sample in states[member['state_id']]['sampling_inputs'].items():
                    expected_sampling[sample_key] += w*sample['sensitivity']
            close(row['potential_kj'],expected,'complete pool partition')
        else:
            if (Chem.GetFormalCharge(Chem.MolFromSmiles(plan['input_smiles']))!=plan['input_charge']
                    or Chem.GetFormalCharge(Chem.MolFromSmiles(plan['reference_smiles']))!=plan['reference_charge']):
                raise ValueError('Reference charges disagree with molecular structures')
            before,after = composition(plan['input_smiles']),composition(plan['reference_smiles'])
            after['H'] += plan['proton_exchange']
            if before!=after or plan['proton_exchange']!=plan['input_charge']-plan['reference_charge']:
                raise ValueError('Molecular reference conversion violates element/charge balance')
            if canonical(plan['reference_smiles'])=='O' and plan['reference_charge']==0:
                solution = ledger['liquid_water_kj_before_standard_state']+std
            else:
                raw_id = digest({'charge':plan['reference_charge'],
                                 'smiles':canonical(plan['reference_smiles'])})
                solution = raw[raw_id]['G']+std
                expected_sampling[raw_id] = 1.
            hydration = row['aqueous_diagnostics'].get('hydration_diagnostics',{}).get('state')
            if hydration:
                if (hydration['status']!='complete' or len(hydration['states'])!=hydration['n_unique_states']
                        or any(m['status']!='estimated' for m in hydration['states'])):
                    raise ValueError('Incomplete declared hydration basis')
                solution += partition([0.]+[m['relative_g_kj'] for m in hydration['states']],hydration_kt)
                shifts = [0.]+[m['relative_g_kj'] for m in hydration['states']]
                low = min(shifts)
                weights = [math.exp(-(v-low)/hydration_kt) for v in shifts]
                norm = math.fsum(weights)
                probabilities = [w/norm for w in weights]
                # Fixed v2 hydration calibration: d(delta_i)/dG_base=-a,
                # d(delta_i)/dG_hydrate=a. The water reference is deterministic.
                slope = .639
                expected_sampling[raw_id] = 1.-slope*(1.-probabilities[0])
                for hydrated,population in zip(hydration['states'],probabilities[1:]):
                    sample_key = digest({'charge':plan['reference_charge'],
                                         'smiles':canonical(hydrated['smiles'])})
                    expected_sampling[sample_key] += slope*population
            close(row['solution_reference_kj'],solution,'solution reference partition')
            phi = independent_phi(plan['pka_sites'],ph,rt_ln10)
            coupled = row.get('coupled_protonation')
            if coupled:
                for model_row in coupled['species'].values():
                    model = model_row['model']
                    if (canonical(model['reference_smiles'])!=canonical(plan['reference_smiles'])
                            or model['conditions']['T_K']!=payload['conditions']['T_K']
                            or model['conditions']['ionic_strength_M']!=payload['conditions']['ionic_strength_M']):
                        raise ValueError('Coupled model reference/conditions mismatch')
                    phi -= independent_phi(model_row['removed_independent_terms'],ph,rt_ln10)
                    phi += macro_phi(model,ph,rt_ln10)
            close(row['protonation_kj'],phi,'owned protonation transform')
            close(row['reference_proton_kj'],plan['proton_exchange']*q_h,'reference proton reservoir')
            close(row['potential_kj'],solution+phi+row['reference_proton_kj'],'molecular state potential')
        variance = []
        if set(expected_sampling)!=set(row['sampling_inputs']):
            raise ValueError('Retained sampling basis has missing or extra raw states')
        for sample_key,sample in row['sampling_inputs'].items():
            if sample_key not in raw:
                raise ValueError('Unknown molecular sampling identity')
            close(sample['sigma_kj'],raw[sample_key]['sigma'],'sampling input identity')
            if sample['charge']!=raw[sample_key]['charge'] or canonical(sample['smiles'])!=raw[sample_key]['smiles']:
                raise ValueError('Retained sampling identity changed structure or charge')
            close(sample['sensitivity'],expected_sampling[sample_key],'sampling partition derivative',1e-10)
            variance.append((sample['sensitivity']*sample['sigma_kj'])**2)
        close(row['sampling_sigma_kj'],math.sqrt(math.fsum(variance)),'state sampling propagation')
    available = {}
    for rid,row in payload['results'].items():
        if row['dG'] is None:
            if row['status']=='assembled':
                raise ValueError('Unavailable result declared assembled')
            continue
        if row['status']!='assembled':
            raise ValueError('Available result lacks assembled status')
        reaction = payload['inputs'][rid]
        if set(row['assignments'])!=set(reaction['species']):
            raise ValueError('Reaction has missing molecular assignments')
        terms = [reaction.get('n_Hplus',0)*q_h]
        sampling = Counter()
        element_balance = Counter()
        charge_balance = Fraction(str(reaction.get('n_Hplus',0)))
        for name,(coefficient,charge,smiles) in reaction['species'].items():
            for element,count in composition(smiles).items():
                element_balance[element] += Fraction(str(coefficient))*count
            charge_balance += Fraction(str(coefficient))*charge
            state = states[row['assignments'][name]]
            plan = state['plan']
            if plan['input_charge']!=charge:
                raise ValueError('Reaction assignment changed input charge')
            if name in reaction.get('species_state_pools',{}):
                if plan.get('kind')!='reducing_sugar_equilibrium':
                    raise ValueError('Declared pool replaced by one drawing')
                if composition(plan['source_states'][0])!=composition(smiles):
                    raise ValueError('Declared pool changed the input elemental formula')
            elif canonical(plan.get('input_smiles',''))!=canonical(smiles):
                raise ValueError('Reaction assignment changed molecular identity')
            terms.append(coefficient*state['potential_kj'])
            for sample_key,sample in state['sampling_inputs'].items():
                sampling[sample_key] += coefficient*sample['sensitivity']
        element_balance['H'] += Fraction(str(reaction.get('n_Hplus',0)))
        if charge_balance or any(element_balance.values()):
            raise ValueError('Assembled input reaction violates element/charge balance')
        expected = math.fsum(terms)
        close(row['dG_unrounded'],expected,'reaction assembly '+rid)
        close(row['dG'],round(expected,1),'reported reaction rounding '+rid)
        if row.get('sampling_sigma_kj') is not None:
            expected_sigma = math.sqrt(math.fsum((v*raw[k]['sigma'])**2 for k,v in sampling.items()))
            close(row['sampling_sigma_kj'],expected_sigma,'reaction sampling propagation '+rid)
        available[rid] = expected
    residuals = []
    if cycles is not None:
        for cycle in cycles['independent_cycles']:
            net = Counter()
            for rid,weight in cycle.items():
                reaction = payload['inputs'][rid]
                for name,(coefficient,charge,smiles) in reaction['species'].items():
                    identity = json.dumps({'charge':charge,'smiles':canonical(smiles),
                        'pool':reaction.get('species_state_pools',{}).get(name),
                        'protonation':reaction.get('species_protonation_models',{}).get(name)},sort_keys=True)
                    net[identity] += Fraction(str(weight))*Fraction(str(coefficient))
                net['explicit-proton'] += Fraction(str(weight))*Fraction(str(reaction.get('n_Hplus',0)))
            if any(net.values()):
                raise ValueError('Declared input cycle does not close')
            if set(cycle)<=set(available):
                residuals.append(math.fsum(w*available[r] for r,w in cycle.items()))
        if any(abs(v)>1e-7 for v in residuals):
            raise ValueError('Physical state ledger has a nonclosing cycle')
    return {'scope':'Independent retained-component arithmetic and declared-state bookkeeping; not experimental accuracy, complete aqueous sampling, or calibrated uncertainty.',
            'arm':payload['arm'],'n_raw_molecules':len(raw),'n_raw_records_authenticated':authenticated,
            'n_states_checked':len(states),'n_pools_checked':n_pools,'n_reactions_checked':len(available),
            'n_evaluable_cycles':len(residuals),'max_cycle_residual_kj':max(map(abs,residuals),default=None)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('ledger',type=Path)
    parser.add_argument('--raw-root',type=Path)
    parser.add_argument('--raw-index',type=Path)
    parser.add_argument('--protocol',type=Path)
    parser.add_argument('--cycles',type=Path)
    parser.add_argument('--output',type=Path)
    args = parser.parse_args()
    read = lambda p:json.loads(p.read_text()) if p else None
    payload = read(args.ledger)
    if args.protocol and hashlib.sha256(args.protocol.read_bytes()).hexdigest()!=payload['protocol_sha256']:
        raise ValueError('Protocol file identity mismatch')
    protocol = read(args.protocol)
    if args.raw_index and protocol and hashlib.sha256(args.raw_index.read_bytes()).hexdigest()!=protocol['raw_index_sha256']:
        raise ValueError('Original raw index identity mismatch')
    report = verify(payload,raw_root=args.raw_root,raw_index=read(args.raw_index),
                    protocol=protocol,cycles=read(args.cycles))
    if args.output:
        with args.output.open('x') as stream:
            json.dump(report,stream,indent=2);stream.write('\n')
    print(json.dumps(report,indent=2))


if __name__=='__main__':
    main()
