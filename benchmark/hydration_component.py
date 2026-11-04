"""Offline, unfiltered analysis of a source-defined hydration component panel.

These are fixed external geometries and solvent increments, not an end-to-end
MetaG reaction benchmark. Only reference-eligible compounds enter error scoring.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import random
import xml.etree.ElementTree as ET


def verify_manifest(packet):
    packet = Path(packet)
    manifest = json.loads((packet/'manifest.json').read_text())
    actual = {str(p.relative_to(packet)) for p in packet.rglob('*')
              if p.is_file() and '__pycache__' not in p.parts}
    if actual != set(manifest['files'])|{'manifest.json'}:
        raise ValueError('Evidence inventory changed')
    for name, entry in manifest['files'].items():
        path = packet/name
        if Path(name).is_absolute() or '..' in Path(name).parts or path.is_symlink():
            raise ValueError('Invalid evidence path')
        data = path.read_bytes()
        if len(data)!=entry['bytes'] or hashlib.sha256(data).hexdigest()!=entry['sha256']:
            raise ValueError(f'Changed evidence bytes: {name}')
    return len(manifest['files'])


def aggregate(values, method, rt):
    """Stable finite-conformer aggregation of energies in kJ/mol."""
    if not values or not all(math.isfinite(v) for v in values):
        raise ValueError('A complete finite energy list is required')
    minimum = min(values)
    shifted = [v - minimum for v in values]
    weights = [math.exp(-v / rt) for v in shifted]
    if method == 'minimum':
        return minimum
    if method == 'boltzmann_mean':
        return minimum + math.fsum(v * w for v, w in zip(shifted, weights)) / math.fsum(weights)
    if method == 'log_partition':
        return minimum - rt * math.log(math.fsum(weights))
    raise ValueError(f'Unknown aggregation: {method}')


def metrics(errors):
    return dict(n=len(errors), MAE_kj=math.fsum(map(abs, errors))/len(errors),
        RMSE_kj=math.sqrt(math.fsum(v*v for v in errors)/len(errors)),
        bias_kj=math.fsum(errors)/len(errors), max_abs_error_kj=max(map(abs, errors)),
        n_abs_error_gt10=sum(abs(v)>10 for v in errors),
        n_abs_error_gt20=sum(abs(v)>20 for v in errors))


def percentile(values, p):
    values = sorted(values)
    index = (len(values)-1)*p
    lo = math.floor(index)
    hi = math.ceil(index)
    return values[lo] + (index-lo)*(values[hi]-values[lo])


def analyze(packet):
    packet = Path(packet)
    if (packet/'manifest.json').exists():
        verify_manifest(packet)
    read = lambda p: json.loads((packet/p).read_text())
    sha = lambda p: hashlib.sha256((packet/p).read_bytes()).hexdigest()
    protocol = read('protocol.json')
    for name, digest in protocol['source_sha256'].items():
        if sha(name) != digest:
            raise ValueError(f'Changed source: {name}')
    references = read('references.json')
    # Reconstruct every target from the original measurement table, including
    # acid-recipe footnotes. Do not substitute a curated secondary-table value.
    table2 = ET.parse(packet/'Endo2023.xml').find('.//table-wrap[@id="tbl2"]')
    source_rows = {}
    for tr in table2.findall('.//tbody/tr'):
        cells = [' '.join(''.join(cell.itertext()).split()) for cell in tr]
        source_rows[cells[0]] = (float(cells[7].replace('–','-').replace('−','-')),
            any(x.get('rid')=='t2fn3' for x in tr.iter('xref')))
    if len(source_rows) != 21 or len(references) != 21:
        raise ValueError('Original source cohort changed')
    for ref in references:
        logk, acid = source_rows[ref['source_label']]
        if logk != ref['source_log10_Kaw']:
            raise ValueError('Source equilibrium constant mismatch')
        target = protocol['R_kj']*protocol['temperature_K']*math.log(10)*logk
        if abs(target-ref['target_molar_hydration_kj'])>1e-12:
            raise ValueError('Molar target convention mismatch')
        if acid != (ref['aqueous_recipe']=='1 mM HCl'):
            raise ValueError('Source recipe mismatch')
        eligible = ref['covalent_identity_matches_source'] and not acid
        if eligible != ref['component_accuracy_eligible']:
            raise ValueError('Eligibility differs from source condition/identity rule')
    admitted = {r['name']:r for r in references if r['component_accuracy_eligible']}
    if sorted(admitted) != protocol['expected_compounds']:
        raise ValueError('Eligible reference set changed')
    inputs = read('inputs.json')
    expected = set(protocol['expected_keys'])
    selected = [r for r in inputs if r['key'] in expected]
    if len(selected) != len(expected):
        raise ValueError('Input coverage changed')
    predictions = {}
    batch_checks = []
    for path in sorted((packet/'predictions').glob('*.json')):
        doc = json.loads(path.read_text())
        if doc['protocol_sha256'] != sha('protocol.json'):
            raise ValueError('Prediction protocol mismatch')
        for r in doc['records']:
            if r['key'] in predictions or r['key'] not in expected:
                raise ValueError('Duplicate/unexpected single point')
            if not math.isfinite(r['energy_eV']):
                raise ValueError('Nonfinite single point')
            if doc['producer_sha256'] != sha('predict_component.py'):
                raise ValueError('Single-point producer mismatch')
            predictions[r['key']] = r
            if 'individual_minus_batch_kj' in r:
                batch_checks.append(r['individual_minus_batch_kj'])
    if set(predictions) != expected:
        raise ValueError('Incomplete predictions; do not score a silently reduced cohort')
    if len(batch_checks)!=3 or max(map(abs,batch_checks))>protocol['inference']['max_individual_batch_difference_kj']:
        raise ValueError('Batch/single numerical gate failed or missing')
    grouped = {}
    for row in selected:
        pred = predictions[row['key']]
        if pred['coordinate_sha256'] != row['coordinate_sha256']:
            raise ValueError('Coordinate identity changed')
        if sha(row['coordinate_file']) != row['coordinate_sha256']:
            raise ValueError('Changed coordinate bytes')
        xyz = (packet/row['coordinate_file']).read_text().splitlines()
        parsed = [line.split() for line in xyz[2:] if line.strip()]
        if int(xyz[0])!=len(row['symbols']) or [r[0] for r in parsed]!=row['symbols']:
            raise ValueError('Ordered source atoms changed')
        if [[float(x) for x in r[1:]] for r in parsed]!=row['positions_A']:
            raise ValueError('Source coordinates changed')
        if len(pred['forces_eV_A'])!=len(row['symbols']) or any(
                len(v)!=3 or not all(math.isfinite(x) for x in v) for v in pred['forces_eV_A']):
            raise ValueError('Invalid atomic-force single-point output')
        values = {k:v*protocol['Eh_to_kj'] for k,v in row['published_energies_Eh'].items()}
        values['UMA'] = pred['energy_eV']*protocol['eV_to_kj']
        grouped.setdefault((row['name'],row['structuremode']), []).append(values)
    names = sorted(admitted)
    rt = protocol['R_kj']*protocol['temperature_K']
    result = dict(schema='metag-hydration-component-results-v1', scope=protocol['scope'],
        experimental_validation_regime='development_exposed_component_ablation',
        source_exposure_addendum=read('source_exposure_addendum.json'),
        protocol_sha256=sha('protocol.json'), n_source_compounds=len(references),
        n_reference_eligible=len(admitted), n_predicted=len(admitted),
        reference_eligible_coverage=1.0, original_source_panel_scored_fraction=len(admitted)/len(references),
        holds={r['name']:r['reason'] for r in references if not r['component_accuracy_eligible']},
        n_single_points=len(predictions), individual_batch_checks_kj=batch_checks,
        maximum_individual_batch_difference_kj=max(map(abs,batch_checks)),
        per_compound=[],arms={},primary_paired_contrasts={},
        zero_baseline=metrics([-admitted[n]['target_molar_hydration_kj'] for n in names]),
        no_fitted_offsets=True, primary_metabolic_references_admitted=0,
        uncertainty_scope=protocol['primary_specification']['uncertainty'])
    for electronic in protocol['electronic_arms']:
        for solvent in protocol['solvent_arms']:
            for geometry in protocol['geometry_arms']:
                for method in protocol['aggregation_arms']:
                    arm = '|'.join([electronic,solvent,geometry,method])
                    errors = []
                    for name in names:
                        gp = 'solv' if geometry == 'solv' else 'gas'
                        sp = 'gas' if geometry == 'gas' else 'solv'
                        gas = grouped[(name,gp)]
                        solv = grouped[(name,sp)]
                        value = aggregate([r[electronic]+r[solvent] for r in solv],method,rt)
                        value -= aggregate([r[electronic] for r in gas],method,rt)
                        target = admitted[name]['target_molar_hydration_kj']
                        error = value-target
                        errors.append(error)
                        result['per_compound'].append(dict(name=name,arm=arm,
                            hydration_kj=value,target_kj=target,error_kj=error,
                            gas_conformers=len(gas),solv_conformers=len(solv)))
                    result['arms'][arm] = metrics(errors)
    primary = 'UMA|alpb|full|minimum'
    table = {(r['arm'],r['name']):r for r in result['per_compound']}
    rng = random.Random(20261104)
    resamples = [[rng.randrange(len(names)) for _ in names] for _ in range(10000)]
    absolute = [abs(table[primary,n]['error_kj']) for n in names]
    boot = [math.fsum(absolute[i] for i in sample)/len(names) for sample in resamples]
    result['primary_arm'] = primary
    result['primary_MAE_descriptive_compound_bootstrap95_kj'] = [percentile(boot,.025),percentile(boot,.975)]
    for compare in ['UMA|smd|full|minimum','el_wb97m-v_def2-tzvppd|alpb|full|minimum']:
        differences = [abs(table[compare,n]['error_kj'])-abs(table[primary,n]['error_kj']) for n in names]
        boot = [math.fsum(differences[i] for i in sample)/len(names) for sample in resamples]
        result['primary_paired_contrasts'][compare] = dict(
            comparison_minus_primary_MAE_kj=math.fsum(differences)/len(names),
            descriptive_compound_bootstrap95_kj=[percentile(boot,.025),percentile(boot,.975)],
            n_improved=sum(d<0 for d in differences),n_worsened=sum(d>0 for d in differences))
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--packet', type=Path, default=Path(__file__).with_name('hydration_component_validation'))
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    result = analyze(args.packet)
    if args.output:
        args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in {'per_compound','arms'}},indent=2))
    for key,value in result['arms'].items():
        if key.endswith('|full|minimum'):
            print(key,json.dumps(value))


if __name__ == '__main__':
    main()
