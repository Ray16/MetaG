#!/usr/bin/env python3
"""Audit a frozen nicotinamide fixed-geometry UMA/DFT × ALPB/SMD control.

Reconstructs every arm from electronic, solvent and retained RRHO terms. Does
not run a backend, fit a correction, or admit an experimental reaction label.
"""
from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
import math
from pathlib import Path

WEIGHTS = {'full_red': 1, 'full_ox': -1, 'core_red': -1, 'core_ox': 1}

def require(condition, message):
    if not condition:
        raise ValueError(message)

def read(path):
    return json.loads(Path(path).read_text())

def file_sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def close(value, expected, tolerance, context):
    require(math.isfinite(value) and math.isfinite(expected) and
            abs(value-expected) <= tolerance, context)

def audit_point(item, point, protocol):
    """Reject incomplete, basis-truncated or inconsistently summed SCF outputs."""
    require(point['passed'] is True, 'Point did not pass')
    require(point['geometry'] == item['geometry'], 'Point geometry/state mismatch')
    require(point['n_AO'] == item['n_AO'], 'AO count mismatch')
    require(point['n_electrons'] == item['n_electrons'], 'Electron count mismatch')
    require(point['overlap_minimum_eigenvalue'] >= protocol['overlap_minimum_eigenvalue'],
            'Overlap gate failed')
    for name in ['gas', 'water']:
        stage = point[name]
        require(stage['converged'] is True, f'{name}: SCF not converged')
        require(stage['n_molecular_orbitals'] == item['n_AO'], f'{name}: AO directions discarded')
        require(0 <= stage['orbital_gradient_norm'] <= protocol['stationary_gradient_maximum'],
                f'{name}: orbital gradient failed')
        summary = stage['scf_summary']
        close(summary['e2'], summary['coul']+summary['exc'], 1e-8, f'{name}: electronic sum')
        terms = [summary[k] for k in ['e1', 'e2', 'nuc']]
        if name == 'water':
            terms += [summary['e_solvent'], summary['e_cds']]
            close(summary['e_cds'], point['CDS_Eh'], 1e-12, 'CDS mismatch')
        close(stage['energy_Eh'], math.fsum(terms), 1e-8, f'{name}: total sum')
    close(point['SMD_delta_Eh'], point['water']['energy_Eh']-point['gas']['energy_Eh'],
          1e-10, 'SMD total difference mismatch (CDS is included once)')
    require(point['surface_discretization_method'] == 'SWIG', 'Unexpected SMD surface')

def pool(rows, kt):
    """Equal counting of the retained minima, with explicit configurational term."""
    totals = [math.fsum(row.values()) for row in rows]
    minimum = min(totals)
    shifted = [math.exp(-(value-minimum)/kt) for value in totals]
    partition = math.fsum(shifted)
    probabilities = [value/partition for value in shifted]
    components = {key: math.fsum(p*row[key] for p,row in zip(probabilities,rows))
                  for key in rows[0]}
    # p log(p) is evaluated only for nonzero probabilities.
    components['configurational_kj'] = kt*math.fsum(p*math.log(p) for p in probabilities if p)
    free_energy = minimum-kt*math.log(partition)
    close(math.fsum(components.values()), free_energy, 1e-7, 'Pool reconstruction')
    return dict(G_kj=free_energy, components_kj=components, probabilities=probabilities,
                n_retained_minima=len(rows))

def audit(root, raw_directory=None):
    root = Path(root)
    protocol = read(root/'fixed_geometry_protocol.json')
    inputs = read(root/'geometry_inputs.json')
    manifest = read(root/'launch_manifest.json')
    for relative, expected in manifest['owned_file_sha256'].items():
        require(file_sha(root/relative) == expected, f'Owned file changed: {relative}')
    require(protocol['no_empirical_correction_or_offset'] is True, 'Unfrozen correction')
    require(set(inputs['source_records']) == set(WEIGHTS), 'Unexpected chemical families')
    raw = {}
    for family, source in inputs['source_records'].items():
        path = Path(source['path']) if raw_directory is None else Path(raw_directory)/Path(source['path']).name
        require(file_sha(path) == source['sha256'], f'Raw record changed: {family}')
        raw[family] = read(path)
    replay = read(root/'uma_replay.json')
    require(replay['passed'] is True and set(replay['rows']) == set(inputs['inputs']),
            'Incomplete UMA replay')
    for value in raw.values():
        require(replay['producer'] == value['settings']['producer']['model'], 'UMA producer mismatch')
    smoke = read(root/'points/smoke.json')
    cpu = read(root/'cpu_water_smoke.json')
    require(smoke['passed'] is True, 'GPU numerical fixture failed')
    require(smoke['source_sha256'] == file_sha(root/'reference_point.py'), 'Smoke code changed')
    smoke_differences = {}
    for stage,key in [('gas','gas_Eh'),('water','water_total_Eh')]:
        difference = abs(smoke[stage]['energy_Eh']-cpu[key])*protocol['Hartree_to_kj_per_mol']
        close(difference, smoke['CPU_GPU_absolute_difference_kj'][key], 1e-12, 'Smoke difference mismatch')
        require(difference <= protocol['numerical_smoke_tolerance_kj'], 'CPU/GPU disagreement')
        smoke_differences[key] = difference
    family_items = {key: [] for key in WEIGHTS}
    points = {}
    residual = Counter()
    charge = 0
    for slot,item in inputs['inputs'].items():
        family = item['family']
        minimum = raw[family]['meta']['thermal']['ensemble']['validated_minima'][item['minimum_index']]
        require(item['geometry'] == minimum['geometry'], f'Raw geometry mismatch: {slot}')
        require(item['reference_smiles'] == minimum['reference_smiles'], f'Structure mismatch: {slot}')
        for key,value in item['original_components_kj'].items():
            close(value,minimum[key],1e-9,f'Original component mismatch: {slot}/{key}')
        original = item['original_components_kj']
        close(math.fsum(original[k] for k in ['electronic_kj','solvation_kj','thermal_kj']),
              original['total_kj'],1e-7,f'Original total mismatch: {slot}')
        row = replay['rows'][slot]
        require(row['geometry_sha256'] == item['geometry']['sha256'], f'Replay geometry mismatch: {slot}')
        close(row['original_electronic_kj'],original['electronic_kj'],1e-9,'Replay source energy mismatch')
        close(row['difference_kj'],row['replay_electronic_kj']-original['electronic_kj'],1e-9,
              'Replay difference mismatch')
        require(abs(row['difference_kj']) <= protocol['UMA_replay_max_difference_kj'], 'UMA replay failed')
        point = read(root/'points'/(slot+'.json'))
        audit_point(item,point,protocol)
        require(point['source_sha256'] == file_sha(root/'reference_point.py'), 'Point code changed')
        require(point['protocol_sha256'] == file_sha(root/'fixed_geometry_protocol.json'), 'Point protocol changed')
        require(point['input_sha256'] == file_sha(root/'geometry_inputs.json'), 'Point inputs changed')
        points[slot] = point
        family_items[family].append((slot,item))
    require([len(family_items[key]) for key in ['core_ox','core_red','full_ox','full_red']] == [3,4,1,1],
            'Retained geometry rule changed')
    for family, items in family_items.items():
        minima = raw[family]['meta']['thermal']['ensemble']['validated_minima']
        if family.startswith('core'):
            require({item['minimum_index'] for _,item in items} == set(range(len(minima))),
                    'Core pool omits a retained minimum')
        else:
            expected = min(range(len(minima)), key=lambda i:(minima[i]['total_kj'],
                                                           minima[i]['geometry']['sha256']))
            require(items[0][1]['minimum_index'] == expected,
                    'Full-cofactor representative is not the predeclared minimum')
        first = items[0][1]['geometry']
        counts = Counter(first['symbols'])
        for _,item in items:
            require(Counter(item['geometry']['symbols']) == counts and item['geometry']['charge'] == first['charge'],
                    'Within-family formula/charge changed')
        for element,count in counts.items():residual[element] += WEIGHTS[family]*count
        charge += WEIGHTS[family]*first['charge']
    require(not any(residual.values()) and charge == 0 and sum(WEIGHTS.values()) == 0,
            'Closed exchange is not atom/charge/molecule balanced')
    arms = {}
    for electronic in ['UMA','DFT']:
        for solvent in ['ALPB','SMD']:
            families = {}
            for family,items in family_items.items():
                rows = []
                for slot,item in items:
                    original = item['original_components_kj']; point = points[slot]
                    e = original['electronic_kj'] if electronic == 'UMA' else point['gas']['energy_Eh']*protocol['Hartree_to_kj_per_mol']
                    s = original['solvation_kj'] if solvent == 'ALPB' else point['SMD_delta_Eh']*protocol['Hartree_to_kj_per_mol']
                    rows.append(dict(electronic_kj=e,solvation_kj=s,thermal_kj=original['thermal_kj']))
                families[family] = pool(rows, protocol['pool_KT_kj'])
            contrast = math.fsum(WEIGHTS[f]*families[f]['G_kj'] for f in WEIGHTS)
            components = {component:math.fsum(WEIGHTS[f]*families[f]['components_kj'][component] for f in WEIGHTS)
                          for component in families['full_red']['components_kj']}
            close(math.fsum(components.values()),contrast,1e-7,'Exchange component sum')
            arms[electronic+'_'+solvent] = dict(exchange_kj=contrast,components_kj=components,families=families)
    original_pools={}
    for family,value in raw.items():
        rows=[{k:m[k] for k in ['electronic_kj','solvation_kj','thermal_kj']}
              for m in value['meta']['thermal']['ensemble']['validated_minima']]
        original_pools[family]=pool(rows,protocol['pool_KT_kj'])
        close(original_pools[family]['G_kj'],value['G'],1e-7,'Original pool reconstruction')
    full_exchange=math.fsum(WEIGHTS[f]*original_pools[f]['G_kj'] for f in WEIGHTS)
    factorial = dict(
        DFT_minus_UMA_at_ALPB_kj=arms['DFT_ALPB']['exchange_kj']-arms['UMA_ALPB']['exchange_kj'],
        SMD_minus_ALPB_at_UMA_kj=arms['UMA_SMD']['exchange_kj']-arms['UMA_ALPB']['exchange_kj'],
        DFT_minus_UMA_at_SMD_kj=arms['DFT_SMD']['exchange_kj']-arms['UMA_SMD']['exchange_kj'],
        SMD_minus_ALPB_at_DFT_kj=arms['DFT_SMD']['exchange_kj']-arms['DFT_ALPB']['exchange_kj'])
    factorial['interaction_kj']=factorial['DFT_minus_UMA_at_SMD_kj']-factorial['DFT_minus_UMA_at_ALPB_kj']
    return dict(schema='metag-nicotinamide-fixed-geometry-audit-v1',passed=True,
        scope='Fixed retained geometries/RRHO; method sensitivity, not experimental accuracy or a production replacement.',
        n_slots=len(points),weights=WEIGHTS,element_residual={},charge_residual=charge,molecule_residual=sum(WEIGHTS.values()),
        CPU_GPU_difference_kj=smoke_differences,max_UMA_replay_difference_kj=max(abs(r['difference_kj']) for r in replay['rows'].values()),
        arms=arms,factorial_contrasts=factorial,original_all_minima_exchange_kj=full_exchange,
        restricted_geometry_effect_kj=arms['UMA_ALPB']['exchange_kj']-full_exchange,
        independent_experimental_labels_admitted=0,empirical_offsets_fitted=0,
        limitations=['One selected retained minimum per full cofactor; original41/55 ensembles preserved.',
                     'SMD is an alternative continuum model, not measured solvation or an exact reference.',
                     'Original RRHO terms retained at fixed geometries; no new stationary nuclei or Hessians.',
                     'Numerical fixture does not bound large-molecule density-fitting, grid, functional, sampling or state errors.',
                     'No orbital-stability proof or total-error interval.'])

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('root',type=Path)
    parser.add_argument('--raw-directory',type=Path)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    result=audit(args.root,args.raw_directory)
    args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:result[k] for k in ['passed','n_slots','factorial_contrasts',
                                         'original_all_minima_exchange_kj','restricted_geometry_effect_kj']},indent=2))

if __name__=='__main__':main()
