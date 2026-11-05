"""Replay the complete fixed-geometry cofactor study without molecular backends.

This checks retained numerical evidence and partition arithmetic. It does not
rerun SCF, validate experimental accuracy, or adopt a production reference policy.
"""
from collections import Counter
from decimal import Decimal, localcontext
from pathlib import Path
import argparse
import ast
import hashlib
import json
import math
import zipfile

WEIGHTS = dict(full_red=1, full_ox=-1, core_red=-1, core_ox=1)


def audit(root):
    root = Path(root)
    read = lambda name: json.loads((root / name).read_text())
    sha = lambda value: hashlib.sha256(value).hexdigest()
    manifest = read('manifest.json')
    for name, digest in manifest['files'].items():
        path = Path(name)
        assert not path.is_absolute() and '..' not in path.parts
        assert not (root/path).is_symlink()
        assert sha((root/path).read_bytes()) == digest, name
    with zipfile.ZipFile(root/'calculation_records.zip') as archive:
        names = archive.namelist()
        assert len(names) == len(set(names))
        assert all(not Path(n).is_absolute() and '..' not in Path(n).parts for n in names)
        raw = {name: archive.read(name) for name in names}
    obj = lambda name: json.loads(raw[name])
    frozen = obj('study_manifest.json')
    for name, digest in frozen['files_sha256'].items():
        assert sha(raw[name]) == digest, name
    analysis_freeze = obj('analysis_freeze.json')
    assert sha(raw['analyze.py']) == analysis_freeze['analyze_sha256']
    assert sha(raw['fixed_geometry_protocol.json']) == analysis_freeze['protocol_sha256']
    inputs = obj('geometry_inputs.json')
    protocol = obj('fixed_geometry_protocol.json')
    numerical_keys = next(ast.literal_eval(node.value) for node in ast.parse(raw['control_common.py']).body
                          if isinstance(node,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='NUMERICAL_KEYS' for t in node.targets))
    previous_protocol = obj('reused_protocol.json')
    assert all(protocol[key] == previous_protocol[key] for key in numerical_keys)
    reuse = obj('reuse_receipt.json')
    expected = read('ensemble_results.json')
    families = {f: sorted(k for k,v in inputs['inputs'].items() if v['family']==f) for f in WEIGHTS}
    assert [len(families[f]) for f in ['core_ox','core_red','full_ox','full_red']] == [3,4,41,55]
    assert len(inputs['inputs']) == 103 and len(reuse) == 9
    residual = Counter(); charge = 0
    for family, keys in families.items():
        source = inputs['source_records'][family]
        data = raw['raw_records/'+Path(source['path']).name]
        assert sha(data) == source['sha256']
        minima = json.loads(data)['meta']['thermal']['ensemble']['validated_minima']
        assert {inputs['inputs'][key]['minimum_index'] for key in keys} == set(range(len(minima)))
        first = inputs['inputs'][keys[0]]['geometry']
        for key in keys:
            item = inputs['inputs'][key]; minimum = minima[item['minimum_index']]
            assert item['geometry'] == minimum['geometry']
            assert item['reference_smiles'] == minimum['reference_smiles']
            assert all(minimum[k] == v for k,v in item['original_components_kj'].items())
            assert Counter(item['geometry']['symbols']) == Counter(first['symbols'])
            assert item['geometry']['charge'] == first['charge']
        for symbol, count in Counter(first['symbols']).items():
            residual[symbol] += count*WEIGHTS[family]
        charge += first['charge']*WEIGHTS[family]
    assert not any(residual.values()) and charge == 0 and sum(WEIGHTS.values()) == 0
    points = {}
    versions = None
    for key, item in inputs['inputs'].items():
        filename = reuse[key]['point'] if key in reuse else 'points/'+key+'.json'
        point = obj(filename)
        assert point['passed'] is True and point['slot'] == key
        assert point['geometry'] == item['geometry']
        assert point['n_AO'] == item['n_AO'] and point['n_electrons'] == item['n_electrons']
        assert point['overlap_minimum_eigenvalue'] >= protocol['overlap_minimum_eigenvalue']
        if versions is None: versions = point['versions']
        assert point['versions'] == versions
        if key in reuse:
            receipt = reuse[key]
            assert sha(raw[filename]) == receipt['sha256']
            for field, witness in [('source_sha256','original_point_source_sha256'),
                                   ('protocol_sha256','original_protocol_sha256'),
                                   ('input_sha256','original_input_sha256')]:
                assert point[field] == receipt[witness]
        else:
            for field, filename in [('source_sha256','reference_point.py'),
                                    ('protocol_sha256','fixed_geometry_protocol.json'),
                                    ('input_sha256','geometry_inputs.json')]:
                assert point[field] == sha(raw[filename])
        for phase in ['gas','water']:
            state = point[phase]; terms = state['scf_summary']
            assert state['converged'] is True and math.isfinite(state['energy_Eh'])
            assert state['n_molecular_orbitals'] == item['n_AO']
            assert 0 <= state['orbital_gradient_norm'] <= protocol['stationary_gradient_maximum']
            assert abs(terms['e2']-terms['coul']-terms['exc']) < 1e-8
            components = [terms[k] for k in ['e1','e2','nuc']]
            if phase == 'water':
                components += [terms['e_solvent'],terms['e_cds']]
                assert abs(terms['e_cds']-point['CDS_Eh']) < 1e-12
            assert abs(math.fsum(components)-state['energy_Eh']) < 1e-8
        assert abs(point['SMD_delta_Eh']-point['water']['energy_Eh']+point['gas']['energy_Eh']) < 1e-10
        assert point['surface_discretization_method'] == 'SWIG'
        points[key] = point
    replay = obj('uma_replay.json')
    assert replay['passed'] is True and set(replay['rows']) == set(points)
    for key, row in replay['rows'].items():
        original = inputs['inputs'][key]['original_components_kj']['electronic_kj']
        assert row['geometry_sha256'] == inputs['inputs'][key]['geometry']['sha256']
        assert row['original_electronic_kj'] == original
        assert abs(row['replay_electronic_kj']-original-row['difference_kj']) < 1e-9
        assert abs(row['difference_kj']) <= protocol['UMA_replay_max_difference_kj']
    exchanges = {}
    with localcontext() as context:
        context.prec = 50
        D = lambda value: Decimal(str(value))
        kt = D(protocol['pool_KT_kj']); hartree = D(protocol['Hartree_to_kj_per_mol'])
        for subset in ['full','restricted']:
            for electronic in ['UMA','DFT']:
                for solvent in ['ALPB','SMD']:
                    arm = f'{subset}|{electronic}|{solvent}'; pools = {}
                    for family, keys in families.items():
                        selected = keys if subset == 'full' else [key for key in keys if key in reuse]
                        rows = []; components = []
                        for key in selected:
                            original = inputs['inputs'][key]['original_components_kj']; point = points[key]
                            # Construct DFT/SMD directly from the water total. This
                            # avoids the producer analysis' gas-plus-difference path.
                            e = D(original['electronic_kj']) if electronic == 'UMA' else D(point['gas']['energy_Eh'])*hartree
                            s = D(original['solvation_kj']) if solvent == 'ALPB' else (D(point['water']['energy_Eh'])-D(point['gas']['energy_Eh']))*hartree
                            terms = dict(electronic_kj=e,solvation_kj=s,thermal_kj=D(original['thermal_kj']))
                            g = sum(terms.values())
                            if electronic == 'DFT' and solvent == 'SMD':
                                g = D(point['water']['energy_Eh'])*hartree+D(original['thermal_kj'])
                            rows.append(g)
                            components.append(terms)
                        minimum = min(rows); shifted = [(-(g-minimum)/kt).exp() for g in rows]
                        partition = sum(shifted); pools[family] = minimum-kt*partition.ln()
                        saved = expected['arms'][arm]['families'][family]
                        assert selected == saved['slots']
                        assert abs(float(pools[family])-saved['G_kj']) < 1e-7
                        assert max(abs(float(p/partition)-q) for p,q in zip(shifted,saved['probabilities'])) < 1e-7
                        probabilities = [p/partition for p in shifted]
                        reconstructed = {key:sum(p*row[key] for p,row in zip(probabilities,components))
                                         for key in components[0]}
                        reconstructed['configurational_kj'] = kt*sum(p*p.ln() for p in probabilities if p)
                        assert all(abs(float(v)-saved['components_kj'][key]) < 1e-7 for key,v in reconstructed.items())
                        assert abs(math.fsum(saved['components_kj'].values())-saved['G_kj']) < 1e-7
                    exchange = sum(WEIGHTS[f]*pools[f] for f in WEIGHTS)
                    assert abs(float(exchange)-expected['arms'][arm]['exchange_kj']) < 1e-7
                    exchanges[arm] = float(exchange)
    contrasts = {f'{e}|{s}':exchanges[f'full|{e}|{s}']-exchanges[f'restricted|{e}|{s}']
                 for e in ['UMA','DFT'] for s in ['ALPB','SMD']}
    for key, value in contrasts.items():
        assert abs(value-expected['full_minus_restricted_kj'][key]) < 1e-7
    for subset in ['full','restricted']:
        arm = lambda e,s:exchanges[f'{subset}|{e}|{s}']
        effects = dict(electronic_at_ALPB_kj=arm('DFT','ALPB')-arm('UMA','ALPB'),
                       solvent_at_UMA_kj=arm('UMA','SMD')-arm('UMA','ALPB'),
                       electronic_at_SMD_kj=arm('DFT','SMD')-arm('UMA','SMD'),
                       solvent_at_DFT_kj=arm('DFT','SMD')-arm('DFT','ALPB'),
                       interaction_kj=arm('DFT','SMD')-arm('UMA','SMD')-arm('DFT','ALPB')+arm('UMA','ALPB'))
        assert all(abs(value-expected['factorial_contrasts'][subset][key]) < 1e-7 for key,value in effects.items())
    assert abs(contrasts['DFT|SMD']-expected['primary_DFT_SMD_ensemble_contrast_kj']) < 1e-7
    assert expected['new_primary_experimental_references_admitted'] == 0
    assert expected['production_reference_policy_adopted'] is False
    assert expected['exceeds_descriptive_sensitivity_threshold'] == (abs(contrasts['DFT|SMD']) >= protocol['descriptive_sensitivity_threshold_kj'])
    return dict(status='passed',verified_pairs=len(points),reused_pairs=len(reuse),
                independent_decimal_exchange_kj=exchanges,full_minus_restricted_kj=contrasts,
                new_primary_experimental_references_admitted=0,production_reference_policy_adopted=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('root', nargs='?', type=Path,
                        default=Path(__file__).resolve().parent/'nicotinamide_ensemble_validation')
    args = parser.parse_args()
    print(json.dumps(audit(args.root),indent=2,allow_nan=False))
