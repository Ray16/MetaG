"""Outcome-free review of formation-derived stereo pairs; no reference enrollment."""
from __future__ import annotations
from collections import Counter, defaultdict
from itertools import combinations
from pathlib import Path
import argparse
import csv
import hashlib
import json
from rdkit import Chem, DataStructs, rdBase
from rdkit.Chem import rdFingerprintGenerator
from benchmark.benchmark_feasibility import ALLOWED_COLUMNS, concrete, development_panel, equation, family

sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()


def structural_rows(path):
    with Path(path).open(newline='') as stream:
        reader = csv.DictReader(stream)
        if set(reader.fieldnames or []) != set(ALLOWED_COLUMNS):
            raise ValueError('Requires the outcome-free structural projection, with exactly its declared columns')
        return list(reader)


def stereo_candidates(rows):
    groups = defaultdict(list)
    unresolved = []
    for index, row in enumerate(rows):
        if row['type'] != 'formation dg':
            continue
        try:
            mol = concrete(row['reaction'])
            exact = Chem.MolToSmiles(mol, isomericSmiles=True)
            mirror = Chem.Mol(mol)
            for atom in mirror.GetAtoms():
                atom.InvertChirality()
            Chem.AssignStereochemistry(mirror, cleanIt=True, force=True)
            flat = Chem.Mol(mol)
            Chem.RemoveStereochemistry(flat)
            groups[Chem.MolToSmiles(flat, isomericSmiles=False)].append(dict(
                row=index, name=row['description'], smiles=exact,
                mirror=Chem.MolToSmiles(mirror, isomericSmiles=True),
                unspecified_atom_stereocentres=[i for i, tag in Chem.FindMolChiralCenters(mol, includeUnassigned=True) if tag == '?'],
                condition_metadata={k: row[k] for k in ['T', 'I', 'pH', 'pMg']}))
        except (ValueError, RuntimeError) as error:
            unresolved.append(dict(row=index, reason=str(error)))
    pairs = []
    for graph, members in sorted(groups.items()):
        for left, right in combinations(members, 2):
            if right['smiles'] in {left['smiles'], left['mirror']}:
                continue
            pairs.append(dict(nonstereo_graph=graph, left=left, right=right))
    return pairs, unresolved


def family_holdout(rows, smiles):
    """Conservative structural preview; excludes formation AND reaction labels."""
    infos = [family(s) for s in smiles]
    if any(info['broad'] is None or info['pool'] is None or info['reservoir'] for info in infos):
        return dict(status='unresolved_candidate_family', launch_allowed=False)
    excluded_families = {info['broad'] for info in infos}
    retained, excluded, pools = [], [], set()
    for index, row in enumerate(rows):
        reasons, participants = [], []
        try:
            if row['type'] not in {'reaction dg', 'formation dg'} or (
                row['type'] == 'reaction dg') != (' = ' in row['reaction']):
                raise ValueError('Unresolved target type/equation')
            participants = [family(s) for _, s in equation(row['reaction'])]
            if any(p['broad'] is None or p['pool'] is None for p in participants):
                reasons.append('unresolved training family')
            if any(p['broad'] in excluded_families for p in participants):
                reasons.append('held-out molecular family')
        except (ValueError, RuntimeError) as error:
            reasons.append('unresolved training row: ' + str(error))
        if reasons:
            excluded.append(dict(row=index, target_type=row['type'], reasons=reasons))
        else:
            retained.append(index)
            pools.update(p['pool'] for p in participants if not p['reservoir'])
    if not pools:
        raise ValueError('No resolved comparison domain remains')
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048, includeChirality=True)
    keys = sorted(pools)
    fingerprints = [generator.GetFingerprint(concrete(s)) for s in keys]
    neighbours = []
    for supplied, info in zip(smiles, infos, strict=True):
        similarities = DataStructs.BulkTanimotoSimilarity(generator.GetFingerprint(concrete(info['pool'])), fingerprints)
        best = max(range(len(keys)), key=lambda i: similarities[i])
        neighbours.append(dict(supplied_smiles=supplied, pool_key=info['pool'],
            nearest_retained_pool=keys[best], nearest_similarity=similarities[best],
            below_declared_0_5_threshold=similarities[best] < 0.5))
    return dict(status='structural_preview_only', heldout_broad_families=sorted(excluded_families),
        retained_rows=retained, excluded_rows=excluded, neighbours=neighbours,
        counts=dict(retained=len(retained), excluded=len(excluded),
                    retained_by_target_type=dict(Counter(rows[i]['type'] for i in retained))),
        any_participant_below_0_5=any(n['below_declared_0_5_threshold'] for n in neighbours),
        launch_allowed=False, training_performed=False)


def audit(training, historical):
    rows = structural_rows(training)
    pairs, unresolved = stereo_candidates(rows)
    development, _molecules, history_unresolved = development_panel(historical)
    witnesses = {tier: defaultdict(list) for tier in ['exact', 'pool', 'broad']}
    historical_data = json.loads(Path(historical).read_text())
    for rid, row in historical_data.items():
        for name, item in row.get('species', {}).items():
            try:
                info = family(item[2])
                for tier in witnesses:
                    if info[tier] is not None:
                        witnesses[tier][info[tier]].append(dict(input_id=rid, species=name))
            except (ValueError, RuntimeError, IndexError, TypeError):
                continue
    for pair in pairs:
        for key in ['left', 'right']:
            member = pair[key]
            info = family(member['smiles'])
            member['state_aliases'] = info
            member['historical_structure_witnesses'] = {
                tier: witnesses[tier].get(info[tier], []) if info[tier] is not None else [] for tier in witnesses}
            member['native_formation_training_row'] = member['row']
        pair['family_holdout'] = family_holdout(rows, [pair[k]['smiles'] for k in ['left', 'right']])
        pair['released_model_training_independent'] = False
        pair['reference_status'] = 'not_enrolled; source ancestry, evaluated quantity, uncertainty and condition transfer unverified'
    counts = dict(formation_rows=sum(r['type'] == 'formation dg' for r in rows),
        stereo_distinct_pairs=len(pairs), nonstereo_graph_groups=len({p['nonstereo_graph'] for p in pairs}),
        pairs_with_both_exact_historical_structure_witnesses=sum(all(p[k]['historical_structure_witnesses']['exact'] for k in ['left', 'right']) for p in pairs),
        pairs_with_any_exact_historical_structure_witness=sum(any(p[k]['historical_structure_witnesses']['exact'] for k in ['left', 'right']) for p in pairs),
        pairs_with_any_pool_historical_structure_witness=sum(any(p[k]['historical_structure_witnesses']['pool'] for k in ['left', 'right']) for p in pairs),
        pairs_with_any_broad_historical_structure_witness=sum(any(p[k]['historical_structure_witnesses']['broad'] for k in ['left', 'right']) for p in pairs),
        pairs_with_unresolved_candidate_family=sum(p['family_holdout']['status'] != 'structural_preview_only' for p in pairs),
        pairs_with_unassigned_atom_stereocentres=sum(any(p[k]['unspecified_atom_stereocentres'] for k in ['left', 'right']) for p in pairs),
        pair_previews_below_0_5=sum(p['family_holdout'].get('any_participant_below_0_5', False) for p in pairs),
        released_training_independent_pairs=0, independent_experimental_references_admitted=0)
    return dict(schema='metag-stereo-reference-feasibility-v1',
        input_sha256=dict(training=sha(training), historical=sha(historical)),
        source_sha256=sha(__file__), family_helper_sha256=sha(Path(__file__).with_name('benchmark_feasibility.py')), rdkit_version=rdBase.rdkitVersion,
        selection='All distinct stereochemical encodings of concrete formation rows sharing a nonstereo graph, excluding exact mirror encodings. Incomplete stereochemistry is retained and flagged, not asserted to define a distinct experimental stereoisomer. No numerical targets or predictions.',
        counts=counts, pairs=pairs, unresolved_formation_rows=unresolved,
        historical_inputs=len(development), historical_unresolved_terms=history_unresolved,
        outcome_values_exported=False, predictions_performed=False, model_training_performed=False,
        experimental_references_admitted=0, production_policy_adopted=False, launch_allowed=False,
        limits='Bounded structural history is not an outcome-exposure audit. Missing witnesses do not establish independence. Holdout previews are not native released-model novelty, source provenance, physical state admission or a trained comparison. Shared compounds/sources prevent treating all pairs as independent clusters.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--training', type=Path, default=Path(__file__).with_name('feasibility_validation')/'training671_structural.csv')
    parser.add_argument('--historical', type=Path, default=Path(__file__).with_name('feasibility_validation')/'historical_inputs.json')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    result = audit(args.training, args.historical)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    print(json.dumps(result['counts'], indent=2))


if __name__ == '__main__':
    main()
