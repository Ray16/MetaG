"""Explore paired errors by structural overlap with dGbyG training compounds.

This is a retrospective, development-exposed analysis; novelty is not proof of
independent experimental data. It uses released-model predictions, not OOF CV.
"""
import argparse
import csv
import hashlib
import json
import math
import re
from collections import defaultdict
from functools import lru_cache
from pathlib import Path

from rdkit import Chem, DataStructs, rdBase
from rdkit.Chem import rdFingerprintGenerator
from rdkit.Chem.MolStandardize import rdMolStandardize

from score_checked import load_predictions
from score import metrics

HERE = Path(__file__).resolve().parent
UNCHARGER = rdMolStandardize.Uncharger()
FP = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
FP_CHIRAL = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048, includeChirality=True)
TERM = re.compile(r'^(\d+(?:\.\d+)?)\s+(.+)$')


@lru_cache(maxsize=None)
def structure(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None or any(atom.GetAtomicNum() == 0 for atom in mol.GetAtoms()):
        raise ValueError(f'unsupported structure: {smiles}')
    mol = UNCHARGER.uncharge(mol)
    return (Chem.MolToSmiles(mol, isomericSmiles=True),
            Chem.MolToSmiles(mol, isomericSmiles=False), mol)


def reaction_signature(species, stereo=True):
    counts = defaultdict(float)
    for coefficient, smiles in species:
        iso, plain, _ = structure(smiles)
        if iso in ('O', '[H+]', '[H]'):
            continue
        counts[iso if stereo else plain] += coefficient
    values = tuple(sorted((s, round(c, 6)) for s, c in counts.items() if abs(c) > 1e-9))
    if not values:
        raise ValueError('empty reaction signature')
    return min(values, tuple((s, -c) for s, c in values))


@lru_cache(maxsize=None)
def count_fingerprint(smiles, chiral=False):
    generator = FP_CHIRAL if chiral else FP
    return generator.GetSparseCountFingerprint(structure(smiles)[2]).GetNonzeroElements()


def reaction_change(signature, chiral=False):
    diff = defaultdict(float)
    for smiles, coefficient in signature:
        for bit, count in count_fingerprint(smiles, chiral).items():
            diff[bit] += coefficient * count
    return {(bit, value > 0): abs(value) for bit, value in diff.items() if abs(value) > 1e-9}


def generalized_jaccard(a, b):
    keys = a.keys() | b.keys()
    denominator = sum(max(a.get(k, 0), b.get(k, 0)) for k in keys)
    return sum(min(a.get(k, 0), b.get(k, 0)) for k in keys) / denominator if denominator else 1.0


def nearest_change(signature, training_changes, chiral=False):
    change = reaction_change(signature, chiral)
    if not change or not training_changes:
        return None
    reverse = {(bit, not sign): value for (bit, sign), value in change.items()}
    return max(max(generalized_jaccard(change, t), generalized_jaccard(reverse, t)) for t in training_changes)


def parse_training_equation(equation):
    if ' = ' not in equation:
        return [(1.0, equation.strip())]
    lhs, rhs = equation.split(' = ')
    result = []
    for side, sign in ((lhs, -1), (rhs, 1)):
        for term in side.split(' + '):
            match = TERM.match(term.strip())
            result.append((sign * (float(match.group(1)) if match else 1.0),
                           match.group(2) if match else term.strip()))
    return result


def training_domain(rows):
    compounds, reactions, reactions_no_stereo, skipped = {}, set(), set(), []
    for i, row in enumerate(rows):
        try:
            parts = parse_training_equation(row['reaction'])
            parsed = [structure(s) for _, s in parts]
            for iso, _, mol in parsed:
                if iso not in ('O', '[H+]', '[H]'):
                    compounds[iso] = FP.GetFingerprint(mol)
            if ' = ' in row['reaction']:
                reactions.add(reaction_signature(parts))
                reactions_no_stereo.add(reaction_signature(parts, stereo=False))
        except (ValueError, KeyError) as exc:
            skipped.append({'row': i, 'reason': str(exc)})
    if not compounds:
        raise ValueError('no training compounds parsed')
    changes = [change for sig in reactions if (change := reaction_change(sig))]
    chiral_changes = [change for sig in reactions if (change := reaction_change(sig, chiral=True))]
    return compounds, reactions, reactions_no_stereo, skipped, changes, chiral_changes


def classify(reaction, domain):
    compounds, rx, rx_no_stereo, _, changes, chiral_changes = domain
    parts = [(float(v[0]), v[2]) for v in reaction['species'].values()]
    try:
        molecules = [structure(s) for _, s in parts]
        molecules = [x for x in molecules if x[0] not in ('O', '[H+]', '[H]')]
        if not molecules:
            raise ValueError('reaction has no non-water compounds')
        unseen = sum(iso not in compounds for iso, _, _ in molecules)
        fps = list(compounds.values())
        nearest = min(max(DataStructs.BulkTanimotoSimilarity(FP.GetFingerprint(mol), fps))
                      for _, _, mol in molecules)
        signature = reaction_signature(parts)
        exact = signature in rx
        try:
            stereofree = reaction_signature(parts, stereo=False) in rx_no_stereo
        except ValueError as exc:
            if str(exc) != 'empty reaction signature':
                raise
            stereofree = None  # Stereo removal erased the transformation; overlap is indeterminate.
        tier = ('reaction_seen' if exact or stereofree else 'stereochemistry_collapsed' if stereofree is None
                else 'all_compounds_seen' if unseen == 0
                else 'unseen_compound_sim_ge_0.5' if nearest >= 0.5 else 'unseen_compound_sim_lt_0.5')
        return {'tier': tier, 'n_unseen_compounds': unseen, 'minimum_nearest_similarity': nearest,
                'reaction_exact': exact, 'reaction_stereofree': stereofree,
                'nearest_reaction_change_similarity': nearest_change(signature, changes),
                'nearest_chiral_change_similarity': nearest_change(signature, chiral_changes, chiral=True)}
    except (ValueError, KeyError, TypeError) as exc:
        return {'tier': 'unclassifiable', 'reason': str(exc)}


def read_comparator(path):
    out = {}
    with path.open(newline='') as stream:
        for row in csv.DictReader(stream, delimiter='\t'):
            if row['method'] != 'dGbyG_eQ_bridge' or row['status'] != 'ok':
                continue
            if row['reaction'] in out:
                raise ValueError(f"duplicate dGbyG prediction {row['reaction']}")
            conditions = tuple(float(row[k]) for k in ('T_K', 'pH', 'I_M', 'pMg'))
            if conditions != (298.15, 7.0, 0.0, 14.0):
                raise ValueError(f"unexpected conditions for {row['reaction']}: {conditions}")
            value = float(row['dG_kJ'])
            if not math.isfinite(value):
                raise ValueError('nonfinite dGbyG prediction')
            out[row['reaction']] = value
    return out


def summarize(group):
    return {'MetaG': metrics([(p['MetaG'], p['reference']) for p in group]),
            'dGbyG': metrics([(p['dGbyG'], p['reference']) for p in group]),
            'MetaG_wins': sum(abs(p['MetaG'] - p['reference']) < abs(p['dGbyG'] - p['reference']) for p in group)}


def has_named_nad(reaction):
    """Exploratory name-based flag; not a comprehensive cofactor classifier."""
    return bool(set(reaction['species']) & {'NAD', 'NADH', 'NADP', 'NADPH', 'NAD+', 'NADP+'})


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--training', type=Path, required=True)
    parser.add_argument('--comparator', type=Path, required=True)
    parser.add_argument('--reactions', type=Path, default=HERE.parent / 'src/metag/data/reactions_opentecr_std.json')
    parser.add_argument('--reference', type=Path, default=HERE / 'reference_mg_bounded.json')
    parser.add_argument('--predictions', type=Path, default=HERE / 'metag_corrected_source_predictions.json')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args(argv)
    with args.training.open(newline='') as stream:
        domain = training_domain(list(csv.DictReader(stream)))
    reactions = json.loads(args.reactions.read_text())
    labels = json.loads(args.reference.read_text())['2.0']['labels']
    metag, comparator = load_predictions(args.predictions), read_comparator(args.comparator)
    paired, groups = [], defaultdict(list)
    for rid, label in labels.items():
        if metag.get(rid) is None or rid not in comparator or rid not in reactions:
            continue
        info = classify(reactions[rid], domain)
        row = {'reaction': rid, 'reference': label['exp'], 'MetaG': metag[rid],
               'dGbyG': comparator[rid], **info}
        paired.append(row)
        groups[info['tier']].append(row)
    summaries = {tier: summarize(group) for tier, group in sorted(groups.items())}
    novel = [p for p in paired if p['tier'] not in ('reaction_seen', 'unclassifiable', 'stereochemistry_collapsed')
             and p['nearest_reaction_change_similarity'] is not None]
    change_bins = {f'{low}-{high}': summarize([p for p in novel
                   if low <= p['nearest_reaction_change_similarity'] < high])
                   for low, high in ((0, 0.6), (0.6, 0.8), (0.8, 1.000001))}
    nad_groups = {
        'contains_named_NAD': summarize([p for p in paired if has_named_nad(reactions[p['reaction']])]),
        'no_named_NAD': summarize([p for p in paired if not has_named_nad(reactions[p['reaction']])]),
        'unseen_compound_ge_0.5_and_named_NAD': summarize([p for p in paired
            if p['tier'] == 'unseen_compound_sim_ge_0.5' and has_named_nad(reactions[p['reaction']])]),
    }
    paths = ('training', 'comparator', 'reactions', 'reference', 'predictions')
    output = {'scope': __doc__, 'training_compounds': len(domain[0]), 'training_reaction_signatures': len(domain[1]),
              'training_rows_skipped': domain[3], 'rdkit_version': rdBase.rdkitVersion,
              'input_sha256': {k: hashlib.sha256(getattr(args, k).read_bytes()).hexdigest() for k in paths},
              'summaries': summaries, 'exploratory_novel_reaction_change_bins': change_bins,
              'exploratory_name_based_nad_groups': nad_groups, 'pairs': paired}
    args.output.write_text(json.dumps(output, indent=2, allow_nan=False) + '\n')
    print(json.dumps({k: v for k, v in output.items() if k != 'pairs'}, indent=2))


if __name__ == '__main__':
    main()
