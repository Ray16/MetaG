"""Independently audit ModelSEED preflight accounting and literal graph balances.

This verifier imports no MetaG planner or thermochemistry implementation. It
uses the recorded molecular admission decisions as inputs, independently checks
their basic graph/reference bookkeeping, and reconstructs every reaction gate.
Neither replay agreement nor graph balance establishes energetic accuracy.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
from fractions import Fraction
from functools import lru_cache
import hashlib
import json
import math
import os
from pathlib import Path
import re

from rdkit import Chem, RDLogger


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read_json(path):
    return json.loads(Path(path).read_text())


def records(path, key):
    result = {}
    with Path(path).open() as stream:
        for line in stream:
            row = json.loads(line)
            if row[key] in result:
                raise ValueError('Duplicate retained identity: '+row[key])
            result[row[key]] = row
    return result


def source_tables(root, kind):
    rows, files = {}, {}
    for path in sorted((root/'Biochemistry').glob(kind+'_*.tsv')):
        if not re.fullmatch(kind+r'_\d+\.tsv', path.name):
            continue
        with path.open(newline='') as stream:
            for row in csv.DictReader(stream, delimiter='\t'):
                if row['id'] in rows:
                    raise ValueError('Duplicate source identity: '+row['id'])
                rows[row['id']] = row
        files[str(path.relative_to(root))] = sha(path)
    if not files:
        raise ValueError('No source chunks found')
    return rows, files


def terms(text):
    """Independently aggregate exact coefficients, preserving compartments."""
    values = defaultdict(Fraction)
    for token in text.split(';'):
        if not token.strip():
            continue
        fields = token.split(':')
        if len(fields) < 3 or not fields[1] or not fields[2]:
            raise ValueError('Incomplete stoichiometric identity')
        values[(fields[1], fields[2])] += Fraction(fields[0])
    return {k: v for k, v in values.items() if v}


@lru_cache(maxsize=None)
def graph(smiles):
    if not isinstance(smiles, str) or not smiles.strip():
        return None
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None or not molecule.GetNumAtoms() or any(a.GetAtomicNum() == 0 for a in molecule.GetAtoms()):
        return None
    explicit = Chem.AddHs(molecule)
    formula = Counter((str(a.GetIsotope()) if a.GetIsotope() else '')+a.GetSymbol()
                      for a in explicit.GetAtoms())
    return dict(formula=dict(formula), charge=Chem.GetFormalCharge(molecule),
                canonical=Chem.MolToSmiles(molecule, isomericSmiles=True),
                n_atoms=explicit.GetNumAtoms(), n_heavy=molecule.GetNumHeavyAtoms(),
                fragments=len(Chem.GetMolFrags(molecule)),
                has_isotopes=any(a.GetIsotope() for a in molecule.GetAtoms()),
                has_atom_maps=any(a.GetAtomMapNum() for a in molecule.GetAtoms()))


def declared_charge(source):
    try:
        value = Fraction(source['charge'])
    except (ValueError, TypeError, KeyError, ZeroDivisionError):
        return None
    return int(value) if value.denominator == 1 else None


def balance(column, compounds):
    elements, formal, declared = defaultdict(Fraction), Fraction(), Fraction()
    known, charge_known, signature = True, True, defaultdict(Fraction)
    for (cid, compartment), coefficient in column.items():
        source = compounds.get(cid)
        if source is None:
            known = charge_known = False
            continue
        q = declared_charge(source)
        if q is None:
            charge_known = False
        else:
            declared += coefficient*q
        state = graph(source.get('smiles') or '')
        if state is None:
            known = False
            continue
        for element, count in state['formula'].items():
            elements[element] += coefficient*count
        formal += coefficient*state['charge']
        signature['explicit_proton' if cid == 'cpd00067' else (state['charge'], state['canonical'])] += coefficient
    return dict(all_formulas_known=known,
                element_residual={k: str(v) for k, v in sorted(elements.items()) if v} if known else None,
                formal_charge_residual=str(formal) if known else None,
                declared_charge_residual=str(declared) if charge_known else None,
                graph_net_empty=not any(signature.values()) if known else None)


def reference_bookkeeping(source, recorded):
    q = declared_charge(source)
    original, reference = graph(source['smiles']), graph(recorded['reference_smiles'])
    if original is None or reference is None:
        raise ValueError('Admitted state has no concrete graph')
    assert original['charge'] == q == recorded['charge']
    assert recorded['canonical_input'] == original['canonical']
    assert original['fragments'] == reference['fragments'] == 1
    assert not original['has_isotopes'] and not original['has_atom_maps']
    assert recorded['n_atoms'] == original['n_atoms'] and recorded['n_heavy_atoms'] == original['n_heavy']
    exchange = Fraction(str(recorded['proton_exchange']))
    assert exchange.denominator == 1
    for element in original['formula'].keys() | reference['formula'].keys():
        expected = reference['formula'].get(element, 0)+(exchange if element == 'H' else 0)
        assert expected == original['formula'].get(element, 0), ('reference element conversion', source['id'], element)
    assert reference['charge'] == recorded['reference_charge']
    assert reference['charge']+exchange == original['charge']


def reservoir_completion(checked):
    """Unique literal balance under a declared water/proton basis; no identity inference."""
    if not checked['all_formulas_known']:
        return None
    residual = {k: Fraction(v) for k, v in checked['element_residual'].items()}
    if any(v for k, v in residual.items() if k not in {'H', 'O'}):
        return None
    h, oxygen = residual.get('H', Fraction()), residual.get('O', Fraction())
    q = Fraction(checked['formal_charge_residual'])
    if h-2*oxygen-q:
        return None
    if not h and not oxygen and not q:
        return None
    return dict(delta_water=str(-oxygen), delta_proton=str(-q),
                kind='proton_only' if not oxygen else 'water_and_proton_basis')


def audit(modelseed, retained):
    compounds, compound_files = source_tables(modelseed, 'compound')
    reactions, reaction_files = source_tables(modelseed, 'reaction')
    protocol, summary = read_json(retained/'protocol.json'), read_json(retained/'summary.json')
    assert dict(compound_files, **reaction_files) == protocol['source_files_sha256']
    for name, digest in summary['output_sha256'].items():
        assert sha(retained/name) == digest
    assert sha(retained/'protocol.json') == summary['protocol_sha256']
    states = records(retained/'compounds.jsonl', 'compound')
    outcomes = records(retained/'reactions.jsonl', 'reaction')
    assert set(states) == set(compounds)
    assert set(r['is_obsolete'] for r in reactions.values()) <= {'0', '1'}
    assert set(r['is_transport'] for r in reactions.values()) <= {'0', '1'}
    counts, state_warnings = Counter(), Counter()
    for cid, source in compounds.items():
        row = states[cid]
        assert row['name'] == source['name'] and row['source_charge'] == source['charge']
        assert row['smiles'] == source['smiles']
        if source['is_obsolete'] == '1':
            assert row['status'] == 'obsolete_compound'
        elif source.get('abstract_compound') == '1':
            assert row['status'] == 'abstract_compound'
        elif (source.get('smiles') or '').lower() in {'', 'null', 'none', 'nan'}:
            assert row['status'] == 'missing_structure'
        elif declared_charge(source) is None:
            assert row['status'] == 'missing_charge'
        elif cid == 'cpd00067':
            assert row['status'] == 'explicit_proton_reservoir'
        counts[row['status']] += 1
        if row['status'] == 'passes_cpu_state_gates':
            reference_bookkeeping(source, row)
            state_warnings.update(set(row['protonation_issues']))
    assert dict(counts) == summary['compound_statuses']
    requested = {rid for rid, source in reactions.items()
                 if source['is_obsolete'] == '0' and source['is_transport'] == '0'}
    assert set(outcomes) == requested
    excluded = Counter('obsolete' if source['is_obsolete'] == '1' else 'transport'
                       for source in reactions.values() if source['id'] not in requested)
    assert dict(excluded) == summary['excluded_from_universe']
    status_counts, first_counts, all_categories, sole_categories = Counter(), Counter(), Counter(), Counter()
    pass_flags, flags, warning_reactions = Counter(), Counter(), Counter()
    known_balance, failure_category_counts, completion_counts = Counter(), Counter(), Counter()
    detailed = []
    for rid in sorted(requested):
        source, row = reactions[rid], outcomes[rid]
        assert source['status'] == row['source_status']
        flags[source['status']] += 1
        expected, complete = [], True
        if source.get('abstract_reaction') == '1':
            expected.append(('abstract_reaction', None))
        try:
            column = terms(source['stoichiometry'])
        except (ValueError, ZeroDivisionError):
            column = {}; expected.append(('unsupported_stoichiometry', None))
        if len({compartment for cid, compartment in column}) > 1:
            expected.append(('multiple_compartments_need_conditions', None))
        assert row['compounds'] == sorted({cid for cid, compartment in column})
        for (cid, compartment), coefficient in column.items():
            if cid == 'cpd00067':
                assert compounds[cid]['smiles'] == '[H+]' and declared_charge(compounds[cid]) == 1
                continue
            if cid not in states:
                expected.append(('missing_compound', cid)); complete = False
            elif states[cid]['status'] != 'passes_cpu_state_gates':
                expected.append((states[cid]['category'], cid)); complete = False
        checked = balance(column, compounds)
        if complete and column:
            assert checked['all_formulas_known']
            if checked['element_residual']:
                expected.append(('element_imbalance', None))
            if Fraction(checked['formal_charge_residual']):
                expected.append(('charge_imbalance', None))
            if checked['graph_net_empty']:
                expected.append(('identity_or_graph_collapsed_reaction', None))
        if not column and not expected:
            expected.append(('empty_net_equation', None))
        actual = [(c['category'], c.get('compound')) for c in row['causes']]
        assert expected == actual, (rid, expected, actual)
        passed = not expected
        assert (row['status'] == 'passes_cpu_input_gates') == passed
        status_counts[row['status']] += 1
        categories = sorted({category for category, cid in expected})
        if expected:
            first_counts[expected[0][0]] += 1
            all_categories.update(categories)
            failure_category_counts[len(categories)] += 1
            if len(categories) == 1:
                sole_categories[categories[0]] += 1
        else:
            pass_flags[source['status']] += 1
        warning_set = sorted({issue for cid in row['compounds'] for issue in states[cid].get('protonation_issues', [])})
        if passed:
            warning_reactions.update(warning_set)
            known_balance['passing_with_protonation_warnings' if warning_set else 'passing_without_recorded_protonation_warnings'] += 1
            assert checked['all_formulas_known'] and not checked['element_residual']
            assert Fraction(checked['formal_charge_residual']) == Fraction(checked['declared_charge_residual']) == 0
        if column:
            if checked['all_formulas_known']:
                imbalance = bool(checked['element_residual']) or Fraction(checked['formal_charge_residual']) != 0
                known_balance['literal_formula_known'] += 1
                known_balance['literal_unbalanced' if imbalance else 'literal_balanced'] += 1
                if not complete and imbalance:
                    known_balance['latent_imbalance_with_unavailable_molecular_state'] += 1
            else:
                known_balance['literal_formula_unknown'] += 1
        completion = reservoir_completion(checked) if column else None
        if completion is not None:
            completion_counts[completion['kind']] += 1
            if complete and set(categories) <= {'element_imbalance', 'charge_imbalance'}:
                completion_counts['otherwise_passes_recorded_state_gates'] += 1
                completion['otherwise_passes_recorded_state_gates'] = True
            else:
                completion['otherwise_passes_recorded_state_gates'] = False
        detailed.append(dict(reaction=rid, status=row['status'], recorded_failure_categories=categories,
                             passing_protonation_issues=warning_set if passed else None,
                             reservoir_completion_candidate=completion,
                             **checked))
    assert dict(status_counts) == summary['reaction_input_statuses']
    assert dict(first_counts) == summary['primary_unavailable_categories']
    assert dict(flags) == summary['all_source_statuses'] and dict(pass_flags) == summary['passed_source_statuses']
    assert len(compounds) == summary['n_compound_rows'] and len(reactions) == summary['n_reaction_rows']
    assert len(requested) == summary['n_requested_live_nontransport']
    report = dict(schema='metag-modelseed-accounting-verification-v1', status='passed',
        scope=__doc__, n_source_files=len(protocol['source_files_sha256']),
        n_compounds=len(compounds), n_source_reactions=len(reactions), n_requested=len(requested),
        excluded=dict(excluded), input_statuses=dict(status_counts), first_recorded_failure=dict(first_counts),
        all_recorded_failure_categories=dict(all_categories), sole_recorded_failure_category=dict(sole_categories),
        n_failure_categories_per_unavailable_reaction={str(k): v for k, v in failure_category_counts.items()},
        structural_balance_and_warnings=dict(known_balance), admitted_compound_warning_counts=dict(state_warnings),
        passing_reaction_warning_counts=dict(warning_reactions),
        reservoir_completion_candidates=dict(completion_counts),
        molecular_admission_rules_independently_reimplemented=False,
        literal_graph_and_reference_charge_bookkeeping_independently_checked=True,
        source_sha256=sha(__file__), retained_protocol_sha256=sha(retained/'protocol.json'),
        retained_output_sha256=summary['output_sha256'],
        interpretation='Sole recorded category is not an estimate of restored coverage. Missing structures/states can hide later gates. No energetic inference, prediction coverage or accuracy is established.')
    return report, detailed


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--modelseed', type=Path, required=True)
    parser.add_argument('--retained', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Output directory already exists')
    if hasattr(os, 'sched_getaffinity'):
        os.sched_setaffinity(0, {max(os.sched_getaffinity(0))})
    RDLogger.DisableLog('rdApp.*')
    report, rows = audit(args.modelseed, args.retained)
    args.output.mkdir()
    with (args.output/'reaction_accounting.jsonl').open('x') as stream:
        for row in rows:
            stream.write(json.dumps(row, allow_nan=False)+'\n')
    report['reaction_accounting_sha256'] = sha(args.output/'reaction_accounting.jsonl')
    (args.output/'verification.json').write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    print(json.dumps({k: report[k] for k in ['status', 'n_compounds', 'n_source_reactions', 'n_requested',
                                           'input_statuses', 'structural_balance_and_warnings']}, indent=2))


if __name__ == '__main__':
    main()
