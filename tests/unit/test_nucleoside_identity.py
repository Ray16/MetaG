"""Preserve source-specific sugar identities through database normalization."""
import hashlib
import importlib.util
import json
from pathlib import Path
import shutil

import pytest
from rdkit import Chem

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT/'benchmark/nucleoside_identity_validation'
SPEC = importlib.util.spec_from_file_location('nucleoside_identity', ROOT/'benchmark/audit_nucleoside_identity.py')
AUDIT = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(AUDIT)


def molecule(code, filename=None):
    record = json.loads((DATA/'properties'/(filename or code+'.json')).read_text())['PropertyTable']['Properties'][0]
    return Chem.MolFromSmiles(record.get('SMILES', record.get('IsomericSMILES')))


def test_complete_source_panel_preserves_covalent_scope_and_no_admission():
    result = AUDIT.audit(DATA)
    assert result['n_compounds'] == 24
    assert result['n_absolute_sugar_checks'] == 24
    assert result['n_relational_checks'] == 20
    assert result['experimental_references_admitted'] == 0
    assert not result['predictions_performed'] and not result['primary_benchmark_launch_allowed']
    assert result['source_specific_resolution']['initial_generic_name_CID'] == 3363
    assert result['source_specific_resolution']['source_specific_CID'] == 5790


def test_generic_name_entry_does_not_supply_missing_stereochemistry():
    with pytest.raises(ValueError, match='unspecified'):
        AUDIT.check_sugar(molecule('5FdUrd'), True)
    AUDIT.check_sugar(molecule('5FdUrd', '5FdUrd_source_stereo.json'), True)


@pytest.mark.parametrize('role', [0, 1, 2, 3])
def test_single_sugar_epimer_is_refused(role):
    structure = molecule('Urd')
    _, indices = AUDIT.sugar_roles(structure)
    structure.GetAtomWithIdx(indices[role]).InvertChirality()
    with pytest.raises(ValueError, match='configuration mismatch'):
        AUDIT.check_sugar(structure, False)


def test_whole_sugar_enantiomer_is_refused():
    structure = molecule('Ado')
    _, indices = AUDIT.sugar_roles(structure)
    for index in indices:
        structure.GetAtomWithIdx(index).InvertChirality()
    with pytest.raises(ValueError, match='configuration mismatch'):
        AUDIT.check_sugar(structure, False)


def test_ribose_cannot_be_relabelled_as_deoxyribose():
    with pytest.raises(ValueError, match='configuration mismatch'):
        AUDIT.check_sugar(molecule('Urd'), True)
    assert AUDIT.canon(AUDIT.remove_2prime_OH(molecule('Urd'))) == AUDIT.canon(molecule('dUrd'))


def test_valid_database_identity_with_wrong_substituent_still_fails_source_check(tmp_path):
    destination = tmp_path/'source'
    shutil.copytree(DATA, destination)
    inputs = json.loads((destination/'source_compounds.json').read_text())
    row = next(r for r in inputs['compounds'] if r['code'] == '5FUrd')
    row['selected_property_file'] = 'properties/5BrUrd.json'
    row['selected_property_sha256'] = hashlib.sha256((destination/row['selected_property_file']).read_bytes()).hexdigest()
    (destination/'source_compounds.json').write_text(json.dumps(inputs))
    with pytest.raises(ValueError, match='substituent or sugar identity mismatch'):
        AUDIT.audit(destination)
