"""Authenticate individual quantum points, not just a cancelling contrast."""
import copy
import importlib.util
import json
from pathlib import Path
import shutil

import pytest

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT/'benchmark/aqueous_geometry_validation'
SPEC = importlib.util.spec_from_file_location('geometry_parent_audit', ROOT/'benchmark/audit_aqueous_geometry_parent.py')
AUDIT = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(AUDIT)


def point_inputs(slot='oxidized_gas'):
    item = AUDIT.read(DATA/'geometry_inputs.json')['inputs'][slot]
    point = AUDIT.read(DATA/'points'/(slot+'.json'))
    log = (DATA/'points'/(slot+'.scf.log')).read_text()
    return item, point, log, AUDIT.read(DATA/'fixed_geometry_protocol.json')


def test_complete_control_reproduces_retained_physics_and_scope():
    report = AUDIT.audit(DATA)
    assert report == AUDIT.read(DATA/'audit_summary.json')
    assert len(report['rows']) == 4
    assert report['experimental_labels_admitted'] == 0
    assert not report['production_promotion']
    assert max(report['water_fixture_GPU_CPU_difference_kj'].values()) < 0.05


def test_equal_offsets_can_cancel_the_contrast_but_cannot_fake_raw_points():
    first = point_inputs('oxidized_gas')
    second = point_inputs('reduced_gas')
    weights = [1, -1]
    assert sum(w*0.01 for w in weights) == 0
    for item, point, log, numerical in [first, second]:
        point['gas']['energy_Eh'] += 0.01
        point['gas']['scf_summary']['e1'] += 0.01
        with pytest.raises(ValueError, match='Native SCF log differs'):
            AUDIT.validate_point(item, point, log, numerical)


def test_compensated_nuclear_error_fails_even_with_unchanged_total():
    item, point, log, numerical = point_inputs()
    point['gas']['scf_summary']['nuc'] += 0.01
    point['gas']['scf_summary']['e1'] -= 0.01
    with pytest.raises(ValueError, match='Nuclear term'):
        AUDIT.validate_point(item, point, log, numerical)


def test_basis_truncation_cannot_be_reported_as_the_full_basis():
    item, point, log, numerical = point_inputs()
    point['water']['n_molecular_orbitals'] -= 1
    with pytest.raises(ValueError, match='AO directions dropped'):
        AUDIT.validate_point(item, point, log, numerical)


def test_nonfinite_solvent_metadata_fails_instead_of_bypassing_comparisons():
    item, point, log, numerical = point_inputs()
    point['CDS_Eh'] = float('nan')
    with pytest.raises(ValueError, match='Nonfinite'):
        AUDIT.validate_point(item, point, log, numerical)


def test_wrong_ordered_geometry_is_not_accepted_as_an_owned_point():
    item, point, log, numerical = point_inputs()
    point['geometry']['positions_angstrom'][0][0] += 0.1
    with pytest.raises(ValueError, match='Point geometry'):
        AUDIT.validate_point(item, point, log, numerical)


@pytest.mark.parametrize('gradient', [0.1, -1.0])
def test_orbital_stationarity_cannot_be_fabricated(gradient):
    item, point, log, numerical = point_inputs()
    point['gas']['orbital_gradient_norm'] = gradient
    with pytest.raises(ValueError, match='orbital-gradient'):
        AUDIT.validate_point(item, point, log, numerical)


def test_nonfinite_json_is_refused(tmp_path):
    path = tmp_path/'bad.json'
    path.write_text('{"energy": NaN}')
    with pytest.raises(ValueError, match='Nonfinite JSON'):
        AUDIT.read(path)


def test_missing_endpoint_never_becomes_zero(tmp_path):
    destination = tmp_path/'control'
    shutil.copytree(DATA, destination)
    (destination/'points/reduced_aqueous.json').unlink()
    with pytest.raises(FileNotFoundError):
        AUDIT.audit(destination)


def test_mixed_package_versions_fail_even_with_unchanged_energies(tmp_path):
    destination = tmp_path/'control'
    shutil.copytree(DATA, destination)
    path = destination/'points/oxidized_gas.json'
    point = AUDIT.read(path)
    point['versions']['pyscf'] = 'unowned-version'
    path.write_text(json.dumps(point))
    with pytest.raises(ValueError, match='dependency versions differ'):
        AUDIT.audit(destination)
