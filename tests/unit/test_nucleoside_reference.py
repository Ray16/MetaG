"""Protect source filtering, stoichiometric units and observation ownership."""
import csv
import importlib.util
import io
import math
from pathlib import Path
from xml.etree import ElementTree as ET
import zipfile

import pytest

ROOT = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location('nucleoside_audit', ROOT/'benchmark/audit_nucleoside_reference.py')
AUDIT = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(AUDIT)
SOURCE = ROOT/'benchmark/nucleoside_reference_review/source_data.zip'


def archive_with_change(tmp_path, member, transform):
    destination = tmp_path/'changed.zip'
    with zipfile.ZipFile(SOURCE) as original, zipfile.ZipFile(destination, 'w') as output:
        for info in original.infolist():
            content = original.read(info)
            output.writestr(info, transform(content) if info.filename == member else content)
    return destination


def raw_change(content, mutate):
    rows = list(csv.DictReader(io.StringIO(content.decode('utf-8-sig'))))
    mutate(rows)
    stream = io.StringIO()
    writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
    writer.writeheader(); writer.writerows(rows)
    return stream.getvalue().encode()


def test_concentration_unit_cancels_for_two_to_two_stoichiometry():
    assert AUDIT.apparent_constant(0.5, 5) == pytest.approx(1/9)
    for initial in (0.002, 2, 2000):  # M, mM, micromolar; same physical mixture.
        product = initial*0.5
        quotient = product*product / ((initial-product)*(5*initial-product))
        assert AUDIT.apparent_constant(0.5, 5) == pytest.approx(quotient)


@pytest.mark.parametrize('conversion,excess', [
    (0, 5), (1, 5), (-0.1, 5), (0.5, 0.5), (math.nan, 5), (0.5, math.inf),
])
def test_invalid_or_boundary_measurements_never_become_finite_constants(conversion, excess):
    with pytest.raises(ValueError, match='physical domain'):
        AUDIT.apparent_constant(conversion, excess)


def test_original_source_replays_all_included_and_excluded_observations():
    report = AUDIT.audit(SOURCE)
    assert report['counts']['raw_observations'] == 576
    assert report['counts']['source_compilation_included'] == 519
    assert report['counts']['source_compilation_excluded'] == 57
    assert len(report['records']) == 576
    assert report['counts']['reaction_replicates_per_condition'] == 2
    assert report['counts']['equilibrium_aliquots_per_reaction'] == 3
    assert report['maximum_reconstructed_cached_constant_difference'] < 1e-10
    assert report['source_worksheet_error_cells'][0]['cell'] == 'R67'
    assert report['identity_issues'][0]['cell'] == 'K145'
    assert not report['predictions_performed'] and report['experimental_references_admitted'] == 0
    assert not report['primary_benchmark_launch_allowed']
    assert report['source_conditions']['ionic_strength_M'] is None


def test_changed_raw_measurement_is_detected_against_compilation(tmp_path):
    def transform(content):
        def mutate(rows):
            rows[0]['Product'] = '0.5'; rows[0]['Substrate'] = '0.5'
        return raw_change(content, mutate)
    changed = archive_with_change(tmp_path, AUDIT.RAW_MEMBER, transform)
    with pytest.raises(ValueError, match='Raw-to-compiled product mismatch'):
        AUDIT.audit(changed)


def test_reaction_replicate_identifier_cannot_be_silently_merged(tmp_path):
    def transform(content):
        def mutate(rows):
            rows[12]['Condition_Name'] = rows[0]['Condition_Name']
        return raw_change(content, mutate)
    changed = archive_with_change(tmp_path, AUDIT.RAW_MEMBER, transform)
    with pytest.raises(ValueError, match='Ambiguous raw condition'):
        AUDIT.audit(changed)


def test_missing_observation_is_an_error_instead_of_an_exclusion(tmp_path):
    changed = archive_with_change(tmp_path, AUDIT.RAW_MEMBER,
        lambda content:raw_change(content, lambda rows:rows.pop()))
    with pytest.raises(ValueError, match='complete 576-row'):
        AUDIT.audit(changed)


def test_changed_cached_equilibrium_constant_is_rejected(tmp_path):
    def transform(workbook):
        original = zipfile.ZipFile(io.BytesIO(workbook))
        stream = io.BytesIO()
        with original, zipfile.ZipFile(stream, 'w') as output:
            for info in original.infolist():
                content = original.read(info)
                if info.filename == 'xl/worksheets/sheet1.xml':
                    root = ET.fromstring(content)
                    cell = next(c for c in root.findall('.//m:c', AUDIT.NS) if c.get('r') == 'E3')
                    cell.find('m:v', AUDIT.NS).text = '99'
                    content = ET.tostring(root)
                output.writestr(info, content)
        return stream.getvalue()
    changed = archive_with_change(tmp_path, AUDIT.WORKBOOK_MEMBER, transform)
    with pytest.raises(ValueError, match='constant arithmetic mismatch'):
        AUDIT.audit(changed)
