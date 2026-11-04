"""Authenticate native apparent-equilibrium arithmetic, without model prediction.

Kaspar et al., doi:10.1002/adsc.201901230; source data doi:10.5281/zenodo.3530526.
The buffered concentration ratios audited here are not activity-based standard
reaction Gibbs energies. Source identities and condition transfer need separate
adjudication before independent benchmark enrollment.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import math
from pathlib import Path
import xml.etree.ElementTree as ET
import zipfile

RAW_MEMBER = 'Data_Figure_4/Results_Figure_4.csv'
WORKBOOK_MEMBER = 'Data_Figure_4/Compiled_Raw_Data_Figure_4.xlsx'
NS = {'m': 'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
# Source SI Table S1 and workbook block order. Aliases are explicit, including
# oxoThd/5MUrd and Ino2; suffix stripping must not invent replicate assignments.
COMPOUNDS = (
    ('Urd', 1, 'Urd', 'Urd_2', 'Urd'),
    ('dUrd', 2, 'dUrd', 'dUrd_2', 'dUrd'),
    ('5MUrd', 9, 'oxoThd', 'oxoThd_2', '5MUrd'),
    ('Thd', 10, 'dThd', 'dThd_2', 'dThd'),
    ('5FUrd', 3, '5FUrd', '5FUrd_2', '5-F-Urd'),
    ('5FdUrd', 4, '5FdUrd', '5FdUrd_2', '5-F-dUrd'),
    ('5BrUrd', 11, '5BUd', '5BUd_2', '5-Br-Urd'),
    ('5BrdUrd', 12, 'd5BUrd', 'd5BUrd_2', '5-Br-dUrd'),
    ('5IUrd', 13, '5IUrd', '5IUrd_2', '5-I-Urd'),
    ('5IdUrd', 14, '5IdUrd', '5IdUrd_2', '5-I-dUrd'),
    ('5EUrd', 15, '5EUrd', '5EUrd_2', '5-Ethynyl-Urd'),
    ('5EdUrd', 16, '5EdUrd', '5EdUrd_2', '5-Ethynyl-dUrd'),
    ('Ado', 5, 'Ado', 'Ado_2', 'Ado'),
    ('dAdo', 6, 'dAdo', 'dAdo_2', 'dAdo'),
    ('2FAdo', 17, '2FAdo', '2FAdo_2', '2-F-Ado'),
    ('2FdAdo', 18, '2FdAdo', '2FdAdo_2', '2-F-dAdo'),
    ('2ClAdo', 7, '2ClAdo', '2ClAdo_2', '2-Cl-Ado'),
    ('2CldAdo', 8, '2CldAdo', '2CldAdo_2', '2-Cl-dAdo'),
    ('2AAdo', 19, '2AAdo', '2AAdo_2', '2-NH2-Ado'),
    # The original right-hand heading duplicates 2-NH2-Ado. This is not silently
    # rewritten: Table S1 and raw identifiers distinguish compound20 explicitly.
    ('2AdAdo', 20, '2AdAdo', '2AdAdo_2', '2-NH2-Ado'),
    ('Guo', 21, 'Guo', 'Guo_2', 'Guo'),
    ('dGuo', 22, 'dGuo', 'dGuo_2', 'dGuo'),
    ('Ino', 23, 'Ino', 'Ino2', 'Ino'),
    ('dIno', 24, 'dIno', 'dIno_2', 'dIno'),
)


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def apparent_constant(conversion: float, phosphate_excess: float) -> float:
    """Nucleoside + Pi = base + pentose1P; initially no products.

    x is fractional conversion; r = initial Pi / initial nucleoside. Because
    there are two solute molecules on each side, molarity units cancel exactly.
    This does not remove activity coefficients or proton/Mg state dependence.
    """
    if (not math.isfinite(conversion) or not math.isfinite(phosphate_excess)
            or not 0 < conversion < 1 or not phosphate_excess > conversion):
        raise ValueError('Conversion or phosphate excess outside physical domain')
    return conversion**2 / ((1-conversion)*(phosphate_excess-conversion))


def worksheet_cells(workbook: bytes) -> dict:
    """Read cached source cells, preserving formulas and errors; never execute them."""
    with zipfile.ZipFile(io.BytesIO(workbook)) as archive:
        shared = []
        if 'xl/sharedStrings.xml' in archive.namelist():
            root = ET.fromstring(archive.read('xl/sharedStrings.xml'))
            shared = [''.join(t.text or '' for t in s.findall('.//m:t', NS))
                      for s in root.findall('m:si', NS)]
        root = ET.fromstring(archive.read('xl/worksheets/sheet1.xml'))
    result = {}
    for cell in root.findall('m:sheetData/m:row/m:c', NS):
        address = cell.attrib['r']
        if address in result:
            raise ValueError('Duplicate worksheet cell: '+address)
        node = cell.find('m:v', NS)
        value = node.text if node is not None else None
        kind = cell.get('t', 'n')
        if kind == 's' and value is not None:
            value = shared[int(value)]
        elif kind == 'inlineStr':
            value = ''.join(t.text or '' for t in cell.findall('.//m:t', NS))
        formula = cell.find('m:f', NS)
        result[address] = dict(value=value, kind=kind,
            formula=formula.text if formula is not None else None,
            formula_attributes=dict(formula.attrib) if formula is not None else None)
    return result


def numeric(cell: dict | None, address: str) -> float | None:
    if cell is None or cell['value'] is None:
        return None
    if cell['kind'] != 'n':
        raise ValueError('Expected numeric source cell: '+address)
    value = float(cell['value'])
    if not math.isfinite(value):
        raise ValueError('Nonfinite source cell: '+address)
    return value


def audit(source_zip: Path) -> dict:
    data = source_zip.read_bytes()
    with zipfile.ZipFile(io.BytesIO(data)) as archive:
        if len(set(archive.namelist())) != len(archive.namelist()):
            raise ValueError('Ambiguous duplicate archive members')
        raw_bytes = archive.read(RAW_MEMBER)
        workbook = archive.read(WORKBOOK_MEMBER)
    raw = list(csv.DictReader(io.StringIO(raw_bytes.decode('utf-8-sig'))))
    if len(raw) != 576:
        raise ValueError('Expected complete 576-row source table')
    cells = worksheet_cells(workbook)
    index = {}
    for position, row in enumerate(raw):
        if int(row['']) != position:
            raise ValueError('Source row identity/order changed')
        key = row['Condition_Name'], row['Time_Point']
        if key in index:
            raise ValueError('Ambiguous raw condition/time key')
        index[key] = position, row
    records, issues = [], []
    maximum_concentration_difference = maximum_constant_difference = 0.0
    for ordinal, (name, compound, first, second, heading) in enumerate(COMPOUNDS):
        start = 1+16*(ordinal//2)
        left = ordinal % 2 == 0
        header_cell = ('B' if left else 'K')+str(start)
        if cells[header_cell]['value'] != heading:
            raise ValueError('Source compound heading changed: '+header_cell)
        if name == '2AdAdo':
            issues.append(dict(kind='duplicated_compound_heading', cell=header_cell,
                source_heading=heading, source_raw_identifier=first, SI_compound_number=compound,
                qualification='Distinct deoxy identity follows raw identifiers and SI TableS1, not this duplicated heading; structures need separate source review.'))
        for replicate, alias in enumerate((first, second), 1):
            product_column = ('C', 'F')[replicate-1] if left else ('L', 'O')[replicate-1]
            constant_column = ('E', 'H')[replicate-1] if left else ('N', 'Q')[replicate-1]
            for temperature_index, temperature_C in enumerate((40, 50, 60, 70)):
                for aliquot in (1, 2, 3):
                    position, row = index.pop((alias, f'{temperature_C}.{aliquot}'))
                    conversion = float(row['Product'])
                    substrate = float(row['Substrate'])
                    if (not math.isfinite(conversion) or not math.isfinite(substrate)
                            or not math.isclose(conversion+substrate, 1.0, abs_tol=1e-12, rel_tol=0)):
                        raise ValueError('Invalid source spectral fractions')
                    line = start+2+3*temperature_index+aliquot-1
                    product_address, constant_address = product_column+str(line), constant_column+str(line)
                    compiled_product = numeric(cells.get(product_address), product_address)
                    included = compiled_product is not None
                    try:
                        constant = apparent_constant(conversion, 5.0)
                    except ValueError:
                        if included:
                            raise ValueError('Compiled point has unphysical conversion') from None
                        constant = None
                    if included:
                        delta = abs(compiled_product-2.0*conversion)
                        maximum_concentration_difference = max(maximum_concentration_difference, delta)
                        if delta > 1e-11:
                            raise ValueError('Raw-to-compiled product mismatch: '+product_address)
                        cached = numeric(cells.get(constant_address), constant_address)
                        if cached is None or not math.isclose(cached, constant, rel_tol=1e-10, abs_tol=1e-12):
                            raise ValueError('Source constant arithmetic mismatch: '+constant_address)
                        maximum_constant_difference = max(maximum_constant_difference, abs(cached-constant))
                    records.append(dict(source_row=position, compound=name, SI_compound_number=compound,
                        source_condition_name=alias, source_time_label=row['Time_Point'],
                        nominal_temperature_K=temperature_C+273.15,
                        reaction_replicate=replicate, equilibrium_aliquot=aliquot,
                        conversion_fraction=conversion, source_compilation_includes=included,
                        compiled_product_cell=product_address, compiled_constant_cell=constant_address,
                        apparent_concentration_equilibrium_ratio=constant,
                        experimental_reference_admitted=False))
    if index:
        raise ValueError('Unmapped source observations remain')
    errors = [dict(cell=address, error=c['value'], formula=c['formula'])
              for address, c in cells.items() if c['kind'] == 'e']
    included = sum(r['source_compilation_includes'] for r in records)
    if included != 519:
        raise ValueError('Source compilation mask differs from reported519points')
    return dict(schema='metag-nucleoside-reference-arithmetic-v1',
        source_sha256=sha(data), raw_member_sha256=sha(raw_bytes), workbook_member_sha256=sha(workbook),
        auditor_sha256=sha(Path(__file__).read_bytes()),
        counts=dict(compounds=24, temperatures=4, reaction_replicates_per_condition=2,
                    equilibrium_aliquots_per_reaction=3, raw_observations=576,
                    source_compilation_included=included, source_compilation_excluded=576-included),
        source_conditions=dict(initial_nucleoside_mM=2.0, initial_phosphate_mM=10.0,
            buffer='50mM glycine, nominal pH9', nominal_pH=9.0,
            actual_pH_at_each_temperature='not established by this arithmetic audit',
            ionic_strength_M=None, free_magnesium_M=None),
        maximum_raw_compiled_product_difference_mM=maximum_concentration_difference,
        maximum_reconstructed_cached_constant_difference=maximum_constant_difference,
        source_worksheet_error_cells=errors, identity_issues=issues,
        records=sorted(records, key=lambda r:r['source_row']),
        interpretation='Source arithmetic and author compilation mask only. Three aliquots from one reaction are not three independent equilibrations. Excluded observations are retained; the physical rationale for every exclusion is not authenticated by a missing cell.',
        scope='Buffered apparent concentration ratios, not activity-based standard Gibbs energies. No temperature, pH, ionic-strength or Mg bridge is admitted. Source/measurement/development/training ancestry and exact chemical-state identities remain to be adjudicated.',
        passed=True, predictions_performed=False, experimental_references_admitted=0,
        primary_benchmark_launch_allowed=False)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source-zip', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    report = audit(args.source_zip)
    args.output.write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k != 'records'}, indent=2))


if __name__ == '__main__':
    main()
