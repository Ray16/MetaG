"""Prepare source-defined PFAS component inputs without inspecting model errors."""
from pathlib import Path
import collections
import csv
import hashlib
import io
import json
import re
import tarfile
import xml.etree.ElementTree as ET

from rdkit import Chem
from rdkit.Chem import rdDetermineBonds

HERE = Path(__file__).resolve().parent
UP = HERE / 'upstream'
sha = lambda data: hashlib.sha256(data).hexdigest()
def write(name, value):
    (HERE / name).write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')
def graph(mol):
    mol = Chem.Mol(mol)
    Chem.RemoveStereochemistry(mol)
    for atom in mol.GetAtoms():
        atom.SetNoImplicit(True)
        atom.SetIsAromatic(False)
    for bond in mol.GetBonds():
        bond.SetBondType(Chem.BondType.SINGLE)
        bond.SetIsAromatic(False)
    return Chem.MolToSmiles(mol, canonical=True, allHsExplicit=True)

all_refs = list(csv.DictReader((UP / 'data/references/dgsolv-references.csv').open()))
refs = [r for r in all_refs if r['Solvent'] == 'water'
        and r['Ref.'].strip() == '10.1021/acs.est.3c02545']
assert len(refs) == len({r['FlexiSol Name'] for r in refs}) == 21
root = ET.parse(HERE / 'Endo2023.xml')
table2 = root.find('.//table-wrap[@id="tbl2"]')
source_rows = table2.findall('.//tbody/tr')
assert len(source_rows) == 21
si = (HERE / 'Endo2023_SI.txt').read_text()
si = si.split('Table S2. SMILES strings for PFAS used in this study. ')[-1].split('S12')[0]
smiles = {}
for line in si.splitlines():
    match = re.fullmatch(r'(.+?)\s+([A-Za-z0-9@+\-=()\[\]#]+)\s*', line)
    if match and Chem.MolFromSmiles(match[2]) is not None:
        smiles[match[1]] = match[2]
assert len(smiles) == 21, smiles.keys()
aliases = {'6:2 FTAC': '6-2fta', '4:2 FTMAC': '4-2ftma'}
reference_records = []
lookup = {r['FlexiSol Name']: r for r in refs}
acid_names = {'PFBSA', 'PFHxSA', 'PFOSA', 'MeFBSA', 'MeFHxSA', 'EtFHxSA'}
for tr in source_rows:
    cells = [' '.join(''.join(x.itertext()).split()) for x in tr]
    assert len(cells) == 9
    label = cells[0]
    name = aliases.get(label, label.lower().replace(':', '-').replace(' ', ''))
    ref = lookup[name]
    mol = Chem.MolFromSmiles(smiles[label])
    cached = Chem.MolFromSmiles(ref['SMILES'])
    identity_ok = Chem.MolToSmiles(mol) == Chem.MolToSmiles(cached)
    logk = float(cells[7].replace('–', '-').replace('−', '-'))
    target = 0.00831446261815324 * 298.15 * __import__('math').log(10) * logk
    refvalue = float(next(v for k, v in ref.items() if k.startswith('Value')))
    # A bounded rounding/unit check, never an empirical offset adjustment.
    delta = 4.184 * refvalue - target
    assert abs(delta) < 0.03, (name, delta)
    reference_records.append(dict(name=name, source_label=label, source_SMILES=smiles[label],
        source_log10_Kaw=logk, target_molar_hydration_kj=target,
        FlexiSol_reference_kcal=refvalue, source_conversion_difference_kj=delta,
        aqueous_recipe='1 mM HCl' if label in acid_names else 'pure water',
        upstream_SMILES=ref['SMILES'], covalent_identity_matches_source=identity_ok,
        component_accuracy_eligible=identity_ok and label not in acid_names,
        reason='different covalent molecule in upstream structures' if not identity_ok else
        ('native pure-water neutral-solute observable' if label not in acid_names
         else 'acid medium; intrinsic neutral activity/speciation bridge not established'),
        commercial_isomer_mixture=label in {'PFHxSA', 'PFOSA'}))
reference_records.sort(key=lambda r: r['name'])
assert sum(r['component_accuracy_eligible'] for r in reference_records) == 13
assert {r['name'] for r in reference_records if not r['covalent_identity_matches_source']} == {'6-2fta', '4-2ftma'}
names = set(lookup)
structures = [r for r in csv.DictReader((UP / 'data/raw_energies/structures.csv').open())
              if r['solvent'] == 'water' and r['name'] in names]
assert len(structures) == len({r['path'] for r in structures}) == 648
assert {r['charge'] for r in structures} == {'0'}
assert {r['tautomer'] for r in structures} == {'0'}
structures.sort(key=lambda r: r['path'])
selected = {r['path']: r for r in structures}
raw_energies = [r for r in csv.DictReader((UP / 'data/raw_energies/energies.csv').open())
                if r['path'] in selected]
assert len(raw_energies) == 648
energy_keys = ['el_gfn2', 'el_r2scan-3c', 'el_wb97m-v_def2-tzvppd', 'alpb', 'smd']
energy_lookup = {r['path']: {k: float(r[k]) for k in energy_keys} for r in raw_energies}
for row in energy_lookup.values():
    assert all(__import__('math').isfinite(x) for x in row.values())
xyz_data = {}
with tarfile.open(HERE / 'flexisol_pinned.tar.gz', 'r:gz') as tar:
    for item in tar:
        rel = '/'.join(item.name.split('/')[1:])
        if rel.endswith('/coord.xyz') and rel[:-10] in selected:
            assert item.isfile()
            data = tar.extractfile(item).read()
            xyz_data[rel[:-10]] = data
assert len(xyz_data) == 648
expected = {r['name']: Chem.AddHs(Chem.MolFromSmiles(r['upstream_SMILES']))
            for r in reference_records}
inputs = []
(HERE / 'coordinates').mkdir(exist_ok=True)
for row in structures:
    data = xyz_data[row['path']]
    block = data.decode()
    mol = Chem.MolFromXYZBlock(block)
    assert mol is not None
    rdDetermineBonds.DetermineConnectivity(mol, charge=0)
    assert graph(mol) == graph(expected[row['name']]), row['path']
    lines = block.splitlines()
    n = int(lines[0])
    assert len(lines[2:]) == n
    symbols, xyz = [], []
    for line in lines[2:]:
        symbol, *coords = line.split()
        assert len(coords) == 3
        symbols.append(symbol)
        xyz.append(list(map(float, coords)))
    key = sha(row['path'].encode())[:24]
    (HERE / 'coordinates' / (key + '.xyz')).write_bytes(data)
    inputs.append(dict(key=key, **row, symbols=symbols, positions_A=xyz,
        coordinate_sha256=sha(data), coordinate_file=f'coordinates/{key}.xyz',
        published_energies_Eh=energy_lookup[row['path']], upstream_connectivity_passed=True,
        covalent_identity_matches_source=next(r['covalent_identity_matches_source']
            for r in reference_records if r['name'] == row['name'])))
write('references.json', reference_records)
write('inputs.json', inputs)
with (HERE / 'selected_upstream_energy_rows.csv').open('w', newline='') as stream:
    writer = csv.DictWriter(stream, fieldnames=list(raw_energies[0]))
    writer.writeheader()
    writer.writerows(sorted(raw_energies, key=lambda r: r['path']))
write('source_preparation.json', dict(structures=len(inputs), compounds=len(names),
    pure_water_accuracy_compounds=13, acid_medium_held_compounds=6,
    covalent_identity_held_compounds=2,
    phases=dict(collections.Counter(r['structuremode'] for r in inputs)),
    ordered_coordinate_and_reference_graphs_checked=648,
    no_prediction_or_error_selection=True,
    source_conventions='Kaw=c_air/c_water; molar gas-to-water hydration = RT ln(Kaw).',
    script_sha256=sha(Path(__file__).read_bytes()),
    inputs_sha256=sha((HERE/'inputs.json').read_bytes()),
    references_sha256=sha((HERE/'references.json').read_bytes())))
print('Prepared all648 source structures for21 PFAS;13 pure-water targets,6 acid-medium and2 identity holds.')
