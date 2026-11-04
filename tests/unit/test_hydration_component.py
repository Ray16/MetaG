"""Thermodynamic aggregation and independent primary-source component arithmetic."""
from decimal import Decimal, localcontext
import importlib.util
import json
import math
from pathlib import Path
import shutil

import pytest

ROOT=Path(__file__).resolve().parents[2]
spec=importlib.util.spec_from_file_location('hydration_component',ROOT/'benchmark/hydration_component.py')
module=importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
PACKET=ROOT/'benchmark/hydration_component_validation'


def test_mean_energy_and_partition_free_energy_are_distinct():
    rt=2.478957
    assert module.aggregate([0.,0.],'minimum',rt)==0
    assert module.aggregate([0.,0.],'boltzmann_mean',rt)==0
    assert module.aggregate([0.,0.],'log_partition',rt)==pytest.approx(-rt*math.log(2))
    values=[-18.,-17.,-5.]
    for method in ['minimum','boltzmann_mean','log_partition']:
        assert module.aggregate([x+91.5 for x in values],method,rt)-91.5==pytest.approx(module.aggregate(values,method,rt))


def test_public_panel_replay_retains_holds_and_coverage():
    assert module.verify_manifest(PACKET)>700
    result=module.analyze(PACKET)
    assert result==json.loads((PACKET/'results.json').read_text())
    assert (result['n_source_compounds'],result['n_reference_eligible'],result['n_single_points'])==(21,13,461)
    assert len(result['holds'])==8
    assert result['experimental_validation_regime']=='development_exposed_component_ablation'
    assert result['primary_metabolic_references_admitted']==0
    assert result['primary_paired_contrasts']['UMA|smd|full|minimum']['n_improved']==13


def test_independent_decimal_primary_minima():
    """No call to the production aggregation or metric code for these sums."""
    p=json.loads((PACKET/'protocol.json').read_text())
    inputs=json.loads((PACKET/'inputs.json').read_text())
    refs=[r for r in json.loads((PACKET/'references.json').read_text()) if r['component_accuracy_eligible']]
    preds={r['key']:r for path in (PACKET/'predictions').glob('*.json')
           for r in json.loads(path.read_text())['records']}
    results=json.loads((PACKET/'results.json').read_text())
    D=lambda value:Decimal(str(value))
    with localcontext() as context:
        context.prec=50
        for electronic in ['UMA','el_wb97m-v_def2-tzvppd']:
            for solvent in ['alpb','smd']:
                errors=[]
                for ref in refs:
                    gas=[];aqueous=[]
                    for row in inputs:
                        if row['name']!=ref['name']:continue
                        e=(D(preds[row['key']]['energy_eV'])*D(p['eV_to_kj']) if electronic=='UMA'
                           else D(row['published_energies_Eh'][electronic])*D(p['Eh_to_kj']))
                        if row['structuremode']=='gas':gas.append(e)
                        else:aqueous.append(e+D(row['published_energies_Eh'][solvent])*D(p['Eh_to_kj']))
                    target=D(p['R_kj'])*D(p['temperature_K'])*D(10).ln()*D(ref['source_log10_Kaw'])
                    errors.append(abs(min(aqueous)-min(gas)-target))
                mae=sum(errors)/D(len(errors))
                assert abs(mae-D(results['arms'][f'{electronic}|{solvent}|full|minimum']['MAE_kj']))<D('1e-7')


def test_missing_geometry_cannot_silently_improve_error(tmp_path):
    # Use a copied real packet; removing even one conformer must invalidate it.
    packet=tmp_path/'packet'
    shutil.copytree(PACKET,packet)
    first=next((packet/'predictions').glob('*.json'))
    first.unlink()
    with pytest.raises(ValueError,match='inventory'):
        module.analyze(packet)


def test_source_identity_holds_are_not_prediction_failures():
    refs=json.loads((PACKET/'references.json').read_text())
    mismatch={r['name'] for r in refs if not r['covalent_identity_matches_source']}
    assert mismatch=={'4-2ftma','6-2fta'}
    from rdkit import Chem
    from rdkit.Chem import rdMolDescriptors
    for ref in refs:
        source=Chem.MolFromSmiles(ref['source_SMILES'])
        upstream=Chem.MolFromSmiles(ref['upstream_SMILES'])
        if ref['name'] in mismatch:
            assert upstream.GetNumAtoms()==source.GetNumAtoms()+1
            assert rdMolDescriptors.CalcMolFormula(upstream)!=rdMolDescriptors.CalcMolFormula(source)
        else:
            assert Chem.MolToSmiles(source)==Chem.MolToSmiles(upstream)
