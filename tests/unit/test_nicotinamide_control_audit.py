"""Adversarial checks for physical-control bookkeeping and rank ownership."""
from copy import deepcopy
import importlib.util
from pathlib import Path

import pytest

SOURCE=Path(__file__).resolve().parents[2]/'benchmark/audit_nicotinamide_control.py'
spec=importlib.util.spec_from_file_location('nicotinamide_audit',SOURCE)
audit=importlib.util.module_from_spec(spec)
spec.loader.exec_module(audit)

def fixture():
    geometry=dict(symbols=['H','H'],positions_angstrom=[[0,0,0],[0,0,.74]],charge=0)
    item=dict(geometry=geometry,n_AO=4,n_electrons=2)
    summary=dict(e1=-7.,e2=2.,nuc=1.,coul=3.,exc=-1.)
    gas=dict(energy_Eh=-4.,converged=True,n_molecular_orbitals=4,
             orbital_gradient_norm=1e-7,scf_summary=summary)
    water=deepcopy(gas)
    water['energy_Eh']=-4.17
    water['scf_summary'].update(e_solvent=-.2,e_cds=.03)
    point=dict(passed=True,geometry=deepcopy(geometry),n_AO=4,n_electrons=2,
               overlap_minimum_eigenvalue=1e-5,gas=gas,water=water,
               SMD_delta_Eh=-.17,CDS_Eh=.03,surface_discretization_method='SWIG')
    protocol=dict(overlap_minimum_eigenvalue=1e-8,stationary_gradient_maximum=5e-5)
    return item,point,protocol

def test_correct_total_includes_cds_once():
    audit.audit_point(*fixture())

@pytest.mark.parametrize('mutation',[
    lambda p:p['water'].update(energy_Eh=-4.14),  # CDS added twice
    lambda p:p['water'].update(energy_Eh=-4.2),   # CDS omitted
    lambda p:p['gas'].update(n_molecular_orbitals=3),
    lambda p:p['water'].update(converged=False),
    lambda p:p['gas'].update(orbital_gradient_norm=.1),
    lambda p:p['gas'].update(energy_Eh=float('nan')),
    lambda p:p.update(n_electrons=4),
    lambda p:p['geometry'].update(charge=1),
    lambda p:p.update(SMD_delta_Eh=-.2),
    lambda p:p.update(CDS_Eh=.04),
    lambda p:p.update(overlap_minimum_eigenvalue=1e-10),
])
def test_rejects_invalid_physical_control(mutation):
    item,point,protocol=fixture()
    mutation(point)
    with pytest.raises(ValueError):audit.audit_point(item,point,protocol)

def test_pool_changes_reference_without_changing_probabilities():
    rows=[dict(electronic_kj=-100.,solvation_kj=-5.,thermal_kj=1.),
          dict(electronic_kj=-98.,solvation_kj=-4.,thermal_kj=1.5)]
    original=audit.pool(rows,2.4789)
    shifted=audit.pool([dict(row,electronic_kj=row['electronic_kj']+1234.) for row in rows],2.4789)
    assert shifted['G_kj']-original['G_kj']==pytest.approx(1234.,abs=1e-10)
    assert shifted['probabilities']==pytest.approx(original['probabilities'],abs=1e-12)

def test_published_packet_reconstructs_all_four_arms():
    packet=SOURCE.parent/'nicotinamide_validation'
    report=audit.audit(packet,packet/'raw_records')
    assert report['n_slots']==9
    assert report['independent_experimental_labels_admitted']==0
    # Independent decimal reconstruction in the manuscript checker supplies
    # these values; no backend or external workspace is needed by this test.
    expected={'UMA_ALPB':-41.79569172164487,'UMA_SMD':-43.52053462832513,
              'DFT_ALPB':-41.91735383775112,'DFT_SMD':-43.68302454934846}
    for arm,value in expected.items():
        assert report['arms'][arm]['exchange_kj']==pytest.approx(value,abs=1e-7,rel=0)
