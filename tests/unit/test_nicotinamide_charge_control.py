"""Reject invalid SCF bookkeeping without treating coupling as cycle closure."""
import importlib.util
from pathlib import Path
import pytest

SOURCE=Path(__file__).resolve().parents[2]/'benchmark/audit_nicotinamide_charge_control.py'
spec=importlib.util.spec_from_file_location('charge_control_audit',SOURCE)
audit=importlib.util.module_from_spec(spec);spec.loader.exec_module(audit)

def fixture():
    return dict(passed=True,converged=True,n_AO=4,n_molecular_orbitals=4,n_electrons=2,
        overlap_minimum_eigenvalue=1e-5,orbital_gradient_norm=1e-7,
        energy_Eh=-4.17,CDS_Eh=.03,surface_discretization_method='SWIG',
        HOMO_Eh=-.2,LUMO_Eh=.1,negative_HOMO_diagnostic=True,
        scf_summary=dict(e1=-7.,e2=2.,nuc=1.,coul=3.,exc=-1.,e_solvent=-.2,e_cds=.03))
NUMERIC=dict(overlap_minimum_eigenvalue=1e-8,stationary_gradient_maximum=5e-5)

def test_aqueous_total_accounts_for_cds_once():
    audit.point_checks(fixture(),NUMERIC,4,2)

@pytest.mark.parametrize('mutation',[
    lambda p:p.update(energy_Eh=-4.14),
    lambda p:p.update(energy_Eh=-4.2),
    lambda p:p.update(n_molecular_orbitals=3),
    lambda p:p.update(n_electrons=4),
    lambda p:p.update(converged=False),
    lambda p:p.update(orbital_gradient_norm=float('nan')),
    lambda p:p.update(overlap_minimum_eigenvalue=1e-10),
    lambda p:p.update(negative_HOMO_diagnostic=False),
    lambda p:p.update(HOMO_Eh=float('-inf')),
    lambda p:p.update(CDS_Eh=.06),
])
def test_rejects_invalid_aqueous_point(mutation):
    point=fixture();mutation(point)
    with pytest.raises(ValueError):audit.point_checks(point,NUMERIC,4,2)

def test_homo_sign_is_recorded_without_claiming_stability():
    point=fixture();point.update(HOMO_Eh=.05,negative_HOMO_diagnostic=False)
    audit.point_checks(point,NUMERIC,4,2)

def test_published_charge_exchange_is_nonzero_and_balanced():
    packet=SOURCE.parent/'nicotinamide_charge_validation'
    report=audit.audit(packet)
    assert report['coupling_kj']==pytest.approx(-5.799882878988578,abs=1e-7,rel=0)
    assert report['element_residual']=={} and report['charge_residual']==0
    assert report['n_new_points']==2 and report['n_parent_points']==2
    assert report['empirical_offsets_fitted']==0
