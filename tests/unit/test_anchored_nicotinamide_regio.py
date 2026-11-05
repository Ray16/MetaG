"""Two regioisomers require distinct cores and one consistent scaffold reference."""
import math
import pytest
from rdkit import Chem

from benchmark.anchored_nicotinamide_regio import anchor_plan,RegioAnchoredEstimator,CORE_12,CORE_14,OX_CORE
from metag.state_function import _combine_sampling
from metag.settings import Settings

RIBOSE='[C@@H]2O[C@H](CO)[C@@H](O)[C@H]2O'
RED14='NC(=O)C1=CN('+RIBOSE+')C=CC1'
RED12='NC(=O)C1=CC=CN('+RIBOSE+')C1'
canonical=lambda s:Chem.MolToSmiles(Chem.MolFromSmiles(s),isomericSmiles=True)


def test_distinct_cores_share_exact_oxidized_reference_and_preserve_atom_order_invariance():
    a,b=anchor_plan(0,RED14),anchor_plan(0,RED12)
    assert a['raw_terms'][1]['smiles']==CORE_14
    assert b['raw_terms'][1]['smiles']==CORE_12 and CORE_14!=CORE_12
    assert a['oxidized_sibling_smiles']==b['oxidized_sibling_smiles']
    molecule=Chem.MolFromSmiles(RED12)
    reordered=Chem.RenumberAtoms(molecule,list(reversed(range(molecule.GetNumAtoms()))))
    assert anchor_plan(0,Chem.MolToSmiles(reordered,canonical=False))==b
    phosphorylated=RED12.replace('CO)','COP(=O)(O)O)')
    assert anchor_plan(0,phosphorylated)['oxidized_sibling_smiles']!=b['oxidized_sibling_smiles']


def test_unknown_isomer_cannot_silently_use_a_declared_core():
    other='N1('+RIBOSE+')CC=CC(C(N)=O)=C1'
    with pytest.raises(ValueError,match='no declared core'):
        anchor_plan(0,other)
    assert anchor_plan(0,CORE_12) is None
    with pytest.raises(ValueError,match='charge'):anchor_plan(1,RED12)
    with pytest.raises(ValueError,match='Isotope'):anchor_plan(0,RED12.replace('NC','[15NH2]C',1))


class Backend:
    settings=Settings()
    model_name='synthetic-regio-transfer'
    producer={'cacheable':False}
    def species(self,q,smi,label,warnings,diagnostics):
        ox=anchor_plan(0,RED14)['oxidized_sibling_smiles']
        return {ox:(100.,7.),CORE_14:(10.,3.),CORE_12:(22.,5.),OX_CORE:(2.,11.)}[smi]
    def liquid_water(self):return 0.


def test_isomer_sampling_cancels_common_full_reference_and_oxidized_core():
    estimator=RegioAnchoredEstimator(Backend())
    a,b=estimator.state(0,RED14),estimator.state(0,RED12)
    assert b['solution_reference_kj']-a['solution_reference_kj']==pytest.approx(12.)
    sigma,inputs=_combine_sampling([(1.,b),(-1.,a)])
    assert sigma==pytest.approx(math.sqrt(34.))
    assert all(p['sensitivity']==0 for p in inputs.values() if p['smiles'] not in {CORE_14,CORE_12})
    other=RegioAnchoredEstimator(Backend())
    assert other.state(0,RED12)==b
    assert other.state(0,RED14)==a
    assert estimator.ledger()['policy'].endswith('v2')
