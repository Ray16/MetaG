"""Failed estimates keep auditable diagnostics; strict thermal validation cannot be bypassed."""
from types import SimpleNamespace

import pytest
from metag import pipeline as p
from metag.energetics import species as sp
from metag.settings import Settings


def _ctx(log=lambda *_: None, **settings):
    return sp._Context(None, Settings(**settings), {'cacheable': True}, log, 'test')


def _track(validated=None):
    return SimpleNamespace(validated=validated or {}, _key=lambda uniq, j: j)


class _Backend:
    model_name = 'fake-model'
    producer = {'cacheable': False}
    settings = Settings()

    def __init__(self, species=None):
        self._species = species

    def species(self, q, smiles, label, warnings, diagnostics):
        return self._species(q, smiles, label, warnings, diagnostics)

    def liquid_water(self):
        return 0.0


def test_missing_estimate_retains_auditable_original_species():
    species={'a':[-1,0,'CCO'],'b':[1,0,'CCO']}
    result=p._no_estimate('case',{'species':species},{},'failed',lambda *_:None)
    assert result['original_species']==species
    assert result['energy_units']=='kJ/mol'
    assert result['quantity']=='standard_transformed_reaction_gibbs'


def test_suspect_result_has_no_point_estimate(monkeypatch):
    def species(backend, rx, routes, G, sig, std, log):
        G.update(a=0, b=1000); sig.update(a=0.1, b=0.1)
        rx['n_Hplus'] = 1  # inject a late charge leak after input/routing validation
    monkeypatch.setattr(p, '_species_free_energies', species)
    result = p.score_reaction({'species': {'a': [-1, 0, 'C'], 'b': [1, 0, 'C']}}, _Backend(), log=lambda *_: None)
    assert result['dG'] is None and result['dG_raw'] is None
    assert 'charge imbalance' in result['suspect']


def test_seed_offset_is_in_reproducibility_config():
    backend = _Backend()
    backend.settings = Settings.from_env({'CONFORMER_SEED_OFFSET': '1000'})
    cfg = p.effective_config(backend)
    assert cfg['settings']['seed_offset'] == 1000
    assert '"seed_offset": 1000' in cfg['species_settings']


def test_strict_thermal_result_reused_only_for_exact_geometry(monkeypatch):
    import numpy as np
    from ase import Atoms
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from metag.energetics.thermal_ensemble import ThermalTrack
    mol=Chem.AddHs(Chem.MolFromSmiles('CO')); AllChem.EmbedMolecule(mol,randomSeed=17)
    at=Atoms([a.GetSymbol() for a in mol.GetAtoms()],positions=mol.GetConformer().GetPositions())
    u=sp.UniqueMinima(template=mol); u.add(at,-100.,-110.)
    calls=[]
    def correction(*args,**kwargs):
        calls.append(kwargs)
        return 5.,dict(n_imag=0,max_imag_cm=0.,_imag_vecs=[])
    track=ThermalTrack(None,mol,0,1,correction=correction,combine=sp.boltz,window=15,maximum=10)
    assert track.update(u)==pytest.approx(-105.)
    def forbidden(*args,**kwargs):raise AssertionError('exact validated geometry should reuse strict Hessian')
    monkeypatch.setattr(sp,'_thermal_at_minimum',forbidden)
    out=sp._validated_minima(_ctx(print),u,mol,None,0,1,track=track,diagnostics={})
    assert out[1]==[5.] and out[2]['n_strict_results_reused']==1
    row=out[2]['validated_minima'][0]
    np.testing.assert_array_equal(row['geometry']['positions_angstrom'],
                                  out[0].ref[0].GetConformer().GetPositions())
    assert row['geometry']['symbols']==[a.GetSymbol() for a in mol.GetAtoms()]
    assert row['electronic_kj']+row['solvation_kj']+row['thermal_kj']==row['total_kj']
    assert len(calls)==1
    # A tiny coordinate change must invalidate the strict evidence.
    x=np.array(u.ref[0].GetConformer().GetAtomPosition(0)); x[0]+=1e-6
    u.ref[0].GetConformer().SetAtomPosition(0,x)
    with pytest.raises(AssertionError,match='exact validated'):
        sp._validated_minima(_ctx(print),u,mol,None,0,1,track=track,diagnostics={})


def test_provisional_soft_modes_cannot_supply_strict_validation():
    from types import SimpleNamespace
    import numpy as np
    from metag.energetics.thermal_ensemble import ThermalTrack
    calls=[]
    def correction(*args,**kw):
        calls.append(kw.get('imag_as_soft',False))
        return 4.,dict(n_imag=0 if calls[-1] else 1,max_imag_cm=100.)
    mol=SimpleNamespace(GetAtoms=lambda:[])
    u=SimpleNamespace(G=[-10.],ref=[SimpleNamespace(GetConformer=lambda:SimpleNamespace(GetPositions=lambda:np.zeros((0,3))))])
    tr=ThermalTrack(None,mol,0,1,correction=correction,combine=sp.boltz,window=15,maximum=10)
    tr.update(u)
    assert calls==[False,True]
    assert tr.validated=={}


def test_failed_species_keeps_thermal_diagnostics():
    def fail(q, smi, label, warnings, diagnostics):
        diagnostics['thermal'] = {'status': 'incomplete', 'n_sampled': 3}
        return None, None
    result = p.score_reaction({'species': {'a': [-1, 0, 'CCO'], 'b': [1, 0, 'CC=O'], 'h2': [1, 0, '[HH]']}},
                              _Backend(fail), log=lambda *_: None)
    assert result['dG'] is None
    assert 'minimum validation failed' in result['suspect']
    assert result['routes']['species_diagnostics']['a']['thermal']['n_sampled'] == 3
    assert 'a' in result['routes']['failed_species']


def test_thermal_resolution_budget_reaches_optimizer(monkeypatch):
    import numpy as np
    seen=[]
    monkeypatch.setattr(sp,'batched_energies',lambda *a,**k:np.array([-1.]))
    def relax(pu,ats,**kw):
        seen.append(kw['steps']);return ats,np.array([-1.,-1.]),np.array([False,False])
    monkeypatch.setattr(sp,'batched_fire',relax)
    monkeypatch.setattr(sp,'same_connectivity',lambda *a:True)
    verdict,_=sp._follow_imaginary_mode(_ctx(thermal_relax_steps=1200),['H'],np.zeros((1,3)),np.zeros((1,3)),0,1,set())
    assert seen==[1200] and verdict=='inconclusive'


def test_validation_budget_failure_is_distinct_from_saddle(monkeypatch):
    from types import SimpleNamespace
    info={}
    assert sp._validated_minima(_ctx(thermal_validation_max=1),SimpleNamespace(G=[0.,1.]),None,None,0,1,
                                track=_track(),diagnostics=info) is None
    assert info['reason']=='candidate validation budget exceeded'
    assert info['n_sampled']==2


def test_unresolved_candidate_retried_without_dropping_ensemble(monkeypatch):
    from ase import Atoms
    from rdkit import Chem
    from rdkit.Chem import AllChem
    mol=Chem.AddHs(Chem.MolFromSmiles('CO'));AllChem.EmbedMolecule(mol,randomSeed=17)
    at=Atoms([a.GetSymbol() for a in mol.GetAtoms()],positions=mol.GetConformer().GetPositions())
    u=sp.UniqueMinima(template=mol);u.add(at,-100.,-110.)
    seen=[]
    def resolve(*args,**kwargs):
        steps=kwargs.get('optimization_steps',400);seen.append(steps)
        return (5. if steps==1200 else None),{},dict(index=0,kind='as_is',pos=at.positions,rejected=[])
    monkeypatch.setattr(sp,'_thermal_at_minimum',resolve)
    out=sp._validated_minima(_ctx(thermal_relax_steps=400,thermal_retry_steps=1200),u,mol,None,0,1,
                             track=_track(),diagnostics={})
    assert seen==[400,1200]
    assert out[0].G==[-110.] and out[1]==[5.]
    assert out[2]['n_candidates_retried']==1


@pytest.mark.parametrize('second,accepted',[('artefact',True),('inconclusive',False)])
def test_strongest_artefact_does_not_validate_other_negative_modes(monkeypatch,second,accepted):
    import numpy as np
    from ase import Atoms
    from rdkit import Chem
    from rdkit.Chem import AllChem
    mol=Chem.AddHs(Chem.MolFromSmiles('CO'));AllChem.EmbedMolecule(mol,randomSeed=17)
    at=Atoms([a.GetSymbol() for a in mol.GetAtoms()],positions=mol.GetConformer().GetPositions())
    u=sp.UniqueMinima(template=mol);u.add(at,-100.,-110.)
    soft_calls=[];modes=[np.zeros_like(at.positions),np.ones_like(at.positions)]
    def correction(*args,**kwargs):
        soft=kwargs.get('imag_as_soft',False)
        if soft:soft_calls.append(True)
        return 5.,dict(n_imag=0 if soft else 2,max_imag_cm=100.,_imag_vecs=list(modes))
    verdicts=iter(['artefact',second]); probes=[]
    def follow(*args,**kwargs):
        probes.append(True);return next(verdicts),None
    monkeypatch.setattr(sp,'uma_gibbs_corr',correction)
    monkeypatch.setattr(sp,'_follow_imaginary_mode',follow)
    monkeypatch.setattr(sp,'batched_fire',lambda pu,ats,**kw:(ats,np.array([-1.]),np.array([True])))
    monkeypatch.setattr(sp,'same_connectivity',lambda *a:True)
    monkeypatch.setattr(sp,'same_stereochemistry',lambda *a:True)
    g,_,_=sp._thermal_at_minimum(_ctx(),u,mol,None,0,1)
    assert len(probes)==2
    assert (g is not None)==accepted
    assert bool(soft_calls)==accepted
