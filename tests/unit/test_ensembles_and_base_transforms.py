"""Regression tests for ensemble and speciation errors; no external QM calls."""
import math

import numpy as np
import pytest
from ase import Atoms
from rdkit import Chem
from rdkit.Chem import AllChem

from types import SimpleNamespace

from metag import pipeline as P
from metag.constants import EV2KJ
from metag.energetics import species as SP
from metag.settings import Settings
from metag.energetics.microstates import combine
from metag.routing import aldehyde_hydration as hydration
from metag.routing import pka_transform as pk


def methane():
    mol = Chem.AddHs(Chem.MolFromSmiles("C"))
    AllChem.EmbedMolecule(mol, randomSeed=1)
    atoms = Atoms([a.GetSymbol() for a in mol.GetAtoms()], mol.GetConformer().GetPositions())
    return mol, atoms


def run_methane(settings=None, warnings=None, diagnostics=None):
    """species_free_energy for methane with the mocked QM of `species_backend` (3 seed batches max)."""
    return SP.species_free_energy(None, 0, "C", "methane", settings=settings or Settings(conv_max=3),
                                  model_name="mock-uma", producer={"cacheable": True}, log=lambda *_: None,
                                  warnings=warnings, diagnostics=diagnostics)


@pytest.fixture
def species_backend(monkeypatch):
    _, atoms = methane()
    monkeypatch.setattr(SP.cache, "get", lambda *a, **k: None)
    writes = []
    monkeypatch.setattr(SP.cache, "put", lambda *a, **k: writes.append((a, k)))
    monkeypatch.setattr(SP, "pool_confs", lambda *a, **k: [atoms.copy()])
    monkeypatch.setattr(SP, "batched_energies", lambda *a, **k: np.array([-100.0 / EV2KJ]))
    monkeypatch.setattr(SP, "batched_fire", lambda pu, ats, **k:
                        (ats, np.full(len(ats), -100.0 / EV2KJ), np.ones(len(ats), dtype=bool)))
    monkeypatch.setattr(SP, "dgsolv", lambda *a, **k: -10.0)
    monkeypatch.setattr(SP, "same_connectivity", lambda *a: True)
    monkeypatch.setattr(SP, "uma_gibbs_corr", lambda *a, **k:
                        (5.0, {"n_imag": 0, "max_imag_cm": 0.0, "_imag_vecs": []}))
    return writes


def test_vertical_ensemble_combines_electronic_solvation_and_rrho_once(species_backend):
    value, _ = run_methane()
    assert value == pytest.approx(-105.0)  # electronic -100, solvation -10, RRHO +5
    assert species_backend[-1][0][4] == pytest.approx(-105.0)


@pytest.mark.parametrize("gap", [-3000.0, 3000.0])
def test_hydration_partition_is_stable_for_large_state_gaps(gap):
    result = hydration.mixture_G(-1e6, -1e6 + gap, 0.0)
    assert math.isfinite(result)
    assert result == pytest.approx(-1e6 + min(gap, 0))
    assert result == pytest.approx(hydration.mixture_G_states(-1e6, [(1, -1e6 + gap)], 0))


def test_duplicate_microstate_labels_are_rejected():
    with pytest.raises(ValueError, match="duplicate microstate"):
        combine([{"name": "a", "G": 0.0}, {"name": "a", "G": 10.0}])


@pytest.mark.parametrize("neutral,charged", [
    ("N", "[NH4+]"), ("CN", "C[NH3+]"), ("CNC", "C[NH2+]C"),
    ("CN(C)C", "C[NH+](C)C"), ("NCC(=O)O", "[NH3+]CC(=O)[O-]"),
    ("CN(C)CC(=O)O", "C[NH+](C)CC(=O)[O-]"),
])
def test_base_transform_independent_of_input_protonation_on_neutral_route(neutral, charged):
    assert pk._neutralize_v2(neutral) == pk._neutralize_v2(charged)


@pytest.mark.parametrize("smiles", ["Nc1ccccc1", "CC(=O)N", "C=N", "NO", "[NH2]", "C[N+](C)(C)C"])
def test_neutral_base_extension_excludes_other_nitrogen_chemistry(smiles):
    assert pk._neutral_base_pkas(Chem.MolFromSmiles(smiles)) == []


def test_fractional_ammonia_transform_cancels_and_reverses():
    species = {"ammonia": [-0.5, 0, "N"], "ammonium": [0.5, 1, "[NH4+]"]}
    routed, sites, nh = pk.build_ph0_reaction(species, n_Hplus=-0.5, force_base=True)
    assert nh == 0
    assert routed["ammonia"][2] == routed["ammonium"][2] == "N"
    assert sites == [["react", pk.AMMONIA_PKA, "base", 0.5], ["prod", pk.AMMONIA_PKA, "base", 0.5]]


@pytest.mark.parametrize("constant", ["AMMONIA_PKA", "AAA_AMINE_PKA", "PRIMARY_AMINE_PKA", "IMIDAZOLE_PKA", "GUANIDINIUM_PKA"])
def test_base_constants_participate_in_calibration_fingerprint(monkeypatch, constant):
    backend = SimpleNamespace(settings=Settings(), model_name="mock-uma", producer={"cacheable": True})
    before = P.effective_config(backend)["pka_constants"]
    monkeypatch.setattr(pk, constant, getattr(pk, constant) + 0.1)
    assert P.effective_config(backend)["pka_constants"] != before


def test_independent_seed_offsets_reach_sampler_and_cache(monkeypatch, species_backend):
    _, atoms=methane()
    seen=[]
    def pool(smi,q,seed,*args,**kw):
        seen.append(seed)
        return [atoms.copy()]
    monkeypatch.setattr(SP,'pool_confs',pool)
    diagnostics={}
    value,_=run_methane(Settings(conv_max=3,seed_offset=1000),diagnostics=diagnostics)
    assert value==pytest.approx(-105.)
    assert seen==[1001,1002,1003]
    assert species_backend[-1][0][3]['seed_offset']==1000
    assert diagnostics['sampling']['seed_offset']==1000
