"""Prevent stale energies and falsely identical estimator configurations (CPU only)."""

from types import SimpleNamespace

import pytest

from metag.energetics import provenance as pr
from metag.energetics import species_cache as cache


def trusted_model(tmp_path):
    checkpoint = tmp_path / "model.pt"
    checkpoint.write_bytes(b"first weights")
    identity = pr.model_identity(checkpoint, atom_refs={"H": 1.2}, form_elem_refs=None,
                                 inference_settings={"precision": "float32"})
    return SimpleNamespace(metag_provenance=identity)


@pytest.fixture
def producer_env(tmp_path, monkeypatch):
    xtb = tmp_path / "xtb"
    xtb.write_text("#!/bin/sh\nexit 0\n")
    xtb.chmod(0o755)
    monkeypatch.setattr(pr, "software_versions", lambda: {"fairchem-core": "test-v1"})
    monkeypatch.setattr(pr, "species_source_identity", lambda: "source-v1")
    monkeypatch.setattr(cache, "CACHE_DIR", str(tmp_path / "cache"))
    monkeypatch.setattr(cache, "_ENABLED", True)
    return xtb


@pytest.mark.parametrize("change", ["software", "xtb", "model", "references", "inference", "source"])
def test_changed_producer_cannot_reuse_old_energy(tmp_path, monkeypatch, producer_env, change):
    predictor = trusted_model(tmp_path)
    old = {"physics": "same", "producer": pr.producer_identity(predictor, producer_env)}
    assert old["producer"]["cacheable"]
    cache.put("CCO", 0, "implicit", old, -123, 0.5)
    assert cache.get("OCC", 0, "implicit", old) == (-123, 0.5)
    if change == "software":
        monkeypatch.setattr(pr, "software_versions", lambda: {"fairchem-core": "test-v2"})
    elif change == "xtb":
        producer_env.write_text("#!/bin/sh\nexit 1\n")
    elif change == "source":
        monkeypatch.setattr(pr, "species_source_identity", lambda: "source-v2")
    else:
        field = {"model": "checkpoint_sha256", "references": "references_sha256",
                 "inference": "inference_sha256"}[change]
        predictor.metag_provenance[field] = "f" * 64
    new = {"physics": "same", "producer": pr.producer_identity(predictor, producer_env)}
    assert cache.get("CCO", 0, "implicit", new) is None
    assert cache.get("CCO", 0, "implicit", old) == (-123, 0.5)


def test_unknown_producer_never_persists_or_reads_energies(tmp_path, producer_env):
    settings = {"producer": pr.producer_identity(SimpleNamespace(metag_model="named"), producer_env)}
    assert settings["producer"]["cacheable"] is False
    cache.put("CCO", 0, "implicit", settings, -123, 0.5)
    assert cache.get("CCO", 0, "implicit", settings) is None
    assert not (tmp_path / "cache").exists()


def test_missing_executable_disables_persistent_cache(tmp_path, producer_env):
    assert not pr.producer_identity(trusted_model(tmp_path), tmp_path / "missing")["cacheable"]


def test_additional_solvent_binary_is_part_of_producer_identity(tmp_path, producer_env):
    predictor = trusted_model(tmp_path)
    cpcmx = tmp_path/'xtbcpx'
    cpcmx.write_text('#!/bin/sh\nexit 0\n'); cpcmx.chmod(0o755)
    first = pr.producer_identity(predictor, producer_env, additional_executables={'cpcmx':cpcmx})
    assert first['cacheable']
    cpcmx.write_text('#!/bin/sh\nexit 1\n')
    second = pr.producer_identity(predictor, producer_env, additional_executables={'cpcmx':cpcmx})
    assert first != second
    cpcmx.unlink()
    assert not pr.producer_identity(predictor, producer_env,
                                    additional_executables={'cpcmx':cpcmx})['cacheable']


def test_file_identity_changes_when_replaced_and_survives_relocation(tmp_path):
    path = tmp_path / "model"
    path.write_bytes(b"weights-1")
    before = pr.file_hash(path)
    elsewhere = tmp_path / "elsewhere"
    elsewhere.write_bytes(path.read_bytes())
    assert pr.file_hash(elsewhere) == before
    path.write_bytes(b"weights-2")
    assert pr.file_hash(path) != before


def test_fairchem_reference_containers_have_content_based_identity():
    from omegaconf import OmegaConf
    plain = {'refs': {'omol': [1.0, 2.0]}, 'metadata': {'description': 'atomic references'}}
    configured = OmegaConf.create(plain)
    assert pr.content_hash(configured) == pr.content_hash(plain)
    configured.refs.omol[0] = 3.
    assert pr.content_hash(configured) != pr.content_hash(plain)


def test_reference_only_edits_reuse_raw_energies_but_change_estimator_identity(tmp_path):
    from metag.validation.provenance import source_tree_fingerprint
    root = tmp_path/'metag'
    (root/'energetics').mkdir(parents=True)
    (root/'validation').mkdir()
    for name in ('pipeline.py', 'chem.py', 'symmetry.py', 'settings.py', 'constants.py', 'timing.py', 'profile.py'):
        (root/name).write_text('# numerical dependency\n')
    (root/'energetics/thermal.py').write_text('temperature = 298.15\n')
    reference = root/'validation/chemical_reference.py'
    reference.write_text('policy = 1\n')
    species_before = pr.species_source_fingerprint(root)
    estimator_before = source_tree_fingerprint(root)
    reference.write_text('policy = 2\n')
    assert pr.species_source_fingerprint(root) == species_before
    assert source_tree_fingerprint(root) != estimator_before
    (root/'constants.py').write_text('T = 310.15\n')
    assert pr.species_source_fingerprint(root) != species_before
    (root/'constants.py').write_text('# numerical dependency\n')
    assert pr.species_source_fingerprint(root) == species_before
    (root/'energetics/thermal.py').write_text('temperature = 310.15\n')
    assert pr.species_source_fingerprint(root) != species_before


class _Backend:
    model_name = "fake-uma"
    producer = {"cacheable": False}

    def __init__(self, settings=None):
        from metag.settings import Settings
        self.settings = settings or Settings()

    def species(self, *args):
        return -50, 0.2

    def liquid_water(self):
        return -10.


def test_per_call_policy_has_distinct_identities():
    from metag import pipeline as p
    from metag.settings import Settings
    backend = _Backend()
    automatic = p.effective_config(backend)
    full = p.effective_config(backend, allow_truncate=False)
    neutral = p.effective_config(backend, force_neutral=True)
    assert automatic != full and automatic != neutral
    assert automatic == p.effective_config(_Backend())
    assert automatic != p.effective_config(_Backend(Settings(seed_offset=1)))


@pytest.mark.parametrize("invalid", [False, True])
def test_result_retains_call_policy_even_on_failed_estimate(invalid):
    from metag import pipeline as p
    backend = _Backend()
    reaction = {"species": {"ethanol": [-1, 0, "CCO"], "ether": [1, 0, "COC"]}}
    if invalid:
        reaction["n_Hplus"] = 1
    result = p.score_reaction(reaction, backend, allow_truncate=False, log=lambda *_: None)
    assert result["config"] == p.effective_config(backend, allow_truncate=False)
    assert (result["dG"] is None) == invalid


def test_equivalent_smiles_produce_identical_fresh_sampling_inputs(monkeypatch):
    import numpy as np
    from metag.energetics import species as sp
    from metag.energetics.conformers import pool_confs
    from metag.settings import Settings
    class SamplingObserved(Exception):
        pass
    observed = []
    def sample(smiles, q, seed, pool, spin):
        atoms = pool_confs(smiles, q, seed, pool, spin=spin)
        observed.append((smiles, [(a.get_chemical_symbols(), a.positions.copy()) for a in atoms]))
        raise SamplingObserved
    monkeypatch.setattr(sp.cache, 'get', lambda *a, **kw: None)
    monkeypatch.setattr(sp, 'sampling_budget', lambda smi, **kw: (1, 1, 3))
    monkeypatch.setattr(sp, 'pool_confs', sample)
    for smi in ('CCO', 'OCC'):
        with pytest.raises(SamplingObserved):
            sp.species_free_energy(None, 0, smi, 'ethanol', settings=Settings(), model_name='test',
                                   producer={'cacheable': False}, log=lambda *_: None)
    assert observed[0][0] == observed[1][0]
    assert len(observed[0][1]) == len(observed[1][1])
    for (sa, xa), (sb, xb) in zip(observed[0][1], observed[1][1]):
        assert sa == sb
        np.testing.assert_array_equal(xa, xb)


def test_loaded_module_versions_are_distinct_from_distribution_metadata(monkeypatch):
    import sys
    from types import SimpleNamespace
    from metag.energetics import provenance as pr
    monkeypatch.setitem(sys.modules,'numpy',SimpleNamespace(__version__='actual-numerical-build'))
    monkeypatch.setattr(pr,'software_versions',lambda:{'numpy':'metadata-build'})
    identity=pr.runtime_identity('unavailable-test-xtb')
    assert identity['software']['numpy']=='metadata-build'
    assert identity['imported_module_versions']['numpy']=='actual-numerical-build'
