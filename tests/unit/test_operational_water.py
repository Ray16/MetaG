from types import SimpleNamespace

import numpy as np
import pytest

from benchmark.operational_records import authenticate
from benchmark.operational_water import water_record
from metag.energetics.backend import UMABackend
from metag.energetics import species, thermal, uma
from metag.settings import Settings


def fixture(monkeypatch, imaginary=0):
    symbols = ['O', 'H', 'H']
    coords = np.array([[0., 0., 0.], [.757, .586, 0.], [-.757, .586, 0.]])
    monkeypatch.setattr(species, 'water_geometry', lambda _: (symbols, coords))
    monkeypatch.setattr(uma, 'batched_energies', lambda _, atoms: np.array([-10.]))

    def correction(*args, **kwargs):
        assert kwargs['hessian'] == 'autograd'
        if kwargs.get('return_info'):
            return -3., dict(n_imag=imaginary, thermodynamics=dict(
                temperature_K=298.15, pressure_Pa=101325.,
                enthalpy_correction_kj=2., minus_T_entropy_kj=-5.))
        return -3.

    monkeypatch.setattr(thermal, 'uma_gibbs_corr', correction)
    model = {k: '0' * 64 for k in ['checkpoint_sha256', 'references_sha256', 'inference_sha256']}
    return SimpleNamespace(pu=object(), settings=Settings(hessian='autograd'),
                           producer=dict(cacheable=False, model=model),
                           _water=None, log=lambda _: None)


def test_recording_adapter_executes_water_imports_and_matches_backend_recipe(monkeypatch, tmp_path):
    backend = fixture(monkeypatch)
    record = water_record(backend)
    assert record['G'] == UMABackend.liquid_water(backend)
    assert record['symmetrynumber'] == 2 and record['geometry_kind'] == 'nonlinear'
    # Exercise the real raw-record authenticator on the adapter's emitted schema.
    import hashlib
    import json
    task = dict(slot='water', kind='liquid_water', source_input_key='liquid_water',
                charge=0, spin=1, n_atoms=3)
    path = tmp_path / 'cache/water/water_reference.json'
    path.parent.mkdir(parents=True)
    path.write_text(json.dumps(record))
    result = dict(status='molecular_record_complete', energy_kj=record['G'],
                  sampling_sigma_kj=0., producer=record['producer'], cache_path=str(path),
                  cache_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                  controller_execution=dict(exit_code=0, timed_out=False))
    _, receipt = authenticate(task, result, tmp_path, {}, backend.producer['model'])
    assert receipt['n_validated_minima'] == 1


def test_imaginary_water_reference_is_refused(monkeypatch):
    backend = fixture(monkeypatch, imaginary=1)
    with pytest.raises(ValueError, match='imaginary internal modes'):
        water_record(backend)
