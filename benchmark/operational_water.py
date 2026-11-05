"""Retain the components of the existing backend's liquid-water reference."""
from dataclasses import asdict
import hashlib
import json
import math

import numpy as np
from ase import Atoms


def water_record(backend):
    """Same geometry, gas RRHO, hydration and concentration recipe as UMABackend.

    This recording adapter adds no solute-standard-state term. Reaction assembly
    applies its existing conversion. It supplies no experimental validation.
    """
    from metag.constants import EV2KJ, R_KJ, LIQUID_WATER_CONC, WATER_DGSOLV_KJ
    from metag.energetics.species import water_geometry
    from metag.energetics.thermal import uma_gibbs_corr
    from metag.energetics.uma import batched_energies
    from metag.symmetry import geometry_and_sigma

    symbols, coords = water_geometry(backend.pu)
    atoms = Atoms(symbols=list(symbols), positions=coords, info=dict(charge=0, spin=1))
    electronic = float(batched_energies(backend.pu, [atoms])[0]) * EV2KJ
    correction, info = uma_gibbs_corr(backend.pu, list(symbols), coords, 0,
                                    hessian=backend.settings.hessian, return_info=True)
    if info['n_imag'] != 0:
        raise ValueError('Water reference has imaginary internal modes')
    geometry = dict(symbols=list(symbols), positions_angstrom=np.asarray(coords).tolist(),
                    charge=0, spin_multiplicity=1)
    geometry['sha256'] = hashlib.sha256(
        json.dumps(geometry, sort_keys=True, allow_nan=False).encode()).hexdigest()
    kind, symmetry, _ = geometry_and_sigma(list(symbols), coords)
    concentration = float(R_KJ * 298.15 * np.log(LIQUID_WATER_CONC))
    components = [electronic, WATER_DGSOLV_KJ, correction, concentration]
    if not all(math.isfinite(x) for x in components):
        raise ValueError('Nonfinite water-reference component')
    return dict(kind='liquid_water_reference', G=float(sum(components)), sigma=0.,
                producer=backend.producer, settings=asdict(backend.settings), geometry=geometry,
                electronic_kj=electronic, experimental_water_hydration_kj=WATER_DGSOLV_KJ,
                thermal_kj=correction, concentration_kj=concentration,
                thermodynamics=info['thermodynamics'], n_imag=info['n_imag'],
                geometry_kind=kind, symmetrynumber=symmetry)
