"""Species-energy backend: binds a loaded UMA model, the runtime settings and the producer identity.

:func:`metag.pipeline.score_reaction` needs only two operations from a backend,

    species(q, smiles, label, warnings, diagnostics) -> (G, sigma)    # kJ/mol, 1 atm gas standard state
    liquid_water() -> G                                                # one liquid-water molecule

plus ``model_name`` and ``producer`` attributes for the configuration fingerprint. Any object with this
interface (for example a replay of recorded energies) can drive the reaction-level pipeline on CPU.
"""

import numpy as np
from ase import Atoms

from metag.constants import EV2KJ, R_KJ, LIQUID_WATER_CONC, WATER_DGSOLV_KJ
from metag.settings import Settings


class UMABackend:
    """UMA electronic energies + per-minimum RRHO + GFN2-xTB/ALPB solvation."""

    def __init__(self, pu, settings=None, log=print):
        from metag.energetics import thermal
        from metag.energetics.provenance import producer_identity
        self.pu = pu
        self.settings = Settings.from_env() if settings is None else settings
        loaded = getattr(pu, "metag_model", None)
        if loaded is not None and loaded != self.settings.model:
            # The model name keys the species cache: scoring with one model under another's key poisons it.
            raise ValueError(f"loaded UMA model {loaded!r} != configured model {self.settings.model!r}")
        self.model_name = loaded or self.settings.model
        self.producer = producer_identity(pu, thermal.XTB)
        self.log = log
        self._water = None

    def species(self, q, smiles, label, warnings, diagnostics):
        from metag.energetics.species import species_free_energy
        return species_free_energy(self.pu, q, smiles, label, settings=self.settings,
                                   model_name=self.model_name, producer=self.producer, log=self.log,
                                   warnings=warnings, diagnostics=diagnostics)

    def liquid_water(self):
        """Free energy of ONE liquid-water molecule:

            G*_liq(H2O) = E_UMA(H2O) + ΔG_hyd,exp(H2O) + G_corr(H2O, 1 atm) + RT ln 55.34

        The experimental hydration free energy replaces a continuum model of water in water. Like every
        species value this is on the 1 atm gas scale; the caller adds the 1 atm -> 1 M term, after which the
        RT ln 55.34 (1 M -> pure liquid) term is correct."""
        if self._water is not None and self.producer["cacheable"]:
            return self._water
        from metag.energetics.species import water_geometry
        from metag.energetics.thermal import uma_gibbs_corr
        from metag.energetics.uma import batched_energies
        sym, coord = water_geometry(self.pu)
        atoms = Atoms(symbols=list(sym), positions=coord, info={"charge": 0, "spin": 1})
        E = float(batched_energies(self.pu, [atoms])[0]) * EV2KJ
        thermal = uma_gibbs_corr(self.pu, list(sym), coord, 0, hessian=self.settings.hessian)
        conc = R_KJ * 298.15 * float(np.log(LIQUID_WATER_CONC))     # +9.96 kJ/mol, 1 M -> liquid 55.34 M
        G = E + WATER_DGSOLV_KJ + thermal + conc
        if self.producer["cacheable"]:
            self._water = G
        self.log(f"    [water ref] G*_liq(H2O) = E {E:.1f} + solv {WATER_DGSOLV_KJ:.1f} + thermal {thermal:.1f} "
                 f"+ conc {conc:.1f} = {G:.1f}")
        return G

