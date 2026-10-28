"""Thermodynamic conditions and physical constants shared by every MetaG number.

Target: 298.15 K, pH 7, zero ionic strength, no Mg2+ binding, 1 M solutes, pure liquid water. Anything
compared against MetaG -- experiment or another estimator -- must first be transformed to these
conditions. The rounded constants (2.303 for ln 10, R = 8.314e-3 kJ/mol/K) are the deployed
convention of the evaluated release and are kept so results stay reproducible.
"""
import math

T = 298.15                      # K
PH = 7.0
R_KJ = 8.314e-3                 # kJ/(mol K), deployed rounding
RT_LN10 = 2.303 * R_KJ * T      # kJ/mol per pH unit (deployed rounding of ln 10)
EV2KJ = 96.485                  # kJ/mol per eV

CONDITIONS = {
    "T_K": T, "pH": PH, "ionic_strength_M": 0.0, "pMg": None,
    "water_activity": 1.0,
    "standard_state": "solutes 1 M; pure liquid water; H+ at the stated pH (transformed ΔrG'°)",
}

# Aqueous proton at pH 7 (computational hydrogen electrode): G(H+, gas) + Tissandier 1 atm(gas) -> 1 M(aq)
# solvation (-1104.5 kJ/mol = -264.0 kcal/mol) + the pH term. Already on the 1 M scale.
G_HPLUS = -26.3 - 1104.5 - 2.303 * R_KJ * T * 7.0          # ~ -1170.8 kJ/mol

# Gas RRHO terms are at 1 atm while xtb ΔG_solv is 1 M(gas) -> 1 M(aq) (Ben-Naim). Every solute needs
# +RT ln(RT/P°V°) = RT ln 24.46 = +7.93 kJ/mol to reach 1 M; without it a reaction inherits -7.93·Δn.
STD_STATE_KJ = R_KJ * T * math.log(0.082057 * T)

# Liquid water: experimental hydration free energy (-6.3 kcal/mol) and the 1 M -> 55.34 M pure-liquid term.
WATER_DGSOLV_KJ = -26.4
LIQUID_WATER_CONC = 55.34
