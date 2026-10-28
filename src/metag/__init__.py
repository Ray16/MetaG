"""MetaG: standard transformed reaction Gibbs energies (ΔrG'°) of biochemical reactions from structure.

ΔrG'° (298.15 K, pH 7, I = 0, 1 M) = Σν·G_species + n_H+·G(H+, pH 7) + Alberty pKa transforms, with
G_species from UMA electronic energies, per-conformer RRHO thermochemistry and GFN2-xTB/ALPB solvation,
after structure-based routing (cofactor cores, spectator truncation, pH-0 neutral reference species).
No parameter fitted to reaction free energies enters an estimate.

CPU (routing, protonation, preflight)::

    from metag.routing.planner import route_reaction

GPU scoring (needs the ``qm`` extra, UMA weights and an xtb executable)::

    from metag.energetics.uma import load_uma
    from metag.pipeline import score_reaction
    result = score_reaction(reaction, model=load_uma())
"""
__version__ = "0.3.0"
