"""Structure-based reaction routing (pure CPU logic; no QM, no reaction-energy data).

    planner          -- stage order and balance guards (entry point: route_reaction)
    cofactor_cores   -- isodesmic NAD(P) nicotinamide and glutathione thiol cores
    ntp_core         -- methyl-capped spectator nucleoside polyphosphates
    truncate         -- spectator truncation around the reaction centre (C-C cuts only)
    truncation_gate  -- full-molecule gate for compact ring-forming/opening reactions
    pka_transform    -- pH-0 neutral reference species and Alberty pKa transforms
    pka_domain       -- supported phosphoryl environments (unsupported ones fail explicitly)
    aldehyde_hydration -- carbonyl hydration microstate enumeration and K_hyd calibration
    sugar_states     -- declared reducing-sugar solution pools
"""
