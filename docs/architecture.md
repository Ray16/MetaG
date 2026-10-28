# Architecture

```
reaction dict
  │  reactions.validate_reaction        shape, charges vs SMILES, single connected species
  ▼
routing.planner.route_reaction          pure CPU; each step balance-guarded (reverted if it breaks balance)
  │   cofactor_cores   NAD(P) -> 1-methylnicotinamide (+dihydro); GSH/GSSG -> capped cysteine (isodesmic)
  │   ntp_core         spectator nucleoside polyphosphate -> methyl polyphosphate
  │   truncation_gate  compact ring formation/opening -> keep full molecules
  │   truncate         spectator truncation (C-C cuts only; radius 2, radius 3 at anomeric centres)
  │   pka_transform    pH-0: neutral reference species + Alberty pKa sites; bases when the gate passes
  │   zwitterion guard remaining zwitterions -> neutral reference
  ▼
protonation_registry                    source-backed coupled pKa ladders replace independent-site terms
state_pools.prepare_sugar_pools         declared reducing-sugar pools (CPU validation)
  ▼
backend.species / backend.liquid_water  (energetics.backend.UMABackend; any object with this interface)
  │   energetics.species   conformer search -> UMA FIRE relaxation -> xtb ALPB -> strict per-minimum RRHO
  │   energetics.species_cache   content-addressed cache keyed by structure, settings and producer identity
  ▼
energetics.aqueous_states.hydrate_all_sites   carbonyl/gem-diol microstates (K_hyd-calibrated)
state_pools.apply_sugar_pools                 sugar pool mixtures
  ▼
assembly.assemble_reaction              Σν·G + n_H+·G(H+) + pKa transforms; shared-estimate sampling σ
  ▼
pipeline.score_reaction                 guards (charge closure, finiteness) -> result dict (dG or None)
```

**Dependency direction.** `routing/` never imports the QM backend. `energetics/species.py` never imports
the pipeline. The pipeline talks to species physics only through the backend interface
(`species`, `liquid_water`, `model_name`, `producer`). Replaying recorded energies through a fake backend
therefore exercises every reaction-level code path on CPU.

**Configuration.** The scientific model is fixed in code. `metag.settings.Settings` holds only search
effort, the model checkpoint and the Hessian implementation. It is resolved once per backend, and the
same object both drives the computation and builds the species-cache key. Each result's `config`
records the settings, constant-table hashes, source hash, producer identity and software versions.

**Units.** UMA returns eV, converted with `EV2KJ = 96.485`. xtb returns Hartree, converted with
2625.499639. Everything downstream is kJ/mol. Species values are cached on the 1 atm gas standard state.
The pipeline adds +RT ln 24.46 per solute (1 M), and liquid water carries RT ln 55.34.
