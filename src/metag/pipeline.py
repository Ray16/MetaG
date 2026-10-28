"""MetaG scoring pipeline: ``score_reaction`` -- one scheme for every reaction class.

Per reaction:
  1. Structural routing (``metag.routing.planner``): NAD(P)/GSH ring cores, nucleotide caps, spectator
     truncation (radius 3 at anomeric centres), pH-0 neutral reference species with an exact
     Alberty pKa transform, zwitterion guard. Every rewrite is element/charge-balance guarded.
  2. Species free energies from the backend (``metag.energetics.backend.UMABackend``): Boltzmann over
     validated minima of E_UMA + ΔG_solv(xTB-ALPB) + RRHO(UMA Hessian), + 1 atm -> 1 M; liquid water
     from its experimental hydration free energy.
  3. Aqueous states: carbonyl <-> gem-diol hydration microstates (K_hyd-calibrated) and declared
     reducing-sugar pools.
  4. ΔG = Σ ν·G + n_H+·G(H+, pH 7) + Σ pKa-transform terms; charge-closure and finiteness guards fail
     closed (dG = None).

No parameter fitted to reaction free energies enters the estimate.
"""
import functools
import hashlib
import json
import math

from metag.assembly import assemble_reaction
from metag.chem import fmt_num, is_water
from metag.constants import CONDITIONS, G_HPLUS, PH, RT_LN10, STD_STATE_KJ
from metag.energetics.errors import SpeciesRearranged
from metag.reactions import validate_reaction
from metag.settings import PHYSICS_VERSION, SAMPLING_PROFILE_VERSION


def score_reaction(reaction, backend=None, *, model=None, key="rxn", log=print, allow_truncate=True):
    """Standard transformed reaction Gibbs energy ΔrG'° (kJ/mol, 298.15 K, pH 7, I = 0) from structure.

    `reaction` = {"species": {name: [coeff, charge, smiles]}, "n_Hplus": float, "note": str, and
    optionally "exp": [float] (reported as error only), "pka_sites", "species_state_pools",
    "species_protonation_models"}. Products have positive coefficients.

    Pass either a species-energy `backend` or a loaded UMA predict unit as `model` (wrapped in
    :class:`~metag.energetics.backend.UMABackend` with settings from the environment). Returns the result
    dict; ``dG`` is None when no valid estimate exists (the reason is in ``suspect``). With the official UMA
    calculator, ``backend_execution`` records the actual Nitro/FP32 calls made by this invocation."""
    if backend is None:
        if model is None:
            raise TypeError("score_reaction needs a backend or a loaded UMA model")
        from metag.energetics.backend import UMABackend
        backend = UMABackend(model, log=log)
    pu = getattr(backend, "pu", None)
    traced = getattr(pu, "_metag_official_calculator", False) is True
    before = pu.execution_snapshot() if traced else None
    result = _score_reaction(reaction, backend, key=key, log=log, allow_truncate=allow_truncate)
    if traced:
        result["backend_execution"] = pu.execution_since(before)
    return result


def _score_reaction(reaction, backend, *, key, log, allow_truncate, _force_neutral=False):
    validate_reaction(reaction)
    cfg = effective_config(backend, allow_truncate=allow_truncate, force_neutral=_force_neutral)
    if reaction.get('species_state_pools'):
        from metag.energetics.provenance import content_hash
        from metag.routing.sugar_states import POLICY
        cfg['declared_sugar_pool_policy'] = POLICY
        cfg['caller_species_state_pools_sha256'] = content_hash(reaction['species_state_pools'])
    no_estimate = functools.partial(_no_estimate, config=cfg)
    orig_species = {k: list(v) for k, v in reaction["species"].items()}   # pre-routing identity
    from metag.chem import isotope_labels, is_balanced
    isotopic = {name: labels for name, (_, _, smi) in orig_species.items() if (labels := isotope_labels(smi))}
    if isotopic:
        reason = "isotope-specific thermochemistry is unsupported; mass and solvent/reservoir conventions are absent"
        return no_estimate(key, reaction, {"errors": [], "warnings": [reason], "truncated": False,
                           "routing_status": "not_run_unsupported_input",
                           "unsupported_isotopes": isotopic,
                           "input_balanced": is_balanced(orig_species, reaction.get("n_Hplus", 0))}, reason, log)
    from metag.routing.planner import route_reaction
    rx, routes, truncated = route_reaction(reaction, allow_truncate=allow_truncate, log=log,
                                           force_neutral=_force_neutral)
    routes["rearrangement_reroute"] = bool(_force_neutral)
    if routes["errors"]:
        # An unexpected exception inside a routing transform would silently score a different estimator
        # than the configured one -> fail closed, before spending any QM.
        return no_estimate(key, reaction, routes, "routing error: " + "; ".join(routes["errors"]), log)
    if not routes["input_balanced"]:
        # Elemental or charge imbalance makes a reaction free energy undefined.
        return no_estimate(key, reaction, routes, routes["warnings"][-1], log)
    from metag.routing.sugar_states import pool_aware_stereochemistry
    stereo = pool_aware_stereochemistry(orig_species, rx['species'], routes.get('declared_sugar_pools', {}))
    if stereo.get('covered_by_declared_sugar_pools'):
        routes['stereochemistry'] = stereo
    if stereo["input"] or stereo["routed"]:
        routes["stereochemistry"] = stereo
        reason = ("unresolved stereochemical identity: specify each stereogenic atom/bond; "
                  "an unspecified SMILES is not a declared equilibrating stereoisomer pool")
        routes["warnings"].append(reason)
        return no_estimate(key, reaction, routes, reason, log)

    # Source-backed coupled protonation ladders replace the routed species' independent-site terms.
    from metag.protonation_registry import prepare_reaction_protonation
    from metag.energetics.protonation_models import POLICY_VERSION
    try:
        effective_pka_sites, coupled_protonation, model_diagnostics, registry_report, _ = prepare_reaction_protonation(
            reaction, rx['species'], rx.get('pka_sites', []),
            routes.get('protonation_sites_by_species'), CONDITIONS,
            enabled=True, rt_ln10=RT_LN10,
            expected_registry_sha256=cfg.get('source_protonation_registry_sha256'))
    except (ValueError, TypeError, KeyError) as exc:
        reason = 'unsupported source-backed protonation model: ' + str(exc)
        routes['protonation_model_error'] = reason
        return no_estimate(key, reaction, routes, reason, log)
    routes['source_protonation_registry'] = registry_report
    if model_diagnostics is not None:
        routes['coupled_protonation_models'] = model_diagnostics
        cfg['species_protonation_model_policy'] = POLICY_VERSION
        if reaction.get('species_protonation_models'):
            cfg['species_protonation_models_sha256'] = model_diagnostics['models_sha256']
    if reaction.get('pka_sites'):
        from metag.energetics.provenance import content_hash
        cfg['caller_pka_sites_sha256'] = content_hash(reaction['pka_sites'])
    sugar_plans = {}
    if routes.get('declared_sugar_pools'):
        from metag.state_pools import prepare_sugar_pools
        try:
            sugar_plans = prepare_sugar_pools(rx, routes, caller_pka_sites=bool(reaction.get('pka_sites')))
        except (ValueError, TypeError, KeyError) as exc:
            reason = 'unsupported declared sugar ensemble: ' + str(exc)
            routes['sugar_pool_error'] = reason
            return no_estimate(key, reaction, routes, reason, log)
        routes['sugar_state_plan'] = sugar_plans
    log(f"\n=== {key}: {rx['note']}  (n_H+={rx['n_Hplus']}) ===")

    G, sig = {}, {}
    std = STD_STATE_KJ
    try:
        _species_free_energies(backend, rx, routes, G, sig, std, log)
    except SpeciesRearranged as e:
        # Sampled gas relaxations could not retain the requested connectivity. A neutral reference may be
        # retried once, and only for a remaining zwitterion; its second failure is terminal.
        if _force_neutral:
            return no_estimate(key, reaction, routes,
                               f"species identity unresolved on neutral-reference route: {e}", log)
        from metag.routing.pka_transform import has_zwitterion
        permitted = has_zwitterion(rx["species"])
        attempt = {"reason": str(e), "requested_species": dict(rx["species"]),
                   "neutral_reference_permitted": permitted,
                   "requested_pka_sites": rx.get("pka_sites", [])}
        routes["species_rearrangement"] = attempt
        if not permitted or rx.get("pka_sites"):
            return no_estimate(key, reaction, routes,
                               f"species identity unresolved; neutral-reference fallback disabled "
                               f"or incompatible with supplied pKa sites: {e}", log)
        log(f"  [{e} -> re-route through the neutral-microspecies path]")
        result = _score_reaction(reaction, backend, key=key, log=log, allow_truncate=allow_truncate,
                                 _force_neutral=True)
        result["routes"]["species_rearrangement"] = attempt
        return result
    if any(v is None for v in G.values()) or len(G) < len(rx["species"]):
        if truncated:                                        # a truncated core failed -> retry FULL molecules
            log(f"  [retry {key} with truncation OFF (full molecules)]")
            return _score_reaction(reaction, backend, key=key, log=log, allow_truncate=False,
                                   _force_neutral=_force_neutral)
        failed = [name for name in rx["species"] if G.get(name) is None]
        routes["failed_species"] = failed
        return no_estimate(key, reaction, routes,
                           "species energy or minimum validation failed: " + ", ".join(failed), log)

    # Carbonyl hydration: fold each carbonyl's aqueous carbonyl <-> gem-diol microstates into its G.
    from metag.energetics.aqueous_states import hydrate_all_sites
    for name, (coeff, q, smi) in list(rx["species"].items()):
        if G.get(name) is not None:
            hydrate_all_sites(backend, rx, name, q, smi, G, sig, std, log, routes, acid_transform=_acid_transform)
    if routes.get("incomplete_hydration"):
        return no_estimate(key, reaction, routes,
                           "incomplete hydration ensemble: " + ", ".join(sorted(set(routes["incomplete_hydration"]))), log)
    if sugar_plans:
        from metag.state_pools import apply_sugar_pools

        def evaluate_sugar_state(name, charge, smiles):
            detail = {'warnings': [], 'species': {}}
            routes.setdefault('sugar_state_attempts', {}).setdefault(name, []).append(
                {'smiles': smiles, 'diagnostics': detail})
            try:
                value, sigma = backend.species(charge, smiles, name + '(sugar-state)', detail['warnings'],
                                               detail['species'])
            except SpeciesRearranged as exc:
                raise ValueError('Sugar state identity unresolved: ' + str(exc)) from exc
            if value is None or sigma is None:
                return value, sigma, detail
            alt_G, alt_sig = {name: value + std}, {name: sigma}
            alt_routes = {'warnings': detail['warnings']}
            alt_rx = dict(rx, species={name: [rx['species'][name][0], charge, smiles]})
            hydrate_all_sites(backend, alt_rx, name, charge, smiles, alt_G, alt_sig, std, log, alt_routes,
                              acid_transform=_acid_transform)
            detail['aqueous_states'] = alt_routes
            routes['warnings'].extend(detail['warnings'])
            if alt_routes.get('incomplete_hydration'):
                raise ValueError('Incomplete joint sugar/hydration state set: ' + name)
            return alt_G[name], alt_sig[name], detail
        try:
            G, sig, sugar_diagnostics = apply_sugar_pools(
                sugar_plans, rx['species'], G, sig, evaluate_state=evaluate_sugar_state)
        except (ValueError, TypeError, KeyError) as exc:
            reason = 'incomplete declared sugar ensemble: ' + str(exc)
            routes['sugar_pool_error'] = reason
            return no_estimate(key, reaction, routes, reason, log)
        routes['sugar_ensembles'] = sugar_diagnostics

    energy = _assemble(rx, routes, G, sig, effective_pka_sites, coupled_protonation)
    dG, U_samp, pka_total = energy.total, energy.sampling_sigma, energy.protonation
    if rx.get("pka_sites"):
        _protonation_diagnostics(rx, routes, effective_pka_sites, coupled_protonation, model_diagnostics,
                                 pka_total)
        model_label = ('mixed independent and macroscopic models' if routes.get('coupled_protonation_models')
                       else 'independent sites')
        log(f"    pKa transform ({model_label}, {len(effective_pka_sites)} independent terms) += {pka_total:+.1f} kJ/mol")
    stages = {"species_sum": energy.species_sum, "proton": energy.proton}

    # Guards (fail closed). (1) charge must close: Σ coeff·q + n_H+ == 0; a broken transform leaks
    # ~n·G(H+) ≈ ±1170 kJ per unbalanced proton. (2) a NaN/inf is never a result.
    suspect = None
    q_imbalance = sum(coeff * q for name, (coeff, q, smi) in rx["species"].items()) + rx["n_Hplus"]
    if abs(q_imbalance) > 1e-9:
        suspect = f"charge imbalance {fmt_num(q_imbalance)} (Σcoeff·q + n_H+ ≠ 0)"
        log(f"  !! SUSPECT: {suspect} -> ΔG leaks ~{fmt_num(q_imbalance)}·G(H+); result unreliable, do not trust.")
    if not math.isfinite(dG):
        m = "non-finite ΔG"
        suspect = m if suspect is None else f"{suspect}; {m}"
    # (3) magnitude diagnostic -- a warning, not a rejection: ΔG is extensive, so the bound scales with the
    # reaction extent (largest non-water |coefficient|).
    extent = max([abs(c) for (c, q, s_) in rx["species"].values() if not is_water(s_, q)] or [1.0])
    sanity = cfg["settings"]["dg_sanity_kj"]
    if math.isfinite(dG) and abs(dG) > sanity * extent:
        m = f"|ΔG|={abs(dG):.0f} kJ > {sanity:.0f} kJ x extent {fmt_num(extent)}: check (not rejected)"
        routes["warnings"].append(m)
        log(f"  !! {m}")
    # dG_raw = the estimate of the ROUTED reaction (cores, truncation and neutral microspecies already
    # applied). No benchmark-fitted reaction-energy correction follows, so dG == dG_raw.
    dG_raw = dG
    stages["pka_transform"] = pka_total
    stages["unanchored"] = dG_raw
    stages["reported"] = dG
    stages = {k: round(v, 2) for k, v in stages.items()}
    common = dict(reaction=key, stages=stages,
                  original_species={name: list(value) for name, value in orig_species.items()},
                  energy_units="kJ/mol", quantity="standard_transformed_reaction_gibbs", reported_rounding_kj=0.1,
                  conditions=CONDITIONS, U_samp=round(U_samp, 1), exp=rx.get("exp"), suspect=suspect,
                  experimental_reference={k: v for k, v in rx.items() if k.startswith("exp_")},
                  note=rx["note"], routes=routes, std_state_kJ=round(std, 2), config=cfg,
                  species_scored={k: list(v) for k, v in rx["species"].items()})
    if suspect is not None:
        log(f"  ΔG = NONE (suspect: {suspect}); raw value {dG_raw:+.1f} kept as dG_raw_unreliable")
        return dict(common, dG=None, dG_raw=None, dG_raw_unreliable=round(dG_raw, 1), err=[])
    errs = [dG - e for e in rx.get("exp", [])]              # experiment is optional
    log(f"  ΔG = {dG:+.1f} kJ/mol (sampling spread {U_samp:.1f})   vs exp {rx.get('exp')}   "
        f"err {[round(e, 1) for e in errs]}")
    return dict(common, dG=round(dG, 1), dG_raw=round(dG_raw, 1), err=[round(e, 1) for e in errs])


def _species_free_energies(backend, rx, routes, G, sig, std, log):
    """Fill G/sig (kJ/mol, 1 M standard state) for every routed species; stops at the first species that
    fails (G[name] = None). Raises SpeciesRearranged (see energetics.species)."""
    species_info = routes.setdefault("species_diagnostics", {})
    for name, (coeff, q, smi) in rx["species"].items():
        diagnostic = species_info[name] = {"sampling": {"status": "unavailable",
                                                        "reason": "backend did not provide diagnostics"}}
        if is_water(smi, q):                              # liquid water (hydrolysis / hydration)
            diagnostic.update(method="liquid_water_reference", sampling={"status": "not_applicable"})
            G[name] = backend.liquid_water(); sig[name] = 0.0     # deterministic reference: no sampling σ
            log(f"    {name:9s} q+0 [water ref liquid]: {G[name]:.1f}")
        else:
            G[name], sig[name] = backend.species(q, smi, name, routes["warnings"], diagnostic)
        if G[name] is None:
            log(f"    {name}: FAILED")
            return
        G[name] += std                                    # 1 atm -> 1 M, applied outside the species cache
        _release_gpu_memory()


def _release_gpu_memory():
    """Return cached GPU blocks so a later large species does not OOM from fragmentation."""
    import sys
    torch = sys.modules.get("torch")
    if torch is not None and torch.cuda.is_available():
        torch.cuda.empty_cache()


def _assemble(rx, routes, G, sig, effective_pka_sites, coupled_protonation):
    from metag.energetics import species_cache
    sampling_groups = {}
    for name, (_, q, smi) in rx["species"].items():
        if name in routes.get('sugar_ensembles', {}):
            sampling_groups[name] = ('sugar_ensemble', q, routes['sugar_ensembles'][name]['state_set_sha256'])
            continue
        method = routes.get("species_diagnostics", {}).get(name, {}).get("method")
        # Unknown callback provenance remains independent; drawings that merely match share nothing.
        sampling_groups[name] = (method, q, species_cache.canonical(smi)) if method else ("unverified", name)
    routes["sampling_uncertainty"] = {
        "policy": "shared_estimates_combined_before_quadrature",
        "n_scored_labels": len(sampling_groups), "n_sampling_groups": len(set(sampling_groups.values())),
        "scope": "sampling only; distinct species/method estimates independent; no fitted network covariance"}
    return assemble_reaction(rx["species"], G, sig, n_hplus=rx["n_Hplus"],
                             proton_energy=G_HPLUS, pka_sites=effective_pka_sites,
                             ph=PH, rt_ln10=RT_LN10, sampling_groups=sampling_groups,
                             coupled_protonation=coupled_protonation,
                             sampling_offsets={name: row['reference_owned_factor_kj']
                                               for name, row in routes.get('sugar_ensembles', {}).items()})


def _protonation_diagnostics(rx, routes, effective_pka_sites, coupled_protonation, model_diagnostics, pka_total):
    from metag.energetics.protonation import reaction_transform_diagnostics
    from metag.validation.protonation import inspect_species
    routes["protonation_diagnostics"] = reaction_transform_diagnostics(effective_pka_sites, PH, RT_LN10)
    if routes.get('coupled_protonation_models'):
        diagnostic = routes['protonation_diagnostics']
        diagnostic.update(model='mixed_independent_and_sequential_macroscopic',
                          independent_correction_kj=diagnostic['correction_kj'],
                          coupled_correction_kj=coupled_protonation, correction_kj=pka_total)
        diagnostic['dG_d_pH_kj_per_unit'] += model_diagnostics['reaction_dG_d_pH_kj_per_unit']
    coverage = {}
    for name, (_, _, smi) in rx["species"].items():
        if name in (routes.get('coupled_protonation_models') or {}).get('species', {}):
            coverage[name] = {'status': 'source_backed_macroscopic_reference',
                              'issues': ['computed_reference_ensemble_and_source_transfer_require_validation'],
                              'scope': 'Whole fully protonated family; missing-state and solvent/thermal errors remain.'}
            continue
        try:
            report = inspect_species(smi)
        except (ValueError, RuntimeError) as exc:
            report = {"status": "unsupported", "issues": ["structural_diagnostic_unavailable"], "error": str(exc)}
        if report["issues"] or report["status"] == "unsupported":
            coverage[name] = report
    routes["protonation_diagnostics"]["structural_coverage"] = coverage


def _acid_transform(smi_neutral):
    """-RT ln Π(1+10^(pH-pKa)) over the acid sites of a NEUTRAL pH-0 species (its max-anion form
    classified with the pKa table), i.e. G'(pH) - G(neutral) for the acid groups."""
    from metag.energetics.protonation import site_free_energy
    from metag.routing import pka_transform as pk
    _, pkas, _ = pk._neutralize(pk._canonicalize_maxanion(smi_neutral))
    return math.fsum(site_free_energy(p, "acid", PH, RT_LN10) for p in (pkas or []))


def _no_estimate(key, reaction, routes, reason, log, *, config=None):
    """Structured failed estimate, retaining routing and species diagnostics without a point energy."""
    log(f"  ΔG = NONE ({reason})")
    return dict(reaction=key, dG=None, dG_raw=None, dG_raw_unreliable=None, suspect=reason, routes=routes,
                original_species={name: list(value) for name, value in reaction["species"].items()},
                energy_units="kJ/mol", quantity="standard_transformed_reaction_gibbs",
                stages={}, conditions=CONDITIONS, config=config, U_samp=None,
                exp=reaction.get("exp"), err=[], note=reaction.get("note", ""),
                experimental_reference={k: v for k, v in reaction.items() if k.startswith("exp_")},
                std_state_kJ=None, species_scored=None)


def no_estimate(key, reaction, reason, log=print):
    """Public failed-estimate record (used by batch drivers for exceptions outside score_reaction)."""
    return _no_estimate(key, reaction, {}, reason, log)


def _constants_hash(module, names):
    """sha256 (12 hex) of the named module-level constants."""
    payload = json.dumps({n: getattr(module, n, None) for n in names}, sort_keys=True, default=str)
    return hashlib.sha256(payload.encode()).hexdigest()[:12]


def effective_config(backend, *, allow_truncate=True, force_neutral=False):
    """Fingerprint of everything that changes a reported number: model and producer, physics version,
    numerical settings, routing constants, pKa and hydration constant tables, conditions and software."""
    from metag.energetics.minimum_optimizer import policy as minimum_optimizer_policy
    from metag.energetics.protonation_models import POLICY_VERSION as model_policy
    from metag.energetics.provenance import source_identity
    from metag.protonation_registry import POLICY_VERSION as registry_policy, registry_identity
    from metag.routing import aldehyde_hydration, pka_transform as pk
    from metag.routing.planner import TRUNCATION_RADIUS
    from metag.validation.stereochemistry import POLICY_VERSION as stereo_policy
    settings = getattr(backend, "settings", None)
    pu = getattr(backend, "pu", None)
    cfg = {"model": backend.model_name, "physics": PHYSICS_VERSION,
           "estimator_source_sha256": source_identity(),
           "settings": settings.as_dict() if settings is not None else {"dg_sanity_kj": 500.0},
           "species_settings": json.dumps(settings.species_key(), sort_keys=True) if settings is not None else None,
           "sampling_profile_version": SAMPLING_PROFILE_VERSION,
           "producer": backend.producer,
           **({"uma_runtime_profile": pu.metag_runtime_profile}
              if getattr(pu, "_metag_official_calculator", False) is True else {}),
           "solv_model": "alpb", "water_reference": "experimental_hydration",
           "aqueous_sampling_uncertainty": "2026-10-09-partition-derivatives",
           "reaction_sampling_uncertainty": "2026-10-09-shared-estimate-groups",
           "aqueous_state_coverage": "2026-10-07-unique-hydrates",
           "trunc_radius": TRUNCATION_RADIUS, "allow_truncate": bool(allow_truncate),
           "force_neutral": bool(force_neutral),
           **pk.pka_model_config(), "pka_table": pk.PKA_TABLE_VERSION,
           "pka_constants": _constants_hash(pk, ("P_LADDER", "ACYL_P_LADDER", "CARBONATE_LADDER", "PPI_LADDER",
                                                 "SULFATE_LADDER", "CARBOXYL_PKA", "CARBOXYL_PKA_ALPHA",
                                                 "SULFONATE_PKA", "THIOL_PKA", "PHENOL_PKA", "AMMONIA_PKA",
                                                 "AAA_AMINE_PKA", "PRIMARY_AMINE_PKA", "IMIDAZOLE_PKA",
                                                 "GUANIDINIUM_PKA")),
           "hydration_constants": _constants_hash(aldehyde_hydration, ("HYDRATION_CAL", "MAX_HYDRATION_SITES")),
           "source_protonation_registry_sha256": registry_identity(),
           "source_protonation_registry_policy": registry_policy,
           "species_protonation_model_policy": model_policy,
           "conditions": json.dumps(CONDITIONS, sort_keys=True),
           "minimum_optimizer_fallback": minimum_optimizer_policy(),
           "stereochemistry_policy": stereo_policy,
           "software": _software_versions()}
    return cfg


@functools.lru_cache(maxsize=1)
def _software_versions():
    """Versions of the external codes whose numerics enter the result."""
    import importlib.metadata as md
    import subprocess
    out = {}
    for pkg in ("fairchem-core", "torch", "rdkit", "ase", "numpy"):
        try:
            out[pkg] = md.version(pkg)
        except md.PackageNotFoundError:
            out[pkg] = None
    try:
        from metag.energetics.thermal import XTB
        r = subprocess.run([XTB, "--version"], capture_output=True, text=True, timeout=30)
        m = [ln for ln in (r.stdout + r.stderr).splitlines() if "xtb version" in ln]
        out["xtb"] = m[0].split("version")[1].split()[0] if m else None
    except Exception:
        out["xtb"] = None
    return json.dumps(out, sort_keys=True)
