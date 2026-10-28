"""Aqueous free energy of one species: conformer search, UMA relaxation, ALPB solvation, per-minimum RRHO.

    G = -RT ln Σ_i exp(-(E_UMA,i + ΔG_solv,i + G_corr,i)/RT)

over the unique, validated gas-phase minima i. Conformer seed batches are added until the ensemble free
energy and the lowest electronic energy both stop moving (a plateau heuristic, not a completeness proof).
Every sampled minimum must pass strict thermal validation (no imaginary mode after one tight
re-optimisation, mode following or artefact probing); any unresolved minimum fails the species instead of
being dropped. Values are on the 1 atm gas standard state; the 1 M correction is added by the caller.
"""
import time
from concurrent.futures import ThreadPoolExecutor

import numpy as np
from ase import Atoms
from rdkit import Chem

from metag.constants import EV2KJ
from metag.energetics import species_cache as cache
from metag.energetics.conformers import (XTB_WORKERS, UniqueMinima, bond_graph, boltz, pool_confs,
                                         same_connectivity, same_stereochemistry, spin_multiplicity)
from metag.energetics.errors import SpeciesRearranged
from metag.energetics.provenance import ordered_geometry_record
from metag.energetics.sampling import sampling_budget, sampling_diagnostics, species_diagnostics
from metag.energetics.thermal import dgsolv, uma_gibbs_corr
from metag.energetics.thermal_ensemble import ThermalTrack
from metag.energetics.uma import batched_energies, batched_fire

# Fraction of relaxed conformers whose xtb solvation may fail before G is considered degraded: a degraded G
# is used for this reaction but never cached, so a transient overload cannot poison later reactions.
MAX_XTB_FAIL_FRAC = 0.25
# Every relaxed conformer within this gas-phase window of the lowest must have a solvation value. Failures
# are not random (compact / highly charged conformers time out more), so a species with missing low-energy
# members fails instead of being biased.
SOLV_SUPPORT_KJ = 30.0
# Lowest-G unique minima tried in order for a strict thermal minimum.
N_THERMAL_CANDIDATES = 3


def _gibbs_corr(ctx, syms, pos, q, mult, **kwargs):
    return uma_gibbs_corr(ctx.pu, syms, pos, q, spin=mult, return_info=True, hessian=ctx.settings.hessian,
                          **kwargs)


class _Context:
    """Per-species constants threaded through the helpers."""

    def __init__(self, pu, settings, producer, log, name):
        self.pu, self.settings, self.producer, self.log, self.name = pu, settings, producer, log, name


def species_free_energy(pu, q, smi, name, *, settings, model_name, producer, log, warnings=None,
                        diagnostics=None):
    """(G, sigma) in kJ/mol for species `smi` with formal charge `q`, or (None, None) if it fails.

    Raises SpeciesRearranged when every sampled conformer changed bonding during gas relaxation (e.g. a
    zwitterion's proton transfer): the species is not a gas-phase minimum under this SMILES. Species-level
    problems are appended to `warnings` and stored with the cached value, so a cache hit reports them too.
    `diagnostics` receives sampling and thermal metadata."""
    # Cache identity and the finite conformer search use the same atom ordering.
    smi = cache.canonical(smi)
    if diagnostics is not None:
        diagnostics.update(method="implicit", source="fresh", sampling={"status": "incomplete"})
    warnings = warnings if warnings is not None else []
    ctx = _Context(pu, settings, producer, log, name)
    mult = spin_multiplicity(smi, q)                      # ground-state spin (O2 triplet, radicals doublet)
    key = dict(settings.species_key(), model=model_name, producer=producer)
    if mult != 1:
        key["spin"] = mult
    template = Chem.AddHs(Chem.MolFromSmiles(smi))
    cached = cache.get(smi, q, "implicit", key, with_meta=True)
    if cached is not None:
        G_c, s_c, meta_c = cached
        if diagnostics is not None:
            diagnostics.update(species_diagnostics(meta_c, "cache"))
        warnings.extend(meta_c.get("warnings", []))
        log(f"    {name:9s} q{q:+d} [implicit CACHED]: {G_c:.1f}")
        return G_c, s_c
    if mult != 1:
        log(f"    {name:9s} q{q:+d} [open-shell: spin multiplicity {mult}]")
    _, keep, pool = sampling_budget(smi, sample_scale=settings.sample_scale, keep_scale=settings.keep_scale)
    ref_graph = bond_graph(template)
    uniq = UniqueMinima(template=template)
    # CONNECTIVITY GUARD: a conformer whose bonding or stereochemistry changed during relaxation is a
    # different species; such conformers are set aside.
    rearranged = []
    track = ThermalTrack(pu, template, q, mult,
                         correction=lambda pu_, syms, pos, q_, **kw: uma_gibbs_corr(
                             pu_, syms, pos, q_, hessian=settings.hessian, **kw),
                         combine=boltz, window=settings.thermal_window_kj, maximum=settings.thermal_max)
    n_relaxed = n_xtb_fail = 0
    failed_E = []                                         # gas energies of conformers whose solvation failed
    best_E = 1e18
    prev_Gens = prev_best = None
    hits = 0
    seed = 0
    last_dG = float("nan")
    gens_traj, successful_seeds, minimum_traj = [], [], []
    while seed < settings.conv_max:
        seed += 1
        cands = pool_confs(smi, q, seed + settings.seed_offset, pool, spin=mult)
        order = np.argsort(batched_energies(pu, cands))[:keep]
        sel = [cands[i] for i in order]
        rel, E, conv = batched_fire(pu, sel, fmax=0.05, steps=300, stop_frac=0.9,
                                    return_converged=True, label=f"{name}s{seed}")
        sel = [a for a, c in zip(rel, conv) if c]
        Eg = E[conv] * EV2KJ
        with ThreadPoolExecutor(max_workers=XTB_WORKERS) as ex:
            solv = list(ex.map(lambda a: dgsolv(a.get_chemical_symbols(), a.get_positions(), q, mult), sel))
        for i, a in enumerate(sel):                       # one retry: most xtb failures are transient
            if solv[i] is not None and not np.isfinite(solv[i]):
                solv[i] = None                            # nonfinite = failed evaluation, never a population
            if solv[i] is None:
                retry = dgsolv(a.get_chemical_symbols(), a.get_positions(), q, mult)
                solv[i] = retry if retry is not None and np.isfinite(retry) else None
        eligible_in_batch = 0
        for a, e, s in zip(sel, Eg, solv):
            n_relaxed += 1
            if s is None:
                n_xtb_fail += 1
                if np.isfinite(e):
                    failed_E.append(float(e))
            if not np.isfinite(e) or s is None:
                continue
            if not same_connectivity(a, ref_graph) or not same_stereochemistry(a, template):
                rearranged.append(a)
                continue
            uniq.add(a, float(e), float(e) + s)
            # Rediscovering an existing valid minimum is real search evidence; a batch with no eligible
            # conformer is not a repeated measurement and cannot contribute a plateau hit.
            eligible_in_batch += 1
            best_E = min(best_E, float(e))
        if not uniq.G or eligible_in_batch == 0:
            hits = 0
            continue
        Gens = track.update(uniq)                          # convergence on E + ΔGsolv + provisional RRHO
        gens_traj.append(Gens)
        successful_seeds.append(seed)
        minimum_traj.append(best_E)
        if prev_Gens is not None:
            last_dG = abs(Gens - prev_Gens)
            if last_dG < settings.conv_tol and abs(best_E - prev_best) < settings.conv_tol:
                hits += 1
                if hits >= settings.conv_hits:
                    break
            else:
                hits = 0
        prev_Gens, prev_best = Gens, best_E
    sp_warn = []
    if not uniq.G:
        if rearranged:
            raise SpeciesRearranged(name, smi, len(rearranged))
        return None, None
    if rearranged:
        log(f"    {name}: {len(rearranged)} conformer(s) changed bonding during relaxation -> excluded")
    missing = [e for e in failed_E if e - min(uniq.E) < SOLV_SUPPORT_KJ]
    if missing:
        sp_warn.append(f"{name} ({smi}): solvation failed for {len(missing)} conformer(s) within "
                       f"{SOLV_SUPPORT_KJ:.0f} kJ of the minimum -> ensemble incomplete, species failed")
        log(f"    !! {sp_warn[-1]}")
        warnings.extend(sp_warn)
        return None, None
    validation = {}
    validated = _validated_minima(ctx, uniq, template, ref_graph, q, mult, track=track, diagnostics=validation)
    if validated is None:
        message = f"{name}: sampled thermal ensemble failed ({validation.get('reason', 'unresolved')})"
        warnings.append(message)
        if diagnostics is not None:
            diagnostics["thermal"] = {"status": "incomplete", "n_sampled": len(uniq.G),
                                      "validation_maximum": settings.thermal_validation_max, **validation}
        log(f"    !! {message}")
        return None, None
    uniq, corrections, validation_info = validated
    ref_j = min(range(len(uniq.G)), key=lambda j: uniq.G[j])
    G, _, ens_info = track.final(uniq, {}, corrections[ref_j], ref_j, validated=dict(enumerate(corrections)))
    ens_info.update(validation_info)
    t_info = {"ensemble": ens_info}
    # sampling uncertainty: spread of the ensemble G over the last few batches
    tail = gens_traj[-3:]
    sigma = float(np.std(tail)) if len(tail) > 1 else (last_dG if np.isfinite(last_dG) else 3.0)
    tag = "conv" if hits >= settings.conv_hits else "CAPPED"
    log(f"    {name:9s} q{q:+d} [implicit {tag} seeds={seed} minima={len(uniq)}/{uniq.n_seen} σ={sigma:.1f}]: "
        f"Gens {G:.1f} + thermal 0.0 = {G:.1f}")
    sampling_info = sampling_diagnostics(gens_traj, successful_seeds, minimum_traj, seeds=seed, hits=hits,
                                         required_hits=settings.conv_hits, maximum=settings.conv_max,
                                         tolerance=settings.conv_tol)
    sampling_info["ensemble_energy_scope"] = "gas minima plus vertical solvation and per-minimum thermal approximations"
    sampling_info["seed_offset"] = settings.seed_offset
    sampling_info["solution_relaxation_in_trajectory"] = False
    sampling_info["final_validation_in_trajectory"] = False
    sampling_info["final_minus_last_provisional_kj"] = float(G - gens_traj[-1]) if gens_traj else None
    if diagnostics is not None:
        diagnostics.update(species_diagnostics({"sampling": sampling_info, "thermal": t_info}, "fresh"))
    fail_frac = n_xtb_fail / max(n_relaxed, 1)
    if fail_frac > MAX_XTB_FAIL_FRAC:                     # degraded ensemble: use it, but never cache it
        sp_warn.append(f"{name} ({smi}): xtb solvation failed on {n_xtb_fail}/{n_relaxed} conformers "
                       f"-> G not cached")
        log(f"    !! {sp_warn[-1]}")
    else:
        meta = {"warnings": sp_warn, "n_minima": len(uniq), "n_seen": uniq.n_seen,
                "n_xtb_fail": n_xtb_fail, "n_rearranged": len(rearranged), "seeds": seed, "thermal": t_info,
                "sampling": sampling_info}
        cache.put(smi, q, "implicit", key, G, sigma, meta=meta)
    warnings.extend(sp_warn)
    return G, sigma


def _validated_minima(ctx, uniq, template, ref_graph, q, mult, *, track, diagnostics):
    """Resolve every sampled candidate, then deduplicate the resulting geometries.

    Exceeding the candidate budget or an unresolved candidate fails the species rather than silently
    discarding its population. Works on singleton copies so index changes during saddle resolution cannot
    invalidate another candidate's bookkeeping. Returns (minima, corrections, info) or None."""
    s = ctx.settings
    if not uniq.G or len(uniq.G) > s.thermal_validation_max:
        diagnostics.update(reason="candidate validation budget exceeded" if uniq.G else "empty ensemble",
                           n_sampled=len(uniq.G), validation_maximum=s.thermal_validation_max)
        return None
    output = UniqueMinima(template=template)
    corrections, kinds = [], []
    optimizer_trials = []
    diagnostics["optimizer_trials"] = optimizer_trials
    n_reused = n_retried = 0
    syms = [a.GetSymbol() for a in template.GetAtoms()]
    for j in range(len(uniq.G)):
        one = UniqueMinima(template=template)
        at = Atoms(syms, positions=uniq.ref[j].GetConformer().GetPositions())
        one.add(at, uniq.E[j], uniq.G[j])
        strict = track.validated.get(track._key(uniq, j))
        if strict is None:
            g, info, resolution = _thermal_at_minimum(ctx, one, template, ref_graph, q, mult)
        else:
            g, info = strict
            resolution = {"index": 0, "kind": "as_is", "pos": at.positions, "rejected": []}
            n_reused += 1
        if (g is None or not np.isfinite(g)) and s.thermal_retry_steps > s.thermal_relax_steps:
            # Retry the SAME candidate with more optimization work, not a looser threshold.
            n_retried += 1
            ctx.log(f"    {ctx.name}: retry unresolved minimum with {s.thermal_retry_steps} optimization steps")
            g, info, resolution = _thermal_at_minimum(ctx, one, template, ref_graph, q, mult,
                                                      optimization_steps=s.thermal_retry_steps)
        if g is None or not np.isfinite(g):
            diagnostics.update(reason="unresolved minimum or nonfinite thermal correction", candidate_index=j,
                               thermal_attempts=info,
                               failure_geometry={
                                   "scope": "failed sampled gas-state candidate; not an accepted thermochemical minimum",
                                   "symbols": syms, "charge": q, "spin_multiplicity": mult,
                                   "reference_smiles": Chem.MolToSmiles(Chem.RemoveHs(template)),
                                   "initial_positions_angstrom": at.positions.tolist(),
                                   "initial_electronic_kj": float(uniq.E[j]),
                                   "initial_vertical_kj": float(uniq.G[j]),
                                   "thermal_trials": resolution.get("thermal_trials", []),
                                   "optimization_trials": resolution.get("optimization_trials", []),
                                   "seed_offset": s.seed_offset, "producer": ctx.producer})
            return None
        optimizer_trials.extend(dict(trial, candidate_index=j)
                                for trial in resolution.get("optimization_trials", []))
        try:
            idx = _apply_resolution(ctx, one, resolution, template, q, mult)
        except RuntimeError as exc:
            ctx.log(f"    {ctx.name}: thermal geometry update failed: {exc}")
            diagnostics.update(reason="validated geometry energy update failed", candidate_index=j, error=str(exc))
            return None
        if idx is None or not np.isfinite(one.E[idx]) or not np.isfinite(one.G[idx]):
            diagnostics.update(reason="invalid resolved geometry energy", candidate_index=j)
            return None
        at.positions = one.ref[idx].GetConformer().GetPositions()
        added = output.add(at, one.E[idx], one.G[idx])
        k = output.last_index
        if added:
            corrections.append(float(g))
        elif output.last_replaced:
            corrections[k] = float(g)
        kinds.append(resolution["kind"])
    output.n_seen = uniq.n_seen
    return output, corrections, {
        "n_candidates_validated": len(uniq.G), "n_strict_results_reused": n_reused,
        "n_candidates_retried": n_retried, "retry_optimization_steps": s.thermal_retry_steps,
        "n_merged_after_resolution": len(uniq.G) - len(output.G),
        "validation_maximum": s.thermal_validation_max, "resolution_kinds": kinds,
        "optimizer_trials": optimizer_trials,
        "validation_scope": "all sampled candidates; heuristic mode-following artefact test; no unsampled-basin bound",
        "validated_minima": [
            {"electronic_kj": e, "solvation_kj": vertical - e, "thermal_kj": g,
             "total_kj": vertical + g,
             "reference_smiles": Chem.MolToSmiles(Chem.RemoveHs(output.ref[j])),
             "geometry": ordered_geometry_record(
                 syms, output.ref[j].GetConformer().GetPositions(), q, mult),
             "geometry_scope": "Accepted gas minimum; these component energies refer to the exported ordered geometry."}
            for j, (e, vertical, g) in enumerate(zip(output.E, output.G, corrections))]}


def _thermal_at_minimum(ctx, uniq, template, ref_graph, q, mult, *, optimization_steps=None):
    """RRHO correction at the lowest-G unique minimum that is a TRUE minimum. Candidates are the unique
    minima in order of G; one with an imaginary mode above thermal.IMAG_TOL_CM is re-optimised once more
    tightly (fmax 0.01 eV/Å); if still not a minimum, every strong imaginary mode is followed; else the
    next candidate is tried. Returns (Gcorr, info, resolution) or (None, per-candidate info, resolution).

    resolution = {"index": accepted minimum index or None, "kind": "as_is" | "artefact" | "retightened" |
    "mode_followed" | None, "pos": geometry the correction belongs to, "rejected": [indices that are not
    minima]}. _apply_resolution makes E, ΔG_solv, geometry and Gcorr refer to the SAME structure."""
    pu, s, name, log = ctx.pu, ctx.settings, ctx.name, ctx.log
    optimization_steps = s.thermal_relax_steps if optimization_steps is None else optimization_steps
    syms = [a.GetSymbol() for a in template.GetAtoms()]
    order = sorted(range(len(uniq.G)), key=lambda j: uniq.G[j])[:N_THERMAL_CANDIDATES]
    cands = [(j, uniq.ref[j].GetConformer().GetPositions()) for j in order]
    tried, rejected = [], []
    thermal_trials, optimization_trials = [], []

    def record_thermal(positions, correction, info):
        thermal_trials.append({
            "positions_angstrom": np.asarray(positions, dtype=float).tolist(),
            "correction_kj": float(correction) if correction is not None and np.isfinite(correction) else None,
            "n_imag": int(info.get("n_imag", 0)),
            "max_imag_cm": float(info.get("max_imag_cm", 0))})
    for j, pos0 in cands:
        pos = pos0
        for attempt in (0, 1):
            g, info = _gibbs_corr(ctx, syms, pos, q, mult)
            record_thermal(pos, g, info)
            vecs = info.pop("_imag_vecs")
            if info["n_imag"] == 0:
                info["retightened"] = bool(attempt)
                return g, info, {"index": j, "kind": "retightened" if attempt else "as_is", "pos": pos,
                                 "rejected": rejected}
            if attempt == 1:
                # A harmless strongest mode does not establish that every other negative mode is harmless.
                all_artefacts = bool(vecs) and len(vecs) == info["n_imag"]
                for mode in vecs:
                    verdict, new_pos = _follow_imaginary_mode(ctx, syms, pos, mode, q, mult, ref_graph,
                                                              template=template,
                                                              optimization_steps=optimization_steps)
                    if verdict != "artefact":
                        all_artefacts = False
                    if verdict == "saddle":
                        g2, info2 = _gibbs_corr(ctx, syms, new_pos, q, mult)
                        record_thermal(new_pos, g2, info2)
                        info2.pop("_imag_vecs")
                        if info2["n_imag"] == 0:
                            info2.update(retightened=True, mode_followed=True)
                            return g2, info2, {"index": j, "kind": "mode_followed", "pos": new_pos,
                                              "rejected": rejected}
                if all_artefacts:
                    g_soft, info_soft = _gibbs_corr(ctx, syms, pos, q, mult, imag_as_soft=True)
                    info_soft.pop("_imag_vecs")
                    info_soft.update(retightened=True, imag_artefact_cm=info["max_imag_cm"],
                                     n_imag_probed=len(vecs))
                    kind = "retightened" if pos is not pos0 else "artefact"
                    return g_soft, info_soft, {"index": j, "kind": kind, "pos": pos, "rejected": rejected}
            tried.append(info["max_imag_cm"])
            if attempt == 0:                              # re-optimise tightly once, keep only if unrearranged
                at = Atoms(symbols=syms, positions=pos, info={"charge": int(q), "spin": int(mult)})
                fire_started = time.monotonic()
                (at,), _, conv = batched_fire(pu, [at], fmax=0.01, steps=optimization_steps, return_converged=True,
                                              label=f"{name}-tight")
                identity_ok = same_connectivity(at, ref_graph) and same_stereochemistry(at, template)
                optimization_trials.append({"positions_angstrom": at.positions.tolist(),
                                            "converged": bool(conv[0]), "identity_preserved": bool(identity_ok),
                                            "optimizer": "batched_FIRE", "elapsed_seconds": time.monotonic()-fire_started,
                                            "step_budget": optimization_steps, "fmax_eV_A": 0.01})
                if not identity_ok:
                    break
                if (not conv[0] and s.thermal_retry_steps > s.thermal_relax_steps and
                        optimization_steps >= s.thermal_retry_steps):
                    from metag.energetics.minimum_optimizer import recover_force_minimum
                    fallback_start = Atoms(symbols=syms, positions=pos,
                                           info={"charge": int(q), "spin": int(mult)})
                    recovered, trial = recover_force_minimum(
                        pu, fallback_start, template, ref_graph, steps=optimization_steps)
                    optimization_trials.append(trial)
                    if trial["force_accepted"]:
                        g_recovered, info_recovered = _gibbs_corr(ctx, syms, recovered.positions, q, mult)
                        record_thermal(recovered.positions, g_recovered, info_recovered)
                        info_recovered.pop("_imag_vecs", None)
                        trial["thermal_n_imag"] = int(info_recovered["n_imag"])
                        trial["thermal_max_imag_cm"] = float(info_recovered["max_imag_cm"])
                        trial["strict_thermal_accepted"] = bool(
                            info_recovered["n_imag"] == 0 and g_recovered is not None and np.isfinite(g_recovered))
                        if trial["strict_thermal_accepted"]:
                            info_recovered.update(retightened=True, optimizer="ASE_LBFGS")
                            log(f"    {name}: strict minimum recovered by ASE_LBFGS "
                                f"({trial['steps']} steps, {trial['model_evaluations']} model evaluations, "
                                f"{trial['elapsed_seconds']:.2f} s)")
                            return g_recovered, info_recovered, {
                                "index": j, "kind": "retightened", "pos": recovered.positions,
                                "rejected": rejected, "optimization_trials": optimization_trials,
                                "thermal_trials": thermal_trials}
                    log(f"    {name}: LBFGS fallback did not establish a strict minimum: {trial}")
                if conv[0]:                               # an unconverged tight re-opt (flat torsions) still
                    pos = at.get_positions()              # gets the mode-following test from the original pos
        rejected.append(j)                                # not a minimum: excluded from the ensemble
        log(f"    {name}: minimum candidate has an imaginary mode ({tried[-1]} cm-1) -> next candidate")
    return None, tried, {"index": None, "kind": None, "pos": None, "rejected": rejected,
                         "thermal_trials": thermal_trials, "optimization_trials": optimization_trials}


def _apply_resolution(ctx, uniq, res, template, q, mult):
    """Make the ensemble consistent with the thermal validation. Returns the index (after the update) of the
    minimum the validated correction belongs to."""
    syms = [a.GetSymbol() for a in template.GetAtoms()]
    j, kind, pos = res["index"], res["kind"], res["pos"]
    drop = set(res["rejected"])
    if kind in ("retightened", "mode_followed"):
        at = Atoms(symbols=syms, positions=pos, info={"charge": int(q), "spin": int(mult)})
        e = float(batched_energies(ctx.pu, [at])[0]) * EV2KJ
        sd = dgsolv(syms, pos, q, mult)
        if not np.isfinite(e) or sd is None or not np.isfinite(sd):
            raise RuntimeError(f"{ctx.name}: solvation failed on the validated thermal geometry")
        if kind == "retightened":                         # same basin: the representative becomes this geometry
            uniq.E[j] = e; uniq.G[j] = e + sd; uniq.ref[j] = uniq._mol_at(at)
        else:                                             # saddle j removed; the lower structure enters
            drop.add(j)
            uniq.E.append(e); uniq.G.append(e + sd); uniq.ref.append(uniq._mol_at(at))
            j = len(uniq.G) - 1
        ctx.log(f"    {ctx.name}: thermal validation -> {kind} geometry carries E, ΔGsolv and RRHO")
    if drop:
        keep = uniq.drop(drop)
        j = keep.index(j) if j in keep else None
        ctx.log(f"    {ctx.name}: {len(drop)} non-minimum candidate(s) removed from the ensemble")
    return j


def _follow_imaginary_mode(ctx, syms, pos, mode, q, mult, ref_graph, step=0.1, drop_kJ=1.0,
                           template=None, optimization_steps=None):
    """('saddle', lower_pos) if relaxing from pos +- step·mode reaches an energy > drop_kJ below pos (same
    connectivity); ('artefact', None) if both probes return; otherwise ('inconclusive', None)."""
    pu, name = ctx.pu, ctx.name
    at0 = Atoms(symbols=syms, positions=pos, info={"charge": int(q), "spin": int(mult)})
    e0 = float(batched_energies(pu, [at0])[0]) * EV2KJ
    disp = [Atoms(symbols=syms, positions=pos + sgn * step * mode, info={"charge": int(q), "spin": int(mult)})
            for sgn in (+1.0, -1.0)]
    steps = ctx.settings.thermal_relax_steps if optimization_steps is None else optimization_steps
    rel, E, conv = batched_fire(pu, disp, fmax=0.01, steps=steps, return_converged=True, label=f"{name}-mode")
    best = None
    n_valid = 0
    probes = []
    for a, e, c in zip(rel, E, conv):
        e = float(e) * EV2KJ
        identity = same_connectivity(a, ref_graph) and (template is None or same_stereochemistry(a, template))
        probes.append({"converged": bool(c), "identity": bool(identity),
                       "delta_kj": float(e-e0) if np.isfinite(e) and np.isfinite(e0) else None})
        if c and np.isfinite(e) and np.isfinite(e0) and identity:
            n_valid += 1
            if e < e0 - drop_kJ and (best is None or e < best[0]):
                best = (e, a.get_positions())
    ctx.log(f"    {name}: imaginary-mode probes {probes}")
    if best is not None:
        return "saddle", best[1]
    # A failed/rearranged/nonfinite probe cannot establish that negative curvature is numerical noise.
    return ("artefact", None) if n_valid == 2 else ("inconclusive", None)


def water_geometry(pu):
    """Lowest converged gas-phase water geometry from a small conformer pool."""
    cands = pool_confs("O", 0, 1, 48)
    sel = [cands[i] for i in np.argsort(batched_energies(pu, cands))[:8]]
    rel, E, conv = batched_fire(pu, sel, fmax=0.05, steps=300, stop_frac=0.9, return_converged=True)
    template = Chem.AddHs(Chem.MolFromSmiles("O"))
    graph = bond_graph(template)
    valid = [(a, float(e)) for a, e, c in zip(rel, E, conv)
             if c and np.isfinite(e) and same_connectivity(a, graph) and same_stereochemistry(a, template)]
    if not valid:
        raise ValueError("no converged water geometry preserving molecular identity")
    a, _ = min(valid, key=lambda item: item[1])
    return a.get_chemical_symbols(), a.get_positions()
