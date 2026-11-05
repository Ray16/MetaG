"""Bounded same-Hamiltonian fallback for force-unconverged gas minima.

This helper does not validate vibrations, assemble energies, or select states.
The caller must validate the returned geometry's Hessian and recompute its
electronic, solvation and thermal terms before admitting it to an ensemble.
"""
import hashlib
import time

import numpy as np
from ase.calculators.calculator import Calculator, all_changes
from ase.optimize import LBFGS

POLICY_VERSION = '2026-11-04-lbfgs-force-and-curvature-recovery-v2'
MAX_STEPS = 1200
MAX_SECONDS = 300.0
FMAX = 0.01
MAXSTEP = 0.2
CURVATURE_MAX_STEPS = 400
CURVATURE_MAX_SECONDS = 180.0
CURVATURE_FMAX = 0.001


def curvature_policy():
    return {'version': POLICY_VERSION, 'method': 'ASE_LBFGS',
            'trigger': 'strong imaginary mode persists after final FIRE/mode-following retry',
            'starting_geometry': 'same identity-preserving force-converged candidate',
            'maximum_steps': CURVATURE_MAX_STEPS, 'maximum_seconds': CURVATURE_MAX_SECONDS,
            'fmax_eV_A': CURVATURE_FMAX, 'maxstep_A': MAXSTEP,
            'acceptance': 'fresh endpoint force, identity and subsequent unchanged strict Hessian gate'}


def policy():
    return {'version': POLICY_VERSION, 'method': 'ASE_LBFGS',
            'trigger': 'identity-preserving force-unconverged FIRE at final long thermal retry',
            'starting_geometry': 'same pre-FIRE candidate', 'maximum_steps': MAX_STEPS,
            'maximum_seconds': MAX_SECONDS, 'fmax_eV_A': FMAX, 'maxstep_A': MAXSTEP,
            'acceptance': 'endpoint force, identity and subsequent strict Hessian; all states retained',
            'persistent_curvature': curvature_policy()}


def _energy_forces(predictor, atoms):
    from metag.energetics.uma import _predict
    energy, forces, _ = _predict(predictor, [atoms])
    return float(energy.detach().cpu().numpy()[0]), forces.detach().cpu().numpy()


def recover_force_minimum(predictor, atoms, template, ref_graph, *, steps=MAX_STEPS):
    """Bounded recovery of a force-unconverged endpoint; vibrational QC follows."""
    return _recover_lbfgs(predictor, atoms, template, ref_graph, steps=steps,
                         maximum_steps=MAX_STEPS, maximum_seconds=MAX_SECONDS,
                         force_target=FMAX, policy_record=policy())


def recover_curvature_minimum(predictor, atoms, template, ref_graph, *, steps=CURVATURE_MAX_STEPS):
    """Search more tightly when force convergence did not establish a minimum.

    This does not accept a negative mode or supply a thermal correction. The
    caller must recompute and pass the endpoint's strict Hessian and energies.
    """
    return _recover_lbfgs(predictor, atoms, template, ref_graph, steps=steps,
                         maximum_steps=CURVATURE_MAX_STEPS, maximum_seconds=CURVATURE_MAX_SECONDS,
                         force_target=CURVATURE_FMAX, policy_record=curvature_policy())


def _recover_lbfgs(predictor, atoms, template, ref_graph, *, steps, maximum_steps,
                   maximum_seconds, force_target, policy_record):
    """Return a copied endpoint and JSON-safe diagnostic; never mutate the input.

    Time is checked between model calls; the caller's hard job timeout must
    still bound a backend that hangs inside an individual evaluation.
    """
    from metag.energetics.conformers import same_connectivity, same_stereochemistry
    if isinstance(steps, bool) or int(steps) != steps or steps < 1:
        raise ValueError('positive integer optimizer step budget required')
    limit = min(int(steps), maximum_steps)
    started = time.monotonic()
    endpoint = atoms.copy()
    record = dict(policy_record, status='started', steps=0, model_evaluations=0,
                  step_budget=limit, converged=False, identity_preserved=False,
                  force_accepted=False, thermal_validation='required_by_caller',
                  starting_positions_sha256=hashlib.sha256(
                      np.asarray(endpoint.positions, dtype='<f8').tobytes()).hexdigest(),
                  charge=int(endpoint.info['charge']), spin_multiplicity=int(endpoint.info['spin']))
    optimizer = None

    class BudgetExceeded(RuntimeError):
        pass

    def evaluate(current):
        if time.monotonic() - started > maximum_seconds:
            raise BudgetExceeded('LBFGS wall-time budget exhausted')
        record['model_evaluations'] += 1
        energy, forces = _energy_forces(predictor, current)
        forces = np.asarray(forces, dtype=float)
        if not np.isfinite(energy) or forces.shape != (len(current), 3) or not np.isfinite(forces).all():
            raise ValueError('nonfinite or malformed UMA energy/forces')
        return float(energy), forces

    class SameHamiltonian(Calculator):
        implemented_properties = ['energy', 'forces']

        def calculate(self, atoms=None, properties=('energy', 'forces'), system_changes=all_changes):
            super().calculate(atoms, properties, system_changes)
            energy, forces = evaluate(atoms)
            self.results = {'energy': energy, 'forces': forces}

    try:
        if (not np.isfinite(endpoint.positions).all() or
                not same_connectivity(endpoint, ref_graph) or
                not same_stereochemistry(endpoint, template)):
            record['status'] = 'rejected_initial_identity'
            return endpoint, record
        endpoint.calc = SameHamiltonian()
        optimizer = LBFGS(endpoint, maxstep=MAXSTEP, logfile=None)
        converged = bool(optimizer.run(fmax=force_target, steps=limit))
        # Reevaluate explicitly: acceptance must describe the returned endpoint,
        # independent of optimizer/cached force or convergence bookkeeping.
        energy, forces = evaluate(endpoint)
        maximum = float(np.linalg.norm(forces, axis=1).max())
        identity = same_connectivity(endpoint, ref_graph) and same_stereochemistry(endpoint, template)
        record.update(status='complete', converged=converged,
                      energy_eV=energy, max_force_eV_A=maximum,
                      identity_preserved=bool(identity),
                      force_accepted=bool(converged and maximum <= force_target and identity))
    except Exception as exc:
        record.update(status='budget_exit' if isinstance(exc, BudgetExceeded) else 'exception',
                      exception=repr(exc))
    finally:
        positions_finite = bool(np.isfinite(endpoint.positions).all())
        record.update(elapsed_seconds=time.monotonic() - started,
                      steps=int(optimizer.nsteps) if optimizer is not None else 0,
                      positions_finite=positions_finite,
                      positions_angstrom=endpoint.positions.tolist() if positions_finite else None)
        endpoint.calc = None
    return endpoint, record
