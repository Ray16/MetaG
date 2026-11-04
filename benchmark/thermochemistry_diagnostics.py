"""Portable CPU diagnostics for retained-basis transfer and reaction-energy consistency.

Run ``python benchmark/thermochemistry_diagnostics.py transfer PLAN.json -o REPORT.json``
or ``... cycles NETWORK.json -o REPORT.json``. Needs NumPy, SymPy and SciPy.
No model, energy cache, reaction-label fitting or deployed prediction is modified.
The inequalities/optimization are standard mathematics, not a new thermodynamic law.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path


def finite(value, label):
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ValueError(label + ' must be finite numerical data')
    return float(value)


def _partition(values, rt):
    lo = min(values)
    return lo - rt * math.log(math.fsum(math.exp(-(v-lo)/rt) for v in values))


def transfer_diagnostic(plan):
    """Exact log partition and envelope on a declared common neutral-state basis.

    A center gap is not a significance statement. An optional interval is supplied
    by the caller, never inferred from data spread or assumed Gaussian uncertainty.
    """
    if plan.get('charge') != 0 or isinstance(plan.get('charge'), bool):
        raise ValueError('Only neutral transfer is supported; ions require a single-ion convention')
    if plan.get('standard_state') != '1M_gas_to_1M_solution':
        raise ValueError('Explicit common 1M gas/solution standard states required')
    rt = finite(plan['rt_kj'], 'RT')
    if rt <= 0:
        raise ValueError('RT must be positive')
    rows = plan['states']
    if not isinstance(rows, list) or not rows:
        raise ValueError('A nonempty retained state basis is required')
    ids, gas, solvent = set(), [], []
    for r in rows:
        name = r['id']
        if not isinstance(name, str) or not name or name in ids:
            raise ValueError('Each physical retained state must have one unique nonempty ID')
        ids.add(name)
        gas.append(finite(r['electronic_kj'], 'electronic energy') + finite(r['thermal_kj'], 'thermal energy'))
        solvent.append(finite(r['solvation_kj'], 'solvation shift'))
    if not all(math.isfinite(x) for x in gas):
        raise ValueError('Gas energy component addition overflowed')
    aq = [a+b for a,b in zip(gas,solvent)]
    if not all(math.isfinite(x) for x in aq):
        raise ValueError('Aqueous energy component addition overflowed')
    value = _partition(aq,rt)-_partition(gas,rt)
    lo, hi = min(solvent), max(solvent)
    center = plan.get('reference_center_kj')
    gap = None
    if center is not None:
        center = finite(center,'reference center')
        gap = max(lo-center, center-hi, 0.)
    interval = plan.get('reference_interval_kj')
    interval_gap = None
    if interval is not None:
        if not isinstance(interval,list) or len(interval)!=2:
            raise ValueError('Reference interval must contain two endpoints')
        left,right = [finite(x,'reference interval') for x in interval]
        if left>right or (center is not None and not left<=center<=right):
            raise ValueError('Ordered interval containing the supplied center required')
        interval_gap = max(lo-right,left-hi,0.)
    return {'n_declared_states':len(ids),'rt_kj':rt,'transfer_kj':value,
            'solvent_envelope_kj':[lo,hi],'reference_center_kj':center,
            'minimum_center_gap_under_reweighting_kj':gap,
            'supplied_interval_gap_kj':interval_gap,
            'center_excluded_on_retained_basis':None if gap is None else gap>0,
            'supplied_interval_excluded':None if interval_gap is None else interval_gap>0,
            'scope':'Common retained neutral-state basis and fixed solvent shifts. Molecular identity, state completeness and minimum validity are supplied assumptions, not certified here. Missing states/relaxation and reference uncertainty are not inferred. This is neither a total-error interval nor proof that a different solvent model is accurate.'}


def radius_from_cycles(cycles, energies, *, numerical_tolerance=1e-7):
    """Minimum max per-reaction perturbation required to satisfy supplied cycle constraints.

    min t subject to C e = -C g and |e_j| <= t. The primal solution is a
    mathematical witness only and is not returned as a new physical prediction.
    Equality dual multipliers supply a lower-bound certificate: ||C^T lambda||_1<=1.
    Basis invariance refers to the same row space and fixed reaction extent/units.
    """
    import numpy as np
    from scipy.optimize import linprog
    import warnings
    g = np.asarray([finite(x,'reaction energy') for x in energies],dtype=float)
    if not len(g):
        raise ValueError('At least one complete reaction energy is required')
    if len(cycles)==0:
        return {'radius_kj':0.,'n_constraints':0,'primal_dual_gap_kj':0.,
                'scope':'No cycle constraints: consistency is untested, not validated.',
                'diagnosis':'no_testable_cycles'}
    import sympy as sp
    if any(len(row)!=len(g) for row in cycles):
        raise ValueError('Cycle rows need one coefficient per supplied reaction')
    for row in cycles:
        for v in row:finite(v,'cycle coefficient')
    # Canonicalize the exact rational ROW SPACE before floating point optimization.
    # Directly solving an ill-conditioned but invertible basis transform can drop
    # nearly dependent constraints and produce a spurious basis-dependent radius.
    exact=sp.Matrix([[sp.Rational(str(v)) for v in row] for row in cycles])
    reduced,pivots=exact.rref()
    if not pivots:
        raise ValueError('Zero constraints carry no cycle test')
    C=np.asarray(reduced[:len(pivots),:].tolist(),dtype=float)
    if not np.isfinite(C).all():
        raise ValueError('Cycle row-space coefficients exceeded finite precision')
    scales = np.max(np.abs(C),axis=1)
    if (scales==0).any():
        raise ValueError('Zero constraint rows carry no test and must be omitted')
    C=C/scales[:,None]
    # fsum avoids losing small residuals in sums of extensive contributions.
    residual=np.array([math.fsum(float(c)*float(x) for c,x in zip(row,g)) for row in C])
    n=len(g);objective=np.r_[np.zeros(n),1.]
    bounds_matrix=np.vstack([np.c_[np.eye(n),-np.ones(n)],np.c_[-np.eye(n),-np.ones(n)]])
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore',message='Unrecognized options detected.*threads')
        answer=linprog(objective,A_ub=bounds_matrix,b_ub=np.zeros(2*n),
                       A_eq=np.c_[C,np.zeros(len(C))],b_eq=-residual,
                       bounds=[(None,None)]*n+[(0,None)],method='highs',
                       options={'primal_feasibility_tolerance':1e-9,
                                'dual_feasibility_tolerance':1e-9,'threads':1})
    if not answer.success:
        raise RuntimeError('Consistency-radius LP failed: '+answer.message)
    e=answer.x[:-1];radius=float(answer.x[-1]);dual=answer.eqlin.marginals
    witness=C.T@dual;norm=float(np.linalg.norm(witness,1));lower=float(-residual@dual)
    violation=float(np.max(np.abs(C@e+residual)));gap=radius-lower
    verification_tol=max(1e-7,abs(radius)*1e-8)
    if (violation>verification_tol or norm>1+1e-7 or abs(gap)>verification_tol
            or np.max(np.abs(e))>radius+verification_tol):
        raise RuntimeError('Primal/dual consistency certificate did not verify')
    individual=max(abs(float(r))/float(np.linalg.norm(c,1)) for r,c in zip(residual,C))
    return {'radius_kj':radius,'n_constraints':len(C),'individual_cycle_bound_kj':individual,
            'dual_lower_bound_kj':lower,'primal_dual_gap_kj':gap,
            'primal_constraint_violation_kj':violation,'dual_witness_l1_norm':norm,
            'dual_reaction_witness':witness.tolist(),'primal_perturbation_witness_kj':e.tolist(),'row_space_canonicalization':'exact_rational_rref',
            'diagnosis':'nonclosing' if radius>numerical_tolerance else 'numerically_consistent_on_supplied_space',
            'scope':'Minimum uniform change per supplied reaction at its declared extent, conditional on this constraint space. Not a calibrated physical error, accuracy estimate, prevalence, or corrected deployed prediction. Closure does not imply accuracy.'}


def cycle_diagnostic(plan):
    """Build the complete exact input-state nullspace, rather than inspect a chosen cycle list."""
    import sympy as sp
    reactions=plan['stoichiometry'];energies=plan['energies_kj']
    if not isinstance(reactions,dict) or not reactions or set(reactions)!=set(energies):
        raise ValueError('Stoichiometry and finite energies must identify the same nonempty complete reaction set')
    if plan.get('energy_units')!='kJ/mol':
        raise ValueError('Explicit kJ/mol reaction energy units required')
    if not isinstance(plan.get('conditions'),dict) or not plan['conditions']:
        raise ValueError('One explicit common thermodynamic condition record required')
    rids=sorted(reactions);states=sorted({s for col in reactions.values() for s in col})
    if any(not isinstance(s,str) or not s for s in states):
        raise ValueError('Nonempty explicit state identities required')
    state_index={s:i for i,s in enumerate(states)};entries={}
    for j,r in enumerate(rids):
        if not isinstance(r,str) or not r or not isinstance(reactions[r],dict):
            raise ValueError('Reaction IDs and stoichiometric columns must be well formed')
        finite(energies[r],'reaction energy')
        for s,v in reactions[r].items():
            finite(v,'stoichiometric coefficient')
            entries[(state_index[s],j)]=sp.Rational(str(v))
    matrix=sp.MutableSparseMatrix(len(states),len(rids),entries)
    vectors=matrix.nullspace();cycles=[]
    for vector in vectors:
        den=sp.ilcm(*[x.q for x in vector]) if len(vector)>1 else vector[0].q
        row=[int(x*den) for x in vector];divisor=math.gcd(*row);row=[v//divisor for v in row]
        if next(v for v in row if v)<0:row=[-v for v in row]
        assert matrix*sp.Matrix(row)==sp.zeros(len(states),1)
        cycles.append(row)
    values=[energies[r] for r in rids]
    report=radius_from_cycles(cycles,values)
    participating=[r for j,r in enumerate(rids) if any(row[j] for row in cycles)]
    report.update(n_reactions=len(rids),n_states=len(states),nullity=len(cycles),
                  n_cycle_participating_reactions=len(participating),
                  cycle_participating_reaction_ids=participating,
                  reaction_ids=rids,independent_cycles=[{r:w for r,w in zip(rids,row) if w} for row in cycles],
                  independent_cycle_residuals_kj=[math.fsum(w*x for w,x in zip(row,values)) for row in cycles],
                  conditions=plan['conditions'])
    return report


def _unique_pairs(pairs):
    out={}
    for k,v in pairs:
        if k in out:raise ValueError('Duplicate JSON key: '+k)
        out[k]=v
    return out


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('mode',choices=['transfer','cycles']);p.add_argument('input',type=Path)
    p.add_argument('-o','--output',type=Path,required=True);args=p.parse_args(argv)
    if args.output.exists():p.error('Output already exists; preserve it or choose a new path')
    raw=args.input.read_bytes();plan=json.loads(raw,object_pairs_hook=_unique_pairs)
    if args.mode=='transfer' and 'systems' in plan:
        result={name:transfer_diagnostic(v) for name,v in plan['systems'].items()}
    else:result=(transfer_diagnostic if args.mode=='transfer' else cycle_diagnostic)(plan)
    import importlib.metadata as md
    versions={name:md.version(name) for name in ['numpy','sympy','scipy']}
    output={'software_versions':versions,'input_sha256':hashlib.sha256(raw).hexdigest(),
            'tool_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'result':result,'method':'Retained-basis transfer envelope / complete input-state cycle nullspace with primal-dual minimax certificate'}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    with args.output.open('x') as f:json.dump(output,f,indent=2,allow_nan=False);f.write('\n')
    print('Wrote',args.output)


if __name__=='__main__':main()
