# Experimental state function profile

`StateFunctionEstimator` assigns each admitted molecular input one retained
potential for a calculation session. Reaction class, direction, extent and
neighbouring metabolites cannot change that assignment. This experimental API
addresses reaction-dependent reference selection and molecular reduction. It
does not change `metag.pipeline.score_reaction` or inherit that estimator's
accuracy and coverage results.

```python
from metag.state_function import StateFunctionEstimator

# backend implements the same species/liquid_water interface as UMABackend.
estimator = StateFunctionEstimator(backend)
results = {key: estimator.score_reaction(rx, key=key)
           for key, rx in reactions.items()}
ledger = estimator.ledger()
```

Use the same instance across the network. Raw molecular energies and the liquid
water reference are retained once, including when the backend is stochastic or
does not normally cache. Returned records and ledger exports are copies.
Backend producer or settings changes invalidate the session. The exported
ledger retains provenance; it is not accepted as a resumable molecular cache.
Execution is serial.

The thermodynamic target remains 298.15 K, pH 7, zero ionic strength, no Mg
binding, 1 M solutes and pure liquid water. An optional reaction-level
`conditions` record constrains that target. Unsupported requirements return no
estimate before molecular work; they do not change constants or reweight an
existing ledger. For example, a request at pH 7.5 and ionic strength 0.12 M is
unavailable. Reference measurements at those conditions need an independently
justified transfer before comparison with this profile.

## Reference treatment

The profile uses full molecules without cofactor replacement, nucleotide caps,
generic truncation or a whole-reaction neutral-reference fallback. Each molecule
uses the existing global acid/base reference construction. Its conversion must
preserve every element and charge when the explicit proton exchange is included.
Source-backed acid ladders are resolved through the existing exact-reference
registry; other admitted sites retain the existing class-based approximations.

For an input state with charge `q` and computed reference charge `q_ref`, the
retained input-gauge potential is

```text
mu = G_solution_reference + G_protonation + (q - q_ref) * G_Hplus
dG = sum(coefficient * mu) + input_n_Hplus * G_Hplus
```

These potentials are not absolute experimental formation energies. Elemental
reference gauges cancel in balanced reactions. The explicit proton terms make
the reference conversion auditable. The chemical accuracy of the class pKas,
aqueous proton convention and solvent model requires separate validation.

Hydration is evaluated in a molecule-local context using the same retained raw
energies and water value. Incomplete hydration fails explicitly. No adjustment
to a reaction result is made to force cycle closure.

## Supported scope

The profile supports specified individual molecular inputs, automatic
source-registry acid models and explicitly declared reducing-sugar pools.
For a pool, each required ring/open-chain state receives its own protonation
and hydration treatment before the complete partition is formed. Its identity
depends on the declared chemical state set rather than the starting drawing or
the source label. The source declaration remains attached to each reaction.
An unavailable constituent refuses the whole pool; retained states are never
renormalized to rescue an incomplete estimate.

The profile refuses flat caller pKa terms and caller protonation
models. It also refuses unresolved stereochemistry outside an admitted pool,
isotopes, unbalanced reactions and failed reference or molecular calculations.
It must not silently substitute a single sugar drawing for a declared solution
pool. Full-molecule calculations can be substantially more expensive than the
established reduced-molecule workflow.

`dG_unrounded` is retained for closure and reversal checks. `dG` is rounded to
0.1 kJ/mol. Sampling derivatives propagate through both hydration and ring/chain
partitions. Shared raw molecular estimates are combined before quadrature, even
when a fixed state is also a constituent of a pool. Distinct raw estimates are
assumed independent. Physical
model error, missing states and parameter uncertainty are not calibrated by
`U_samp`; it is not a total prediction interval.

## Verification and promotion

The first retained-record development check gives MAE 16.30 kJ/mol on 111
available labels, versus 12.12 for zero. On 109 labels shared with the
reduction-based estimator, the corresponding MAEs are 16.14 and 7.59 kJ/mol.
These conditional reassemblies do not establish fresh-run coverage or independent
accuracy. Internal consistency is insufficient grounds to promote this profile
over the established scorer.

Run the focused tests with
`python -m pytest tests/unit/test_state_function.py tests/unit/test_state_function_pools.py`.
The all-input synthetic audit accepts a separately retained input-state cycle
basis:

```bash
python benchmark/audit_state_function.py --cycles input_nullspace.json --output audit.json
```

The audit verifies original-state cancellation before evaluating each cycle. Its
synthetic values include structural terms and changing repeated-call values to
exercise retained-state behaviour. These checks establish bookkeeping properties,
not experimental accuracy or physical calculation coverage.

`benchmark/verify_state_function_ledger.py` independently checks v2 retained
component arithmetic, chemical/proton bookkeeping, complete declared partitions,
shared sampling derivatives and reaction sums. With the original raw directory,
index and intervention protocol it also authenticates molecular records and
reconstructs the declared component changes. This is stronger than cycle closure
alone: its tests reject opposite reaction offsets that preserve a cycle while
disagreeing with the molecular ledger.

```bash
python benchmark/verify_state_function_ledger.py STUDY/full.json \
  --raw-root ORIGINAL/shared_record_ledger \
  --raw-index ORIGINAL/shared_record_index.json \
  --protocol STUDY/protocol.json --cycles ORIGINAL/input_nullspace.json
```

Promotion requires owned physical checks of the reference treatment and
declared state-pool populations, full-molecule resource and
coverage measurements, and evaluation with the frozen experimental protocol.
The pre-enrollment design and candidate inventory are in
`benchmark/state_function_validation/`. They currently admit no experimental
reactions and are explicitly not a frozen prediction benchmark.
