# Twelve source-curated ModelSEED boundary reactions

This numerical pilot evaluates three fixed state-function estimators: full molecules,
the earlier 1,4-only nicotinamide anchor, and separate 1,2/1,4 methyl cores sharing one
exact full oxidized reference per scaffold. Molecular structures were chosen from
the ChEBI cross-references already present in the source database. All 18 participating
source identities and exact proton rebasing were fixed before new calculations.

The initial 14 molecular attempts yielded 9 complete wrapper records, 4 quality refusals,
and one wrapper cache-key error. Triplet oxygen's existing complete record was
independently verified: the wrapper had omitted its spin=3 cache-key field. No
oxygen calculation or original outcome was overwritten. The physical first-pass
count is therefore 10 accepted records and 4 quality refusals.

All 4 refusals reached a 400-step curvature-recovery limit without satisfying its
unchanged 0.001 eV/Angstrom force criterion, with identity preserved. A separately
frozen numerical profile raised only that limit to 1200 steps and recovered 3 records;
full 1,2-NADP remained unavailable. The model, solvent, Hessian/force acceptance
criteria, stochastic seeds and 96-candidate ceiling were unchanged. Accepted
legacy 64-ceiling records retain their owning source/settings; none is relabelled.

Primary availability is 5/12, 5/12, 8/12; after the separate numerical extension it is
11/12, 11/12, 12/12. The complete available stoichiometric nullspace has one cycle per
arm in the augmented panel. All 136 augmented reversal/extent/atom-order/species-order
checks pass. These are targeted pilot results, not genome-scale coverage.

For NADH 1,4-to-1,2 isomerization, the modeled energies are +12.7246 kJ/mol (full),
-29.1694 (earlier mixed approximation),+14.4842 (two-isomer anchor). The latter two
approximations remain cycle consistent. Agreement within 1.76 kJ/mol of the full
model is an internal model comparison, not experimental accuracy. The local
approximation also imposes equal NAD/NADP isomer-core differences, and can miss
scaffold-specific effects. Sampling spreads exclude these physical-model errors.
No experimental labels are added, and no production policy is selected.

From the repository root, using its audit environment:

```sh
python benchmark/boundary_reference.py --cold-replay
```

The data audit reconstructs all raw partitions, declared protonation terms,
reaction sums and correlated sampling propagation with 50-digit decimal arithmetic.
A cold process replays all 12 requests under all 3 arms and both availability profiles,
with quantum imports/network forbidden. Four stored fixed-geometry Hessians support
an independent QR/RRHO numerical check; this is not an all-minimum Hessian replay,
proof of exact stationarity, complete aqueous sampling or experimental validation.
Retained-geometry force residuals are reported and are not uniformly 0.001 eV/Angstrom.

All 25 accepted primitive records, initial/secondary outcomes, source code and logs
are included. Raw per-minimum thermal values are model outputs; independently
regenerating every quantum/Hessian calculation additionally requires the pinned
FairChem SDK, licensed UMA weights and xTB. Original component references, their
experimental-hydration water value and each producer identity remain explicit.
No model weights are redistributed. The original 8 GPU-hour common deadline was
retained across numerical profiles. Molecular-stage elapsed times exclude model
startup and are not a whole-method speedup benchmark.
