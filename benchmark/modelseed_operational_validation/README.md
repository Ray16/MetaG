# ModelSEED operational probability sample: pre-computation freeze

This is a label-free operational study of numerical availability and cost. It is
separate from the independent experimental benchmark, whose admission and success
rules remain unchanged. No experimental Gibbs-energy labels, predictions or
accuracy comparisons select the sample.

The source denominator is all 43,650 live nontransport source records at pinned
ModelSEED revision b539adb69a955dac15cc30dbeb3debbd87739018. The sampling frame is
the 9,719 passes under the existing **four-cofactor** curation profile. Four IDs
were drawn without replacement from each of eight nonempty EC strata, using a
recorded random seed fixed before sampling. All 32 selected IDs remain in the
denominator. Source IDs, including duplicate chemical equations, are the units;
these are not 32 independent experimental or transformation-family clusters.

The full-molecule state-function and two-regioisomer nicotinamide approximations
are fixed arms. They target 298.15 K, pH 7, I=0, no Mg binding, 1 M solutes and unit
water activity. This does not authenticate the completeness or physical accuracy
of each source annotation or protonation model. Original protonation warnings
remain. No production reference policy is selected.

Structural prerequisite discovery exports no synthetic energies. It finds 117
primitive inputs. One reaction, rxn23598, exceeds the predeclared 100-total-atom
capacity bound and remains unavailable under both arms. The remaining arms need
115 fresh raw tasks, including one liquid-water reference. No earlier molecular
energy is reused. Each unique primitive is calculated once and shared across the
batch and arms; whole-molecule source references remain distinct.

Numerical settings retain the model/solvent/force/Hessian/identity criteria and
use the previously tested 1200-step curvature-recovery limit and 96-candidate
validation ceiling. Each raw input receives one execution opportunity of at most
1800 seconds. Two H100 workers share an eight-GPU-hour ceiling and one common
deadline, fixed at launch. All failures, capacity refusals and budget-censored
inputs remain visible. No resampling or post-outcome rescue replaces them.

Availability estimates require known stratum weights. Intervals invert exact
finite-population hypergeometric tests with a Bonferroni allocation across strata.
The full-source denominator includes known structural failures. **If the common
budget truncates any required opportunity, population numerical-availability
inference is withheld.** The small sample may give broad intervals; no predefined
minimum-coverage success or experimental reliability claim is made here.

Initialization, per-primitive and batch walltimes are reported separately. Shared
cache work is not standalone reaction runtime, and no genome-scale walltime is
extrapolated. Archive/ledger checks do not establish experimental accuracy,
complete aqueous ensembles or calibrated total prediction intervals.

`design.json` preserves the pre-selection rule and seed. `protocol.json`,
`freeze.json`, the complete sample and exact raw prerequisites were recorded
before target quantum outputs. The later code guard only enforces a boolean
coverage-completeness flag; sampling is unchanged and replays identically.
Frozen source, worker adapter tests, fixture coordinates and upstream notices are
included. Model weights and the separately obtained FairChem SDK are not
redistributed. Provenance and hardware checks are separate from study completion.

This packet initially contains **no physical outputs**. Subsequent outcome and
verification records must retain this pre-computation manifest and the original
sample; the design alone is not a completed numerical benchmark.
