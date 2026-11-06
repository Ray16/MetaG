# Aconitate raw-source reanalysis

This is reference construction and diagnostic evidence. **No independent primary
Gibbs-energy reference is admitted**, no MetaG/dGbyG energy is predicted, and no
production policy or primary benchmark gate changes.

The source is S. N. Walford, *A study of the isomerisation of aconitic acid with
reference to sugarcane processing* (University of Natal, 2000),
[primary thesis](https://researchspace.ukzn.ac.za/bitstreams/3e8dad8b-266e-4195-9df6-cc6bf4c43ed0/download).
Appendix 6, Tables A6.12–A6.16, printed pp.125–126/PDF pp.140–141, contains
HPLC peak areas and run-specific external standards. The response-factor method
is described on printed p.45/PDF p.60. The complete primary PDF is owned locally
and SHA256-pinned in source_provenance.json; it is not redistributed here.

All 20 temperature/buffer requests are retained: 19 measured runs and the
source-declared missing 25 mM/110 C cell. The 314 timestamps contain 309 cis and
314 trans peak areas. Five blank initial cis fields remain missing, not zero or
inferred detection bounds. There are 1,018 intended numeric field locations,
including these blanks. Hidden PDF OCR and a separate Tesseract image pass agree
at 913 locations. All 105 disagreements were image-reviewed; ten retain explicit
digit alternatives. The correction ledger and small numeric-cell image panels
preserve those decisions. Repeated conditions and timepoints are one reaction
and one experimental source, not hundreds of independent reference labels.

Calibrated concentrations use each run's own standard area and concentration.
The main model fits closed, reversible two-state relaxation to both channels.
Separate-channel and zero-initial-cis fits are sensitivities. Native measured
concentration totals decrease by about 7–24% in the three 110 C runs, compared
with much smaller drift below 110 C. This is an apparent concentration/recovery
change; its chemical or analytical cause is not established by these data.

A postbaseline exploratory sensitivity fits cis/(cis+trans) at all 309 paired
timestamps; the five unpaired trans observations remain in the main data/fits.
Common multiplicative recovery factors cancel in this ratio; selective loss does
not. At 110 C, closed-model K values are about 0.19–0.24 and ratio-fit values about
0.16–0.17. Neither arm is selected as the reference. A synthetic selective-loss
control has true underlying k_forward/k_reverse=0.15 but gives ratio-fit K≈0.072,
showing why better fit or stable ratios alone cannot establish the reaction K.

All 81 fixed-K profiles per measured condition, their nuisance parameters,
transcription variants and adaptive threshold crossings are retained. Fixed-rate
linear nuisances use a nonnegative active-set solution followed by scalar rate
optimization. These working-IID Gaussian likelihood-ratio diagnostics are **not
authenticated experimental confidence intervals**: calibration, unequal/correlated
errors, physical model adequacy, native pH(T), activity/speciation and transfer
from 70–110 C to 298.15 K remain unresolved. Lower-bound calibration flags record
values below the source's lowest tested calibration point; they are not detection
limits or automatic exclusions.

## Replay

Run `python audit_aconitate_raw.py --packet .` from this directory, or
`python benchmark/audit_aconitate_raw.py` from the MetaG checkout. The standard-
library checker authenticates the frozen files and reconstructs concentrations,
missing-data bookkeeping, all stored curve/profile sums of squares, linear
nuisance stationarity, alternative-digit fits and ratio summaries. It does not
re-optimize the source experiment or certify reference eligibility.

The exact fitting, profile and synthetic-control generators are retained.
Their original CPU-affinity/thread settings and library versions are in the
protocol and results. Work on a scratch copy when regenerating their outputs;
compressed JSON files contain the exact original JSON bytes. The OCR assets are
recognition aids, not thermodynamic model training. No model weights or quantum
calculations are needed for the arithmetic replay. The supplied numeric-cell
panels are inspection evidence; hashes and arithmetic alone cannot establish
that an experiment or transcription is correct.

`python refit_aconitate_source.py --packet . --profiles` additionally reruns all19
main fits and all1,539profile optimizations using NumPy/SciPy, on one available
CPU with native thread caps of1. It reads the frozen inputs and does not overwrite
source records. This is algorithm reexecution, not independent experiment validation.

Three coarse-grid supports are empty because the 81-point grid misses the narrow
working-model support region. Adaptive threshold crossings are finite in all19
cases; an empty coarse grid is not evidence of non-identifiability. Those
working-noise support endpoints still do not establish experimental uncertainty.
