# Nucleoside reference-source review

This review authenticates source arithmetic and records eligibility problems.
It admits **no independent standard Gibbs-energy reference**, runs no MetaG or
dGbyG prediction, and leaves the existing validation design unchanged.

## Retained data and reconstruction

The unchanged `source_data.zip` comes from Kaspar et al.,
[Thermodynamic Reaction Control of Nucleoside Phosphorolysis](https://doi.org/10.1002/adsc.201901230),
with data at [Zenodo record3530526](https://doi.org/10.5281/zenodo.3530526).
The authors released the archive under
[CC BY4.0](https://creativecommons.org/licenses/by/4.0/); attribution, upstream
checksum and original/projection identities are in `source_ownership.json`.

All 576 raw observations are retained: 24 compounds, four temperatures, two
reaction mixtures per condition and three equilibrium aliquots per mixture.
The compiled worksheet includes 519 points and omits 57. Every included product
concentration agrees with its original spectral fraction within
5.2e-15mM; every included cached equilibrium constant reconstructs within
1.9e-14. The source calculation is

$$
K_{c,\mathrm{app}}=\frac{x^2}{(1-x)(r-x)},\qquad r=P_0/N_0=5,
$$

where `x` is fractional conversion and both products initially vanish. The
two-solute-to-two-solute stoichiometry makes concentration units cancel. It does
not remove activity coefficients or protonation and metal-binding effects.

Three aliquots from a reaction do not constitute three independent equilibrations.
The retained compilation mask authenticates inclusion, not the physical reason
for every exclusion. Raw observations omitted by the authors remain visible.

The worksheet duplicates a compound heading at `K145`; the raw identifier and
supplementary TableS1 distinguish compound20, the deoxy analogue. Exact molecular
structures still need separate authentication. A source summary cell, `R67`,
contains `#DIV/0!` for an omitted condition. That error is preserved and does not
enter the 519 authenticated constants. The upstream workbook is unchanged.

## Remaining reference requirements

These are buffered apparent concentration ratios. They have not been converted
to activity-based standard Gibbs energies. Actual pH at each measurement
temperature, ionic strength, magnesium convention, chemical-state pools and
reference uncertainties remain unresolved here. An extrapolation of the authors'
apparent temperature fit would not by itself establish the primary benchmark's
298.15K, pH7, zero-ionic-strength, Mg-free reference. Repeated observations also
do not supply additional experimental-source clusters.

`additional_source_review.json` records this dataset and a second lead,
[Westarp et al., Nucleoside Phosphorylases make N7-xanthosine](https://doi.org/10.1038/s41467-024-47287-4).
The latter reports an approximately21C native equilibrium. Its temperature,
speciation and measurement basis must be owned before constructing a matched
N7/N9 comparison. Neither source is enrolled in the primary benchmark.

`legacy_candidate_adjudication.json` records review of the four distant-case
leads from the earlier feasibility audit. It identifies repeated imidazole
measurements, reaction orientation, an explicit xanthosine-product anomer, and
missing condition bridges. It preserves all original training targets.

## Replay

From the MetaG checkout, using only Python's standard library:

```bash
python benchmark/audit_nucleoside_reference.py \
  --source-zip benchmark/nucleoside_reference_review/source_data.zip \
  --output /tmp/nucleoside_source_audit.json
```

The result should reproduce `arithmetic_audit.json` with the recorded auditor
source. `manifest.json` authenticates all retained review files except itself.
No spreadsheet formulas are executed. Missing or inconsistent included source
values cause failure, and missing observations never become zero energies.
