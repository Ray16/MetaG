# Global amino-acid reference candidate and minimum-recovery control

This packet preserves two separate results: a bounded numerical minimum-recovery
fix, and a source-coupled global-reference prototype that was **not adopted**.
The reference comparisons are development-exposed and use condition-modeled
historical labels. No independent experimental label or genome-scale numerical
coverage claim is established.

The candidate selects one fully protonated NH3+/COOH reference per exact molecule
for glycine, alanine, aspartate and glutamate. Seven explicit stereochemical
references use the previously reviewed Miller/Smith-Magowan ladders and their
declared conditional concentration-scale conversions. Source parameters are fixed
before this candidate's error analysis. No reaction name, direction, neighbouring
compound or residual chooses a reference, and no reaction-level hook is patched.
Unmatched molecules retain the full-molecule reference policy. This differs from
the earlier coupled/reduced 173-label study, which remains a separate result.

The original 589-record full-molecule assembly reproduces all 140 earlier
available estimates and the 111-label MAE of 16.30374 kJ/mol. Five existing raw
extensions and two fresh source-reference parents then provide an explicit
596-record bank. All records and their provenance are retained in
`raw_records.zip`; none is chosen by a lower energy or experimental error.

One of the new parents, S-alanine, initially failed the original strict thermal
workflow. A force-converged geometry retained a 78.5 cm-1 imaginary mode. Three
tighter FIRE controls removed the strong mode but did not reach their 0.001 eV/A
force target. Separate LBFGS controls from two saved starts reached forces of
0.00061 and 0.00074 eV/A and passed the unchanged 50 cm-1 strong-mode gate.
Their full Hessians are retained; independent QR projection and ideal-gas thermal
replay reproduce the corrections within 0.0001 kJ/mol.

The bounded curvature recovery waits until the existing final FIRE/mode-following
path fails. It caps LBFGS at 400 steps, checks a 180-second budget between model
evaluations, and uses a 0.001 eV/A force
criterion, identity checks and a freshly evaluated strict Hessian. Electronic,
ALPB and thermal values must all describe the returned geometry. No strong
negative mode is softened to obtain an estimate. A complete fresh S-alanine
calculation with this separately fingerprinted source passes all sampled-minimum
checks; its original refusal remains unchanged in the control archive.

The full retained-record comparison gives:

| Arm | Available assemblies / 364 requests | Labels available | MAE on those labels, kJ/mol |
|---|---:|---:|---:|
| Full reference, original 589 records | 140 | 111 | 16.30 |
| Full reference, extended bank | 143 | 112 | 16.68 |
| Global amino reference, completed bank | 148 | 117 | 16.50 |

The last two aggregate MAEs use different targets. On their **112 identical
labels**, MAE increases from **16.68 to 17.11**, compared with 12.15 for zero.
On the seven common targets containing a matched source amino family, it
increases from 15.10 to 21.97. The other 105 common predictions are unchanged.
No source-reference policy is promoted from the changed availability or apparent
aggregate improvement. Missing archived calculations are distinguished from
workflow refusals; these bank-limited counts are not physical applicability rates.

All 592 reversal/scaling/atom-order/species-order checks on 148 estimates pass.
All 44 modeled protonation drawings across seven references have consistent
potentials after their explicit proton exchanges. Complete available-subnetwork
nullspaces contain 11 baseline and 12 candidate cycles, all closing numerically.
These are state-function checks, not accuracy evidence. The independently sampled
alanine enantiomer ensembles give a 2.88 kJ/mol racemization residual despite
only 0.017 kJ/mol reported sampling spread. No numerical enantiomer equality was
imposed; this illustrates why that spread is not total predictive uncertainty.

All 32 previously declared source-scale scenarios are retained without selecting
one by error. Their largest reaction range is 0.029 kJ/mol, conditional on those
scale alternatives; it does not bound missing source or physical uncertainties.
Changing molecular reference and protonation polynomial together also prevents
assigning the accuracy change uniquely to a single component.

## Offline reproduction

From the repository root:

```bash
python benchmark/global_amino_reference.py
python benchmark/global_amino_reference.py --cold-replay --physical-replay
```

The basic audit uses the standard library and authenticates every raw record,
reconstructs every molecular partition and recalculates common-target MAEs.
Cold reassembly additionally uses NumPy, RDKit and ASE with the archived assembly
source. A fresh isolated process forbids imports of Torch, FairChem and other
quantum runtimes. Physical replay independently reconstructs source-family
potentials with decimal arithmetic and projected Hessian thermochemistry.
No GPU or model weights are needed for either replay.

`assembly_source.zip` is the frozen study assembly implementation;
`recovery_source.zip` owns the distinct numerical-recovery producer. The recorded
raw profiles differ in source identity while retaining the same model/checkpoint,
settings and component definitions. Original accepted records are not relabelled.
Historical cluster launchers and all failed/partial control outcomes are retained
for provenance; their paths are not a portable GPU environment. This packet
contains numerical model outputs and source code, not UMA weights.
