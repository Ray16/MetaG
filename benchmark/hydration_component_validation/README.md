# Fixed-geometry hydration component ablation

This development-exposed experiment tests an explicit component decision against
experimental air/water partition data. It does not validate the complete MetaG
workflow or establish independent metabolic reaction accuracy.

The original measurement source is Endo, Hammer and Matsuzawa (2023),
[DOI 10.1021/acs.est.3c02545](https://doi.org/10.1021/acs.est.3c02545).
External conformer geometries and calculated energy components come from
Wittmann, Selzer and Grimme's FlexiSol study,
[DOI 10.1039/D5SC06406F](https://doi.org/10.1039/D5SC06406F), at repository
commit `7b44798f26c888ef541faa6143a813136921483f`.
FlexiSol already studies geometry, ensemble and electronic-level effects; this
experiment extends that component comparison using UMA single points.

All 21 aqueous entries associated with the original study were considered.
The frozen source rule holds two different covalent structures and six acid-medium
entries, leaving 13 nominal-water compounds. The unchanged external 6:2 FTAC and
4:2 FTMAC structures each contain one additional CH2 group compared with the
primary SI. Their geometries and the complete eligibility ledger remain here;
their experimental labels are not scored against a different molecule.
The six 1 mM HCl samples require an additional intrinsic-neutral activity/speciation
bridge. PFHxSA and PFOSA also have commercial isomer-mixture qualifications.
The singly stereogenic secondary FTOH is included because enantiomers have the
same intrinsic hydration free energy in an achiral continuum. No diastereomer
transfer is assumed.

The 13-compound cohort and all 461 gas/solution conformers were fixed before the
new UMA calculations and numerical error analysis. Every requested single point
completed. UMA uses the official FP32 OMol inference path, charge zero and singlet
multiplicity. No geometries were optimized and no Hessians were computed.
Three single-versus-batch checks differ by at most 0.000075 kJ/mol. A separate
native xTB replay reproduces one published molar ALPB increment within
0.00000002 kJ/mol. This checks that increment's numerical convention, not accuracy.

Targets are reconstructed from the original Table 2 cycle column as
`R*T*ln(10)*log10(Kaw)`, with `Kaw=c_air/c_water`, 298.15 K and a common molar
concentration reference. No gas-pressure conversion or fitted offset is added.
Nominal-water experiments can retain a small acetone fraction, and source purity,
finite concentration and cycle-covariance qualifications remain. These are not
new strict activity-standard reference admissions.

The primary comparison uses separate minima of the phase-specific energy lists.
Published GFN2 ALPB and r2SCAN-3c SMD increments are each combined with UMA,
r2SCAN-3c or omegaB97M-V/def2-TZVPPD electronic energies. These SMD increments
were not recomputed with UMA or the hybrid density. Gas-only and solution-only
geometry lists, the source's Boltzmann energy mean, and a finite-state
log-partition sensitivity are all retained. There are no extra RRHO terms.
The Boltzmann energy mean is not an ensemble Gibbs free energy. The log-partition
arm has unit conformer weights and incomplete sampling/degeneracy information;
it is not promoted as a complete thermodynamic ensemble.

| Primary fixed-geometry arm | MAE (kJ/mol) |
|---|---:|
| UMA + ALPB | 13.24 |
| UMA + SMD | 4.32 |
| omegaB97M-V/def2-TZVPPD + ALPB | 13.18 |
| omegaB97M-V/def2-TZVPPD + SMD | 4.41 |
| Zero hydration baseline | 9.09 |

SMD improves all 13 paired UMA errors. Its mean absolute-error change is
-8.92 kJ/mol; the descriptive compound-bootstrap interval is [-10.73,-7.29].
Replacing UMA with the published hybrid DFT values while retaining ALPB changes
MAE by -0.06 kJ/mol, with interval [-0.17,0.06]. All 54 specified scoring arms,
per-compound errors and the bootstrap seed are supplied without error clipping.
Compounds share one source and related chemistry: these intervals do not include
experimental cycle covariance or demonstrate source-independent superiority.

The original experiment was previously consulted during MetaG development and
used in a different, source-corrected 14-compound transfer-range panel.
`source_exposure_addendum.json` records that discovery without altering the frozen
protocol, selection or predictions. This new external-geometry panel is not a
replacement for that panel and supplies zero independent primary metabolic labels.
Training-chemistry overlap is not established. No production method is changed.

## Offline reproduction

From the repository root, with Python 3.11 or later:

```bash
python benchmark/hydration_component.py --output /tmp/hydration_component.json
python -m pytest -q tests/unit/test_hydration_component.py
```

The first command uses only the standard library. The tests also use RDKit and
pytest, and independently reconstruct the primary minima with 50-digit decimal
arithmetic. The packet's copy of `hydration_component.py` can run without the
repository by specifying `--packet /path/to/this/folder`.
The manifest authenticates every retained file. The original producer/launcher
files record historical cluster paths; they do not provide a portable GPU
environment or redistribute UMA weights.

`selected_upstream_energy_rows.csv` is an unchanged-value projection of the full
upstream CSV; `inputs.json` retains ordered coordinates and source components.
The full upstream archive and energy CSV can be retrieved at the pinned URLs in
the retrieval receipts to repeat source preparation. They are not duplicated in
this packet. `Endo2023_SI.txt` is a text extraction; the unchanged PDF and XML are
retained as the authoritative sources. Primary-source content is CC BY 4.0;
the upstream repository's MIT license is retained in `upstream/LICENSE.txt`.
The FlexiSol paper itself is not redistributed here.
