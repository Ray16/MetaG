# Additional equilibrium-source review

These reviews update the independent-enrollment source inventory.
**No experimental reference is admitted.** Existing development labels and
predictions are unchanged. Eligibility was not selected using model errors.

McCoy et al. (2006), [DOI 10.1021/bi052232m](https://doi.org/10.1021/bi052232m),
measured ADP-glucose phosphorolysis equilibrium. The two constants, 5 and 40,
appear as their reciprocals in the NIST compilation and reproduce both underlying
observations in row 331 of the released 671-row dGbyG training table. The same
reaction, rxn12158, also occurs in MetaG's historical development inputs.
It is excluded from an independent cohort. The primary article reports 1 mM
added MgCl2 for the second observation, whereas
[NIST Part 7, section 4.46, page 1367](https://srd.nist.gov/jpcrdreprint/1.2789450.pdf)
reports 0.1 mM. This discrepancy remains unresolved. Added salt is not a measured
free-magnesium activity, and buffer concentration does not specify ionic strength.

Ajalla Aleixo et al. (2019),
[DOI 10.1111/febs.14782](https://doi.org/10.1111/febs.14782), report a fumarase
equilibrium ratio of 4.3 ± 0.2 from a reversible initial-rate fit. Rounded Table 2
parameters reconstruct 4.25 using the Haldane relation. This authenticates the
quantity's equilibrium interpretation, but the assay temperature is described
only as room temperature. The native pH is 8.5 in 50 mM Tris-HCl and 150 mM KCl;
a transfer to MetaG's fixed target and measurement/training ancestry remain
unresolved. The companion thesis describes the same fit and is not counted again.

The JSON files retain source locators, condition gaps, overlap witnesses and
arithmetic. `arithmetic_receipt.json` checks the stated numbers against the owned
training snapshot; it does not independently authenticate the linked articles.
Neither review establishes an unfamiliar-chemistry advantage or changes the
number of independently admitted primary references from zero.

[The later source follow-ups](source_followups_20261105.json) retain the newly
owned Walford (2000) aconitate and Westarp (2025) nucleoside theses, their exact
source identities and the reviewed assay conditions. Aconitate ratios derived
at 70–110 C are not 298.15 K labels. Its opposite geometric drawing is absent
from the concrete training tables and has Morgan similarity 0.44, but shares the
conservative stereo-erased family and identical default comparator features
with the training cis drawing. Chemical distance and scored feature novelty
must therefore be reported separately.

The cumulative nucleoside thesis repeats several published experiments and
does not supply independent replications of them. It distinguishes the 2 mM
phosphate equilibrium-estimation recipe from a 0.5 mM preparative recipe;
post-quench methanol is not a reaction cosolvent. Some product constants use
borrowed donor calibrations. The linked raw-data repository explicitly requires
author access, and condition transfer, covariance and ancestry remain unresolved.
The source records also retain an unsuccessful hydantoin full-text access attempt
and flag a new aldolase lead's eQuilibrator-derived constant as a model result.
No source hold, benchmark gate or reference-admission count changes.

The [neutral thermochemical lineage review](neutral_thermochemical_lineage.json)
checks a separate aqueous formation-energy route. The source myo-inositol has
exactly the same full stereochemical InChIKey as both a comparator training
participant and the historical MetaG input named `L-Inositol`. A newer crystal
measurement therefore cannot establish unfamiliar-molecule performance. Its
conversion to an aqueous reference also requires matched dissolution activity,
crystal identity and uncertainty; no aqueous label is constructed here.

The review also traces two long-chain alkylbenzene formation energies to the
specific homologous-series correlation cited in Shock and Helgeson's Table 6.
Those entries cannot supply independent experimental labels for the proposed
four-hydrocarbon exchange. Other entries in that table have different ancestry;
the use of an equation of state alone does not disqualify every reference-state
value. All source decisions were made without running thermodynamic predictors
or inspecting their errors. The primary reference count remains zero.

The [raw aconitate follow-up](../aconitate_raw_source_review/README.md) retains
all 19 measured condition runs, 314 timestamps, five blank initial cis fields,
run-specific calibrants and explicit OCR alternatives. Closed-concentration and
common-recovery ratio fits differ most at 110 C, where apparent total
concentration declines by 7–24%. Arithmetic replay and single-CPU reoptimization
are available; neither establishes a chemical loss mechanism, a calibrated
reference interval or a valid transfer to 298.15 K. No primary label is admitted.

`additional_sources.json` records the six remaining original candidates and two
new leads, with explicit distinctions between owned primary evidence, indexed
publisher methods, and abstract/metadata-only inspection. These access levels
are not interchangeable. The original validation design and launch prohibition
remain unchanged.

| Source | What is established | What prevents admission |
|---|---|---|
| [Kratzer and Nidetzky, 2005](https://pmc.ncbi.nlm.nih.gov/articles/PMC1175129/) | A Haldane reconstruction gives 420.94, rounding to the reported wild-type 421. The mutant estimate is 153. | Consistency/uncertainty, condition transfer and measurement ancestry remain unresolved; mutant rows are not independent reactions. |
| [Nasution et al., 2006](https://doi.org/10.1002/bit.20842) | The abstract describes intracellular metabolite quantification. | Full methods and the equilibrium meaning/ancestry of the indexed constants are unverified. |
| [Barber et al., 2006](https://doi.org/10.1074/jbc.M512727200) | Publisher methods describe an endpoint ratio after repeated fresh-enzyme additions, at 25 C and pH 7.4. | Numerical results/uncertainty, condition conversion and measurement ancestry are unverified. Derived Michaelis constants do not imply that the equilibrium ratio was derived from them. |
| [Genda et al., 2006](https://www.jstage.jst.go.jp/article/bbb/70/5/70_5_1102/_pdf) | All six Table 2 Haldane reconstructions round to the reported values. Native assays used 20–23 C. | Temperature/activity transfer and ancestry remain unresolved; buffer concentrations are not ionic strengths. |
| [Schumacher et al., 2006](https://doi.org/10.1002/biot.200600039) | The abstract includes mixed-solvent kinetics and a buffer control. | Full methods, equilibrium interpretation and row-specific conditions are unverified. |
| [Meng et al., 2008](https://doi.org/10.1016/j.bbapap.2008.02.021) | The abstract reports an approximate ratio of 0.3 for two isoforms. | Definition, method, uncertainty and Mg/activity conditions need primary verification; indexed conditions are 21 C and pH 7.5. |
| [Meersseman Arango et al., 2025](https://doi.org/10.1021/cbe.4c00186) | Owned primary XML specifies 35 C, pH 8 and 10% methanol for the standard reaction. | Mixed-solvent/temperature conversion and the native ratio's definition and uncertainty remain unresolved. Product-crystallization yields are not aqueous equilibrium constants. |
| [Wu et al., 2023](https://doi.org/10.1002/ijch.202200065) | The abstract confirms interconversion of azido-labelled UDP sugars. | Full results, Supporting Information, chemical states, conditions and ancestry remain unverified; azide presence alone does not establish unfamiliarity. |

Post-quench detection conditions in the epimerase study must not be assigned to
the preceding reaction. Likewise, the 23–25 C stability incubation in Genda et al.
must not replace its 20–23 C kinetic assay. Named reaction overlap with a training
table is recorded separately from proven measurement duplication. A different
paper, enzyme variant or database entry does not by itself establish independence.

No predictions were run or consulted in these decisions. All eleven inventory
entries are excluded, held, or pending further source verification; none is an
independently admitted standard-Gibbs reference. The primary hypothesis, novelty
rules, condition requirements and minimum cluster gate have not been relaxed.
