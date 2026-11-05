# Additional equilibrium-source review

These two reviews resolve entries in the independent-enrollment source inventory.
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
