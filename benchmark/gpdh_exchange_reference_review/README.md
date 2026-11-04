# GPDH exchange reference review

The candidate equation is (R)-glycerol 3-phosphate + glycolaldehyde -> DHAP +
ethylene glycol. Common enzyme-bound cofactor symbols cancel between two
reported half reactions. A genuine thermodynamic interpretation additionally
requires common conditions and compatible chemical states; cancellation alone
does not establish it. **No primary or secondary thermodynamic label is admitted.**

The sources are Cristobal et al. (2024),
[DOI 10.1021/acs.biochem.3c00702](https://doi.org/10.1021/acs.biochem.3c00702),
and the reverse-kinetic measurements of Reyes et al. (2015),
[DOI 10.1021/ja5123842](https://doi.org/10.1021/ja5123842).
The original Scheme 6 image, Table 1 and rate arithmetic support K1=0.019.
The conflicting prose value 0.19 is preserved in the ledger and not used.

Both reverse constants are reported using free-carbonyl concentrations.
Using the stated GA and DHAP normalization fractions 0.06 and 0.45 gives:

| Nominal source reconstruction | Ratio | -RT ln(ratio), kJ/mol |
|---|---:|---:|
| Free-carbonyl basis, reported rate constants | 27.63975 | -8.22829 |
| Restored total-substrate concentration basis | 3.68530 | -3.23343 |
| Alternative DHAP scaling of 0.55 | 3.01525 | -2.73598 |

The 4.99486 kJ/mol difference between the first two rows is a state-definition
change, not an empirical correction to a prediction. Undoing a concentration-axis
normalization does not independently measure the true hydration population.
Rate and normalization uncertainties share dependencies and cannot be propagated
as independent measurements.

DHAP's normalization history remains ambiguous. The 2015 article states 0.45;
Tsang et al. (2008), [DOI 10.1021/bi8001743](https://doi.org/10.1021/bi8001743),
used 0.55. The cited primary DHAP study, Reynolds et al. (1971),
[DOI 10.1042/bj1220285](https://doi.org/10.1042/bj1220285), reports 0.55 keto at
20 C and 0.83 at 37 C. Neither number establishes a measured 25 C fraction for
the later assay. The two decoding scenarios are not a confidence interval.

Forward assays used 20 mM TEA and reverse assays 10 mM TEA, despite matching
25 C, pH 7.5 and ionic strength 0.12 M. Differential buffer effects are not bounded.
The complete rate covariance and a source-owned transfer to MetaG's pH 7,
zero-ionic-strength target are also unavailable. The numerical ratios above
therefore remain conditional reconstructions, not benchmark Gibbs-energy labels.
No prediction or residual was consulted to choose between them.

`training_overlap.json` is an outcome-free check against the owned **671-row**
structural projection. Three compound families have resolved training matches.
Ethylene glycol has no resolved pool match and nearest known-pool Morgan
similarity 4/9, but unresolved training aliases and the separate 675-row snapshot
prevent a complete novelty claim. The audit does not establish absence from UMA.

`candidate_request.json` preserves the chemically balanced equation, the primary
article's explicit R stereochemistry and its native pH/I request. Both calculation
APIs now withhold this request before molecular work because those conditions
are unsupported. They do not return a pH 7 estimate as a pH 7.5 prediction.
`compatibility_audit.json` records the source comparison: the packaged 364 inputs
omit condition requests, and removing only the new guards and docstrings leaves
the preceding runtime ASTs unchanged. This is not a new molecular sweep or a
validation of physical accuracy.

Reproduce the arithmetic with standard-library Python:

```bash
python benchmark/gpdh_exchange_reference.py
```

For authenticated primary-source extraction, download the two unchanged XML files
at the URLs in `source_facts.json`, save them with the recorded filenames, then run:

```bash
python benchmark/gpdh_exchange_reference.py --primary-dir /path/to/primary_xml
```

The original articles are linked rather than redistributed. Default replay checks
the curated numerical facts and arithmetic; it does not claim to have authenticated
unavailable original files. The local study did perform the authenticated extraction.
