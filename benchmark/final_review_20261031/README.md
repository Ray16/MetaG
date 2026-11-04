# Final review: source identity and bounded input correction

The molecular sweep and final package have separate full source identities. `independent_arithmetic.json` records both and the unchanged raw-molecular identity. `evaluated/` preserves the two original files changed in this pass. The sole runtime change rejects atom-map annotations before routing; the pipeline module documentation was also corrected. No evaluated input contains these annotations.

Use `verify_source_change.py` as documented in its module header after extracting the code and current-study archives together. The verification reads the actual full source trees and all 364 inputs; it neither recalculates energies nor promotes an old cache to a new physical model. The mapped-water probe uses artificial backend energies only to demonstrate dispatch, not physical numerical error.

The complete final review, commands, issue ledger and additional read-only arithmetic are supplied under `analysis/metag_final_review_20261031_01/` in the code/evidence archive. The original molecular study remains unchanged under `analysis/metag_corrected_validation_20261030_01/` in the current-study archive.
