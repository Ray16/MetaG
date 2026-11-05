# Formation-derived stereochemical reference feasibility

This outcome-free screen retains all 17 distinct paired stereochemical encodings
from the 222 formation rows in the pinned 671-row structural projection. Exact
mirror encodings are filtered. Nine pairs contain unspecified atom stereocentres;
their drawings are not certified as distinct experimental solution states.

The pairs span 11 nonstereochemical graph groups. Both formation components of
every pair are already in the released training table. Thirteen pairs have at
least one pool-alias witness in the bounded MetaG development history; 16 have a
broad-family witness. A witness establishes structural exposure, not which labels
or errors were read. Missing witnesses cannot certify independence.

Every molecular-family holdout preview removes reaction **and formation** labels
involving that family and every unresolved training row. Just one pair—cis/trans
1,2-dichloroethene—has a participant below the predeclared Morgan similarity 0.5
threshold after exclusion (both 0.375). It is a structural control lead, not an
admitted independent metabolic benchmark. No comparator was retrained.

The three neutral-alditol pairs have exact historical structure witnesses for
both partners. Their nearest retained similarities range from 0.667 to 0.833.
The formation-table condition fields (298.15 K, I=0.25 M, pH 7, pMg 14) do not
establish primary experimental provenance, uncertainties or a transfer to the
fixed zero-ionic-strength, no-Mg target. Unknown anomers and solution pools need
source-owned definitions rather than arbitrary stereochemical completion.

From the repository root, with MetaG's CPU audit dependencies installed:

```sh
python -m benchmark.stereo_reference_feasibility --output /tmp/stereo_review.json
```

Inputs are the existing outcome-free projection and bounded historical inventory
in `../feasibility_validation/`; their hashes are in `manifest.json`. The exact
structural source/runtime is also archived. No model weights or primary article
PDFs are redistributed. `registration_metadata.json` records verified titles
and primary registration provenance without republishing abstracts.

`source_review.json` retains three explicit holds. The neutral-alditol source's
primary derivation is unavailable. The earlier Wu azidosugar lead confirms
interconversion, without a thermodynamic constant. A new glycosynthase lead has
two dependent concentration quotients at 343.15 K; source arithmetic reconstructs
the rounded values, but equilibration, activity/proton conventions, uncertainty
and temperature transfer remain unverified. It is not a Gibbs-energy label.
Primary article and ESI locators are retained; local PDF retrieval failed.

This study admits **zero** independent experimental references, computes no MetaG
prediction, changes no production policy and does not relax the primary design.
It resolves these leads' structural feasibility; source or family holdouts still
require independent experimental ownership and an appropriately retrained comparator.
