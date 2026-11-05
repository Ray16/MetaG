# Explicit ModelSEED cofactor curation: separate research profile

Four ChEBI identities were selected from the cross-references already present in
ModelSEED revision b539adb69a955dac15cc30dbeb3debbd87739018: 57540, 57945, 57783,
58349 for NAD, NADH, NADPH, NADP. Their owned ChEBI molfiles and PubChem charged
records agree with the curated development-panel amides. They preserve all eight
declared carbon stereocenters. ModelSEED's imidate-like representations differ
in bond order/proton placement; the reduced structures also differ by one H+
from those referenced molecular species. This is explicit source-identity
curation, not an assertion that the calculated microstates have equal free energy.

The literal source data and existing stereo gates are unchanged. The separate
profile replaces exactly these four source rows, pins their original identities,
and rebases the explicit H+ coefficient in each compartment by minus the sum of
coefficient times added molecular hydrogens. It never adds water or forces an
unbalanced equation to balance. This basis conversion applies globally to the
metabolite identity, not by reaction class or prediction error.

Of 43,650 live nontransport inputs, 7,511 touch the four cofactors. Recomputing all
affected gates admits 3,147 more inputs, changing 6,572 literal passes to 9,719
separate-profile passes, with no losses. All 5,854 known formula residuals and
all 7,511 declared-charge residuals are preserved; 7,407 equations require a
proton-coefficient change. The other 36,139 input records remain unchanged.
Every newly admitted input retains protonation warnings. Of the 3,147 admitted
affected inputs, 3,135 cancel full cofactor-family coefficients and 12 retain
them. Those 12 are prospective development controls, not validated energies.
No quantum calculations, independent Gibbs-energy references, numerical
prediction coverage or production-policy promotion are supplied by this study.

Run the independent offline accounting audit from the repository root:

```sh
python benchmark/audit_modelseed_cofactor_curation.py --cold-targets
# Additionally reauthenticate the owned upstream tables and their raw graphs:
python benchmark/audit_modelseed_cofactor_curation.py --modelseed /path/to/ModelSEEDDatabase
```

The offline audit independently reads all primary structures and the declared
atom mapping, reconstructs proton rebasing, balance residuals, input gates,
gained/lost inputs and scaffold-boundary membership. It uses derived atom counts
and recorded unchanged admission decisions for noncofactor molecules. The raw
mode additionally verifies all source-table hashes, recomputes every participating
compound graph/formula and checks all affected equations against the upstream
tables. Neither mode validates pKa accuracy, complete solution speciation or
thermodynamic prediction quality. Counterfactual gate timings are not a quantum
cost benchmark. The original frozen code is retained in runtime_source.zip.

Primary source URLs and byte hashes are in registry.json/retrieval.json. Owned
payloads include the four ChEBI molfiles and four PubChem property records;
ancillary original HTML retrieval hashes are recorded but those website pages
are not redistributed. Original ModelSEED source tables remain caller-owned;
derived atom accounting and equations use source IDs and rational coefficients.
Source licensing/attribution notices remain in upstream_LICENSE.txt.

Sources: ModelSEEDDatabase, https://github.com/ModelSEED/ModelSEEDDatabase/tree/b539adb69a955dac15cc30dbeb3debbd87739018;
ChEBI, https://www.ebi.ac.uk/chebi/; PubChem, https://pubchem.ncbi.nlm.nih.gov/.

The owned atom mapping addresses the parsed source-literal SMILES and parsed target canonical SMILES, rather than Molfile atom numbering. `--cold-targets` independently checks all four original state-reference plans in a fresh process using the archived runtime, while forbidding quantum imports and network access. It does not rerun every noncofactor state gate.
