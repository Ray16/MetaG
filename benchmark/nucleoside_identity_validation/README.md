# Covalent identities for the nucleoside reference panel

This inventory supplies neutral reference drawings for the complete substrate
panel in [Kaspar et al.](https://doi.org/10.1002/adsc.201901230), using the
published supplementary Figure S1 and Table S1. It supplements the
[source-arithmetic review](../nucleoside_reference_review/README.md).
Experimental standard-state references and solution-state pools remain
unadmitted; no thermodynamic prediction is performed.

Each of the 24 structures passes an explicit beta-D sugar check at the identified
C1′–C4′ positions. Four ribose/deoxyribose pairs differ only by removal of the
2′ hydroxyl, and 16 substituted analogues retain their parent sugar configuration
and match substitution at pyrimidine C5 or purine C2 in the source figure. Those checks authenticate
covalent drawings, not their protonation, tautomer or metal-bound populations.

The original name lookup for 5-fluoro-2′-deoxyuridine returned PubChem CID 3363,
which omits all three sugar stereocentres. Figure S1 specifies the beta-D isomer.
An exact InChIKey lookup of that source-derived structure returns
[CID 5790](https://pubchem.ncbi.nlm.nih.gov/compound/5790). Both unaltered responses
are retained. The generic entry is not silently assigned the specific isomer's
stereochemistry. Tests reject the unspecified entry, individual sugar epimers,
the opposite enantiomer and a substituted-compound swap.

The source figure labels deoxyinosine 64, while Table S1 and the main experimental
panel identify compound 24. The numbering discrepancy and source locations
are recorded; no thermodynamic observation or original source file is changed.

All 25 PubChem responses are preserved under `properties/`. The 24 name-query
receipts and the additional exact-stereo lookup retain URLs and response hashes.
The source-compound file records which response supplies each reference drawing.
No test case is selected using a model prediction or residual.

From the MetaG checkout:

```bash
python benchmark/audit_nucleoside_identity.py \
  benchmark/nucleoside_identity_validation \
  --output /tmp/nucleoside_identity_audit.json
```

The retained audit records its RDKit version and source/input identities.
`manifest.json` authenticates all retained files except itself. The replay uses
cached responses and does not contact PubChem. Primary experimental eligibility,
training/development overlap, novelty and actually scored molecular structures
still require their own checks.
