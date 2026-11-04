# Rare-sugar equilibrium candidates: initial screen

The 1994 study [Itoh et al., DOI 10.1271/bbb.58.2168](https://doi.org/10.1271/bbb.58.2168)
reports equilibria for four ketose epimerizations, using purified D-tagatose 3-epimerase. The paper
reports 30 °C, 25 mM Tris-HCl pH 7.5, 1% substrate, and equilibrium assessed after four days. Its
four substrate/product pairs are represented in [the structure-only input](itoh_1994_structures.json).
These are *candidates*, not an admitted benchmark or an independent-validation claim.

The report identifies total sugar parents, so the inputs declare reducing-sugar solution pools. Each
SMILES is an open-chain representation with stereochemistry derived from PubChem identifiers:
D-tagatose CID 439312, D-sorbose 107428, D-psicose 441036, D-fructose 2723872, D-ribulose 151261,
D-xylulose 5289590, L-ribulose 644111 and L-xylulose 22253. PubChem's cyclic records for tagatose,
psicose and fructose were opened with MetaG's sugar-state generator; this does not establish that the
finite modeled pool has accurate state populations. The full paper and stereochemical assignments
should be rechecked before any calculation is interpreted.

`metag.tools.preflight` reports `no_known_gap_detected` for all four reactions. The declared ketose
pool generator constructs five ring/open states per hexose and three per pentose, before hydration.
Preflight checks structural admissibility; it is not a successful energy calculation or evidence of
accuracy. The input file contains no reference values.

Several issues must be resolved before these candidates are scored against either method:

1. The measurements are at 303.15 K; MetaG targets 298.15 K. Temperature correction needs reaction
   enthalpy or a sensitivity analysis, and buffer/ionic effects and water-solution conventions need
   assessment. Source uncertainty and evidence that equilibrium was reached must be reviewed.
2. [Kim et al. 2006](https://pmc.ncbi.nlm.nih.gov/articles/PMC1392948/) report a different
   fructose/psicose ratio at 30 °C from the 1994 study. This is an unresolved reference discrepancy,
   not a reason to select the value more favorable to a prediction.
3. These species and related sugar chemistry appear in earlier MetaG investigations and dGbyG's
   training corpus. An exact reaction match was not found in the parsed released training table, but
   tautomer and source-level overlap remain to be checked. Their development exposure has not been
   fully reconstructed. Studies using the same chemical transformation are not independent examples.
4. The original domain proxy does not encode chirality. All four reactions collapse when
   stereochemistry is stripped and have an empty nonchiral reaction-change fingerprint. A separate
   chirality-aware proxy gives nearest similarities of 0.33, 0.33, 0.60 and 0.33 respectively. Its
   nearest training match for each is D-ribulose-5-phosphate ↔ D-xylulose-5-phosphate. This is a
   heuristic similarity, not a validated extrapolation test; the source compounds and conditions differ.
5. The literature results are approximate equilibrium ratios, not definitive 298.15 K reference
   labels. Prediction inputs and reference curation must remain separate; no results are admitted
   until reaction direction, full solution state and experimental provenance are reviewed.

The 2026 [Cho and Ahn qNMR study](https://doi.org/10.1002/bkcs.70170) is a possible *component*
validation of fructose and allulose ring populations. It is not a direct reaction-energy benchmark;
the reported solvent and experimental conditions require review before use.

## Published dGbyG graph limitation: a concrete hypothesis

The published dGbyG inference path (`api/utils.py` calling `mol_to_graph_data`) uses atomic number,
hybridization, aromaticity and charge, with bond type as its bond feature. Its default graph features
do not encode chirality. Other files in the package contain chirality-related code; this conclusion
is specific to the published ensemble inference path checked here.

For each of the four **open-chain input-state pairs**, the normalized graphs are isomorphic with
identical model features. The official dGbyG API returned exactly zero difference for each of the 100
ensemble members on each pair. The reproducible [audit](itoh_1994_dgbyg_stereo_audit.json) records
the normalized inputs and source/weight hashes. This is a structural limitation of the tested
single-state inference: identical input graphs cannot encode different graph-based energies.
It is **not** yet a comparison of solution-pool energies against experiment or proof MetaG is more
accurate. Both methods must be evaluated under comparable state and condition conventions; dGbyG's
condition correction could add additional terms in a broader use case.

Reproduce locally with the archived dGbyG installation and weights:

```
../tools/dGbyG/.venv/bin/python benchmark/audit_dgbyg_stereo.py \
  --model-root ../tools/dGbyG --reactions benchmark/candidates/itoh_1994_structures.json \
  --output /tmp/dgbyg_stereo_audit.json --infer
```

The `--infer` option uses the upstream API with CPU selection forced during import. The tested
open-chain representations are not asserted to be the appropriate experimental solution-state
representation for dGbyG; testing that choice and the reference discrepancy remains necessary.
