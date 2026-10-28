# Validation against the evaluated release

MetaG 0.3.0 is a production-only refactor of the development release 0.2.0.dev24 plus four corrections.
The refactor is checked for equivalence first; the corrections are then measured reaction by reaction. Switches for options that were tested and rejected were removed;
the default physics was kept. Equivalence of the refactor (before the corrections) was checked on 2026-10-27 at three levels.

**Reaction level (CPU, exact).** Every reaction of the packaged 364-reaction panel plus the examples
(368 inputs) was scored by both implementations with the same deterministic replay backend. The backend
returns a fixed pseudo-energy for each (charge, canonical SMILES). This exercises every reaction-level
code path: routing, source-backed protonation, hydration microstates, sugar pools, assembly and the fail-closed guards. Three variants
were run: the reactions as written, every reaction reversed, and a variant in which some charged species
fail gas-phase identity checks, which exercises the neutral-reference re-route. All 3 × 368 outcomes are
identical, including the unrounded assembled energies, which match bit for bit. The replay detects a
0.01 kJ/mol change in one constant.

**Species level (GPU).** Free energies of methanol, acetic acid, phosphoric acid, tetramethylammonium (+1),
triplet O₂, glycerol and liquid water were computed from empty caches on a V100. They agree with
0.2.0.dev24 to at most 1.2 × 10⁻⁴ kJ/mol. Two dev24 runs of the same species differ by up to
1.7 × 10⁻⁴ kJ/mol through GPU floating-point nondeterminism. Water and O₂ are bit-identical.

**Whole reactions (GPU).** Four panel reactions were scored end to end from cold caches by both
implementations: pyrophosphate hydrolysis (rxn00001), a hydrated carbonyl (rxn00541), hydration plus
truncation (rxn01117), and truncation with a net proton (rxn45677). They give the same ΔrG'°
(−31.5, 13.0, −10.1, 13.1 kJ/mol), identical stage decompositions to 0.01 kJ/mol, the same scored species
and the same wall time.

**Output schema.** The following fields were removed: the historical class-σ interval block (`sigma_pred`,
`ci95`, `ci_info`, `sigma_breakdown`, `sign_status`, `unresolved`; never coverage-calibrated in
production), the always-null `anchor`/`smd`/`water_ref` records, the `explicit` flag, and the route
flags of removed options. `dG`, `dG_raw`, `stages`, `U_samp`, `routes`, `species_scored` and `suspect`
are unchanged. Species-cache keys differ from dev24 by construction (new source identity), so earlier
caches are never reused.

## Corrections in 0.3.0

The corrections are listed in the [CHANGELOG](../CHANGELOG.md). Their footprint was determined on CPU with the
same replay as above:

- direction-invariant routing changes the route of 13 reactions as written; before it, 16 of 368
  reactions violated ΔG(reverse) = −ΔG(forward), afterwards none do;
- the redox-proton guard changes 1 reaction;
- the deduplication fallback changes 14 species, which occur in 10 reactions;
- failing degraded species changes no reaction of the panel.

These 24 reactions were recomputed with 0.3.0 from empty caches on V100 GPUs. Three unchanged control
reactions recomputed the same way reproduce their dev24 values to within 0.1 kJ/mol (rxn00001 -31.5 → -31.5, rxn01329 +8.8 → +8.8, rxn44278 -6.5 → -6.4).
rxn00824 exceeded the 5400 s per-reaction budget in dev24 and remains a time-budget failure. Values are
in kJ/mol:

| Reaction | Correction | dev24 | 0.3.0 |
|---|---|---|---|
| rxn00191 | orientation | +0.8 | +1.6 |
| rxn00260 | orientation | +10.0 | +10.9 |
| rxn00493 | orientation | +9.9 | +10.8 |
| rxn00527 | orientation | +9.9 | +10.7 |
| rxn00753 | dedup | -23.0 | -21.8 |
| rxn00754 | dedup | -22.0 | -23.2 |
| rxn00806 | orientation | +1.2 | +0.9 |
| rxn00824 | dedup | — | — |
| rxn00849 | orientation | -4.6 | -4.1 |
| rxn00906 | orientation | +7.4 | +9.5 |
| rxn00907 | orientation | +7.4 | +9.5 |
| rxn01201 | orientation | +22.0 | +24.1 |
| rxn01336 | dedup | +3.7 | +3.9 |
| rxn01380 | dedup | +3.6 | +3.6 |
| rxn01409 | dedup | -24.1 | -21.4 |
| rxn01532 | dedup | +21.8 | +21.1 |
| rxn01702 | dedup | -40.5 | -40.5 |
| rxn01726 | dedup | -23.8 | -20.4 |
| rxn01999 | orientation | +8.0 | +8.8 |
| rxn02320 | orientation | +3.8 | +4.7 |
| rxn02902 | redox guard | +26.0 | +40.7 |
| rxn02914 | orientation | +8.4 | +9.3 |
| rxn03281 | dedup | -5.8 | -7.1 |
| rxn16141 | orientation | +22.0 | +24.1 |

On the Mg-bounded benchmark (primary policy) the MAE changes from 8.44 to 8.42 kJ/mol on the same 239
labels. `benchmark/metag_0.3.0_predictions.json` holds the combined panel predictions.
