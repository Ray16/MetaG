# Equivalence with the evaluated release

MetaG 0.3.0 is a production-only refactor of the evaluated development release 0.2.0.dev24, whose
results are reported in the manuscript. Switches for options that were tested and rejected were removed;
the default physics was kept. Equivalence was checked on 2026-10-27 at three levels.

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
