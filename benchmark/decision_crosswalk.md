# Historical decisions versus current-source audit controls

The historical exposure and unknown dates are recorded in the
[historical decision log](history/dev24_selection/DECISION_LOG_CURRENT.md); this table
does not establish dates or claim that a current counterfactual exactly
recreates an earlier implementation. References below to implemented audit
controls describe provisional analysis scripts in the separate development
workspace; they are not included in this repository or available as production
settings. The full current-source audit remains incomplete; an
availability-limited secondary analysis is described in
[development exposure](development_exposure.md). These controls have not been
used to select a production policy. A finite 901-policy catalog is frozen, but
it does not reconstruct all historical choices or interactions. Safeguards
cannot be removed silently to increase matrix coverage.

| Historical decision | Current-source behavior or boundary | Current audit treatment |
| --- | --- | --- |
| `PKA_ENV` | Environment-specific carboxyl constants in `routing/pka_transform.py` | Per-environment flattening to the generic value is implemented; this is a counterfactual approximation, not proof of a historical equivalent. |
| `FREE_PPI_PKA` | Free pyrophosphate uses `PPI_LADDER` | An alternative ladder is not yet implemented or validated. |
| `ANHYDRIDE_PKA` | Phosphoryl-ladder and domain logic in `pka_transform.py` and `pka_domain.py` | Historical alternative and its physical domain remain to be mapped. |
| `CARBOXYL_PAIRS` | Current carboxyl environment classifier treats adjacent carboxyl groups separately from alpha classes | Historical paired-acid alternative is not reconstructed. |
| `POLYACID_PKA` | Groupwise ladders and optional coupled models are distinct mechanisms | Historical polyacid alternative is not reconstructed. |
| `ARYLAMINE_NONBASIC` | Aromatic amines are excluded from `_AMINE_N`; `_base_gate` uses structural counts | Historical gating alternative is not reconstructed; preserve charge and proton consistency. |
| `PH0_REDOX_PROTON` | Anion-only pH-0 build carries a residual proton only when no unneutralized anion remains | No safety guard is removed; an older redox policy is not reconstructed. |
| `SOURCE_PKA_REGISTRY` | Exact source-backed protonation models can replace owned independent-site terms | Analysis-only disabling of automatic registry is implemented; explicit caller models remain active. |
| `RING_PROTON_ENSEMBLE` | No comparable production switch was identified | Historical rejected candidate remains unimplemented, not counted as tested. |
| `COFACTOR_RING` | `_route_cores` applies cofactor ring rewriting with balance guard | Analysis-only bypass of the rewrite is implemented; other routing remains. |
| `COA_CORE` | No matching current switch identified | Historical alternative remains unimplemented. |
| `NTP_CORE` | `_route_cores` applies nucleotide coring with balance guard | Analysis-only bypass is implemented; generic truncation remains a separate decision. |
| `AUTO_TRUNCATE` | Public `allow_truncate` controls generic spectator truncation, not cofactor cores | Public full-route ablation is available. |
| `ROUTE_FULL` | `truncation_gate.prefer_full` can prevent generic truncation | Analysis-only bypass of this preference is implemented, retaining truncation safety checks. |
| `TRUNC_V2` | Current `truncate.py` contains the refactored structural algorithm | Historical alternative not reconstructed. |
| `TRUNC_VALIDATE` | Current default does not automatically compare two truncation radii | Separate radius alternatives are available; no automatic uncertainty or acceptance bound is implied. |
| `PH0_AUTO` | Primary pH-0 routing is followed by a neutral-reference guard | Analysis-only skip of the primary route is implemented; the later guard remains. |
| `PH0_BASES` | `_base_gate` and `build_ph0_reaction(..., base=True)` determine primary base handling | Anion-only primary build is implemented; forced neutral-reference behavior remains. |
| `ZWITTERION_PH0` | Neutral-reference guard and rearrangement handling protect species identity | Guard is not disabled; no valid counterfactual is yet established. |
| `NEUTRAL_QM` | Forced neutral routing is an internal recovery path | Not reconstructed as a standalone scoring policy. |
| `PH0_ISOMERASE` | Primary pH-0 builder operates without a separate exposed isomer toggle | Isomer-only primary-route skip is implemented; later guards remain. |
| `CARBONYL_HYDRATION_ALL` | Complete supported hydrate microstates are folded into the species ensemble | Analysis-only omission is implemented as an incomplete-state ablation, not a physically complete estimator. |
| `RING_CHAIN` | Reducing-sugar pools require explicit declarations and have completeness checks | No undeclared or incomplete pool is fabricated; historical alternative not reconstructed. |
| `HYDRATION_CAL` | Independent-property calibration maps raw hydrate energies to the mixture | Analysis-only raw `(1, 0)` map is implemented; this is not a reaction-fitted parameter search. |
| `STD_STATE_1M` | The same 1 M convention is applied to solutes | Inconsistent standard-state alternatives are quarantined. |
| `ALDEHYDE_HYDRATION` | Production uses general hydration enumeration rather than the old activated-aldehyde helper | Omission and raw calibration are available. The separately hashed v2 overlay adds an activated-aldehyde-only site control; it is a current-code proxy, not a verified exact reconstruction of historical behavior. |
| `ANCHOR_CORRECT` | No reaction-fitted offset in current scorer | Reaction-label-fitted correction is excluded from this project. |
| `SMD_SOLV` | Current production species backend uses ALPB | Alternative SMD species calculations require a separate, validated producer; no values are invented from ALPB caches. |
| `WATER_REF_HYDROLYASE` | Current backend uses a common experimental-hydration liquid-water reference | Historical class-specific water alternative is not implemented. |
| `WATER_REF_EXP` | Common liquid-water reference requires its unrounded backend value | Two producer-matched captures agreed at full precision and were consistent with bounds inferred from rounded baseline stages. The partial default-replay check is still ongoing; historical alternate water values are not imported. |
| `TRUNC_SPECTATOR_CATIONS` | A current guard rejects truncations that mangle retained cation centers | No safety guard is disabled. |
| `TRUNC_MAXANION_RETRY` | Current planner has bounded fallback and neutrality logic | Historical retry alternative is not reconstructed. |
| `TRUNC_FG_CUTS` | Functional-group-preserving `_cut_close` is used with ring closure | The structural safety guard is not disabled; the historical unsafe cuts are not treated as valid predictions. |
| `TRUNC_ANOMERIC_RADIUS` | `truncation_radius` raises the radius at anomeric centers | Analysis-only disabling of the adaptive radius is implemented with remaining guards intact. |

Known residual-informed choices without faithful alternatives remain an
explicit gap. A finite grouped selection test over the implemented union must
be labelled with its exact candidate set and cannot quantify the effect of the
unreconstructed rules or all prior development.
