# Changelog

## Unreleased (0.3.1.dev0)

- Add the released dGbyG ensemble comparison on the fixed 69 development
  reactions, with all 100 member sums and molecular covariance retained.
  Independent decimal/statistical replay checks the complete packet without
  Torch. This is separate from OOF-repeat spread and supplies no independent
  reliability claim; 68 normalized equations match the shipped training table.

- Add a bounded, same-Hamiltonian LBFGS recovery when a force-converged
  candidate still has a strong imaginary mode after the final FIRE and
  mode-following retry. Recovery requires forces <=0.001 eV/A, preserved
  identity and a fresh Hessian passing the unchanged 50 cm-1 gate. Electronic,
  solvation and thermal terms are recomputed at the accepted endpoint. The
  raw-species source/cache identity changes; this is not a new full benchmark.

- Refuse explicit calculation-condition requests outside the implemented
  298.15 K, pH 7, zero-ionic-strength, Mg-free target before molecular work.
  Both scoring APIs previously ignored this optional input metadata. Omitted
  conditions retain the documented fixed target; the check performs no
  temperature, pH, ionic-strength or Mg transformation. Historical evaluated
  results are not relabelled as predictions at native experimental conditions.

- Reject a timed-out RDKit maximum-common-substructure search instead of using
  its partial atom map. The reaction is withheld through the routing error
  guard. This changes routing behavior and the estimator fingerprint; the
  previously reported frozen-source benchmark predates this correction.

- Make the archived dGbyG comparison default to the manuscript's exact-family
  identity cohort, checking paper-row, signed extent and condition-bridge
  metadata. Retain the broader unresolved-identity count separately and report
  descriptive within-class correlations.
- Reject non-finite or invalid sampling settings before scoring, and reject
  scales that produce an empty conformer pool or no retained conformers. Default
  settings and the physical energy model are unchanged. The source fingerprint
  changes, so existing species caches are not silently relabelled.
- Withhold a reaction when a backend returns a missing, non-finite or invalid
  species energy or sampling diagnostic, including required sugar states and
  the liquid-water reference used in hydration.
- Avoid claiming that all converged conformers rearranged when a species has
  mixed identity and solvation failures; such a case is withheld rather than
  triggering the neutral-reference reroute.
- Withhold a species if its heavy-atom automorphisms also exceed the bounded
  deduplication mapping search, instead of treating a truncated search as exact.
  A CPU screen of 829 distinct input and routed structures from the 364-reaction
  panel found none above that limit; auxiliary enumerated states were not screened.
- Validate optional experimental reference values before loading or scoring
  species, so invalid metadata cannot discard an otherwise completed calculation.

## 0.3.0 — 2026-10-27

First public release: a production-only refactor of the development release 0.2.0.dev24, followed by
four corrections found in a code review. The refactor alone reproduces dev24 exactly; the corrections
change 24 of the 364 panel reactions. See [docs/validation.md](docs/validation.md).

Corrections:

- **Direction-invariant routing.** MCS pairing, atom maps and the reaction-centre ring gate depended on
  which side was written first, and 16 of 368 replayed reactions violated ΔG(rev) = −ΔG(fwd). Every
  reaction is now routed in a canonical orientation and mapped back. 13 panel reactions change route.
- **Redox-proton guard.** The pH-0 route carried an H residual as a redox proton after a charge-closure
  test that always passed. It now requires that no neutralised species keeps a negative charge
  (rxn02902 changes).
- **Conformer deduplication.** Species with more than 20,000 all-atom automorphisms (methyl/H rotors) are
  compared by heavy-atom RMSD, so copies of one basin can no longer be counted twice (14 species,
  10 reactions).
- **Degraded species fail closed.** A species with more than 25% failed xtb solvations now fails instead
  of being used uncached. No panel reaction is affected.

Refactor:

- Only the default method remains. Options that were tested and rejected were removed: anchor offsets,
  explicit-water clusters, SMD/CPCM-X/COSMO solvation, solvent relaxation, qRRHO, CoA core, ring-chain and
  ring-proton pools, the site-pKa model and the legacy runtime profile.
- Search and scientific settings live in one `Settings` object that is resolved per backend. The same
  object drives the computation and species-cache key. Some execution and cache settings remain
  import-time values; see [docs/architecture.md](docs/architecture.md).
- `score_reaction(reaction, backend)` gets species energies through a small backend interface
  (`UMABackend`), so all reaction-level code can be replayed on CPU.
- New `metag score` CLI options: `--reactions`, `--only` (repeatable), `--out`, `--sampling-profile`.
- Removed output fields: the historical class-σ interval block, the always-null anchor/SMD/water-patch
  records and the explicit-water flag.
- `benchmark/` reproduces the reported Mg-bounded openTECR metrics.
