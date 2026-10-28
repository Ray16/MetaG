# Changelog

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
- Runtime numerics live in one `Settings` object that is resolved per call. The same object drives the
  computation and the species-cache key; there is no import-time state.
- `score_reaction(reaction, backend)` gets species energies through a small backend interface
  (`UMABackend`), so all reaction-level code can be replayed on CPU.
- New `metag score` CLI options: `--reactions`, `--only` (repeatable), `--out`, `--sampling-profile`.
- Removed output fields: the historical class-σ interval block, the always-null anchor/SMD/water-patch
  records and the explicit-water flag.
- `benchmark/` reproduces the reported Mg-bounded openTECR metrics.
