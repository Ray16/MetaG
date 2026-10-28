# Changelog

## 0.3.0 — 2026-10-27

First public release. This is a production-only refactor of the evaluated development release
0.2.0.dev24. Scored species, energies and ΔrG'° are unchanged; see [docs/validation.md](docs/validation.md).

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
