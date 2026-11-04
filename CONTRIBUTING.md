# Contributing to MetaG

Changes should improve the usability or scientific reliability of the workflow while preserving the meaning of existing results. Start with the [README](README.md), [architecture](docs/architecture.md) and [validation record](docs/validation.md).

## Development environment

Use a separate environment. CPU tools need NumPy and RDKit; the test suite also imports ASE and PyTorch. From the repository root:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install 'torch==2.13.0' --index-url https://download.pytorch.org/whl/cpu
python -m pip install -e '.[test]'
python -m pytest -q -m 'not gpu'
```

GPU work needs a separately provisioned, compatible molecular runtime. Use the [pinned profiles](environments/README.md), obtain model access and configure xTB before running it. Do not use the CPU test environment as evidence of GPU compatibility.

## Scientific behavior and regressions

Keep a refactor separate from a change to a physical model. Changes to routing, constants, reference states, conformer sampling, minimum validation, solvation, thermal terms or numerical precision can alter the estimator even when an average benchmark score improves. Describe the affected chemistry, the reason for the change and the evidence supporting it. Preserve adverse results and report coverage on identical targets.

For an engineering refactor, compare behavior before and after on representative inputs. Add tests for a meaningful failure or invariant, rather than tests that only mirror the implementation. Useful checks include reaction reversal and extent, molecular identity, unchanged-species cancellation, explicit state ownership and complete failure reporting. CPU replay can verify reaction assembly on supplied energies; it does not validate fresh molecular sampling or physical accuracy.

Results must distinguish unavailable estimates from valid zero energies. Do not replace missing results with zero, silently renormalize an incomplete declared ensemble, or treat conformer-sampling spread as a calibrated prediction interval. Any new default needs appropriate physical validation and an updated provenance record.

## Source and package boundaries

- `routing/` and CPU help/preflight must not import the optional QM runtime.
- `energetics/species.py` must not depend on the reaction pipeline. The backend interface supplies molecular results to the pipeline.
- Capture scientific settings, model/checkpoint identity and producer provenance. Do not weaken cache invalidation to reuse incompatible records.
- Keep executable examples small and omit credentials, model weights, caches and machine-specific absolute paths. Document necessary external inputs and their provenance.
- Preserve the existing file layout and interfaces unless a change has a clear maintenance benefit. Historical benchmark files are evidence, not disposable scratch output.

## Lightweight checks

```bash
metag tools
metag score --help
metag preflight examples/ester_hydrolysis.json -o output/contribution_preflight.json
python benchmark/score_checked.py
```

The physical diagnostic companion needs `.[audit]` and is included in the source distribution. Its examples are synthetic controls:

```bash
python -m pip install -e '.[audit]'
python benchmark/thermochemistry_diagnostics.py transfer examples/diagnostics/transfer.json -o output/contribution_transfer.json
python benchmark/thermochemistry_diagnostics.py cycles examples/diagnostics/triangle.json -o output/contribution_cycles.json
```

Choose new diagnostic output paths; the companion refuses to overwrite existing reports. Preflight replaces its selected output file. Scoring checkpoints an accumulated batch to the selected directory and does not resume by reading an earlier batch JSON.

## Build and inspect a distribution

Build from a clean, reviewed source snapshot in a separate environment:

```bash
python -m pip install build
python -m build
```

Check the actual wheel and source archive, not only the working tree. The wheel includes the installable `metag` package and bundled JSON assets. The source distribution also includes examples, documentation, tests, environments and benchmark scripts. Install the built wheel in a fresh environment and run a CPU preflight from outside the source tree to detect accidental editable-install dependencies. Reproducing manuscript evidence additionally requires the identified raw study sources; installing the package alone does not supply them.

## Reporting a change or a problem

Provide the MetaG version or source fingerprint, command, relevant configuration, environment, input or minimal reproducer, expected behavior and actual diagnostic. For numerical differences, include model identity, hardware and execution profile. Do not publish private inputs or access tokens. State exactly which checks ran and whether validation involved retained predictions, cached molecular records or fresh computation.

Before a release, review the source diff, build artifacts, runnable examples, scientific claim-to-evidence links, external data/model rights, license and citation metadata. Software tests cannot establish independent biochemical accuracy or guarantee journal acceptance.
