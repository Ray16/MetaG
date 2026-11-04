# MetaG

MetaG estimates **standard transformed reaction Gibbs energies (ΔrG'°)** from molecular structures. It combines UMA electronic energies, per-conformer harmonic thermochemistry, GFN2-xTB/ALPB solvation and protonation transforms. It also provides CPU tools for checking reaction inputs and diagnosing physical approximations.

The target conditions are **298.15 K, pH 7, zero ionic strength, no Mg²⁺ binding, 1 M solutes and pure liquid water**. Standard-state energy signs alone do not determine physiological reaction direction; concentrations and activities are also needed.

MetaG is approximate research software. No numerical offset fitted to the evaluation's biochemical reaction-energy labels enters a prediction. The electronic potential is trained, hydration uses experimental calibration, and reaction-benchmark residuals informed method development. The workflow does not provide calibrated total-error intervals or an established global network potential.

**Start here:** [CPU quick start](#cpu-quick-start) · [GPU installation](#gpu-installation) · [Inputs and outputs](#inputs-and-outputs) · [Reproduction](#reproducing-results-and-checks) · [Contributing](CONTRIBUTING.md)

## CPU quick start

Clone the repository and install its CPU dependencies in a separate environment:

```bash
git clone https://github.com/Ray16/MetaG.git
cd MetaG
python -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -e .
metag tools
metag preflight examples/ester_hydrolysis.json -o output/ester_preflight.json
```

The package requires Python ≥ 3.10. CPU routing and preflight need NumPy and RDKit; they do not require PyTorch, FairChem, UMA weights, xTB or a GPU. The molecular runtime has stricter requirements, described below.

The bundled example describes ethyl acetate hydrolysis. The command writes a JSON report containing input identity, routed structures, protonation terms and structural admission diagnostics. It receives a `review` flag for class-based protonation assumptions. A preflight report is a structural check, not a prediction or evidence of accurate aqueous speciation. Existing output at the chosen preflight path is replaced; use a new path to retain earlier reports.

See [examples/README.md](examples/README.md) for additional inputs, including stereochemical controls and a source-backed protonation declaration.

## How the workflow works

1. **Validate and route the reaction on CPU.** Check structures, charge, balance and declared stereochemistry. Apply structure-based cofactor replacements, spectator truncation and reference-state transforms where admitted. Every rewrite is element- and charge-balance guarded.
2. **Compute molecular free energies.** Generate conformers, optimize them with UMA, validate each retained minimum and calculate its harmonic RRHO correction. Evaluate vertical GFN2-xTB/ALPB water solvation at the corresponding geometry.
3. **Assemble solution states and the reaction.** Combine retained minima, carbonyl hydration states and declared reducing-sugar pools. Add the solute standard-state correction, explicit proton reservoir and independent-site or source-backed coupled protonation terms.
4. **Record results and diagnostics.** Preserve routed identities, component energies, admission failures, configuration fingerprints and molecular-producer provenance. An unavailable estimate is `null`, never zero.

```text
ΔrG'° = Σ ν_i G_i + n_H+ · G(H+, pH 7) + protonation-transform terms
G_i   = -RT ln Σ_minima exp(-(E_UMA + ΔG_solv,ALPB + G_RRHO)/RT) + RT ln 24.46
```

The electronic model is `uma-s-1p2p1` (OMol25). Routing includes NAD(P)/glutathione core replacements, methyl-capped nucleotide spectators and generic spectator truncation. Balance alone does not establish energetic cancellation. The [architecture guide](docs/architecture.md) describes module responsibilities and dependency boundaries.

## GPU installation

Scoring additionally requires:

- A compatible CUDA GPU and PyTorch wheel.
- The official FairChem SDK pinned to commit `645af7d03e91138ef07e51de387316c9e2e152dc`.
- Access to the UMA weights on Hugging Face (`facebook/UMA`) under the FAIR Chemistry License.
- An [xTB](https://github.com/grimme-lab/xtb) executable, version ≥ 6.7, available through `XTB_BIN` or `PATH`.

For the V100/CUDA 12.6 profile, install the matching PyTorch wheel before the `qm` extra:

```bash
python -m pip install 'torch==2.13.0' --index-url https://download.pytorch.org/whl/cu126
python -m pip install -c environments/qm-official-cu126.txt -e '.[qm]'
export XTB_BIN=/absolute/path/to/xtb
metag score --reactions examples/ester_hydrolysis.json --only ester_hydrolysis --out output/ester
```

The calculation writes `output/ester/metag_ester_hydrolysis.json`. A batch checkpoints the accumulated outcome list after each reaction, including failures. An existing batch file at that path is replaced; the CLI does not resume by reading it. Molecular caches can be reused only when their producer and scientific settings match.

The default `nitro` profile uses UMA Nitro for supported single-system optimization on sm_80–sm_120 GPUs and official FP32 for final energies, convergence checks and Hessians. V100 (sm_70) uses official FP32 throughout. The CUDA 13 wheel has no V100 kernels; use CUDA 12.6 for V100. See [environments/README.md](environments/README.md) for the tested profiles, dependencies and model-access requirements. Compatibility of an untested environment is not established by these instructions.

On a shared cluster, use your scheduler or GPU reservation system. Each worker needs its own GPU allocation and bounded CPU/xTB concurrency. Do not launch competing jobs by setting `CUDA_VISIBLE_DEVICES` manually.

## Inputs and outputs

A reaction JSON maps identifiers to reactions. Each species is `[coefficient, formal_charge, SMILES]`, with **products positive**. Include stoichiometric water explicitly. `n_Hplus` follows the same product-positive convention. The quick-start input is:

```json
{
  "ester_hydrolysis": {
    "note": "ethyl acetate + H2O -> acetic acid + ethanol",
    "n_Hplus": 0,
    "species": {
      "ethyl_acetate": [-1, 0, "CCOC(C)=O"],
      "water": [-1, 0, "O"],
      "acetic_acid": [1, 0, "CC(=O)O"],
      "ethanol": [1, 0, "CCO"]
    }
  }
}
```

Use declared stereochemistry and one connected molecular species per SMILES. Unsupported isotopes, unresolved stereochemistry, unsupported phosphoryl environments and explicit-water cluster requests receive no estimate. Stoichiometric reaction water is supported. Atom-map annotations are rejected before molecular work; remove bookkeeping labels while preserving chemical identity and stereochemistry. This restriction does not remove the correspondences constructed internally by routing.

The Python interface accepts one reaction directly:

```python
import json
from pathlib import Path
from metag.energetics.uma import load_uma
from metag.pipeline import score_reaction

reaction = json.loads(Path("examples/ester_hydrolysis.json").read_text())["ester_hydrolysis"]
result = score_reaction(reaction, model=load_uma())
print(result["dG"], result["U_samp"], result["suspect"])
```

| Result field | Meaning |
|---|---|
| `dG` | ΔrG'° in kJ/mol, rounded to 0.1. `null` means unavailable; zero can be a valid estimate. |
| `suspect` | Reason for a missing estimate, such as routing failure, unresolved identity or failed molecular work. |
| `U_samp` | Conformer-sampling spread in kJ/mol; not total predictive uncertainty. |
| `stages` | `species_sum`, `proton` and `pka_transform` contributions. |
| `species_scored`, `original_species` | Routed chemistry and original input identities. |
| `routes` | Applied routes and per-species sampling/thermal diagnostics. |
| `config` | Model, producer, settings, constant-table and software fingerprints. |

Malformed inputs can raise validation errors in the Python API. The batch CLI records individual scoring exceptions as missing estimates and continues; a failed checkpoint write stops the batch.

### Commands and configuration

Run `metag score --help` or `metag <tool> --help` for the full interface.

| Command | Purpose |
|---|---|
| `metag preflight INPUT -o REPORT` | CPU structural and protonation admission diagnostics. |
| `metag score --reactions INPUT --out DIR` | Molecular scoring; `--only KEY` is repeatable. Without an input file, uses `RXN_FILE` or the packaged 364-reaction panel. |
| `metag evaluate PANEL -o REPORT` | Evaluate predictions against an explicitly declared reference and comparison conditions. |
| `metag build-reference --help` | Reconstruct a reference ledger from external source data; requires the `reference` extra. |

Search settings are resolved once per backend through `metag.settings.Settings.from_env`. Set environment variables **before** starting the CLI or importing MetaG: cache, xTB and some execution options are read at module import. Species caches include result-affecting settings, model identity and the xTB binary hash. Set `METAG_CACHE` to choose the cache directory.

| Variable | Default | Purpose |
|---|---|---|
| `UMA_MODEL` | `uma-s-1p2p1` | Checkpoint identity. |
| `UMA_RUNTIME_PROFILE` | `nitro` | Hardware-dependent Nitro optimization with FP32 final evaluation; `official_fp32` uses FP32 throughout. |
| `UMA_HESSIAN` | `autograd` | Autograd Hessian; `fd` uses central differences. |
| `METAG_SAMPLING_PROFILE` | `current` | `thorough` increases search effort; also available as `--sampling-profile thorough`. It does not guarantee accuracy. |
| `CONV_TOL`, `CONV_HITS`, `CONV_MAX` | 2.5 kJ/mol, 2, 8 | Sampling plateau criterion. |
| `SAMPLE_SCALE`, `KEEP_SCALE` | 1, 1 | Scale generated and retained conformer pools. |
| `CONFORMER_SEED_OFFSET` | 0 | Seed offset for replicate searches. |
| `XTB_BIN` | `xtb` on `PATH` | xTB executable. |
| `METAG_XTB_WORKERS` | 8 | Concurrent single-threaded xTB processes; budget these CPU slots per worker. |
| `UMA_ENERGY_CHUNK`, `UMA_FIRE_CHUNK`, `UMA_HESS_CHUNK` | Automatic | GPU batch sizes controlling memory use. |

Changing execution profiles or chunk sizes can introduce floating-point or sampling variation. Numerical agreement across environments must be checked. Autograd Hessians temporarily wrap FairChem module functions, so concurrent inference within one Python process is unsupported; use separate processes and allocations.

## Physical diagnostic companion

The source distribution includes `benchmark/thermochemistry_diagnostics.py`, a CPU tool for retained-basis solvent bounds and complete cycle-space consistency certificates. It does not alter deployed predictions. These bounds measure conditional model inconsistency, not calibrated physical error. Install the `audit` extra, then run the bundled synthetic controls:

```bash
python -m pip install -e '.[audit]'
python benchmark/thermochemistry_diagnostics.py transfer examples/diagnostics/transfer.json -o output/transfer.json
python benchmark/thermochemistry_diagnostics.py cycles examples/diagnostics/triangle.json -o output/cycles.json
```

Use new diagnostic report paths: this companion refuses to overwrite existing reports. The companion and examples ship in the source distribution, not the minimal scoring wheel. See [docs/physical-diagnostics.md](docs/physical-diagnostics.md) for inputs, assumptions and certificate interpretation.

## Reproducing results and checks

The working source is an unreleased `0.3.1.dev0` development version, distinct from public `0.3.0`. The completed corrected molecular evaluation used source fingerprint `3d42323af063483c2bcc44f96e8b86ba5ae85024369a39ae59d67f20579dcd47`. Later source changes must not be presented as a fresh evaluation of that study. The full raw study archive and final public release remain pending.

On the corrected evaluation, 270 of 364 attempted reactions have estimates; 92 are withheld and two reach the time budget. On 198 of 274 Mg-bounded reference labels, MAE is 8.38 kJ/mol versus 11.96 for zero. The preceding estimator gives 8.37 on those same targets: changed coverage explains the lower aggregate error. These labels informed development.

On 69 admitted out-of-fold dGbyG pairs, MAEs are 7.79 (MetaG), 4.59 (dGbyG) and 7.81 kJ/mol (zero), without fully propagated bridge uncertainty. A shared-record audit evaluates 30 of 49 independent input-state cycle vectors; four fail numerical closure, with a maximum residual of 8.70 kJ/mol. These results support approximate reaction-level thermochemistry and diagnostics. They do not establish independent biochemical accuracy, complete aqueous-state sampling or validated network flux constraints.

Reproduction has three distinct levels:

| Level | Starting evidence | What it verifies |
|---|---|---|
| Recalculate metrics | Retained prediction and reference files in `benchmark/` | Cohorts and reported summary statistics. No molecular calculations. |
| Reassemble reactions | Immutable molecular records, exact source and settings from the study archive | Reaction assembly and controls on that retained evidence. Not new sampling. |
| Recompute molecules | Model access, pinned runtime, xTB and GPU resources | New optimization, Hessians and sampling. Runtime and numerical variation must be assessed. |

Recalculate the current retained benchmark metrics on CPU:

```bash
python benchmark/score_checked.py
```

The [benchmark guide](benchmark/README.md) identifies current and historical cohorts, reference assumptions, required external datasets and comparator reconstruction. The [validation record](docs/validation.md) distinguishes source equivalence tests, cached replay and fresh computation. Do not compare version-level MAEs without matching the target set, and never treat a missing estimate as zero.

### Tests

The CPU unit/regression suite imports PyTorch and ASE even though it does not require FairChem, model weights or xTB. Install a compatible CPU PyTorch wheel first, then the test extra:

```bash
python -m pip install 'torch==2.13.0' --index-url https://download.pytorch.org/whl/cpu
python -m pip install -e '.[test]'
python -m pytest -q -m 'not gpu'
```

GPU integration tests under `tests/integration/` are opt-in and require the molecular runtime. A passing CPU suite does not reproduce molecular energies or validate scientific accuracy. See [CONTRIBUTING.md](CONTRIBUTING.md) for development, regression and package-build checks.

## Repository layout

```text
src/metag/          Installable package and bundled reaction/protonation data
  routing/          CPU structural rewrites and reference-state routing
  energetics/       UMA, conformers, RRHO, solvation, state models and cache
  tools/            CPU preflight, evaluation and reference construction
benchmark/          Retained predictions, references, diagnostics and audit scripts
examples/           Small reaction inputs and synthetic diagnostic controls
environments/       Pinned molecular-runtime profiles and installation notes
docs/               Architecture, validation, protonation and solution-state guides
tests/              CPU unit/regression tests and opt-in molecular integration tests
```

## Citation, license and acknowledgments

Use [CITATION.cff](CITATION.cff) for the available software citation metadata. The associated manuscript is in preparation; no article DOI is claimed here.

MetaG's existing [MIT License](LICENSE) covers the repository software. UMA weights are not included and have separate FAIR Chemistry License terms. MetaG builds on FairChem/UMA, xTB, RDKit, ASE and eQuilibrator. Benchmark labels derive from openTECR and TECRDB; consult their provenance and rights before redistribution.
