# MetaG

MetaG estimates **standard transformed reaction Gibbs energies, ΔrG'°**, of biochemical reactions from
molecular structure alone. It needs no training on reaction thermodynamics.

Target conditions are 298.15 K, pH 7, zero ionic strength, no Mg²⁺ binding, 1 M solutes and pure liquid
water.

```
ΔrG'° = Σ ν_i G_i  +  n_H+ · G(H+, pH 7)  +  Σ Alberty pKa-transform terms
G_i   = -RT ln Σ_minima exp(-(E_UMA + ΔG_solv,ALPB + G_RRHO)/RT)  +  RT ln 24.46
```

- **Electronic energies:** Meta's UMA machine-learned interatomic potential (`uma-s-1p2p1`, OMol25).
- **Thermochemistry:** per-conformer harmonic RRHO from the exact UMA Hessian.
- **Solvation:** vertical GFN2-xTB/ALPB(water).
- **Routing before QM:** reactions are rewritten so the computed species are ones these models describe
  well. Every rewrite is element- and charge-balance guarded:
  - isodesmic NAD(P) and glutathione cores, methyl-capped spectator nucleotides;
  - spectator truncation around the reaction centre;
  - neutral (pH-0) reference species, bridged to pH 7 with exact Alberty transforms.
- **Solution states:** carbonyl hydration microstates and declared reducing-sugar pools (anomers, ring
  sizes, open chain) are folded into the species free energies.
- **Fail-closed results:** invalid or incomplete results return `dG: null` with a reason. They never
  return a number.

On the Mg-bounded openTECR benchmark, MetaG gives MAE 8.4 kJ/mol on 239 of 274 labels; predicting zero
gives 11.7 on the same labels. See [benchmark/](benchmark/) for the caveats. MetaG is intended as a
training-free complement to fitted group-contribution methods for chemistry outside their training data.
It does not replace them in-distribution.

## Installation

CPU only (routing, protonation, preflight; Python ≥ 3.10):

```bash
pip install -e .
metag tools
```

GPU scoring additionally needs:

- the pinned official FairChem SDK and Torch;
- access to the UMA weights on Hugging Face (`facebook/UMA`, FAIR Chemistry License);
- an [xtb](https://github.com/grimme-lab/xtb) ≥ 6.7 executable.

Install the Torch wheel that matches your GPU first, then the `qm` extra with the matching constraints
file (see [environments/](environments/README.md)):

```bash
pip install 'torch==2.13.0' --index-url https://download.pytorch.org/whl/cu126     # V100 / CUDA 12.6
pip install -c environments/qm-official-cu126.txt -e '.[qm]'
export XTB_BIN=/absolute/path/to/xtb
```

## Quick start

A reaction file maps identifiers to balanced reactions. Each species is
`[coefficient, formal_charge, SMILES]`, with products positive. Include water explicitly. `n_Hplus` uses
the same product-positive sign convention.

```json
{
  "ester_hydrolysis": {
    "note": "ethyl acetate + H2O -> acetic acid + ethanol",
    "n_Hplus": 0,
    "species": {
      "ethyl_acetate": [-1, 0, "CCOC(C)=O"],
      "water":         [-1, 0, "O"],
      "acetic_acid":   [ 1, 0, "CC(=O)O"],
      "ethanol":       [ 1, 0, "CCO"]
    }
  }
}
```

Check structures, protonation and routing on CPU before spending GPU time:

```bash
metag preflight reactions.json -o preflight.json
```

Score on a GPU:

```bash
metag score --reactions reactions.json --only ester_hydrolysis --out metag_out
```

From Python:

```python
import json
from metag.energetics.uma import load_uma
from metag.pipeline import score_reaction

reaction = json.load(open("reactions.json"))["ester_hydrolysis"]
result = score_reaction(reaction, model=load_uma())
print(result["dG"], result["U_samp"], result["suspect"])
```

Species free energies are cached on disk, keyed by structure, charge, model identity and every
result-affecting setting. Set `METAG_CACHE` to choose the directory. On a shared cluster, launch through
your scheduler or GPU reservation tool rather than setting `CUDA_VISIBLE_DEVICES` by hand.

More inputs are in [examples/](examples/). The full benchmark panel ships as
`src/metag/data/reactions_opentecr_std.json`.

## Results

| Field | Meaning |
|---|---|
| `dG` | ΔrG'° in kJ/mol, rounded to 0.1. `null` when no valid estimate exists; never zero. |
| `suspect` | Reason for a missing estimate (routing error, unbalanced input, unresolved stereochemistry, failed species, charge leak, ...). |
| `U_samp` | Conformer-sampling spread in kJ/mol. This is not a predictive uncertainty. |
| `stages` | Decomposition: `species_sum`, `proton`, `pka_transform`. |
| `species_scored` | The routed chemistry that was computed, next to the input in `original_species`. |
| `routes` | Which routes fired, and per-species sampling and thermal diagnostics. |
| `config` | Fingerprint of model, producer, settings, constant tables and software versions. |

MetaG does not report calibrated prediction intervals.

## Settings

These are read when a backend is created (`metag.settings.Settings.from_env`), not at import.
Scientific settings, the model identity and the xtb binary hash are part of the species-cache key.
Execution settings (xtb workers, GPU chunk sizes) do not change results.

| Variable | Default | Meaning |
|---|---|---|
| `UMA_MODEL` | `uma-s-1p2p1` | UMA checkpoint |
| `UMA_RUNTIME_PROFILE` | `nitro` | `nitro`: UMA Nitro for single-system optimisation on sm_80–sm_120 GPUs, FP32 for everything else. `official_fp32`: FP32 throughout. |
| `UMA_HESSIAN` | `autograd` | Exact Hessian. `fd` uses central differences. |
| `METAG_SAMPLING_PROFILE` | `current` | `thorough`: stricter plateau, larger pools (also `metag score --sampling-profile`). |
| `CONV_TOL`, `CONV_HITS`, `CONV_MAX` | 2.5 kJ, 2, 8 | Conformer-search plateau criterion. |
| `SAMPLE_SCALE`, `KEEP_SCALE` | 1, 1 | Scale the conformer pool and the number of conformers relaxed and solvated. |
| `CONFORMER_SEED_OFFSET` | 0 | Seed offset for independent replicate searches. |
| `XTB_BIN` | `xtb` on `PATH` | xtb executable |
| `METAG_XTB_WORKERS` | 8 | Concurrent single-threaded xtb processes. This changes speed only. |
| `UMA_ENERGY_CHUNK`, `UMA_FIRE_CHUNK`, `UMA_HESS_CHUNK` | | GPU batch sizes, to avoid running out of memory on large species. They change memory use only. |

## Scope and limitations

- **Development exposure.** Discrete method choices were made against TECRDB residuals, so benchmark
  scores are development-exposed. No parameter fitted to reaction free energies is applied.
- **Fixed conditions.** Results are for the stated conditions only. Compare against data transformed to
  them.
- **Dominant error.** The main residual is aqueous solvation of charged and flexible solutes. Electronic
  and thermal terms were checked against independent references.
- **Refused inputs.** Reactions with unspecified stereocentres, isotopes, unsupported phosphoryl
  environments (for example phosphoramidates) or explicit-water requests get no estimate.

## Repository layout

```
src/metag/
  pipeline.py        score_reaction: routing -> species energies -> solution states -> assembly
  routing/           structure-based routes (pure CPU)
  energetics/        UMA, conformers, thermal (RRHO/Hessian), xtb solvation, species cache, backend
  protonation_registry.py, state_pools.py, assembly.py, settings.py, constants.py
  tools/             preflight, evaluate, build-reference
  data/              benchmark panel and source-backed protonation models
benchmark/           Mg-bounded openTECR labels, archived predictions, scorer
examples/            small input files
tests/               CPU unit tests (GPU tests are marked `gpu`)
```

Run the CPU tests with `pip install -e '.[test]' && pytest`; they need a (CPU) Torch install but no
FairChem, UMA weights or xtb.

## License and acknowledgements

MetaG is released under the [MIT License](LICENSE). The UMA weights are not included; they are governed
by the FAIR Chemistry License. MetaG builds on FairChem/UMA, xtb, RDKit, ASE and eQuilibrator. Benchmark
labels derive from openTECR and TECRDB.

A manuscript describing MetaG is in preparation.
