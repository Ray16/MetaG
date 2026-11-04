# Examples

| File | Content |
|---|---|
| `ester_hydrolysis.json` | ethyl acetate + water -> acetic acid + ethanol; the README quick-start input |
| `amp_hydrolysis.json` | AMP + H2O -> adenosine + phosphate (pH-0 route, phosphate ladders) |
| `citrate_isocitrate.json` | citrate -> D-threo-isocitrate with both stereocentres specified |
| `stereo_controls.json` | E -> Z isomerisation of 2-butene and 1,2-dichloroethene (no experimental labels) |
| `citrate_NIST_2002_1M.json` | a source-backed sequential pKa model; see [docs/source-backed-protonation.md](../docs/source-backed-protonation.md) |

CPU checks, without the QM dependencies:

```bash
metag preflight examples/ester_hydrolysis.json -o output/ester_preflight.json
```

GPU scoring:

```bash
metag score --reactions examples/ester_hydrolysis.json --out output/ester
```

For replicate searches, run separate processes with different `CONFORMER_SEED_OFFSET` values, for example
0 and 1000, each with its own output directory. These are development examples, not an independent
benchmark.

Preflight reports structure/protonation admission issues rather than molecular energies. The ester example receives a `review` status for its class-based protonation assumptions. A successful preflight is not evidence of complete aqueous-state sampling or predictive accuracy. Scoring requires the [molecular runtime](../environments/README.md) and writes an outcome list to `output/ester/metag_ester_hydrolysis.json`; a refused or failed prediction has `dG: null`.

Preflight replaces its selected report, and scoring replaces the accumulated batch JSON at the selected output path. Choose new output paths to preserve earlier runs. Replicate seeds do not remove model or reference errors.

The [diagnostic examples](diagnostics/) are separate synthetic controls for `benchmark/thermochemistry_diagnostics.py`. They are described in [the diagnostic guide](../docs/physical-diagnostics.md).
