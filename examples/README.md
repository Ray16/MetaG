# Examples

| File | Content |
|---|---|
| `amp_hydrolysis.json` | AMP + H2O -> adenosine + phosphate (pH-0 route, phosphate ladders) |
| `citrate_isocitrate.json` | citrate -> D-threo-isocitrate with both stereocentres specified |
| `stereo_controls.json` | E -> Z isomerisation of 2-butene and 1,2-dichloroethene (no experimental labels) |
| `citrate_NIST_2002_1M.json` | a source-backed sequential pKa model; see [docs/source-backed-protonation.md](../docs/source-backed-protonation.md) |

CPU checks, without the QM dependencies:

```bash
metag preflight examples/amp_hydrolysis.json -o amp_preflight.json
```

GPU scoring:

```bash
metag score --reactions examples/amp_hydrolysis.json --out amp_results
```

For replicate searches, run separate processes with different `CONFORMER_SEED_OFFSET` values, for example
0 and 1000, each with its own output directory. These are development examples, not an independent
benchmark.
