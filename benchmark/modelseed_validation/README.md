# ModelSEED input admission

This is a CPU-only preflight of the literal annotations in ModelSEED revision
[`b539adb69a95`](https://github.com/ModelSEED/ModelSEEDDatabase/tree/b539adb69a955dac15cc30dbeb3debbd87739018).
It uses the experimental full-molecule state policy. Its counts are not physical
prediction coverage and do not inherit the reduction-based benchmark's accuracy.

| Stage | Reactions |
|---|---:|
| Source entries | 56,006 |
| Obsolete entries excluded first | 7,622 |
| Live transport entries excluded | 4,734 |
| Requested live nontransport inputs | 43,650 |
| Pass recorded input gates | 6,572 |
| Require curation or are unsupported | 37,078 |
| Passes with protonation warnings | 6,109 |
| Passes without recorded protonation warnings | 463 |

The warnings include unvalidated class-pKa transfers, independent acid-site
approximations and incomplete nitrogen speciation. Absence of a warning is not
validation. All source status flags and overlapping failure causes are retained.
First-cause counts are mutually exclusive; any-cause counts overlap. A sole
recorded failure category does not estimate restored coverage because later
checks can be hidden by an earlier failure.

An independent implementation checks source membership, reaction accounting and
reference element/proton/charge conversion without importing MetaG planning or
thermochemistry. Literal balances are known for 23,742 nonempty equations:
19,968 balance and 3,774 do not. Of those imbalances, 2,226 were hidden behind
unavailable molecular states in the original preflight. This does not establish
that an intended biochemical transformation is wrong; missing reservoir
conventions and source annotations require curation.

The exploratory ordinary-water/proton basis admits 301 uniquely balanced
completion candidates. In 85 cases the other recorded state gates pass; one
completion collapses the net graph. All remain curation candidates. No source
equation is silently changed and none becomes an additional energy prediction.

Use the recorded environment and the pinned, caller-owned source tables:

```bash
python benchmark/audit_modelseed_inputs.py \
  --modelseed /path/to/ModelSEEDDatabase --output /path/to/new_preflight
python benchmark/verify_modelseed_inputs.py \
  --modelseed /path/to/ModelSEEDDatabase \
  --retained /path/to/new_preflight --output /path/to/new_verification
```

The audit limits native execution and affinity to one CPU and does no molecular
energy work. The primary run took 159 seconds; this is not a quantum timing or
speedup benchmark. Published-runtime admission records match the preceding
working-source snapshot byte for byte. That does not prove numerical
thermochemistry equivalence.

The concrete-atom guard rejects wildcard/dummy molecules before molecular work.
The existing 364 benchmark input structures contain none and remain admitted by
that guard. Unrelated pending engineering changes are not included in this
runtime repair. Full source fingerprints, derived records and additional
verification are retained with manuscript SI S16; upstream chemical tables and
model weights are not redistributed here. Source licenses remain with their
owners; see the [upstream notices](https://github.com/ModelSEED/ModelSEEDDatabase/blob/b539adb69a955dac15cc30dbeb3debbd87739018/LICENSE).

The database is described by [Seaver et al.](https://doi.org/10.1093/nar/gkaa746).
Numerical genome-scale estimates, experimental accuracy, total uncertainty and
physiological direction remain unassessed by this preflight.
