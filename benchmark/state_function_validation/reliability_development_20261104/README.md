# Developmental reliability feasibility

This analysis uses the existing 69 exact-family matched targets. MetaG and the
reference labels were already exposed during development. The comparator spread
is the population standard deviation across 100 archived author-selected CV
repeats, normalized by absolute reaction extent. It is not the released
pretrained-model uncertainty, and bridge uncertainty is excluded.

The four ranking policies were fixed before calculating these new ranking
outcomes. They retain the same number of targets, averaging across all exact-size
subsets of a tied boundary tier. The target of the reliability comparison is the
error of the comparator mean. MetaG and zero errors are also evaluated on each
identical selected subset. No predictor, threshold or fitted combination is
adopted from these results.

| Ranking | MAE at 18/69 | MAE at 35/69 | MAE at 52/69 | MAE at 69/69 | AURC |
|---|---:|---:|---:|---:|---:|
| CV repeat spread | 2.12 | 3.28 | 3.66 | 4.59 | 2.85 |
| MetaG–comparator disagreement | 2.14 | 2.08 | 2.61 | 4.59 | 2.38 |
| Equal-rank combination | 1.38 | 2.15 | 3.13 | 4.59 | 2.28 |
| Comparator magnitude control | 2.95 | 3.11 | 3.47 | 4.59 | 3.15 |

Units are kJ/mol. AURC is the mean comparator MAE over every integer retained
count from 1 to 69. Disagreement has lower aggregate selective error, but is
slightly worse than repeat spread at 18 retained targets. Of three comparator
errors above 20 kJ/mol, disagreement retains none at 52 targets while repeat
spread retains one. This does not establish zero future risk or superiority to
the deployed uncertainty measure.

The original 67,500 selected prediction values were authenticated against the
retained vector fingerprints. Reassembled means differ from the preceding
comparison by at most 7.2e-15 kJ/mol. A separate implementation reproduces 2,484
risk entries, all boundary-tie expectations and all four AURCs without importing
the evaluator. These checks establish retained arithmetic, not reference or
experimental independence.

Source/family dependence clusters have not been curated for the aggregate
labels. No reaction-independent confidence intervals or significance claims are
made. The evaluator offers conditional paired cluster intervals only if complete
clusters are provided; fewer than ten clusters yields an inconclusive status.
Ten clusters is a minimum gate, not proof of adequate power or independence.

Run from the source checkout with caller-owned inputs:

```bash
python benchmark/reliability.py input.json --output report.json
python benchmark/plot_reliability.py report.json figure_output_stem
```

Input rows require distinct `reaction` IDs and finite `reference`, `metag`,
`comparator` and nonnegative `spread` values. Include `cohort_exposure` and
`spread_definition` in the input metadata. Optional `cluster` values on every
row and a justified `cluster_definition` enable conditional intervals. Metadata
does not itself prove the supplied clustering or experimental independence.

The current protocol and aggregate report are retained here. Per-target inputs,
all curves and portable arithmetic checks are included with manuscript SI S15.
This development result motivates an independently frozen reliability test;
that requirement remains open.
