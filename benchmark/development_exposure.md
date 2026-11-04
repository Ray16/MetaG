# Development exposure and selection-audit status

## Final guarded finite-policy audit

The completed retrospective exercise used 901 declared policies on 243 guarded
primary-label targets. Of these, 525 had complete outcomes (521 distinct prediction
vectors) and 376 had at least one recorded refusal; no otherwise viable policy
remained unresolved. Across twenty grouped five-fold outer partitions per grouping,
mean apparent-to-outer MAE differences were 0.15 and 0.17 kJ/mol, with split ranges
0–0.88 and 0–1.06, respectively. These are descriptive ranges, not confidence
intervals. No policy was promoted. The exercise covers this finite catalog and
observed cohort, not undocumented human choices or total historical optimism.
It is not independent validation. The original evaluated source and outcomes are
retained in the manuscript's guarded selection evidence and local portable archive;
a reviewed public archive remains pending.

The independent-audit corrections to canonical graph routing, light-species
linearity and strong-imaginary-mode acceptance postdate that evaluated source.
Their changed cache physics identity prevents old molecular records being
advertised as new calculations. Their accuracy, coverage and selection audit have
not yet been rerun; the numbers above must not be attributed to the corrected code.

## Earlier incomplete frozen-source follow-up (historical)

The analysis below predates the MCS timeout guard in the current working code;
that guarded estimator requires a separate full evaluation. A fresh frozen-source
calculation returns 245 estimates on the primary
Mg-bounded reference cohort. The full finite policy audit is still incomplete.
An availability-limited secondary analysis used 521 candidate policies with
complete outcomes on all 245 targets, from a catalog of 901; 352 policies had
at least one observed workflow refusal and 28 otherwise viable policies were
unresolved at the analysis snapshot. Across 20 grouped five-fold outer
partitions with inner selection, the mean apparent-to-outer MAE difference was
0.15 kJ/mol (split range 0.00–0.87); a second grouping gave 0.14 kJ/mol
(range 0.00–0.42). These are descriptive split ranges, not confidence
intervals. The analysis subset was fixed after availability inspection and
before those errors were scored.

This analysis is conditional on the observed outcomes and does not bound
selection optimism for all catalogued policies, undocumented human decisions,
or the full development history. It is not an independent prospective test,
and no policy was adopted from its results. Compact evidence and fold
verification are included in the revised manuscript's supporting evidence;
the full candidate matrix and audit sources still need a separately reviewed
public archive.

## Earlier development history

MetaG's reported 0.3.0 reaction errors are development-exposed. No numerical
reaction-energy offset is fitted to TECRDB in the scorer, but reaction
residuals informed discrete choices in earlier development. The absence of a
fitted offset does not make the benchmark a held-out test.

The earlier dev24 inventory recorded direct benchmark-outcome evidence for
NTP core handling, the pH-0 base and isomerase gates, ring-chain handling,
aldehyde hydration, and two truncation policies. For other choices the
available history was often insufficient to establish when or why a choice was
adopted. The [historical decision log, interpretation and inventory table](history/dev24_selection/)
are retained here with hashes. The full historical JSON inventory and frozen
source remain in the archived development workspace; they are not bundled in
this repository. These are historical policy names, not switches supported by
the released 0.3.0 interface.

The [current-source crosswalk](decision_crosswalk.md) maps those historical
decisions to the refactored code and identifies alternatives that have not been
reconstructed. Its analysis-only prototype controls are not part of the
released workflow and do not constitute a completed selection audit.

The 0.3.0 source history records four later correctness changes: solvation
failure handling (`ba9a18d`), heavy-atom symmetry fallback (`2b9b693`), the
redox-proton guard (`f7e4037`), and direction-invariant routing (`0b48206`).
The commit messages and code state their technical rationale; they do not, by
themselves, establish complete blindness to all previous benchmark results.
The 24 affected panel reactions were subsequently recomputed, as described in
[validation](../docs/validation.md). The panel combines these results with
previously evaluated predictions for the other reactions.

The previously declared 426-arm audit was designed for dev24. It did not
finish, and the production refactor removed several historical alternatives.
Its intermediate results cannot quantify the selection optimism of the 0.3.0
workflow. A current-version audit needs a prespecified, executable policy space,
fixed target cohorts, retained failures and grouped outer evaluation, followed
by a genuinely independent prospective reference set. Even a complete finite
policy audit cannot reconstruct undocumented human choices. The magnitude of total selection optimism remains unknown, including after
the limited current-source sensitivity analysis above.

An outcome-blind routing inventory of the 364 input reactions found 448 unique
non-water structures in the 358 routes without a recorded routing error. The
recent hydrolase study's current-source caches contain structure-and-charge
matches for 57 of those 448. The other 391 were not found in that study's
caches; this is not an inventory of every historical cache, and it excludes
auxiliary hydration and sugar states. A structure match does not by itself
authorize reuse: settings, source, model and producer identities must also
match. Earlier full-panel and 0.3.0 correction caches carry different raw-source
fingerprints. The inventory does not compute selection statistics.
