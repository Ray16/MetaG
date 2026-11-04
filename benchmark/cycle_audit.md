# Reaction-cycle consistency diagnostic

A standard biochemical reaction energy must be additive when the same chemical states, reservoirs
and conventions are used throughout a cycle. The archived panel uses reaction-specific routing;
balanced rewrites do not guarantee the resulting reaction energies come from one shared set of
species energies. The following **development-panel diagnostics** compare rounded archived 0.3.0
predictions. They are physical-consistency tests, not evidence that any one prediction is accurate.

| Cycle (signed sum of reactions) | Archived sum | Maximum rounding contribution |
|---|---:|---:|
| rxn00577 + rxn12158 − rxn14104 | +4.4 kJ/mol | 0.15 kJ/mol |
| rxn00191 + rxn00182 − rxn00278 | +2.9 kJ/mol | 0.15 kJ/mol |

The [cycle checker](check_cycle.py) verifies that the original reaction structures, formal charges,
stoichiometric coefficients and explicit proton counts cancel, and that repeated identities have
the same declared sugar-pool metadata. The tests use the packaged structure and predictions, including
the documented mixed dev24/0.3.0 provenance. Energy sums are larger than the most that 0.1 kJ/mol
rounding could account for.

In the first cycle, the archived route for rxn12158 truncates the glucose-containing species, while
rxn00577 and rxn14104 preserve their full sugar structures. This route difference is a plausible
contributor to the nonclosure, but this diagnostic does not isolate it: independently sampled
species, pKa assignments and other approximations could also contribute. In the second cycle,
rxn00191 was recomputed for 0.3.0 while the other two use dev24 results; their routing and protonation
assignments differ. The hybrid provenance and reaction-specific state models should be tested before
attributing that residual to a specific component.

The cycles were discovered by exploratory analysis of the development panel. They are **not** a
prespecified independent validation set. A fix must not force a zero sum by changing reaction values
without addressing the underlying state and routing inconsistency. A suitable follow-up is a controlled
shared-backend replay and matched full-versus-truncated calculation for all three members, with
identical scientific settings and sampling replicates; report both cycle closure and accuracy changes.

The shared-backend replay has now been performed with deterministic artificial energies for every
canonical species. For the sucrose cycle it closes to numerical precision with truncation disabled;
with normal routing, rxn12158 is truncated while the other two reactions are not, and it does not
close. This demonstrates that the route structure can cause nonadditivity even when the underlying
species callback is consistent. It does not measure the physical size of that effect. The alanine /
glutamate cycle remains nonclosed without truncation: the same amino acids are scored as protonated
in rxn00191 and neutral in the two dehydrogenase reactions. Protonation ownership and reference-state
consistency therefore need a separate fix. Arbitrary replay values are not physical energy estimates.

An experimental replay override that forces neutral-base routing closes the amino-acid cycle to
within its 0.15 kJ/mol rounding bound when truncation is disabled. It makes the relevant amino-acid
reference states match across the three reactions. This establishes a possible consistency remedy,
not its chemical accuracy. The production base gate has exclusions for amide and other nitrogen
chemistry; forcing it universally could assign inappropriate pKa terms. A CPU impact audit with
truncation disabled found altered species, pKa sites or proton counts in 64 of 364 panel reactions.
Both baseline and override produced routing errors for six inputs, with no newly introduced routing
error IDs. This structural check does not validate any altered energy or protonation model.

A narrower experimental screen in `experimental_transamination.py` restricts the override to balanced
four-species exchanges containing one primary alpha-ammonium species and one alpha-oxo acid per side.
On the 364 input structures it changes nine routes and introduces no new routing-error IDs. The
deterministic artificial-energy replay of the alanine/glutamate cycle then sums to 0.1 kJ/mol with
or without truncation, within the 0.15 kJ/mol rounding bound. This is an exploratory bookkeeping
diagnostic, not a physical-energy validation or an independent test; the affected chemistry was
identified in the development panel. The rule remains outside the production workflow pending
matched state, pKa and energy checks that do not select it by benchmark residuals.

Reproduce the structural and rounded-energy checks:

```
python benchmark/check_cycle.py --cycles benchmark/reaction_cycles.json --output /tmp/metag_cycles.json
PYTHONPATH=src python benchmark/replay_cycles.py --output /tmp/metag_cycle_replay.json
PYTHONPATH=src python benchmark/replay_cycles.py --force-base --output /tmp/metag_cycle_forced.json
PYTHONPATH=src python benchmark/compare_protonation_routes.py --output /tmp/metag_route_impact.json
PYTHONPATH=src python benchmark/evaluate_transamination_route.py --output /tmp/metag_transamination.json
```

The checker does not verify the runtime backend, automatically assigned pKa ownership or thermodynamic
equivalence of structural routes. It refuses caller-supplied protonation models instead of silently
equating them. Inconsistent declared sugar pools are reported separately from stoichiometric closure.
