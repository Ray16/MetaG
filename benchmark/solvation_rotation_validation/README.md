# Fixed-pose solvation rotational-invariance control

All 432 requested gas/water pairs were attempted: 426 succeeded and six failed
at gas-phase SCC convergence. Those failures occur for two proper orientations
of the vertically deprotonated oxidized cofactor. No failed energy is zero-filled;
six downstream water calls were not attempted. The original complete-panel
analysis remains withheld and its gate is unchanged.

The frozen panel contains six cofactor/core poses, identity plus 16 seeded
non-cubic proper rotations, a separate inversion control, and four numerical
arms (default/dense SASA grid crossed with default/tighter SCC accuracy).
All original-reference poses reproduce their recorded production solvation
values. A separate attachment audit covers the 103 retained reference geometries.
These controls do not validate protonation populations or Gibbs accuracy.

The worker circuit breaker stopped after three consecutive SCC failures. Its
records remain untouched. After native logs and terminal ownership were verified,
only the 47 previously unattempted tasks continued with identical worker and
protocol bytes under the original deadline; one further SCC refusal was retained.

Reconstruct the complete inventory and finite complete-subcase diagnostics:

    python benchmark/solvation_rotation.py --strict-endpoint-replay

The optional replay executes the original frozen scalar analysis without xTB or
UMA and verifies that the failed complete-panel endpoint remains withheld.
Raw native stdout/stderr, XYZ inputs, charge/phase flags, failures, successful
outputs, task definitions and frozen source are included. The denser grid's
numerical invariance does not establish experimental accuracy, and no production
profile is adopted from this control. Ranges over sampled orientations are not
bounds over all rotations or calibrated total uncertainty intervals.
