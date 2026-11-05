# Completed fresh liquid-water reference

A fresh reference completed on one gate-owned V100 using the unchanged checkpoint,
source, atom references and existing water recipe. It has zero imaginary internal
modes and rotational symmetry number2. All three startup point checks pass the
0.001kJ/mol compatibility gate. The actual V100 runtime and producer are retained;
optimization uses the generalFP32 fallback rather than H100's preferred Nitro.
Final energies and autograd Hessians use generalFP32 on both. Point checks do not
prove identical optimization or Hessians across devices.

The original water-wrapper ImportError and the subsequent H100 resource
interruption remain explicit, separate outcomes. The successful V100 reference
supports only a separately labeled corrected-wrapper assembly. It does not erase
either failed attempt, rescue original population inference, or supply an
experimental accuracy result. The original21:30UTCdeadline/eightGPUhceiling remain.

Run `python -m benchmark.audit_operational_water_reference` from MetaG. The auditor
checks all packet bytes, binds the pre-computation source freeze, authenticates the
raw geometry/components/thermochemistry, reconstructs the partition-independent
water sum in high precision, checks source and hardware-aware producer identities,
and verifies terminal worker/predictor/gate identities. It performs no model
inference or network access. The prior failed attempts' recorded raw-job timing
is retained; startup time and reservation occupancy are separately scoped.
