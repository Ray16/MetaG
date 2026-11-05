# Frozen dGbyG stereo-feature audit

This version-specific representation audit does not enroll experimental labels,
validate a final MetaG policy or establish independent predictive superiority.
The preliminary feature omission and native-vector equalities were already
observed before the structural audit freeze. All 17 formation encoding pairs and
all 12 stereo-distinct same-connectivity pairs in the archived native molecular
bank are retained, including unresolved stereochemistry and enantiomers. Equality
of enantiomer energies is appropriate in an achiral environment.

The exact active API and feature/network files are copied from the previously
frozen native-ensemble experiment; hashes are bound to its protocol. Other dGbyG
versions, optional stereochemical features and other graph architectures are
outside this audit. The upstream MIT license is in owned_source/LICENSE.
source_run.py preserves the original generating script and its historical paths;
it is not the portable entry point. Original SMILES are in inputs/ and retain the
atom ordering needed to reproduce every archived tensor digest.

Run `python audit_stereo_features.py --packet .` from this directory using RDKit.
This independently reconstructs the actual default graph arrays and replays the
archived 100-member contrasts without Torch, models, network access or inference.
With Torch/PyG installed, add `--actual-featurizer` to execute the frozen feature
function on CPU and compare every array to the independent reconstruction. The
original runtime was Torch 2.8.0+cu128, PyG 2.8.0.post1 and RDKit 2026.03.5; CUDA was
not initialized. Neither mode imports the comparator API or loads its weights.

All five graph-cancelled development reactions are retained in the retrospective
error join: historical MetaG MAE 2.61 kJ/mol versus comparator/zero 1.03 kJ/mol.
These are historical corrected-source MetaG values, not the newer state-function
research profile. Existing reference-condition, training and development exposure
qualifications remain. Zero member variation is not a calibrated error interval;
its omission of representation/state/condition errors does not itself establish
an independent reliability advantage. Explicit chemical stereoisomers and a
complete equilibrating state pool are different targets. Any future validation
must match those definitions, use independent experimental references and include
both zero and a stereo-aware graph-model control.
