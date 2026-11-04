"""Synthetic all-input state-function checks. No physical accuracy is measured.

Usage: python benchmark/audit_state_function.py --cycles input_nullspace.json --output audit.json
The cycle file must contain independent_cycles, each mapping reaction IDs to weights.
"""
import argparse
from collections import Counter
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from rdkit import Chem
from metag.settings import Settings
from metag.state_function import StateFunctionEstimator


class SyntheticBackend:
    model_name = "synthetic-state-function-control"
    producer = {"cacheable": False, "kind": "synthetic-not-physical"}

    def __init__(self):
        self.settings = Settings()
        self.calls = Counter()

    def species(self, charge, smiles, label, warnings, diagnostics):
        self.calls[(charge, smiles)] += 1
        mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
        h = int(hashlib.sha256(smiles.encode()).hexdigest()[:6], 16) % 97
        diagnostics["method"] = "synthetic"
        return -100.*sum(a.GetAtomicNum() for a in mol.GetAtoms()) + h + self.calls[(charge, smiles)], .5

    def liquid_water(self):
        return -999.


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reactions", type=Path, default=ROOT/"src/metag/data/reactions_opentecr_std.json")
    parser.add_argument("--cycles", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error("output already exists")
    panel = json.loads(args.reactions.read_text())
    cycles = json.loads(args.cycles.read_text())["independent_cycles"]
    estimator = StateFunctionEstimator(SyntheticBackend())
    results = {rid: estimator.score_reaction(rx, key=rid) for rid, rx in panel.items()}
    reversal, scaling = {}, {}
    for rid, rx in panel.items():
        if results[rid]["dG"] is None:
            continue
        for factor, errors in ((-1, reversal), (2.5, scaling)):
            changed = dict(rx, species={n: [factor*c, q, s] for n, (c, q, s) in rx["species"].items()},
                           n_Hplus=factor*rx.get("n_Hplus", 0))
            alt = estimator.score_reaction(changed, key=rid)
            errors[rid] = (None if alt["dG"] is None else
                           alt["dG_unrounded"] - factor*results[rid]["dG_unrounded"])
    cycle_rows = []
    for cycle in cycles:
        # Verify declared cycles independently in the original input-state basis.
        states, protons = Counter(), Fraction(0)
        for rid, weight in cycle.items():
            for name, (coefficient, charge, smiles) in panel[rid]["species"].items():
                smi = Chem.MolToSmiles(Chem.MolFromSmiles(smiles), isomericSmiles=True)
                # Declared pools and protonation models are part of the input
                # thermodynamic state, not interchangeable with a fixed drawing.
                identity = json.dumps({"charge": charge, "smiles": smi,
                    "pool": panel[rid].get("species_state_pools", {}).get(name),
                    "protonation": panel[rid].get("species_protonation_models", {}).get(name)}, sort_keys=True)
                states[identity] += Fraction(str(weight))*Fraction(str(coefficient))
            protons += Fraction(str(weight))*Fraction(str(panel[rid].get("n_Hplus", 0)))
        if protons or any(states.values()):
            raise ValueError("supplied cycle does not close on original states")
        missing = [rid for rid in cycle if results[rid]["dG"] is None]
        residual = None if missing else math.fsum(w*results[rid]["dG_unrounded"] for rid, w in cycle.items())
        cycle_rows.append({"cycle": cycle, "unavailable": missing, "residual_kj": residual})
    evaluated = [r["residual_kj"] for r in cycle_rows if r["residual_kj"] is not None]
    result = {"scope": "Synthetic bookkeeping checks only; these energies are not physical predictions or coverage validation.",
              "inputs_sha256": hashlib.sha256(args.reactions.read_bytes()).hexdigest(),
              "cycles_sha256": hashlib.sha256(args.cycles.read_bytes()).hexdigest(),
              "n_requested": len(panel), "n_synthetic_estimates": sum(r["dG"] is not None for r in results.values()),
              "unavailable": {rid: r["suspect"] for rid, r in results.items() if r["dG"] is None},
              "n_input_cycles": len(cycles), "n_evaluable_cycles": len(evaluated),
              "max_cycle_residual_kj": max(map(abs, evaluated), default=None), "cycles": cycle_rows,
              "max_reversal_error_kj": max((abs(v) for v in reversal.values() if v is not None), default=None),
              "max_scaling_error_kj": max((abs(v) for v in scaling.values() if v is not None), default=None),
              "reversal_failures": [rid for rid, v in reversal.items() if v is None or abs(v)>1e-7],
              "scaling_failures": [rid for rid, v in scaling.items() if v is None or abs(v)>1e-7],
              "raw_backend_max_calls_per_identity": max(estimator.backend.calls.values(), default=0),
              "ledger_sha256": estimator.ledger()["sha256"],
              "physical_validation": "pending; do not freeze a prediction benchmark from this check"}
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write("\n")
    print(json.dumps({k:v for k,v in result.items() if k not in {"cycles", "unavailable"}}, indent=2))
    if any(abs(v)>1e-7 for v in evaluated) or result["reversal_failures"] or result["scaling_failures"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
