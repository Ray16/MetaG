"""CPU-only reaction/protonation preflight: python -m metag.tools.preflight INPUT -o REPORT."""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path

from rdkit import Chem, rdBase
from metag.reactions import validate_reaction
from metag.routing.pka_transform import pka_model_config
from metag.routing.planner import route_reaction
from metag.validation.protonation import inspect_species
from metag.constants import CONDITIONS, RT_LN10
from metag.protonation_registry import prepare_reaction_protonation
from metag.energetics.protonation_models import reference_assignments


def audit_reactions(reactions):
    from metag.protonation_registry import registry_identity
    from metag.routing.planner import TRUNCATION_RADIUS
    config = dict(pka_model_config(), trunc_radius=TRUNCATION_RADIUS,
                  source_protonation_registry_sha256=registry_identity())
    source_root = Path(__file__).parents[1]
    source_hashes = {str(path.relative_to(source_root)): hashlib.sha256(path.read_bytes()).hexdigest()
                     for path in sorted(source_root.rglob("*.py"))}
    species, records = {}, {}
    for key, reaction in reactions.items():
        try:
            validate_reaction(reaction)
            from metag.chem import isotope_labels
            isotopic = {name: labels for name, (_, _, smi) in reaction["species"].items()
                        if (labels := isotope_labels(smi))}
            if isotopic:
                records[key] = {"status": "blocked", "unsupported_isotopes": isotopic,
                                "error": "isotope-specific mass and solvent/reservoir conventions are not implemented"}
                continue
        except (ValueError, TypeError) as exc:
            records[key] = {"status": "invalid", "error": str(exc)}
            continue
        keys = []
        for name, (_, charge, smi) in reaction["species"].items():
            canonical = Chem.MolToSmiles(Chem.MolFromSmiles(smi))
            keys.append(canonical)
            if canonical not in species:
                species[canonical] = dict(inspect_species(canonical), uses=[])
            species[canonical]["uses"].append({"reaction": key, "name": name, "input_charge": charge})
        try:
            routed, routes, _ = route_reaction(reaction, log=lambda _: None)
            routed_species = {name: inspect_species(smi) for name, (_, _, smi) in routed["species"].items()}
            from metag.routing.sugar_states import pool_aware_stereochemistry
            stereo = pool_aware_stereochemistry(reaction['species'], routed['species'],
                                                routes.get('declared_sugar_pools', {}))
            record = {"routes": routes, "routed_species": routed_species, "stereochemistry": stereo,
                      "species_scored": {name:list(value) for name,value in routed["species"].items()},
                      "n_Hplus_scored": routed.get("n_Hplus", 0),
                      "pka_sites_scored": routed.get("pka_sites")}
            if (routes['input_balanced']
                    and not routes['errors'] and not stereo['input'] and not stereo['routed']):
                try:
                    sites, coupled, diagnostics, registry, references = prepare_reaction_protonation(
                        reaction, routed['species'],
                        routed.get('pka_sites', []), routes.get('protonation_sites_by_species'),
                        CONDITIONS, enabled=True, rt_ln10=RT_LN10,
                        expected_registry_sha256=config.get('source_protonation_registry_sha256'))
                except (ValueError, TypeError, KeyError) as exc:
                    records[key] = dict(record, status='blocked_source_protonation_model', error=str(exc))
                    continue
                routes['source_protonation_registry'] = registry
                if diagnostics is not None:
                    routes['coupled_protonation_models'] = diagnostics
                record.update(pka_sites_scored=sites, coupled_protonation_kj=coupled,
                              protonation_reference_models=references)
            elif not reaction.get('species_protonation_models'):
                record['protonation_reference_models'] = reference_assignments(
                    routed['species'], routes.get('protonation_sites_by_species'))
            if (routes.get('declared_sugar_pools') and routes['input_balanced'] and not routes['errors']
                    and not stereo['input'] and not stereo['routed']):
                from metag.state_pools import prepare_sugar_pools
                try:
                    plans = prepare_sugar_pools(routed, routes, caller_pka_sites=bool(reaction.get('pka_sites')))
                except (ValueError, TypeError, KeyError) as exc:
                    records[key] = dict(record, status='blocked_sugar_ensemble', error=str(exc))
                    continue
                routes['sugar_state_plan'] = plans
            if not routes["input_balanced"] or routes["errors"]:
                status = "blocked"
            elif stereo["input"] or stereo["routed"]:
                status = "blocked_unresolved_stereochemistry"
            elif any(species[s]["status"] == "unsupported" for s in keys) or any(
                    row["status"] == "unsupported" for row in routed_species.values()):
                status = "review_unsupported_species"
            else:
                status = "review" if any(species[s]["issues"] for s in keys) or any(
                    row["issues"] for row in routed_species.values()) else "no_known_gap_detected"
            records[key] = dict(record, status=status)
        except (ValueError, RuntimeError) as exc:
            records[key] = {"status": "blocked", "error": str(exc)}
    return {"scope": "Structural preflight only; no_known_gap_detected is not validated speciation or accuracy.",
            "config": config, "source_sha256": source_hashes, "rdkit_version": rdBase.rdkitVersion,
            "n_reactions": len(reactions), "n_unique_structures": len(species),
            "status_counts": dict(Counter(r["status"] for r in records.values())),
            "issue_counts": dict(Counter(i for r in species.values() for i in r["issues"])),
            "species": species, "reactions": records}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("-o", "--output", type=Path, required=True)
    args = parser.parse_args(argv)
    raw = args.input.read_bytes()
    reactions = json.loads(raw)
    if not isinstance(reactions, dict):
        parser.error("input must be a JSON mapping of reaction identifiers to reactions")
    report = audit_reactions(reactions)
    report.update(input_file=str(args.input.resolve()), input_sha256=hashlib.sha256(raw).hexdigest())
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    print(json.dumps({k: report[k] for k in ("n_reactions", "n_unique_structures", "status_counts", "issue_counts")}, indent=2))


if __name__ == "__main__":
    main()
