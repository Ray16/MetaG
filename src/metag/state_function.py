"""Experimental, full-molecule state-function estimator.

Each molecular input has one retained potential, independent of reaction class,
direction, neighbours and stoichiometric extent. This establishes bookkeeping
consistency; it does not validate the inherited protonation or solvent model.
The established score_reaction entry point and its benchmark are unchanged.
"""
import copy
import json
import math

from rdkit import Chem

from metag.chem import is_balanced, is_water, require_supported_isotopes
from metag.constants import CONDITIONS, G_HPLUS, PH, RT_LN10, STD_STATE_KJ
from metag.energetics.errors import SpeciesRearranged
from metag.energetics.protonation import reaction_transform
from metag.energetics.provenance import content_hash
from metag.reactions import validate_reaction

POLICY = "full-molecule-species-potential-v1"


def _finite(value):
    return (isinstance(value, (int, float)) and not isinstance(value, bool)
            and math.isfinite(value))


def plan_state(charge, smiles):
    """Plan one state without consulting any reaction or experimental label.

    The reference conversion is accepted only if its explicit proton exchange
    preserves both elements and charge. Global basic-site treatment uses the
    existing class constants; its physical accuracy remains to be validated.
    """
    validate_reaction({"species": {"state": [1, charge, smiles]}})
    require_supported_isotopes(smiles)
    from metag.validation.stereochemistry import unresolved_species
    if unresolved_species({"state": [1, charge, smiles]}):
        raise ValueError("unresolved stereochemical identity")
    canonical = Chem.MolToSmiles(Chem.MolFromSmiles(smiles), isomericSmiles=True)
    from metag.routing.pka_transform import _neutralize_v2
    reference, acids, bases, reference_charge = _neutralize_v2(canonical)
    if reference is None:
        raise ValueError("molecular reference construction failed")
    reference = Chem.MolToSmiles(Chem.MolFromSmiles(reference), isomericSmiles=True)
    exchange = int(charge) - int(reference_charge)
    conversion = {"input": [-1, charge, canonical],
                  "reference": [1, reference_charge, reference]}
    if not is_balanced(conversion, exchange):
        raise ValueError("molecular reference conversion does not conserve elements and charge")
    from metag.validation.protonation import inspect_species
    return {"input_charge": int(charge), "input_smiles": canonical,
            "reference_charge": int(reference_charge), "reference_smiles": reference,
            "proton_exchange": exchange,
            "pka_sites": [["prod", p, "acid", 1.] for p in acids]
                         + [["prod", p, "base", 1.] for p in bases],
            "protonation_admission": inspect_species(canonical)}


class StateFunctionEstimator:
    """One serial calculation session with immutable, exportable state records.

    Full molecules are scored; no cofactor replacement, truncation, reaction
    gate or reaction-wide fallback is applied. Unsupported state declarations
    fail explicitly. Reuse this instance across a network so a stochastic
    backend cannot assign different values to repeated states. An exported
    ledger is evidence, not an importable molecular cache.
    """
    def __init__(self, backend, *, log=lambda *_: None):
        from metag.pipeline import effective_config
        self.backend, self.log = backend, log
        self.config = dict(effective_config(backend, allow_truncate=False),
                           state_function_policy=POLICY, base_reference="global",
                           cofactor_replacement=False, nucleotide_replacement=False,
                           physical_validation="pending")
        self._backend_identity = self._identity()
        self._records = {}
        self._raw = {}
        self._water = None

    def _identity(self):
        settings = getattr(self.backend, "settings", None)
        return content_hash({"model": self.backend.model_name, "producer": self.backend.producer,
                             "settings": settings.as_dict() if settings is not None else None})

    def _check_identity(self):
        if self._identity() != self._backend_identity:
            raise ValueError("backend identity/settings changed within the state ledger")

    def _water_energy(self):
        if self._water is None:
            value = self.backend.liquid_water()
            if not _finite(value):
                raise ValueError("invalid liquid-water reference")
            self._water = value
        return self._water

    def _species_energy(self, charge, smiles, label, warnings, diagnostics):
        canonical = Chem.MolToSmiles(Chem.MolFromSmiles(smiles), isomericSmiles=True)
        key = (int(charge), canonical)
        if key not in self._raw:
            detail = {}
            value, sigma = self.backend.species(charge, canonical, "state", warnings, detail)
            if not _finite(value) or not _finite(sigma) or sigma < 0:
                raise ValueError("invalid molecular energy or sampling diagnostic")
            self._raw[key] = json.dumps([value, sigma, detail], sort_keys=True, allow_nan=False)
        value, sigma, detail = json.loads(self._raw[key])
        diagnostics.update(detail)
        return value, sigma

    def state(self, charge, smiles):
        """Return a copy of a retained state record; callers cannot mutate it."""
        self._check_identity()
        plan = plan_state(charge, smiles)
        key = content_hash({"charge": plan["input_charge"], "smiles": plan["input_smiles"]})
        if key not in self._records:
            record = self._evaluate(plan)
            self._check_identity()
            record["state_id"] = key
            record["record_sha256"] = content_hash(record)
            self._records[key] = json.dumps(record, sort_keys=True, allow_nan=False)
        return json.loads(self._records[key])

    def _evaluate(self, plan):
        from metag.energetics.aqueous_states import hydrate_all_sites
        from metag.pipeline import _acid_transform
        from metag.protonation_registry import prepare_reaction_protonation
        q, smi = plan["reference_charge"], plan["reference_smiles"]
        sites = plan["pka_sites"]
        original = {"species": {"state": [1, plan["input_charge"], plan["input_smiles"]]}}
        species = {"state": [1, q, smi]}
        remaining, coupled, model, registry, _ = prepare_reaction_protonation(
            original, species, sites, {"state": sites}, CONDITIONS, enabled=True,
            rt_ln10=RT_LN10,
            expected_registry_sha256=self.config["source_protonation_registry_sha256"])
        routes = {"warnings": []}
        detail = {}
        if is_water(smi, q):
            value, sigma = self._water_energy(), 0.
        else:
            value, sigma = self._species_energy(q, smi, "state", routes["warnings"], detail)
        energies, sigmas = {"state": value + STD_STATE_KJ}, {"state": sigma}
        # This local context contains only this molecule's own reference terms.
        rx = {"species": species, "pka_sites": sites}
        hydrate_all_sites(self, rx, "state", q, smi, energies, sigmas, STD_STATE_KJ,
                          self.log, routes, acid_transform=_acid_transform)
        if routes.get("incomplete_hydration"):
            raise ValueError("incomplete molecular hydration ensemble")
        transform = reaction_transform(remaining, PH, RT_LN10) + coupled
        reservoir = plan["proton_exchange"] * G_HPLUS
        potential = energies["state"] + transform + reservoir
        if not _finite(potential) or not _finite(sigmas["state"]):
            raise ValueError("nonfinite transformed molecular potential")
        return {"plan": plan, "potential_kj": potential, "sampling_sigma_kj": sigmas["state"],
                "solution_reference_kj": energies["state"], "protonation_kj": transform,
                "reference_proton_kj": reservoir, "molecular_diagnostics": detail,
                "aqueous_diagnostics": routes, "source_protonation": registry,
                "coupled_protonation": model}

    # Hydration callbacks also use the same retained raw molecular values.
    def species(self, charge, smiles, label, warnings, diagnostics):
        return self._species_energy(charge, smiles, label, warnings, diagnostics)

    def liquid_water(self):
        return self._water_energy()

    def score_reaction(self, reaction, *, key="rxn"):
        validate_reaction(reaction)
        self._check_identity()
        result = {"reaction": key, "dG": None, "dG_unrounded": None, "U_samp": None,
                  "conditions": dict(CONDITIONS), "config": json.loads(json.dumps(self.config)),
                  "original_species": copy.deepcopy(reaction["species"]), "state_records": {},
                  "suspect": None, "uncertainty_scope": "sampling only; no calibrated total-error interval"}
        for field in ("pka_sites", "species_state_pools", "species_protonation_models"):
            if reaction.get(field):
                result["suspect"] = "unsupported in state-function profile: " + field
                return result
        if not is_balanced(reaction["species"], reaction.get("n_Hplus", 0)):
            result["suspect"] = "input reaction not element/charge balanced"
            return result
        terms = [reaction.get("n_Hplus", 0) * G_HPLUS]
        groups = {}
        for name, (coefficient, charge, smiles) in reaction["species"].items():
            try:
                record = self.state(charge, smiles)
            except (ValueError, SpeciesRearranged) as exc:
                result["suspect"] = f"{name}: {exc}"
                result["failed_state"] = {"name": name, "charge": charge, "smiles": smiles}
                return result
            result["state_records"][name] = record
            terms.append(coefficient * record["potential_kj"])
            # Different proton drawings can share exactly the same sampled
            # reference. Their reservoir/protonation offsets are deterministic.
            reference = record["plan"]
            identity = (reference["reference_charge"], reference["reference_smiles"])
            group = groups.setdefault(identity, [0., record["sampling_sigma_kj"]])
            if not math.isclose(group[1], record["sampling_sigma_kj"], rel_tol=0, abs_tol=1e-12):
                raise ValueError("shared reference has inconsistent sampling diagnostics")
            group[0] += coefficient
        total = math.fsum(terms)
        sigma = math.sqrt(math.fsum((c * s)**2 for c, s in groups.values()))
        if not _finite(total) or not _finite(sigma):
            result["suspect"] = "reaction assembly exceeded finite numeric range"
            return result
        result.update(dG=round(total, 1), dG_unrounded=total, U_samp=sigma)
        return result

    def ledger(self):
        """Export copies of the complete retained raw and transformed records."""
        payload = {"policy": POLICY, "config": self.config,
                   "states": {k: json.loads(v) for k, v in sorted(self._records.items())},
                   "raw_molecules": [{"charge": q, "smiles": s, "record": json.loads(v)}
                                     for (q, s), v in sorted(self._raw.items())],
                   "liquid_water_kj_before_standard_state": self._water}
        return dict(json.loads(json.dumps(payload)), sha256=content_hash(payload))
