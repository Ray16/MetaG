"""Conservative structural coverage diagnostics; no QM or experimental reaction energies.

A built neutral reference is not evidence of complete speciation. This diagnostic
identifies some known gaps, not every possible ionizable group.
"""
from rdkit import Chem
from metag.routing import pka_transform as pk


def inspect_species(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return {"status": "invalid", "issues": ["invalid_smiles"]}
    canonical = pk._canonicalize_maxanion(smiles)
    cmol = Chem.MolFromSmiles(canonical)
    issues = set()
    anions = pk._anion_sites(cmol)
    covered = {i for i, _ in anions}
    if any(a.GetFormalCharge() < 0 and a.GetIdx() not in covered for a in cmol.GetAtoms()):
        issues.add("unmodeled_anionic_site")
    if any(a.GetAtomicNum() not in {1, 6, 7, 8, 9, 15, 16, 17, 35, 53} for a in mol.GetAtoms()):
        issues.add("element_outside_protonation_rules")
    guan_c = {m[1] for m in cmol.GetSubstructMatches(pk._GUAN_C)}
    for a in cmol.GetAtoms():
        if a.GetSymbol() == "P" and any(n.GetSymbol() == "N" for n in a.GetNeighbors()):
            issues.add("phosphoramidate_pka_source_unverified")
        if a.GetSymbol() != "N":
            continue
        if a.GetIsAromatic():
            issues.add("aromatic_n_speciation_incomplete")
        if any(n.GetIdx() in guan_c for n in a.GetNeighbors()):
            issues.add("guanidine_amidine_speciation_incomplete")
        if a.GetFormalCharge() == 1 and a.GetTotalNumHs() > 0 and not a.GetIsAromatic():
            heavy = [n for n in a.GetNeighbors() if n.GetAtomicNum() != 1]
            if heavy and not all(n.GetAtomicNum() == 6 and n.GetHybridization() == Chem.HybridizationType.SP3
                                 for n in heavy):
                issues.add("protonated_n_uses_broad_class_approximation")
    if len(anions) > 1:
        issues.add("multiple_acid_sites_independent_site_approximation")
    try:
        reference, acids, bases, charge = pk._neutralize_v2(smiles)
        if reference is None:
            return {"status": "unsupported", "issues": sorted(issues | {"reference_build_failed"})}
    except (ValueError, RuntimeError) as exc:
        return {"status": "unsupported", "issues": sorted(issues), "error": str(exc)}
    if acids or bases:
        issues.add("class_pka_transfer_unvalidated" if pk.pka_model_config()["pka_model"] == "table"
                   else "site_model_transfer_unvalidated")
    if charge:
        issues.add("reference_retains_charge")
    return {"status": "review" if issues else "no_known_gap_detected", "issues": sorted(issues),
            "reference_smiles": reference, "reference_charge": charge,
            "acid_pkas": acids, "base_pkas": bases}
