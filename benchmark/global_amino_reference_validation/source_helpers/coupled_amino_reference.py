"""Analysis prototype: complete source ladders on fully protonated primary-amino acids.

This module does not supply empirical constants or modify the production package.
Use one registry per reaction. A neutral NH2/carboxylic acid is NOT the reference
for a whole macroscopic ladder: the reference here has NH3+ and neutral COOH.
Source provenance is caller declared and must be reviewed separately; validation
of this schema cannot authenticate a publication or certify aqueous completeness.
"""
from collections import Counter
import copy
import math

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

from metag.energetics.protonation import macroscopic_transform, normalize_site
from metag.routing.pka_transform import _neutralize_v2

SCOPE = 'fully_protonated_primary_amino_carboxylic_acid_family'
ACID = Chem.MolFromSmarts('[CX3](=[OX1])[OX2H1]')


def canonical(smiles):
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None or len(Chem.GetMolFrags(molecule)) != 1:
        raise ValueError('A reference must be one valid connected molecule')
    if any(a.GetIsotope() or a.GetAtomMapNum() or a.GetNumRadicalElectrons()
           for a in molecule.GetAtoms()):
        raise ValueError('Mapped, isotope-specific or radical references are unsupported')
    return Chem.MolToSmiles(molecule, isomericSmiles=True), molecule


def neutral_key(smiles):
    canonical(smiles)
    neutral, _, _, _ = _neutralize_v2(smiles)
    if neutral is None:
        raise ValueError('Cannot identify the acid/amine reference skeleton')
    return canonical(neutral)[0]


def hydrogen_count(smiles):
    return sum(a.GetAtomicNum() == 1 for a in Chem.AddHs(Chem.MolFromSmiles(smiles)).GetAtoms())


def include_enantiomer_counterparts(models, conditions):
    """Derive mirror-family constants explicitly in achiral, metal-free water.

    This transfers equilibrium constants, not calculated molecular energies.
    Each reference remains a separate stereoisomer; no racemic mixing entropy
    or extra degeneracy is added. Mirror-state sampling still needs checks.
    """
    output, seen = [], set()
    for model in models:
        checked = validate_family(model, conditions)
        for mirror in (False, True):
            derived = copy.deepcopy(checked)
            molecule = Chem.MolFromSmiles(checked['reference_smiles'])
            if mirror:
                for atom in molecule.GetAtoms():
                    if atom.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED:
                        atom.InvertChirality()
                derived['source']['conversion'] += (
                    ' Mirror-enantiomer constants derived by achiral-water thermodynamic equivalence '
                    '(primary section2.2); not independent D-isomer measurements. No mixing factor.')
            derived['reference_smiles'] = Chem.MolToSmiles(molecule, isomericSmiles=True)
            if derived['reference_smiles'] not in seen:
                seen.add(derived['reference_smiles'])
                output.append(validate_family(derived, conditions))
    return output


def validate_family(model, conditions):
    required = {'kind', 'reference_smiles', 'reference_deprotonations', 'reference_scope',
                'reference_formula', 'pKas', 'conditions', 'source'}
    if not isinstance(model, dict) or set(model) != required:
        raise ValueError('Complete source-family declaration required')
    if (model['kind'] != 'sequential_macroscopic_acid_pka'
            or model['reference_scope'] != SCOPE
            or type(model['reference_deprotonations']) is not int
            or model['reference_deprotonations'] != 0):
        raise ValueError('Whole ladders require the fully protonated reference, not an intermediate microstate')
    smiles, molecule = canonical(model['reference_smiles'])
    if model['reference_formula'] != rdMolDescriptors.CalcMolFormula(molecule):
        raise ValueError('Reference formula mismatch')
    if any(a.GetAtomicNum() not in (1, 6, 7, 8) for a in molecule.GetAtoms()):
        raise ValueError('Only closed-shell C/H/N/O primary-amino carboxylic families are supported')
    nitrogens = [a for a in molecule.GetAtoms() if a.GetAtomicNum() == 7]
    if len(nitrogens) != 1:
        raise ValueError('Exactly one primary ammonium group is supported')
    nitrogen = nitrogens[0]
    neighbors = [a for a in nitrogen.GetNeighbors() if a.GetAtomicNum() != 1]
    if (nitrogen.GetFormalCharge() != 1 or nitrogen.GetTotalNumHs() != 3
            or len(neighbors) != 1 or neighbors[0].GetAtomicNum() != 6
            or neighbors[0].GetHybridization() != Chem.HybridizationType.SP3
            or nitrogen.IsInRing()
            or any(a.GetFormalCharge() != 0 for a in molecule.GetAtoms() if a.GetIdx() != nitrogen.GetIdx())):
        raise ValueError('Reference must have NH3+ on saturated carbon and neutral remaining groups')
    acids = molecule.GetSubstructMatches(ACID)
    if not acids or len({m[0] for m in acids}) != len(acids):
        raise ValueError('Complete protonated carboxylic acid groups required')
    for carbon, oxo, hydroxyl in acids:
        others = [a for a in molecule.GetAtomWithIdx(carbon).GetNeighbors()
                  if a.GetIdx() not in (oxo, hydroxyl)]
        if len(others) != 1 or others[0].GetAtomicNum() != 6:
            raise ValueError('Carbamic/carbonic acids need a separate reference model')
    acid_oxygen = {i for match in acids for i in match[1:]}
    alcohol_counts = Counter()
    for atom in molecule.GetAtoms():
        if atom.GetAtomicNum() != 8 or atom.GetIdx() in acid_oxygen:
            continue
        neighbors = list(atom.GetNeighbors())
        if (atom.GetTotalNumHs() != 1 or len(neighbors) != 1
                or neighbors[0].GetAtomicNum() != 6
                or neighbors[0].GetHybridization() != Chem.HybridizationType.SP3):
            raise ValueError('Other carbonyl, ether, tautomer and hydration branches are unsupported')
        alcohol_counts[neighbors[0].GetIdx()] += 1
    if any(n > 1 for n in alcohol_counts.values()):
        raise ValueError('Gem-diols need a coupled hydration model')
    values = model['pKas']
    if (not isinstance(values, (list, tuple)) or len(values) != len(acids)+1
            or any(type(v) not in (float, int) or not math.isfinite(v) for v in values)):
        raise ValueError('Sequential ladder must cover all carboxyls and the ammonium step')
    expected_conditions = {'T_K', 'ionic_strength_M', 'solvent', 'standard_state', 'magnesium_binding'}
    source_conditions = model['conditions']
    if not isinstance(source_conditions, dict) or set(source_conditions) != expected_conditions:
        raise ValueError('Explicit source conditions required')
    for key in ('T_K', 'ionic_strength_M'):
        value = source_conditions[key]
        if type(value) not in (int, float) or not math.isfinite(value) or value != conditions[key]:
            raise ValueError('Source condition mismatch: '+key)
    if (source_conditions['solvent'] != 'water' or source_conditions['standard_state'] != 'solute_1M_activity'
            or source_conditions['magnesium_binding'] != 'absent'
            or conditions.get('pMg') is not None or conditions.get('water_activity') != 1):
        raise ValueError('Only pure-water activity, 1M solutes and no Mg binding are supported')
    source = model['source']
    required_source = {'locator', 'location', 'artifact_sha256', 'determination', 'parameter_scale',
                       'conversion', 'uncertainty', 'uncertainty_kind', 'covariance', 'development_exposure'}
    if not isinstance(source, dict) or set(source) != required_source:
        raise ValueError('Complete parameter-source provenance required')
    if source['determination'] not in ('measured', 'evaluated') or source['parameter_scale'] != 'molarity_activity':
        raise ValueError('Explicit measured/evaluated target-scale parameters required')
    if any(not isinstance(source[k], str) or not source[k].strip()
           for k in ('locator', 'location', 'conversion', 'uncertainty_kind')):
        raise ValueError('Nonempty source locations, conversions and uncertainty descriptions required')
    digest = source['artifact_sha256']
    if not isinstance(digest, str) or len(digest) != 64 or any(c not in '0123456789abcdef' for c in digest):
        raise ValueError('Source artifact SHA256 required')
    if type(source['development_exposure']) is not bool or source['covariance'] is not None:
        raise ValueError('Exposure declaration required; covariance propagation is not implemented')
    uncertainty = source['uncertainty']
    if uncertainty is not None and (not isinstance(uncertainty, (list, tuple)) or len(uncertainty) != len(values)
            or any(type(v) not in (int, float) or not math.isfinite(v) or v < 0 for v in uncertainty)):
        raise ValueError('Uncertainties must cover every step, or explicitly be unavailable')
    return dict(copy.deepcopy(model), reference_smiles=smiles)


class FamilyRegistry:
    """One-reaction source ownership and reference replacement; no empirical table."""
    def __init__(self, models, conditions):
        for key in ('T_K', 'ionic_strength_M', 'pH'):
            value = conditions[key]
            if type(value) not in (int, float) or not math.isfinite(value):
                raise ValueError('Finite target conditions required')
        if conditions['T_K'] <= 0 or conditions['ionic_strength_M'] < 0:
            raise ValueError('Invalid target temperature or ionic strength')
        self.conditions = dict(conditions)
        self.models = {}
        self.assigned = {}
        self.replaced_terms = {}
        self.reference_shifts = {}
        self.original_species = None
        for model in models:
            checked = validate_family(model, conditions)
            key = neutral_key(checked['reference_smiles'])
            if key in self.models:
                raise ValueError('Duplicate source model for one chemical family')
            self.models[key] = checked

    def model_for(self, smiles):
        """Look up complete supported families without changing unrelated chemistry."""
        molecule = Chem.MolFromSmiles(smiles)
        if molecule is None:
            raise ValueError('Invalid molecular input')
        if (len(Chem.GetMolFrags(molecule)) != 1
                or sum(a.GetAtomicNum() == 7 for a in molecule.GetAtoms()) != 1
                or any(a.GetIsotope() or a.GetAtomMapNum() or a.GetNumRadicalElectrons()
                       for a in molecule.GetAtoms())):
            return None
        return self.models.get(neutral_key(smiles))

    def bind_input(self, species):
        """Own original full-molecule identities before any cofactor/core rewrite."""
        if self.original_species is not None or self.assigned:
            raise ValueError('Use a fresh registry for each reaction')
        self.original_species = copy.deepcopy(species)

    def _parents(self, name, model):
        candidates = [name]+([name[:-2]] if name.endswith('_t') else [])
        return [parent for parent in candidates if parent in self.original_species
                and self.model_for(self.original_species[parent][2]) == model]

    def preserve_full_reference(self, molecule, keep, original_close):
        """Keep full source-owned molecules; this is not an energetic cancellation test."""
        if self.model_for(Chem.MolToSmiles(molecule, isomericSmiles=True)) is not None:
            return set(range(molecule.GetNumAtoms()))
        return original_close(molecule, keep)

    def replace_built_reference(self, result, ownership):
        """Replace neutral-base references and remove their old transforms exactly once."""
        if result is None:
            return None
        if self.original_species is None:
            raise ValueError('Bind original full-molecule identities before reference replacement')
        species, sites, protons = copy.deepcopy(result)
        remaining = Counter(normalize_site(site) for site in sites)
        covered = set()
        for name, (coefficient, charge, smiles) in species.items():
            model = self.model_for(smiles)
            if model is None or not self._parents(name, model):
                continue
            covered.update(self._parents(name, model))
            if name not in ownership:
                raise ValueError('Source reference has no per-species transform ownership')
            self.replaced_terms[name] = copy.deepcopy(ownership[name])
            for site in ownership[name]:
                term = normalize_site(site)
                if term[0] != ('prod' if coefficient > 0 else 'react') or term[3] != abs(coefficient):
                    raise ValueError('Old transform ownership does not match stoichiometry')
                if remaining[term] <= 0:
                    raise ValueError('Duplicate or missing old transform replacement')
                remaining[term] -= 1
            reference = model['reference_smiles']
            delta_h = hydrogen_count(reference)-hydrogen_count(smiles)
            if delta_h != 1-charge:
                raise ValueError('Reference mapping is not a pure proton transfer')
            self.reference_shifts[name] = {'old_smiles': smiles, 'old_charge': charge,
                                           'reference_smiles': reference, 'reference_charge': 1,
                                           'protons_added_per_molecule': delta_h}
            protons -= coefficient*delta_h
            species[name] = [coefficient, 1, reference]
            ownership[name] = []
            self.assigned[name] = model
        required = {name for name, (c, q, s) in self.original_species.items()
                    if c and self.model_for(s) is not None}
        if not required <= covered:
            raise ValueError('A source-owned full species was removed or changed by structural routing')
        # Preserve order and multiplicity for all unrelated transforms.
        kept = []
        for site in sites:
            term = normalize_site(site)
            if remaining[term] > 0:
                kept.append(site)
                remaining[term] -= 1
        return species, kept, protons

    def contribution(self, scored_species, *, rt_ln10):
        """Signed whole-family transforms on the FINAL oriented/scaled reaction."""
        rows = {}
        if self.original_species is None:
            raise ValueError('Original full-molecule ownership is required')
        for name, (coefficient, charge, smiles) in scored_species.items():
            model = self.model_for(smiles)
            if model is None or not self._parents(name, model):
                continue
            if name not in self.assigned or charge != 1 or canonical(smiles)[0] != model['reference_smiles']:
                raise ValueError('A whole ladder cannot be attached to an unassigned/intermediate reference')
            transform = macroscopic_transform(model['pKas'], self.conditions['pH'], rt_ln10,
                                               reference_deprotonations=0)
            rows[name] = {'model': model, 'coefficient': coefficient, 'transform': transform,
                          'reaction_contribution_kj': coefficient*transform['G_transform']}
        if set(self.assigned) != set(rows):
            raise ValueError('A source-owned reference disappeared after assignment')
        return math.fsum(row['reaction_contribution_kj'] for row in rows.values()), rows
