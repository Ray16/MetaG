"""Explicit source-aware parameter construction, independent of reaction outcomes."""
import copy

from coupled_amino_reference import canonical, include_enantiomer_counterparts

DELTA = 0.00127
FAMILIES = {
    '[NH3+]CC(=O)O': ('glycine', 2),
    '[NH3+][C@@H](C)C(=O)O': ('alanine', 2),
    '[NH3+][C@@H](CC(=O)O)C(=O)O': ('aspartate', 3),
    '[NH3+][C@@H](CCC(=O)O)C(=O)O': ('glutamate', 3),
}
FAMILIES = {canonical(smiles)[0]: value for smiles, value in FAMILIES.items()}
FLAGS = [('aspartate', 0), ('aspartate', 1), ('aspartate', 2), ('glutamate', 0), ('glutamate', 1)]


def candidate_models(packet, conditions, flags=(0, 0, 0, 0, 0)):
    """Construct one declared scale scenario; never infer a probability distribution.

    Gly molality is primary-source verified; Ala is supported by curated source
    annotations with original confirmation pending. The five flags are unresolved
    step-scale interpretations; zero denotes the nominal conditional choice.
    """
    if len(flags) != 5 or any(type(value) is not int or value not in (0, 1) for value in flags):
        raise ValueError('Five explicit zero/one source-scale flags required')
    if packet.get('parameter_state') != 'unconverted_primary_transcription':
        raise ValueError('Expected an unconverted source transcription')
    shifts = {'glycine': [DELTA, DELTA], 'alanine': [DELTA, DELTA],
              'aspartate': [0., 0., 0.], 'glutamate': [0., 0., 0.]}
    for (family, step), flag in zip(FLAGS, flags):
        shifts[family][step] = flag*DELTA
    output, used = [], set()
    for raw in packet['models']:
        key = canonical(raw['reference_smiles'])[0]
        if key not in FAMILIES:
            raise ValueError('Unexpected source family in the fixed candidate')
        family, steps = FAMILIES[key]
        if family in used or len(raw['pKas']) != steps:
            raise ValueError('Duplicate or incomplete primary table transcription')
        used.add(family)
        model = copy.deepcopy(raw)
        if model['source']['artifact_sha256'] != packet['source_pdf_sha256']:
            raise ValueError('Inconsistent primary source identity')
        model['pKas'] = [value+shift for value, shift in zip(raw['pKas'], shifts[family])]
        model['conditions']['standard_state'] = 'solute_1M_activity'
        model['source']['parameter_scale'] = 'molarity_activity'
        model['source']['conversion'] = (
            f"Unprimed primary Table6 transcription; explicit per-step shifts {shifts[family]}. "
            'Gly verified molal; Ala molal ancestry curation-supported, original confirmation pending; '
            'Asp and first two Glu steps remain conditional scale choices. Glu step3 follows primary '
            'molar activity conversion. Delta0.00127 from owned NIST2002 Eq10. '
            'No reaction errors select parameters. Conditional corner envelope is not total uncertainty.')
        output.append(model)
    if used != {'glycine', 'alanine', 'aspartate', 'glutamate'}:
        raise ValueError('All four fixed primary families are required')
    return include_enantiomer_counterparts(output, conditions)
