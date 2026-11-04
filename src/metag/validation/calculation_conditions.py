"""Check explicit calculation requests against the implemented thermodynamic target.

This is a calculation-admission check, separate from reference-condition schema
validation. It performs no condition transformation and cannot normalize data.
"""
import math
from collections.abc import Mapping

from metag.constants import CONDITIONS


def calculation_condition_issues(reaction):
    """Return reasons an optional ``conditions`` request cannot be calculated.

    Omission keeps the documented fixed target. A supplied partial record
    constrains that target; it never changes global constants or backend settings.
    Unknown requirements are refused rather than silently ignored.
    """
    if 'conditions' not in reaction:
        return []
    requested = reaction['conditions']
    if not isinstance(requested, Mapping):
        return ['conditions must be a mapping of calculation requirements']
    numeric = {
        'T_K': CONDITIONS['T_K'],
        'temperature_K': CONDITIONS['T_K'],
        'pH': CONDITIONS['pH'],
        'ionic_strength_M': CONDITIONS['ionic_strength_M'],
        'water_activity': CONDITIONS['water_activity'],
        'solute_standard_state_M': 1.0,
    }
    exact = {
        'pMg': None,
        'magnesium_binding': 'absent',
        'standard_state': CONDITIONS['standard_state'],
        'water_standard_state': 'pure_liquid',
        'energy_units': 'kJ/mol',
        'quantity': 'standard_transformed_reaction_gibbs',
    }
    issues = []
    for key, value in requested.items():
        if key in numeric:
            valid = (isinstance(value, (int, float)) and not isinstance(value, bool)
                     and math.isfinite(value))
            if not valid:
                issues.append(f'{key} must be a finite numeric calculation requirement')
            elif not math.isclose(value, numeric[key], rel_tol=0, abs_tol=1e-9):
                issues.append(f'{key}: requested {value!r}, supported {numeric[key]!r}')
        elif key in exact:
            if value != exact[key] or isinstance(value, bool):
                issues.append(f'{key}: requested {value!r}, supported {exact[key]!r}')
        else:
            issues.append(f'unsupported calculation-condition requirement {key!r}')
    return issues
