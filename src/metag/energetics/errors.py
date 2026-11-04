"""Failures shared by species evaluation and reaction orchestration."""

import math
from numbers import Real


def valid_energy(value):
    """A finite scalar energy; booleans and missing values are not energies."""
    return isinstance(value, Real) and not isinstance(value, bool) and math.isfinite(value)


def valid_sampling_sigma(value):
    """A finite nonnegative sampling diagnostic, not a predictive interval."""
    return valid_energy(value) and value >= 0


class SymmetryEnumerationError(ValueError):
    """The bounded atom-mapping search cannot safely identify duplicate minima."""


class SpeciesRearranged(Exception):
    """Every converged, evaluated conformer failed connectivity or specified-stereo checks.

    This does not establish what happened to unconverged conformers or prove the
    absence of a physical minimum. The requested state is not scored or cached;
    score_reaction may retry a permitted neutral-reference route or withhold it.
    """
    def __init__(self, name, smi, n):
        super().__init__(f"{name} ({smi}): connectivity or specified stereochemistry changed in all "
                         f"{n} converged, evaluated conformers")
        self.smi = smi


