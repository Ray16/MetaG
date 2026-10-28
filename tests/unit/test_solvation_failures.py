"""Partial or failed xtb output is not evidence of a successful solvent calculation (ALPB dgsolv)."""
from types import SimpleNamespace
import pytest
from metag.energetics import thermal


@pytest.mark.parametrize('returncode,solv_energy,expected', [
    (0, '-1.010000', -26.25499639),
    (1, '-1.010000', None),
    (0, None, None),
])
def test_dgsolv_requires_successful_exit_and_complete_result(monkeypatch, returncode, solv_energy, expected):
    def producer(command, *, cwd, **kwargs):
        assert command[command.index('--uhf') + 1] == '2'
        if '--alpb' in command:
            out = '' if solv_energy is None else f'  | TOTAL ENERGY  {solv_energy} Eh  |\n'
            return SimpleNamespace(returncode=returncode, stdout=out, stderr='')
        return SimpleNamespace(returncode=0, stdout='  | TOTAL ENERGY  -1.000000 Eh  |\n', stderr='')
    monkeypatch.setattr(thermal.subprocess, 'run', producer)
    actual = thermal.dgsolv(['O', 'O'], [[0, 0, 0], [0, 0, 1.2]], 0, mult=3)
    assert actual == pytest.approx(expected) if expected is not None else actual is None
