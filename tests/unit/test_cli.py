"""metag.cli: argument defaults, the tool listing and early refusal of unknown reaction keys."""
import json

import pytest

from metag import cli
from metag.tools import COMMANDS


def test_parse_args_defaults(monkeypatch):
    monkeypatch.delenv("RXN_FILE", raising=False)
    monkeypatch.delenv("METAG_OUT", raising=False)
    a = cli.parse_args([])
    assert a.reactions is None and a.only is None and a.sampling_profile is None
    assert a.out == "./metag_out"
    b = cli.parse_args(["--only", "r1", "--only", "r2", "--sampling-profile", "thorough", "--out", "o"])
    assert b.only == ["r1", "r2"] and b.sampling_profile == "thorough" and b.out == "o"


def test_parse_args_rejects_unknown_profile():
    with pytest.raises(SystemExit):
        cli.parse_args(["--sampling-profile", "fast"])


def test_tools_lists_every_command(capsys):
    assert set(COMMANDS) == {"preflight", "evaluate", "build-reference"}
    cli.main(["tools"])
    out = capsys.readouterr().out
    for name in COMMANDS:
        assert "metag " + name in out


def test_unknown_only_key_exits_before_loading_a_model(tmp_path):
    path = tmp_path / "rx.json"
    path.write_text(json.dumps({"r1": {"species": {"A": [-1, 0, "CCO"], "B": [1, 0, "CCO"]},
                                       "n_Hplus": 0, "note": ""}}))
    with pytest.raises(SystemExit) as exc:
        cli.main(["score", "--reactions", str(path), "--only", "nope", "--out", str(tmp_path / "o")])
    assert "nope" in str(exc.value)
    assert not (tmp_path / "o").exists()          # refused before the model load / output directory
