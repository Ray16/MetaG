"""Input-file failures are actionable and occur before molecular initialization."""
import pytest

from metag.cli import main


@pytest.mark.parametrize('contents', [b'{broken', b'\xff\xfe\xfd'])
def test_score_unreadable_json_exits_before_output_creation(tmp_path, contents):
    source = tmp_path/'reactions.json'
    source.write_bytes(contents)
    output = tmp_path/'results'
    with pytest.raises(SystemExit, match='cannot read reactions'):
        main(['score', '--reactions', str(source), '--out', str(output)])
    assert not output.exists()


def test_score_directory_is_an_input_error(tmp_path):
    with pytest.raises(SystemExit, match='cannot read reactions'):
        main(['score', '--reactions', str(tmp_path)])


@pytest.mark.parametrize('kind', ['missing', 'directory', 'invalid_json', 'invalid_encoding'])
def test_preflight_input_error_is_usage_diagnostic(tmp_path, capsys, kind):
    source = tmp_path/'input.json'
    if kind == 'directory':
        source.mkdir()
    elif kind == 'invalid_json':
        source.write_text('{broken')
    elif kind == 'invalid_encoding':
        source.write_bytes(b'\xff\xfe\xfd')
    output = tmp_path/'report.json'
    with pytest.raises(SystemExit) as error:
        main(['preflight', str(source), '-o', str(output)])
    assert error.value.code == 2
    assert not output.exists()
    diagnostic = capsys.readouterr().err
    assert 'cannot read reactions from' in diagnostic and str(source) in diagnostic
    assert 'Traceback' not in diagnostic
