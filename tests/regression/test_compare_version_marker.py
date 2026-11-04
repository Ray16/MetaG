"""Version-only proof must not hide executable or package-data differences."""
import importlib.util
from pathlib import Path
import pytest
SCRIPT=Path(__file__).resolve().parents[2]/'benchmark/compare_version_marker.py'
spec=importlib.util.spec_from_file_location('compare_version_marker',SCRIPT)
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)

def make_tree(path,version):
    path.mkdir(parents=True)
    (path/'__init__.py').write_text(f'__version__ = "{version}"\n')
    for name in ('chem.py','symmetry.py','settings.py','constants.py','timing.py'):(path/name).write_text('value = 1\n')
    (path/'data').mkdir();(path/'data/example.json').write_text('{"value": 1}\n')

def test_accepts_only_literal_version_change(tmp_path):
    a=tmp_path/'a/metag';b=tmp_path/'b/metag';make_tree(a,'0.3.0');make_tree(b,'0.3.1.dev0')
    result=module.compare(a,b)
    assert result['normalized_init_ast_equal'] and result['current_version']=='0.3.1.dev0'
    (b/'__init__.py').write_text('__version__ = "0.3.1.dev0"\nvalue = 2\n')
    with pytest.raises(ValueError,match='Non-version syntax'):
        module.compare(a,b)

def test_rejects_changed_physical_source_or_package_data(tmp_path):
    a=tmp_path/'a/metag';b=tmp_path/'b/metag';make_tree(a,'0.3.0');make_tree(b,'0.3.1.dev0')
    (b/'chem.py').write_text('value = 2\n')
    with pytest.raises(ValueError,match='Changes exceed version marker'):
        module.compare(a,b)
    (b/'chem.py').write_text('value = 1\n');(b/'data/example.json').write_text('{"value": 2}\n')
    with pytest.raises(ValueError,match='Changes exceed version marker'):
        module.compare(a,b)
