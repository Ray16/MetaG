"""The source comparison must not hide executable or data changes."""
import importlib.util
from pathlib import Path

SCRIPT=Path(__file__).resolve().parents[2]/'benchmark/compare_source_trees.py'
spec=importlib.util.spec_from_file_location('compare_source_trees',SCRIPT)
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)

def test_docstrings_and_comments_are_distinguished_from_executable_code():
    a='"""old"""\ndef f(x):\n """old"""\n return x+1\n'
    b='"""new"""\n# note\ndef f(x):\n """new"""\n return x+1\n'
    assert module.semantic_tree(a)==module.semantic_tree(b)
    assert module.semantic_tree(a)!=module.semantic_tree(b.replace('x+1','x+2'))
    assert module.semantic_tree('def f():\n x=1\n "later"\n')!=module.semantic_tree('def f():\n x=1\n')

def test_changed_data_and_missing_files_are_reported(tmp_path):
    left=tmp_path/'left/metag';right=tmp_path/'right/metag';left.mkdir(parents=True);right.mkdir(parents=True)
    (left/'example.py').write_text('x=1\n');(right/'example.py').write_text('x=1\n')
    (left/'data.json').write_text('{"x":1}');(right/'data.json').write_text('{"x":2}')
    (left/'only.py').write_text('x=1\n')
    result=module.compare(left,right)
    assert result['missing_current']==['only.py']
    assert result['semantic_tree_differences']==['data.json']
