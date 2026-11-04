"""Compare two MetaG package trees, preserving their distinct source identities.

Example::

    python benchmark/compare_source_trees.py /path/to/evaluated/metag src/metag --output /tmp/source-comparison.json

AST equality after removing docstrings is a source-level diagnostic. It does
not prove runtime equivalence, and code that inspects changed docstrings or
source bytes can behave differently. No numerical benchmark data are used.
"""
import argparse
import ast
import hashlib
import json
import sys
from pathlib import Path


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def semantic_tree(text):
    class StripDocs(ast.NodeTransformer):
        def visit(self, node):
            node = super().visit(node)
            if isinstance(node, (ast.Module, ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)) and node.body:
                first = node.body[0]
                if isinstance(first, ast.Expr) and isinstance(first.value, ast.Constant) and isinstance(first.value.value, str):
                    node.body.pop(0)
            return node
    return ast.dump(StripDocs().visit(ast.parse(text)), include_attributes=False)


def fingerprint(root):
    """Match MetaG's full Python-source fingerprint convention."""
    result = hashlib.sha256()
    for path in sorted(root.rglob('*.py')):
        result.update(str(path.relative_to(root.parent)).encode())
        result.update(path.read_bytes())
    return result.hexdigest()


def compare(evaluated, current):
    evaluated, current = Path(evaluated).resolve(), Path(current).resolve()
    if evaluated.name != 'metag' or current.name != 'metag' or not evaluated.is_dir() or not current.is_dir():
        raise ValueError('Both arguments must be existing metag package directories')
    def files(root):
        return {str(p.relative_to(root)):p for p in root.rglob('*')
                if p.is_file() and '__pycache__' not in p.parts and p.suffix in ('.py','.json')}
    a,b=files(evaluated),files(current)
    missing_current=sorted(set(a)-set(b));missing_evaluated=sorted(set(b)-set(a))
    changed=[];semantic_differences=[];records={}
    for name in sorted(set(a)&set(b)):
        left,right=digest(a[name]),digest(b[name])
        records[name]={'evaluated_sha256':left,'current_sha256':right}
        if left==right:continue
        changed.append(name)
        if not name.endswith('.py') or semantic_tree(a[name].read_text())!=semantic_tree(b[name].read_text()):
            semantic_differences.append(name)
    return {'evaluated_fingerprint':fingerprint(evaluated),'current_fingerprint':fingerprint(current),
            'n_evaluated_files':len(a),'n_current_files':len(b),'missing_current':missing_current,
            'missing_evaluated':missing_evaluated,'byte_changed_files':changed,
            'semantic_tree_differences':semantic_differences,'files':records,'python_version':sys.version,
            'scope':'Python AST comparison after removing docstrings and byte comparison of packaged JSON. Not a runtime equivalence proof; source-byte or docstring introspection may differ.'}


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('evaluated',type=Path);parser.add_argument('current',type=Path)
    parser.add_argument('--output',type=Path,required=True);args=parser.parse_args(argv)
    result=compare(args.evaluated,args.current)
    with args.output.open('x') as stream:json.dump(result,stream,indent=2,allow_nan=False);stream.write('\n')
    print(json.dumps({k:result[k] for k in ('evaluated_fingerprint','current_fingerprint','n_evaluated_files',
      'n_current_files','missing_current','missing_evaluated','byte_changed_files','semantic_tree_differences')},indent=2))
    return int(bool(result['missing_current'] or result['missing_evaluated'] or result['semantic_tree_differences']))

if __name__=='__main__':raise SystemExit(main())
