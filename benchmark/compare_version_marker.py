"""Compare an evaluated source tree with a version-marker-only working tree.

This is a source-level diagnostic, not runtime-equivalence proof. It refuses
any package file difference other than one literal __version__ assignment in
__init__.py, and verifies the raw-species source fingerprints separately.
"""
import argparse,ast,hashlib,json
from pathlib import Path
from metag.validation.provenance import source_tree_fingerprint
from metag.energetics.provenance import species_source_fingerprint

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def package_files(root):return {str(p.relative_to(root)):p for p in root.rglob('*') if p.is_file() and (p.suffix=='.py' or (p.parent.name=='data' and p.suffix=='.json'))}
def normalized_init(path):
 tree=ast.parse(path.read_text());found=[]
 for node in tree.body:
  if isinstance(node,ast.Assign) and len(node.targets)==1 and isinstance(node.targets[0],ast.Name) and node.targets[0].id=='__version__':
   if not isinstance(node.value,ast.Constant) or not isinstance(node.value.value,str):raise ValueError('Version assignment must be a literal string')
   found.append(node.value.value);node.value=ast.Constant(value='<VERSION>')
 if len(found)!=1:raise ValueError('Expected exactly one version assignment')
 return ast.dump(tree,include_attributes=False),found[0]
def compare(evaluated,current):
 a=package_files(evaluated);b=package_files(current)
 if set(a)!=set(b):raise ValueError('Package file sets differ')
 changed=sorted(k for k in a if sha(a[k])!=sha(b[k]))
 if changed!=['__init__.py']:raise ValueError('Changes exceed version marker: '+repr(changed))
 aa,av=normalized_init(a['__init__.py']);bb,bv=normalized_init(b['__init__.py'])
 if aa!=bb or av==bv:raise ValueError('Non-version syntax differs or version unchanged')
 raw_a=species_source_fingerprint(evaluated);raw_b=species_source_fingerprint(current)
 if raw_a!=raw_b:raise ValueError('Raw-species fingerprint changed')
 return {'evaluated_fingerprint':source_tree_fingerprint(evaluated),'current_fingerprint':source_tree_fingerprint(current),'raw_species_fingerprint':raw_a,'evaluated_version':av,'current_version':bv,'n_files':len(a),'byte_different_files':changed,'normalized_init_ast_equal':True,'scope':'All package files byte-identical except a literal version-marker assignment; raw-species fingerprint identical. Source-level comparison only; package metadata or runtime introspection may differ.'}
def main():
 parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('evaluated',type=Path);parser.add_argument('current',type=Path);parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
 if args.output.exists():raise SystemExit('Refusing overwrite')
 out=compare(args.evaluated,args.current)
 with args.output.open('x') as f:json.dump(out,f,indent=2);f.write('\n')
 print(json.dumps(out,indent=2))
if __name__=='__main__':main()
