"""Verify the final input guard against the actual evaluated package and all inputs.

Example from MetaG_submission/MetaG after extracting code and current-study ZIPs:
python benchmark/final_review_20261031/verify_source_change.py \
 --evaluated-package ../analysis/metag_corrected_validation_20261030_01/source/src/metag \
 --final-package src/metag \
 --inputs ../analysis/metag_corrected_validation_20261030_01/inputs.json
"""
import argparse,ast,json,sys
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--evaluated-package',required=True,type=Path);p.add_argument('--final-package',required=True,type=Path);p.add_argument('--inputs',required=True,type=Path);a=p.parse_args();oldroot=a.evaluated_package.resolve();newroot=a.final_package.resolve();sys.path.insert(0,str(newroot.parent))
from metag.reactions import validate_reaction
from metag.validation.provenance import source_tree_fingerprint
from metag.energetics.provenance import species_source_fingerprint
proof=json.loads(Path(__file__).with_name('independent_arithmetic.json').read_text())['source_change_proof']
assert source_tree_fingerprint(oldroot)==proof['evaluated_source']
assert source_tree_fingerprint(newroot)==proof['final_source']
assert species_source_fingerprint(oldroot)==species_source_fingerprint(newroot)==proof['raw_species_final']
changed=[p.relative_to(newroot).as_posix() for p in sorted(newroot.rglob('*.py')) if p.read_bytes()!=(oldroot/p.relative_to(newroot)).read_bytes()];assert changed==['pipeline.py','reactions.py']
x=ast.parse((oldroot/'pipeline.py').read_text());y=ast.parse((newroot/'pipeline.py').read_text());x.body.pop(0);y.body.pop(0);assert ast.dump(x)==ast.dump(y)
x=ast.parse((oldroot/'reactions.py').read_text());y=ast.parse((newroot/'reactions.py').read_text());fn=next(v for v in y.body if isinstance(v,ast.FunctionDef) and v.name=='validate_reaction');loop=next(v for v in fn.body if isinstance(v,ast.For) and isinstance(v.target,ast.Tuple) and [n.id for n in v.target.elts]==['name','value']);guard=[v for v in loop.body if isinstance(v,ast.If) and 'molAtomMapNumber' in ast.unparse(v.test)];assert len(guard)==1;loop.body.remove(guard[0]);assert ast.dump(x)==ast.dump(y)
inputs=json.loads(a.inputs.read_text());assert len(inputs)==proof['n_panel_inputs_accepted_by_guard']
for r in inputs.values():validate_reaction(r)
print('Verified exact source delta, unchanged raw molecular code and all',len(inputs),'evaluated inputs.')
