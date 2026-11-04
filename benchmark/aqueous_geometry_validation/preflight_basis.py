"""CPU-only four-geometry basis, overlap and installed-source ownership preflight."""
from pathlib import Path
import importlib.metadata as metadata
import hashlib,json,os,sys
import numpy as np
from pyscf import gto,lib
HERE=Path(__file__).resolve().parent
read=lambda p:json.loads(p.read_text())
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
protocol=read(HERE/'fixed_geometry_protocol.json');inputs=read(HERE/'geometry_inputs.json')['inputs']
lib.num_threads(1)
versions={key:metadata.version(key) for key in ['pyscf','gpu4pyscf-cuda12x','gpu4pyscf-libxc-cuda12x','numpy','cupy-cuda12x']}
expected={'pyscf':'2.14.0','gpu4pyscf-cuda12x':'1.8.1','gpu4pyscf-libxc-cuda12x':'0.8.1','numpy':'2.2.6','cupy-cuda12x':'14.1.1'}
assert versions==expected,versions
installed={}
for package,relative in [
 ('pyscf','pyscf/gto/mole.py'),('pyscf','pyscf/gto/basis/def2-tzvpd.dat'),
 ('pyscf','pyscf/gto/basis/def2-universal-jkfit.dat'),('pyscf','pyscf/dft/rks.py'),
 ('gpu4pyscf-cuda12x','gpu4pyscf/scf/hf.py'),('gpu4pyscf-cuda12x','gpu4pyscf/dft/rks.py'),
 ('gpu4pyscf-cuda12x','gpu4pyscf/solvent/smd.py'),('gpu4pyscf-cuda12x','gpu4pyscf/solvent/pcm.py')]:
 path=Path(metadata.distribution(package).locate_file(relative));assert path.is_file(),path
 installed[str(path)]=sha(path)
rows=[]
for slot,item in inputs.items():
 geom=item['geometry']
 mol=gto.M(atom=list(zip(geom['symbols'],geom['positions_angstrom'])),unit='Angstrom',basis=protocol['basis'],charge=geom['charge'],spin=0,verbose=0,max_memory=protocol['max_memory_MB'])
 minimum=float(np.linalg.eigvalsh(mol.intor('int1e_ovlp'))[0])
 assert mol.nao_nr()==item['n_AO'] and mol.nelectron==item['n_electrons']
 assert minimum>=protocol['overlap_minimum_eigenvalue']
 rows.append(dict(slot=slot,n_AO=mol.nao_nr(),n_electrons=mol.nelectron,minimum_overlap_eigenvalue=minimum,passed=True))
report=dict(host=os.uname().nodename,python=sys.executable,cpu_affinity=sorted(os.sched_getaffinity(0)),versions=versions,installed_file_sha256=installed,rows=rows,passed=True,geometry_inputs_sha256=sha(HERE/'geometry_inputs.json'),numeric_settings_sha256=sha(HERE/'fixed_geometry_protocol.json'),source_sha256=sha(Path(__file__)))
(HERE/'dependency_identity.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(dict(rows=rows,versions=versions),indent=2))
