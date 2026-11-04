"""CPU numerical fixtures for the later GPU reference implementation."""
from pathlib import Path
import json,math,os,hashlib,time
HERE=Path(__file__).resolve().parent
os.sched_setaffinity(0,{1})
from pyscf import gto,dft,lib
from pyscf.solvent import smd
lib.num_threads(1)
protocol=json.loads((HERE/'fixed_geometry_protocol.json').read_text());start=time.monotonic()
mol=gto.M(atom=[('O',(0,0,0)),('H',(.7586,0,.5043)),('H',(-.7586,0,.5043))],unit='Angstrom',basis=protocol['basis'],charge=0,spin=0,verbose=4,output=str(HERE/'cpu_water_smoke.log'),max_memory=4000)
def calculator():
 mf=dft.RKS(mol,xc=protocol['functional']).density_fit(auxbasis=protocol['auxbasis']);mf.grids.level=protocol['grid_level'];mf.nlcgrids.level=protocol['nlc_grid_level'];mf.conv_tol=protocol['SCF_tolerance_Eh'];mf.conv_tol_grad=protocol['SCF_gradient_tolerance'];mf.max_cycle=protocol['max_cycle'];mf.chkfile=None;assert mf.do_nlc();return mf
gas=calculator();eg=float(gas.kernel());assert gas.converged and math.isfinite(eg)
water=calculator();water=water.SMD();water.with_solvent.solvent='water';water.with_solvent.sasa_ng=590;water.with_solvent.equilibrium_solvation=True;ea=float(water.kernel(dm0=gas.make_rdm1()));assert water.converged and math.isfinite(ea)
result=dict(geometry=[['O',[0,0,0]],['H',[.7586,0,.5043]],['H',[-.7586,0,.5043]]],gas_Eh=eg,water_total_Eh=ea,SMD_delta_Eh=ea-eg,CDS_Eh=float(water.with_solvent.e_cds),source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),elapsed_seconds=time.monotonic()-start,scope='CPU density-fitted numerical implementation fixture, not experimental accuracy.')
(HERE/'cpu_water_smoke.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
