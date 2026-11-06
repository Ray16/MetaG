"""Independent fixed-K profile check: two nonnegative linear coefficients and one rate."""
from pathlib import Path
import hashlib,importlib.util,json,math,time
import numpy as np
from scipy.optimize import minimize_scalar,brentq
from scipy.stats import chi2
h=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('raw_source_fit',h/'fit_timecourses.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
raw=json.loads((h/'calibrated_runs.json').read_text());fit=json.loads((h/'fit_diagnostics.json').read_text());started=time.monotonic();threshold=float(chi2.ppf(.95,1));out=[]

from profile_math import profile_minimum

for source,row in zip(raw['runs'],fit['runs'],strict=True):
 if source['status']!='measured':continue
 t,y=m.arrays(source);n=int(np.isfinite(y).sum());minimum=row['SSE_mM2'];Khat=row['parameters']['K'];differences=[]
 for p in row['K_profile']:
  sse,params=profile_minimum(t,y,p['K']);differences.append(sse-p['SSE_mM2'])
 def lr(logK):
  sse,params=profile_minimum(t,y,10**float(logK));return n*math.log(max(sse/minimum,1))
 center=math.log10(Khat);ends=[]
 for endpoint in [-3,1]:
  if lr(endpoint)<threshold:ends.append(None);continue
  root=brentq(lambda x:lr(x)-threshold,min(center,endpoint),max(center,endpoint),xtol=1e-8);kk=10**root;sse,params=profile_minimum(t,y,kk);ends.append(dict(K=kk,SSE_mM2=sse,nuisance_parameters=params,working_iid_LR=n*math.log(sse/minimum)))
 check_sse,_=profile_minimum(t,y,Khat)
 out.append(dict(table=row['table'],T_C=row['T_C'],buffer_mM=row['buffer_mM'],maximum_original_profile_SSE_difference_mM2=max(abs(x) for x in differences),largest_original_profile_objective_improvement_mM2=-min(differences),baseline_minimum_SSE_difference_mM2=check_sse-minimum,working_iid_support_endpoints=ends,coarse_grid_was_empty=row['working_iid_profile_grid_support']['lower'] is None))
 print(row['T_C'],row['buffer_mM'],'support',[None if x is None else round(x['K'],4) for x in ends],flush=True)
r=dict(schema='metag-aconitate-profile-refinement-v2',fit_sha256=hashlib.sha256((h/'fit_diagnostics.json').read_bytes()).hexdigest(),source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),elapsed_s=time.monotonic()-started,method='Adaptive root refinement using the generator scalar-rate/linearNNLS profile helper; this is not an independent reoptimization after the profile-optimizer correction.',independent_optimizer=False,scope='Working-iidGaussian likelihood-ratio shape only. These are not authenticated experimental95%confidence intervals or standard-Gibbs reference errors. No source-condition or model-adequacy hold is relaxed.',LR_threshold=threshold,runs=out,experimental_labels_admitted=0)
(h/'profile_refinements.json').write_text(json.dumps(r,indent=2)+'\n')
