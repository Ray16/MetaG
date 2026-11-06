"""Closed-state/source-identifiability diagnostics; never enroll or score Gibbs labels."""
from pathlib import Path
from itertools import product
import datetime,hashlib,json,math,os,time
import numpy as np
import scipy
from scipy.optimize import least_squares
from scipy.stats import spearmanr,chi2
from profile_math import profile_minimum
h=Path(__file__).resolve().parent;began=time.monotonic();limit=1800
data=json.loads((h/'calibrated_runs.json').read_text())

def arrays(run,changes=None):
 changes=changes or {};t=np.array([changes.get(('time',o['source_row']),o['time_min']) for o in run['observations']]);ys=[]
 for name in ['cis','trans']:
  cal=run['calibration'][name];conc=changes.get(('stdconc',name),cal['concentration_mM']);area=changes.get(('stdarea',name),cal['area'])
  ys.append([np.nan if o[name+'_area'] is None else changes.get((name,o['source_row']),o[name+'_area'])/area*conc for o in run['observations']])
 return t,np.array(ys).T

def prediction(p,t,feq=None):
 if feq is None:logtotal,f0,feq,lograte=p
 else:logtotal,f0,lograte=p
 total=np.exp(logtotal);relax=np.exp(-np.exp(lograte)*t);cis=total*(feq+(f0-feq)*relax)
 return np.column_stack([cis,total-cis])

def closed(t,y,feq=None,start=None):
 valid=np.isfinite(y);total=np.nanmedian(np.nansum(y,axis=1));scale=total
 f0=np.nan_to_num(y[0,0],nan=0)/total
 last=np.nan_to_num(y[-1,0],nan=0)/total
 fullstarts=[[math.log(total),min(.99,max(0,f0)),min(.99,max(.01,fr)),math.log(mult/max(t.max(),1))] for fr in [last,.5] for mult in [.2,1,5]]
 if start is not None:fullstarts.insert(0,list(start))
 if feq is not None:starts=[[x[0],x[1],x[3]] if len(x)==4 else x for x in fullstarts];lower=[-14,0,-21];upper=[3,1,0]
 else:starts=fullstarts;lower=[-14,0,0,-21];upper=[3,1,1,0]
 best=None
 for guess in starts:
  if time.monotonic()-began>limit:raise RuntimeError('Prespecified CPU wall budget exhausted')
  fit=least_squares(lambda p:((prediction(p,t,feq)-y)/scale)[valid],np.clip(guess,np.array(lower)+1e-10,np.array(upper)-1e-10),bounds=(lower,upper),xtol=1e-11,ftol=1e-11,gtol=1e-11,max_nfev=1000)
  sse=float(np.sum((prediction(fit.x,t,feq)[valid]-y[valid])**2))
  if best is None or sse<best[0]:best=(sse,fit)
 return best

def single(t,y):
 mask=np.isfinite(y);tt=t[mask];yy=y[mask];scale=float(np.median(yy)) or .01;best=None
 def pred(p):return p[1]+(p[0]-p[1])*np.exp(-np.exp(p[2])*tt)
 for factor in [.2,1,5]:
  fit=least_squares(lambda p:(pred(p)-yy)/scale,[max(yy[0],1e-10),max(yy[-1],1e-10),math.log(factor/max(tt.max(),1))],bounds=([0,0,-21],[3,3,0]),xtol=1e-11,ftol=1e-11,gtol=1e-11,max_nfev=1000)
  sse=float(np.sum((pred(fit.x)-yy)**2))
  if best is None or sse<best[0]:best=(sse,fit)
 fit=best[1]
 return dict(initial_mM=float(fit.x[0]),asymptote_mM=float(fit.x[1]),lambda_per_minute=float(np.exp(fit.x[2])),SSE_mM2=best[0],at_numeric_parameter_boundary=bool(np.any(fit.active_mask)),optimizer_success=bool(fit.success))

def ratio_fit(t,y):
 mask=np.all(np.isfinite(y),axis=1);tt=t[mask];ratio=y[mask,0]/np.sum(y[mask],axis=1)
 def pred(p):return p[1]+(p[0]-p[1])*np.exp(-np.exp(p[2])*tt)
 best=None
 for eq in [ratio[-1],.5]:
  for rate_factor in [.2,1,5]:
   fit=least_squares(lambda p:pred(p)-ratio,[np.clip(ratio[0],0,1),np.clip(eq,0,1),math.log(rate_factor/max(tt.max(),1))],bounds=([0,0,-21],[1,1,0]),xtol=1e-12,ftol=1e-12,gtol=1e-12,max_nfev=1000)
   sse=float(np.sum((pred(fit.x)-ratio)**2))
   if best is None or sse<best[0]:best=(sse,fit)
 fit=best[1];eq=fit.x[1]
 return dict(n_paired_timestamps=int(mask.sum()),n_missing_channel_timestamps=int((~mask).sum()),initial_cis_fraction=float(fit.x[0]),asymptotic_cis_fraction=float(eq),K=float(eq/(1-eq)),lambda_per_minute=float(np.exp(fit.x[2])),SSE_fraction2=best[0],residual_time_spearman=float(spearmanr(tt,pred(fit.x)-ratio).statistic),at_numeric_parameter_boundary=bool(np.any(fit.active_mask)),scope='Exploratory common-multiplicative-recovery sensitivity. Ratios do not cancel selective chemical loss, unequal calibration drift or speciation/condition errors. This fit is not adopted as a thermodynamic reference.')

def main():
 global began
 began=time.monotonic()
 assert sorted(os.sched_getaffinity(0))==[79]
 assert all(os.environ.get(k)=='1' for k in ['OMP_NUM_THREADS','OMP_THREAD_LIMIT','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS','NUMEXPR_NUM_THREADS'])
 results=[]
 for run in data['runs']:
  row={k:run[k] for k in ['page','table','T_C','buffer_mM','status']}
  if run['status']!='measured':results.append(row);continue
  t,y=arrays(run);sse,fit=closed(t,y);logtotal,f0,feq,lr=fit.x;pred=prediction(fit.x,t);residual=pred-y;paired=np.all(np.isfinite(y),axis=1);totals=np.sum(y[paired],axis=1)
  row.update(n_timepoints=len(t),n_observed_areas=int(np.isfinite(y).sum()),parameters=dict(total_mM=float(np.exp(logtotal)),initial_cis_fraction=float(f0),equilibrium_cis_fraction=float(feq),K=float(feq/(1-feq)),lambda_per_minute=float(np.exp(lr))),SSE_mM2=sse,optimizer_success=bool(fit.success),at_numeric_parameter_boundary=bool(np.any(fit.active_mask)),observed_total=dict(first_paired_mM=float(totals[0]),last_paired_mM=float(totals[-1]),relative_first_to_last_change=float((totals[-1]-totals[0])/totals[0]),range_relative_to_mean=float(np.ptp(totals)/np.mean(totals)),time_spearman=float(spearmanr(t[paired],totals).statistic)),residual_channel_correlation=float(np.corrcoef(residual[paired].T)[0,1]),residual_RMSE_mM={name:float(np.sqrt(np.nanmean(residual[:,i]**2))) for i,name in enumerate(['cis','trans'])},residual_time_spearman={name:float(spearmanr(t[np.isfinite(y[:,i])],residual[np.isfinite(y[:,i]),i]).statistic) for i,name in enumerate(['cis','trans'])},separate_channel_fits={name:single(t,y[:,i]) for i,name in enumerate(['cis','trans'])},predicted_concentrations_mM=pred.tolist(),residuals_mM=[[None if not np.isfinite(x) else float(x) for x in rr] for rr in residual],raw_calibration_range_flags=dict(cis_below_lowest_tested=sum(o['cis_below_tested_calibration_range'] is True for o in run['observations']),trans_below_lowest_tested=sum(o['trans_below_tested_calibration_range'] is True for o in run['observations'])))
  row['common_recovery_ratio_sensitivity']=ratio_fit(t,y)
  # Initial-condition sensitivity keeps every observed datum and fixes sourceCis0=0.
  zfit=least_squares(lambda p:((prediction([p[0],0,p[1],p[2]],t)-y)/np.exp(logtotal))[np.isfinite(y)],[logtotal,feq,lr],bounds=([-14,0,-21],[3,1,0]),xtol=1e-11,ftol=1e-11,gtol=1e-11,max_nfev=1000)
  row['fixed_initial_zero_sensitivity']=dict(parameters=zfit.x.tolist(),K=float(zfit.x[1]/(1-zfit.x[1])),SSE_mM2=float(np.nansum((prediction([zfit.x[0],0,zfit.x[1],zfit.x[2]],t)-y)**2)))
  # Enumerate every local combination of retained transcription alternatives.
  dimensions=[]
  for o in run['observations']:
   for name in ['time','cis','trans']:
    vals=o['time_candidates_min'] if name=='time' else o[name+'_area_candidates']
    if len(vals)>1:dimensions.append(((name,o['source_row']),vals))
  for name in ['cis','trans']:
   cal=run['calibration'][name]
   for kind,field in [('stdconc','concentration_candidates_mM'),('stdarea','area_candidates')]:
    if len(cal[field])>1:dimensions.append(((kind,name),cal[field]))
  variants=[]
  if dimensions:
   for choice in product(*(v for k,v in dimensions)):
    changes=dict(zip((k for k,v in dimensions),choice));tt,yy=arrays(run,changes);vsse,vfit=closed(tt,yy,start=fit.x);ff=vfit.x[2]
    variants.append(dict(choices=[dict(field=list(k),value=v) for k,v in changes.items()],K=float(ff/(1-ff)),SSE_mM2=vsse,parameters=vfit.x.tolist()))
  row['transcription_sensitivity']=dict(all_local_combinations=variants,K_range=[min([row['parameters']['K']]+[v['K'] for v in variants]),max([row['parameters']['K']]+[v['K'] for v in variants])],scope='Digit-choice sensitivity only; not measurement uncertainty or a confidence interval.')
  # Profile all declared81Kvalues. Working-iid likelihood-ratio values describe
  # shape; residual/condition/calibration errors prevent a reference CI claim.
  profile=[]
  for logK in np.linspace(-3,1,81):
   kk=10**float(logK);eq=kk/(1+kk);psse,profile_parameters=profile_minimum(t,y,kk)
   profile.append(dict(K=kk,SSE_mM2=psse,nuisance_parameters=profile_parameters,working_iid_LR=float(np.isfinite(y).sum()*np.log(max(psse/sse,1)))))
  support=[v['K'] for v in profile if v['working_iid_LR']<=float(chi2.ppf(.95,1))]
  row['K_profile']=profile;row['working_iid_profile_grid_support']=dict(threshold=float(chi2.ppf(.95,1)),lower=min(support) if support else None,upper=max(support) if support else None,lower_truncated=bool(support and support[0]==profile[0]['K']),upper_truncated=bool(support and support[-1]==profile[-1]['K']),scope='Conditional asymptotic working-iidGaussian diagnostic only; not an authenticated experimental95%interval. Boundary/heteroskedastic/correlated residuals and source uncertainty remain.')
  results.append(row)
  (h/'fits.partial.json').write_text(json.dumps(results,indent=2,allow_nan=False)+'\n')
  print(run['T_C'],run['buffer_mM'],'K',round(row['parameters']['K'],4),'mass_change',round(row['observed_total']['relative_first_to_last_change'],3),flush=True)
 report=dict(schema='metag-aconitate-raw-fit-diagnostics-v1',utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),protocol_sha256=hashlib.sha256((h/'protocol.json').read_bytes()).hexdigest(),data_sha256=hashlib.sha256((h/'calibrated_runs.json').read_bytes()).hexdigest(),analysis_source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),runtime=dict(numpy=np.__version__,scipy=scipy.__version__,CPU_affinity=sorted(os.sched_getaffinity(0)),native_threads=1,elapsed_s=time.monotonic()-began,GPUs=0),runs=results,source_reference_uncertainty_intervals=None,temperature_bridge_constructed=False,model_predictions_performed=False,experimental_labels_admitted=0,scope='Exploratory raw-source reanalysis. Conditional fits/profile shapes and transcription sensitivity do not authenticate a thermodynamic reference or validateMetaG.')
 (h/'fit_diagnostics.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')


if __name__=='__main__':main()
