from pathlib import Path
import importlib.util,json
import numpy as np
from scipy.linalg import expm
h=Path(__file__).resolve().parent;spec=importlib.util.spec_from_file_location('fit_math',h/'fit_timecourses.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
K=.15;f=K/(1+K);t=np.linspace(0,2000,41);true=[np.log(.3),.01,f,np.log(.003)];closed=m.prediction(true,t);masked=closed.copy();masked[0,0]=np.nan
sse,fit=m.closed(t,masked);assert abs(fit.x[2]/(1-fit.x[2])-K)<1e-7 and sse<1e-14
common=closed*np.exp(-.4*t/t.max())[:,None];ratio=m.ratio_fit(t,common);assert abs(ratio['K']-K)<1e-7
kf=.003*f;kr=.003*(1-f);dc=.003;A=np.array([[-kr-dc,kf],[kr,-kf]]);selective=np.array([expm(A*x)@closed[0] for x in t]);biased=m.ratio_fit(t,selective);assert abs(biased['K']-K)>.02
r=dict(scope='Synthetic source-fit/interpretation controls only; no physical predictions or experimental reference validation.',true_K=K,closed_recovery_with_missing_first_cis=dict(K=float(fit.x[2]/(1-fit.x[2])),SSE_mM2=sse),common_multiplicative_loss=dict(K=ratio['K'],expected_K=K),selective_loss_counterexample=dict(fitted_ratio_K=biased['K'],true_underlying_kf_over_kr=K,not_the_same=True),conclusion='A ratio fit can survive common recovery changes yet fail to identify the underlying reactionK under selective loss. It must not be adopted as a thermodynamic label from fit quality alone.',primary_labels_admitted=0)
(h/'synthetic_controls.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r,indent=2))
