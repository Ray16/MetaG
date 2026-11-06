from pathlib import Path
import json,math
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import FormatStrFormatter,NullFormatter
h=Path(__file__).resolve().parent;data=json.loads((h/'calibrated_runs.json').read_text());fits=json.loads((h/'fit_diagnostics.json').read_text());rows=[r for r in fits['runs'] if r['status']=='measured'];raw=[r for r in data['runs'] if r['status']=='measured']
plt.rcParams.update({'font.size':9,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42})
colors={25:'#0072B2',50:'#D55E00',75:'#009E73',150:'#CC79A7'}
with PdfPages(h/'native_timecourse_diagnostics.pdf',metadata={'Title':'Walford2000 aconitate raw-source diagnostics; no reference enrollment'}) as pdf:
 for source,row in zip(raw,rows,strict=True):
  t=np.array([o['time_min'] for o in source['observations']]);y=np.array([[np.nan if o[k+'_mM'] is None else o[k+'_mM'] for k in ['cis','trans']] for o in source['observations']]);tt=np.linspace(0,t.max(),300);p=row['parameters'];z=np.exp(-p['lambda_per_minute']*tt);ci=p['total_mM']*(p['equilibrium_cis_fraction']+(p['initial_cis_fraction']-p['equilibrium_cis_fraction'])*z);pred=np.column_stack([ci,p['total_mM']-ci]);fig,axes=plt.subplots(2,2,figsize=(10.5,7.2),layout='constrained')
  for j,name in enumerate(['cis','trans']):
   ax=axes[0,j];ax.plot(t/60,y[:,j],'o',ms=3,color=colors[source['buffer_mM']],label='Calibrated peak areas');ax.plot(tt/60,pred[:,j],color='black',lw=1,label='Closed two-state fit');ax.set(xlabel='Elapsed time (h)',ylabel=f'{name} aconitate (mM)');ax.legend(fontsize=7)
  paired=np.all(np.isfinite(y),axis=1);ax=axes[1,0];ax.plot(t[paired]/60,y[paired].sum(axis=1),'o-',ms=3,lw=.7,color=colors[source['buffer_mM']]);ax.axhline(p['total_mM'],color='black',ls='--',lw=1,label='Closed-fit total');ax.set(xlabel='Elapsed time (h)',ylabel='Measured cis + trans (mM)');ax.legend(fontsize=7)
  ax=axes[1,1];profile=row['K_profile'];ax.plot([v['K'] for v in profile],[v['working_iid_LR'] for v in profile],color='#0072B2',lw=1);ax.axhline(3.841458820694124,color='grey',ls=':',lw=1,label='Working-IID LR threshold');ax.axvline(p['K'],color='black',lw=1,label=f"Closed K={p['K']:.4f}");ratio=row['common_recovery_ratio_sensitivity']['K'];ax.axvline(ratio,color='#CC79A7',ls='--',lw=1,label=f'Ratio sensitivity K={ratio:.4f}');ax.set(xscale='log',xlabel='K = cis / trans',ylabel='Working-IID likelihood ratio',ylim=(0,20),xlim=(.03,.6));ax.legend(fontsize=7);ax.set_xticks([.05,.1,.2,.4]);ax.xaxis.set_major_formatter(FormatStrFormatter('%.2g'));ax.xaxis.set_minor_formatter(NullFormatter())
  fig.suptitle(f"{source['T_C']} °C; {source['buffer_mM']} mM phosphate buffer; nominal pH 7 at room preparation\n{source['table']} — {len(t)} timestamps; {sum(np.isnan(y[:,0]))} blank cis observations retained",fontsize=11)
  fig.supxlabel('Exploratory source reanalysis. Digit alternatives, correlated/unequal errors, loss mechanisms and pH/activity transfer remain unresolved.\nProfile threshold is a working-noise diagnostic, not an authenticated experimental confidence interval. No 298.15 K reference admitted.',fontsize=8)
  pdf.savefig(fig);plt.close(fig)
fig,axes=plt.subplots(1,2,figsize=(10,4),layout='constrained')
for buffer,color in colors.items():
 group=[r for r in rows if r['buffer_mM']==buffer];T=[r['T_C'] for r in group]
 axes[0].plot(T,[r['parameters']['K'] for r in group],'o-',color=color,label=f'{buffer} mM',lw=1,ms=4)
 axes[0].plot(T,[r['common_recovery_ratio_sensitivity']['K'] for r in group],'^--',color=color,lw=1,ms=4)
 axes[1].plot(T,[100*r['observed_total']['relative_first_to_last_change'] for r in group],'o-',color=color,lw=1,ms=4)
axes[0].set(xlabel='Native assay temperature (°C)',ylabel='Fitted K (cis / trans)');axes[0].legend(title='Phosphate buffer',fontsize=7);axes[0].text(.03,.98,'Solid circles: closed concentrations\nDashed triangles: ratio sensitivity',ha='left',va='top',transform=axes[0].transAxes,fontsize=8)
axes[1].set(xlabel='Native assay temperature (°C)',ylabel='First-to-last paired total change (%)');axes[1].axhline(0,color='grey',lw=.7)
fig.suptitle('Aconitate source diagnostics: 19 measured runs, one chemical equation\nNominal pH 7 at room-temperature buffer preparation',fontsize=11);fig.supxlabel('All runs retained; 25 mM / 110 °C was not measured. No model-error comparison or independent reference admission.',fontsize=8)
fig.savefig(h/'source_diagnostics_overview.png',dpi=180);fig.savefig(h/'source_diagnostics_overview.pdf');plt.close(fig)
print('Wrote all19native diagnostic pages and overview.')
