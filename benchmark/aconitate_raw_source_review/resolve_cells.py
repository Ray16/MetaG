"""Image-reviewed resolutions of the dual-OCR disagreement cells."""
from pathlib import Path
from decimal import Decimal
import hashlib,json
h=Path(__file__).resolve().parent;r=json.loads((h/'cell_comparison.json').read_text());disputed=json.loads((h/'disputed_cells.json').read_text())
panels=[
['370','2875','4305','10110','11570','37223616','40853984','3380','0.2856','10','14696','43025888'],
['42930688','41216896','3453584','38868544','0.2093','0.2856','10',None,'39719520','3288744','0.2083','0.2856'],
['28907536','40302240','5','42725888','42740224','42352096','42285664','41517792','41474688','41215200','40512800','40990368'],
['40348128','2445990','37760992','5740','0.2856','51896','40800480','598291','39798464','4585','36305760','36437312'],
['0.2093','0.2856','40715648','5',None,'44014816','6205309','5426509','0.2093','75000','513481','33357312'],
['5345','0.2093','0.2856','5345','0.2093','0.2856',None,'34948128','34187296','0.2093','0.2856','5533299'],
['4680','35066432','65','38534496','35727104','35552224','35598816','5352045','1745','0.2856','39662016',None],
['370549','41033088','38443936','815','5212656','5365021','5545498','5601485','0.2856','39456288',None,'1020736'],
['39646400','5151315','0.2093','0.2856','31653872','0','1800','35584096','43449600']]
choices=sum(panels,[]);assert len(choices)==len(disputed)==105
# Ambiguity is retained rather than resolved using reaction trends or fit quality.
alternatives={6:['40853984','40859984'],22:['0.2083','0.2093'],54:['6205309','5205309'],61:['0.2093','0.2083'],72:['4680','4660','4860'],74:['65','85'],90:['5545498','5645498'],91:['5601485','5801485'],96:['39646400','39846400'],100:['31653872','31683872']}
key=lambda c:(c['table'],c['buffer_mM'],c['row'],c['column'])
lookup={key(c):(i,choice) for i,(c,choice) in enumerate(zip(disputed,choices,strict=True))};cells=[];ledger=[]
for cell in r['cells']:
 if cell['agrees']:
  c=dict(cell,nominal_value=cell['hidden']['value'],value_candidates=[cell['hidden']['value']],verification='dual_OCR_agreement')
 else:
  i,choice=lookup[key(cell)];values=[float(v) for v in alternatives.get(i,[choice]) if v is not None]
  c=dict(cell,nominal_value=None if choice is None else float(choice),value_candidates=values,verification='source_blank_not_zero' if choice is None else 'image_reviewed_with_digit_alternatives' if i in alternatives else 'image_reviewed',verification_id=i)
  ledger.append(dict(id=i,cell=dict(table=cell['table'],buffer_mM=cell['buffer_mM'],row=cell['row'],role=cell['role'],column=cell['column']),original_hidden_text=cell['hidden']['text'] if cell['hidden'] else None,independent_image_OCR=cell['image_OCR_text'],source_clip=disputed[i]['source_clip'],panel='conflicts_'+str(i//12).zfill(2)+'.png',panel_sha256=hashlib.sha256((h/('conflicts_'+str(i//12).zfill(2)+'.png')).read_bytes()).hexdigest(),nominal_value=choice,value_candidates=alternatives.get(i,[choice]) if choice is not None else [],status=c['verification']))
 cells.append(c)
runs=[]
for layout_run in json.loads((h/'raw_cell_tokens.json').read_text())['runs']:
 base={k:layout_run[k] for k in ['page','table','T_C','buffer_mM','status']}
 if base['status']!='measured':runs.append(base);continue
 cc=[c for c in cells if key(c)[:2]==(base['table'],base['buffer_mM'])];indices=sorted({c['row'] for c in cc if c['role']=='observation'})
 cs={c['column']:c for c in cc if c['role']=='standard_concentration_mM'};ar={c['column']:c for c in cc if c['role']=='standard_area'}
 base['calibration']={name:dict(concentration_mM=cs[name+'_area']['nominal_value'],area=ar[name+'_area']['nominal_value'],concentration_candidates_mM=cs[name+'_area']['value_candidates'],area_candidates=ar[name+'_area']['value_candidates']) for name in ['cis','trans']};base['observations']=[]
 for i in indices:
  fields={c['column']:c for c in cc if c['row']==i};assert len(fields)==3
  time=fields['time_min']['nominal_value'];obs=dict(source_row=i,time_min=time,time_candidates_min=fields['time_min']['value_candidates'])
  for name in ['cis','trans']:
   c=fields[name+'_area'];cal=base['calibration'][name];v=c['nominal_value'];obs[name+'_area']=v;obs[name+'_area_candidates']=c['value_candidates'];obs[name+'_mM']=None if v is None else v/cal['area']*cal['concentration_mM']
   obs[name+'_below_tested_calibration_range']=None if v is None else obs[name+'_mM']<(.025 if name=='cis' else .028)
  base['observations'].append(obs)
 assert all(b['time_min']>a['time_min'] for a,b in zip(base['observations'],base['observations'][1:]))
 runs.append(base)
assert len(runs)==20 and sum(r['status']=='measured' for r in runs)==19
assert sum(len(r.get('observations',[])) for r in runs)==314
assert sum(o['cis_area'] is None for r in runs for o in r.get('observations',[]))==5
assert all(o['trans_area'] is not None for r in runs for o in r.get('observations',[]))
report=dict(schema='metag-aconitate-calibrated-source-runs-v1',source_PDF_sha256=json.loads((h/'protocol.json').read_text())['source_PDF_sha256'],scope='314timestamps in19measured runs;5blank cis readings retained as missing, notzero or authenticated detection bounds. Nominal digit choices are provisional where alternatives remain. These are not enrolled thermodynamic reference labels.',correction_ledger=ledger,unassigned_OCR_symbols_disposition='18decimal glyphs lie below the row matching tolerance; their standards were visually resolved. TwoOCR0glyphs overlap the Conc(mM)label, not an observation. Original symbol coordinates retained in cell_comparison.json.',runs=runs,counts=dict(requested_conditions=20,measured_runs=19,source_no_experiment=1,timestamps=314,measured_cis_areas=309,measured_trans_areas=314,source_blank_cis=5,numeric_cell_locations=1018,dual_OCR_agreements=913,image_reviewed_disagreements=105,digit_ambiguity_cells=len(alternatives)),calibration='Csample_mM = A_sample/A_standard * C_standard_mM; run-specific standards. Sourceprintedp45/PDFp60.',new_primary_references_admitted=0,model_predictions_performed=False)
(h/'calibrated_runs.json').write_text(json.dumps(report,indent=2)+'\n');(h/'correction_ledger.json').write_text(json.dumps(ledger,indent=2)+'\n');print(json.dumps(report['counts'],indent=2))
