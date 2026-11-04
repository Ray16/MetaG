"""Audit ModelSEED input completeness with CPU state-reference gates only.

No Gibbs energies, geometries, model weights, or numerical accuracy are computed.
Missing states and source annotations remain explicit; sugar-pool candidates are
curation leads and are never substituted into the strict as-provided profile.
"""
import argparse
from collections import Counter
import csv
import datetime
from fractions import Fraction
from functools import lru_cache
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import re
import signal
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_chunks(directory,kind):
    files = sorted(p for p in directory.glob(kind+'_*.tsv')
                   if re.fullmatch(kind+r'_\d+\.tsv',p.name))
    if not files:raise ValueError('No ModelSEED '+kind+' data chunks')
    rows = {}
    for path in files:
        with path.open(newline='') as stream:
            for row in csv.DictReader(stream,delimiter='\t'):
                identity = row['id']
                if identity in rows:raise ValueError('Duplicate ModelSEED identity: '+identity)
                rows[identity] = row
    return rows,{str(p.relative_to(directory.parent)):sha(p) for p in files}


def parse_equation(value):
    """Aggregate exact coefficients while preserving compartment identities."""
    result = Counter()
    for term in value.split(';'):
        if not term.strip():continue
        fields = term.split(':',3)
        if len(fields)<3:raise ValueError('Incomplete stoichiometric term')
        coefficient = Fraction(fields[0])
        if not fields[1] or not fields[2]:raise ValueError('Missing compound or compartment')
        result[(fields[1],fields[2])] += coefficient
    return {key:value for key,value in result.items() if value}


def pin_one_cpu():
    allowed = set(os.sched_getaffinity(0));occupied=set();parents={os.getpid()};pid=os.getppid()
    while pid>1:
        parents.add(pid)
        try:pid=int(Path('/proc',str(pid),'stat').read_text().rsplit(')',1)[1].split()[1])
        except OSError:break
    for path in Path('/proc').glob('[0-9]*'):
        try:
            pid=int(path.name)
            if pid in parents:continue
            argv=(path/'cmdline').read_bytes().decode().split('\0')
            if not argv[0]:continue
            name=Path(argv[0]).name;affinity=set(os.sched_getaffinity(pid))
            if ('python' in name or name in {'xtb','crest'}) and len(affinity)<os.cpu_count():
                occupied.update(affinity)
        except (OSError,ValueError,UnicodeError):pass
    free=sorted(allowed-occupied)
    if not free or max(os.getloadavg())+1>os.cpu_count():raise RuntimeError('No verified CPU capacity')
    os.sched_setaffinity(0,{free[-1]})
    return {'host':os.uname().nodename,'pid':os.getpid(),'cpu_affinity':[free[-1]],
            'online_cpus':os.cpu_count(),'load':list(os.getloadavg())}


def failure_category(message):
    text=message.lower()
    for needle,category in [('wildcard/dummy','undefined_atom'),('unparseable','invalid_smiles'),
        ('disconnected','disconnected_species'),('atom-mapped','atom_map_annotation'),
        ('declared charge','charge_structure_mismatch'),('isotope','unsupported_isotope'),
        ('transition-metal','unsupported_metal_spin'),('stereochemical','unresolved_stereochemistry'),
        ('unsupported phosphoryl','unsupported_phosphoryl_reference')]:
        if needle in text:return category
    return 'unsupported_reference_model'


def main():
    started = time.monotonic()
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--modelseed',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--state-budget-seconds',type=float,default=5.)
    args=parser.parse_args()
    if args.output.exists():parser.error('output directory already exists')
    if not 0<args.state_budget_seconds<=60:parser.error('state budget must be in (0,60] seconds')
    resources=pin_one_cpu()
    from rdkit import Chem,RDLogger
    from rdkit.Chem import rdMolDescriptors
    from metag.reactions import validate_reaction
    from metag.state_function import plan_state,POLICY
    from metag.routing.sugar_states import enumerate_sugar_states
    from metag.energetics.conformers import spin_multiplicity
    from metag.chem import formula_charge
    RDLogger.DisableLog('rdApp.*')
    compounds,compound_files=load_chunks(args.modelseed/'Biochemistry','compound')
    reactions,reaction_files=load_chunks(args.modelseed/'Biochemistry','reaction')
    args.output.mkdir()
    def alarm(signum,frame):raise TimeoutError('CPU state-reference preflight time budget exceeded')
    signal.signal(signal.SIGALRM,alarm)

    @lru_cache(maxsize=None)
    def inspect(charge,smiles):
        signal.setitimer(signal.ITIMER_REAL,args.state_budget_seconds)
        try:
            validate_reaction({'species':{'state':[1,charge,smiles]}})
            molecule=Chem.MolFromSmiles(smiles)
            sizes={'n_atoms':Chem.AddHs(molecule).GetNumAtoms(),
                   'n_heavy_atoms':molecule.GetNumHeavyAtoms(),
                   'n_rotatable_bonds':rdMolDescriptors.CalcNumRotatableBonds(molecule)}
            spin_multiplicity(smiles,charge)
            try:plan=plan_state(charge,smiles)
            except ValueError as exc:
                category=failure_category(str(exc));candidate=False;count=None
                if category=='unresolved_stereochemistry':
                    family=enumerate_sugar_states(smiles,charge)
                    candidate=family.applicable and family.complete
                    count=len(family.states) if candidate else None
                return dict(sizes,status='input_state_unavailable',category=category,reason=str(exc),
                            possible_reducing_sugar_parent=candidate,n_candidate_states=count)
            return dict(sizes,status='passes_cpu_state_gates',category=None,
                        canonical_input=plan['input_smiles'],reference_charge=plan['reference_charge'],
                        reference_smiles=plan['reference_smiles'],proton_exchange=plan['proton_exchange'],
                        protonation_issues=plan['protonation_admission']['issues'],
                        n_reference_sites=len(plan['pka_sites']))
        except TimeoutError as exc:
            return {'status':'preflight_resource_limit','category':'preflight_time_budget','reason':str(exc)}
        except (ValueError,RuntimeError) as exc:
            return {'status':'input_state_unavailable','category':failure_category(str(exc)),'reason':str(exc)}
        finally:signal.setitimer(signal.ITIMER_REAL,0)

    try:revision=subprocess.check_output(['git','-C',str(args.modelseed),'rev-parse','HEAD'],text=True).strip()
    except subprocess.CalledProcessError:revision=None
    runtime_files = sorted(p for p in (ROOT/'src/metag').rglob('*') if p.is_file() and p.suffix in {'.py','.json'})
    owned_sources = {str(p.relative_to(ROOT)):sha(p) for p in [Path(__file__),*runtime_files]}
    software = {'python':sys.version}
    for name in ['rdkit','numpy','ase']:
        software[name] = importlib.metadata.version(name)
    protocol={'schema':'metag-modelseed-structural-audit-v2','created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
              'scope':'Strict as-provided input completeness only. No thermodynamic estimates, QM success rate, uncertainty calibration, accuracy or physiological direction is established.',
              'modelseed_revision':revision,'source_files_sha256':dict(compound_files,**reaction_files),
              'state_function_policy':POLICY,'state_budget_seconds':args.state_budget_seconds,
              'reaction_universe':'All nonobsolete reactions explicitly annotated is_transport=0. Abstract entries and annotation/state failures stay in the requested denominator.',
              'sugar_policy':'Unresolved reducing-sugar candidates are reported as curation leads only. No pool, anomer, charge or compound alias is inferred or silently substituted.',
              'source_status_policy':'Database status flags are retained; elemental/charge balance is independently recomputed from supplied concrete structures.',
              'resources':resources,'software':software,'source_sha256':owned_sources,
              'source_scope':'Classifier and all runtime Python/JSON package files; numerical model weights are not used.',
              'source_manifest_sha256':hashlib.sha256(json.dumps(owned_sources,sort_keys=True).encode()).hexdigest()}
    (args.output/'protocol.json').write_text(json.dumps(protocol,indent=2)+'\n')
    state_records={};compound_counts=Counter();sugar_candidates=[]
    nulls={'','null','none','nan'}
    with (args.output/'compounds.jsonl').open('x') as stream:
        for i,(cid,row) in enumerate(sorted(compounds.items())):
            base={'compound':cid,'name':row.get('name'),'source_charge':row.get('charge'),'smiles':row.get('smiles')}
            if row.get('is_obsolete')=='1':record={'status':'obsolete_compound','category':'obsolete_compound'}
            elif row.get('abstract_compound')=='1':record={'status':'abstract_compound','category':'abstract_compound'}
            elif (row.get('smiles') or '').lower() in nulls:record={'status':'missing_structure','category':'missing_smiles'}
            else:
                try:
                    charge=Fraction(row['charge']);assert charge.denominator==1
                    charge=int(charge)
                except (ValueError,TypeError,KeyError,AssertionError,ZeroDivisionError):
                    record={'status':'missing_charge','category':'missing_or_noninteger_charge'}
                else:
                    if cid=='cpd00067':
                        valid=charge==1 and Chem.MolToSmiles(Chem.MolFromSmiles(row['smiles']))=='[H+]'
                        if not valid:raise ValueError('Proton reservoir identity changed')
                        record={'status':'explicit_proton_reservoir','category':None}
                    else:record=inspect(charge,row['smiles'])
                    base['charge']=charge
            item=dict(base,**record);state_records[cid]=item;compound_counts[item['status']]+=1
            if item.get('possible_reducing_sugar_parent') and not re.search(r'alpha|beta|α|β|pyran|furan',row.get('name',''),re.I):
                sugar_candidates.append(cid)
            stream.write(json.dumps(item,allow_nan=False)+'\n')
            if (i+1)%1000==0:print('COMPOUNDS',i+1,'of',len(compounds),flush=True)
    outcomes=Counter();primary=Counter();status_flags=Counter();passed_flags=Counter();excluded=Counter();n_requested=0
    with (args.output/'reactions.jsonl').open('x') as stream:
        for rid,row in sorted(reactions.items()):
            if row.get('is_obsolete')=='1':excluded['obsolete']+=1;continue
            if row.get('is_transport')=='1':excluded['transport']+=1;continue
            if row.get('is_transport')!='0':excluded['unknown_transport_annotation']+=1;continue
            n_requested+=1;status_flags[row.get('status','')]+=1
            item={'reaction':rid,'source_status':row.get('status'),'causes':[]}
            if row.get('abstract_reaction')=='1':item['causes'].append({'category':'abstract_reaction'})
            try:terms=parse_equation(row['stoichiometry'])
            except (ValueError,KeyError,ZeroDivisionError) as exc:
                item['causes'].append({'category':'unsupported_stoichiometry','reason':str(exc)});terms={}
            if len({compartment for cid,compartment in terms})>1:
                item['causes'].append({'category':'multiple_compartments_need_conditions'})
            item['compounds']=sorted({cid for cid,compartment in terms})
            complete=True;balance=Counter();charge_balance=Fraction(0);signature=Counter()
            for (cid,compartment),coefficient in terms.items():
                if cid=='cpd00067':balance['H']+=coefficient;charge_balance+=coefficient;signature['H+']+=coefficient;continue
                if cid not in state_records:
                    item['causes'].append({'category':'missing_compound','compound':cid});complete=False;continue
                compound=state_records[cid]
                if compound['status']!='passes_cpu_state_gates':
                    item['causes'].append({'category':compound['category'],'compound':cid});complete=False;continue
                formula,q=formula_charge(compound['smiles'])
                for element,count in formula.items():balance[element]+=coefficient*count
                charge_balance+=coefficient*q
                signature[(q,compound['canonical_input'])]+=coefficient
            if complete and terms:
                if any(balance.values()):item['causes'].append({'category':'element_imbalance','residual':{k:str(v) for k,v in balance.items() if v}})
                if charge_balance:item['causes'].append({'category':'charge_imbalance','residual':str(charge_balance)})
                if not any(signature.values()):item['causes'].append({'category':'identity_or_graph_collapsed_reaction'})
            if not terms and not item['causes']:item['causes'].append({'category':'empty_net_equation'})
            item['status']='passes_cpu_input_gates' if not item['causes'] else 'requires_curation_or_is_unsupported'
            outcomes[item['status']]+=1
            if item['causes']:primary[item['causes'][0]['category']]+=1
            else:passed_flags[row.get('status','')]+=1
            stream.write(json.dumps(item,allow_nan=False)+'\n')
    summary={'scope':protocol['scope'],'n_compound_rows':len(compounds),'compound_statuses':dict(compound_counts),
             'n_cached_molecular_inspections':inspect.cache_info().currsize,
             'reducing_sugar_parent_curation_leads':sugar_candidates,
             'n_reaction_rows':len(reactions),'n_requested_live_nontransport':n_requested,
             'excluded_from_universe':dict(excluded),'reaction_input_statuses':dict(outcomes),
             'primary_unavailable_categories':dict(primary),'all_source_statuses':dict(status_flags),
             'passed_source_statuses':dict(passed_flags),'physical_calculations_attempted':0,
             'physical_predictions_available':None,'experimental_accuracy_assessed':False,
             'preflight_elapsed_seconds':time.monotonic()-started,
             'protocol_sha256':sha(args.output/'protocol.json'),
             'output_sha256':{n:sha(args.output/n) for n in ['compounds.jsonl','reactions.jsonl']}}
    (args.output/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps({k:summary[k] for k in ['n_compound_rows','n_reaction_rows',
          'n_requested_live_nontransport','excluded_from_universe','reaction_input_statuses',
          'primary_unavailable_categories','preflight_elapsed_seconds','physical_calculations_attempted']},indent=2))


if __name__=='__main__':main()
