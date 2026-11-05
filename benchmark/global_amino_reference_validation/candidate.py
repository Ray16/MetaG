"""Frozen research candidate: molecule-local whole-family amino-acid references."""
import copy
import json

from rdkit import Chem

from metag.chem import is_balanced
from metag.constants import CONDITIONS, G_HPLUS, PH, RT_LN10, STD_STATE_KJ
from metag.energetics.protonation import macroscopic_transform
from metag.energetics.provenance import content_hash
from metag.routing.aldehyde_hydration import n_hydration_sites
from metag.state_function import StateFunctionEstimator, plan_state
from coupled_amino_reference import neutral_key, validate_family


class GlobalAminoEstimator(StateFunctionEstimator):
    """No reaction routing or hook patch; one immutable model set per session."""
    def __init__(self, backend, models):
        super().__init__(backend)
        records={}
        for model in models:
            checked=validate_family(model,CONDITIONS)
            key=neutral_key(checked['reference_smiles'])
            if key in records:
                raise ValueError('Ambiguous models for one molecular family')
            assert n_hydration_sites(checked['reference_smiles'])==0
            records[key]=checked
        self._amino_models=json.dumps(records,sort_keys=True,allow_nan=False)
        self._amino_identity=content_hash(records)
        self.config.update(state_function_policy='research-global-source-amino-references-v1',
            amino_reference_models_sha256=self._amino_identity,
            amino_reference_scope='Conditional source-scale interpretation; no production promotion')

    def _check_identity(self):
        super()._check_identity()
        if hasattr(self,'_amino_models') and content_hash(json.loads(self._amino_models))!=self._amino_identity:
            raise ValueError('Source model set changed within the molecular ledger')

    def state(self,charge,smiles,*,pool=None):
        if pool is not None:
            return super().state(charge,smiles,pool=pool)
        self._check_identity()
        plan=plan_state(charge,smiles)  # retain all input/identity guards
        mol=Chem.MolFromSmiles(plan['input_smiles'])
        if sum(a.GetAtomicNum()==7 for a in mol.GetAtoms())!=1:
            return super().state(charge,smiles)
        model=json.loads(self._amino_models).get(neutral_key(plan['input_smiles']))
        if model is None:
            return super().state(charge,smiles)
        reference=model['reference_smiles']
        exchange=int(charge)-1
        conversion={'input':[-1,charge,plan['input_smiles']], 'reference':[1,1,reference]}
        if not is_balanced(conversion,exchange):
            raise ValueError('Source reference conversion does not conserve elements/charge')
        plan.update(reference_smiles=reference,reference_charge=1,proton_exchange=exchange,
            pka_sites=[],source_reference_model=copy.deepcopy(model),
            source_reference_policy='global-exact-primary-amino-family')
        key=content_hash({'charge':plan['input_charge'],'smiles':plan['input_smiles']})
        if key not in self._records:
            warnings=[];detail={}
            energy,sigma=self._species_energy(1,reference,'state',warnings,detail)
            transform=macroscopic_transform(model['pKas'],PH,RT_LN10,reference_deprotonations=0)
            reservoir=exchange*G_HPLUS
            sampling_inputs=self._sampling_inputs(1,reference,{})
            record=dict(plan=plan,potential_kj=energy+STD_STATE_KJ+transform['G_transform']+reservoir,
                solution_reference_kj=energy+STD_STATE_KJ,protonation_kj=transform['G_transform'],
                reference_proton_kj=reservoir,sampling_sigma_kj=sigma,sampling_inputs=sampling_inputs,
                molecular_diagnostics=detail,aqueous_diagnostics={'warnings':warnings},
                source_protonation={'policy':'global-source-amino-family','model':model},
                coupled_protonation=transform,state_id=key)
            record['record_sha256']=content_hash(record)
            self._records[key]=json.dumps(record,sort_keys=True,allow_nan=False)
        self._check_identity()
        return json.loads(self._records[key])

    def ledger(self):
        result=super().ledger()
        result.pop('sha256')
        result['source_amino_models']=json.loads(self._amino_models)
        return dict(result,sha256=content_hash(result))
