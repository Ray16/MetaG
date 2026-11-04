"""Build auditable dGbyG pairs from the authors' selected 10-fold CV outputs.

Inputs are the 100 selected replicate TSVs reconstructed from the authors'
cross-validation archive, the matching 675-row historical training CSV, and
the retained reaction identity/condition bridge audit. No model is trained.
Rows without an audited extent and bridge are listed as unavailable.
"""
import argparse
import csv
import hashlib
import json
import math
import struct
from pathlib import Path

import numpy as np

from score_dgbyg_oof import load_identity

HERE = Path(__file__).resolve().parent
TRAINING_SHA256 = 'cd8ede797d2f1c088cccf6b6561d258a34e753ff1ea8c43ca8cb36ca28553f3c'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_training(path):
    if sha(path) != TRAINING_SHA256:
        raise ValueError('Historical CV training table hash mismatch')
    with path.open(newline='') as stream:
        rows = list(csv.DictReader(stream))
    if len(rows) != 675:
        raise ValueError('The archived CV training table must contain 675 rows')
    values = np.asarray([float(row['standard_dg_prime']) for row in rows],dtype=float)
    if not np.isfinite(values).all():
        raise ValueError('Nonfinite training target')
    return rows, values


def selected_ensemble(directory, targets, *, manifest=None):
    arrays = []; hashes = {}
    for replicate in range(100):
        path = directory/f'{replicate}.tsv'
        with path.open(newline='') as stream:
            reader = csv.DictReader(stream,delimiter='\t')
            rows = list(reader)
        if len(rows) != len(targets):
            raise ValueError(f'Incorrect target count in replicate {replicate}')
        if [int(row['']) for row in rows] != list(range(len(targets))):
            raise ValueError(f'Target index mismatch in replicate {replicate}')
        observed = np.asarray([float(row['standard_dg_prime']) for row in rows],dtype=float)
        prediction = np.asarray([float(row['prediction']) for row in rows],dtype=float)
        if not np.isfinite(prediction).all() or not np.allclose(observed,targets,rtol=1e-8,atol=1e-8):
            raise ValueError(f'Nonfinite predictions or target mismatch in replicate {replicate}')
        if manifest is not None:
            digest=hashlib.sha256()
            for value in prediction:
                digest.update(struct.pack('>d',float(value)))
            if digest.hexdigest() != manifest['replicates'][str(replicate)]['prediction_vector_sha256']:
                raise ValueError(f'Archived prediction fingerprint mismatch in replicate {replicate}')
        arrays.append(prediction); hashes[path.name] = sha(path)
    return np.nansum(np.asarray(arrays),axis=0)/len(arrays), hashes


def build_rows(identity, training, predictions):
    output=[]; unavailable=[]
    for reaction, record in sorted(identity.items()):
        row = record['author_paper_row']
        if row >= len(training) or training[row]['type'] != 'reaction dg':
            raise ValueError(f'Author row is not a reaction target: {reaction}')
        if 'signed_author_scale' not in record or 'bridge_I0_noMg_minus_native_kj' not in record:
            unavailable.append(reaction); continue
        value = float(predictions[row]); scale = float(record['signed_author_scale'])
        bridge = float(record['bridge_I0_noMg_minus_native_kj'])
        if not all(math.isfinite(x) for x in (value,scale,bridge)) or scale == 0:
            raise ValueError(f'Invalid comparison fields: {reaction}')
        output.append({'reaction':reaction,'paper_row':row,'scale':scale,
                       'bridge_delta_I0_minus_I025':bridge,'dgbyg':value})
    return output, unavailable


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--selected-dir',type=Path,required=True,help='directory containing 0.tsv through 99.tsv')
    parser.add_argument('--training-csv',type=Path,required=True,help='historical 675-row CV training data')
    parser.add_argument('--identity',type=Path,default=HERE/'dgbyg_family_identity_complete.json')
    parser.add_argument('--manifest',type=Path,default=HERE/'dgbyg_cv_manifest.json',
                        help='pinned source and prediction fingerprints')
    parser.add_argument('--output',type=Path,required=True,help='new pairing TSV')
    parser.add_argument('--audit',type=Path,required=True,help='new provenance JSON')
    args=parser.parse_args(argv)
    if args.output.exists() or args.audit.exists() or args.output.resolve()==args.audit.resolve():
        raise ValueError('Output and audit must be distinct new files')
    identity=load_identity(args.identity); training,targets=read_training(args.training_csv)
    manifest=json.loads(args.manifest.read_text())
    if (manifest.get('n_replicates')!=100 or manifest.get('n_targets')!=len(training)
            or manifest.get('training_sha256')!=TRAINING_SHA256 or set(manifest.get('replicates',{}))!={str(i) for i in range(100)}):
        raise ValueError('Invalid archived CV fingerprint manifest')
    prediction,replicate_hashes=selected_ensemble(args.selected_dir,targets,manifest=manifest)
    rows,unavailable=build_rows(identity,training,prediction)
    with args.output.open('x',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=['reaction','paper_row','scale','bridge_delta_I0_minus_I025','dgbyg'],delimiter='\t')
        writer.writeheader(); writer.writerows(rows)
    audit={'n_author_targets':len(training),'n_replicates':len(replicate_hashes),
           'n_identity_rows':len(identity),'n_pair_rows':len(rows),'unavailable_without_bridge_or_scale':unavailable,
           'training_sha256':sha(args.training_csv),'identity_sha256':sha(args.identity),'manifest_sha256':sha(args.manifest),
           'replicate_sha256':replicate_hashes,'pairs_sha256':sha(args.output),
           'scope':'Archived author-selected 10-fold predictions only. The selected replicate files must be generated from the original archive using the authors’ criterion; this script does not independently verify fold assignments or select epochs.'}
    with args.audit.open('x') as stream:
        json.dump(audit,stream,indent=2,allow_nan=False);stream.write('\n')
    print(json.dumps({'n_pair_rows':len(rows),'n_missing_bridge_or_scale':len(unavailable)},indent=2))


if __name__=='__main__':
    main()
