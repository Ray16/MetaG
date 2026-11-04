"""Select the authors' minimum-training-loss epoch from archived dGbyG CV runs.

Input is a local extraction of `cross_validation_results/10_fold_cross_validation`
from Zenodo 13885478. For each of 100 replicates, the matching NPY loss array
and CSV (or CSV.GZ) provide the selected epoch and target predictions. No
models are retrained and no validation-target errors are used to pick an epoch.
"""
import argparse
import csv
import gzip
import hashlib
import json
from pathlib import Path

import numpy as np

HERE=Path(__file__).resolve().parent

def sha(path):
    digest=hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda:stream.read(1024*1024),b''):
            digest.update(block)
    return digest.hexdigest()


def uncompressed_sha(path):
    digest=hashlib.sha256()
    with (gzip.open(path,'rb') if path.suffix=='.gz' else path.open('rb')) as stream:
        for block in iter(lambda:stream.read(1024*1024),b''):
            digest.update(block)
    return digest.hexdigest()


def select_one(loss_path, table_path, output, *, expected_rows=675):
    losses=np.load(loss_path,allow_pickle=False)
    if losses.ndim != 2 or not losses.size or not np.isfinite(losses).all():
        raise ValueError('Training loss array must be nonempty, finite and two-dimensional')
    epoch=int(losses.mean(axis=0).argmin())
    opener=gzip.open if table_path.suffix=='.gz' else open
    with opener(table_path,'rt',newline='') as stream:
        reader=csv.reader(stream);header=next(reader)
        if len(header)!=len(set(header)) or 'standard_dg_prime' not in header or str(epoch) not in header:
            raise ValueError('Prediction table has missing or duplicate columns')
        target_col=header.index('standard_dg_prime'); prediction_col=header.index(str(epoch))
        rows=[]
        for index,row in enumerate(reader):
            if len(row)!=len(header) or int(row[0])!=index:
                raise ValueError('Prediction table row/index mismatch')
            target=float(row[target_col]); prediction=float(row[prediction_col])
            if not np.isfinite(target) or not np.isfinite(prediction):
                raise ValueError('Nonfinite target or prediction')
            rows.append((index,target,prediction))
    if len(rows)!=expected_rows:
        raise ValueError('Unexpected prediction table length')
    with output.open('x',newline='') as stream:
        writer=csv.writer(stream,delimiter='\t',lineterminator='\n')
        writer.writerow(['','standard_dg_prime','prediction']);writer.writerows(rows)
    return {'epoch':epoch,'loss_shape':list(losses.shape),'loss_sha256':sha(loss_path),
            'prediction_table_sha256':sha(table_path),'selected_sha256':sha(output)}


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--raw-dir',type=Path,required=True,help='directory containing 0.npy/0.csv through 99.npy/99.csv; .csv.gz also accepted')
    parser.add_argument('--output-dir',type=Path,required=True,help='new output directory')
    parser.add_argument('--manifest',type=Path,default=HERE/'dgbyg_cv_manifest.json',help='pinned archive-member fingerprints')
    args=parser.parse_args(argv)
    manifest=json.loads(args.manifest.read_text())
    if manifest.get('n_replicates')!=100 or manifest.get('n_targets')!=675 or set(manifest.get('replicates',{}))!={str(i) for i in range(100)}:
        raise ValueError('Invalid archived CV fingerprint manifest')
    args.output_dir.mkdir(parents=True,exist_ok=False)
    result={}
    for index in range(100):
        loss=args.raw_dir/f'{index}.npy'; plain=args.raw_dir/f'{index}.csv'; compressed=args.raw_dir/f'{index}.csv.gz'
        if plain.exists()==compressed.exists():
            raise ValueError(f'Exactly one CSV or CSV.GZ required for replicate {index}')
        table=plain if plain.exists() else compressed; expected=manifest['replicates'][str(index)]
        if sha(loss)!=expected['loss_sha256'] or uncompressed_sha(table)!=expected['uncompressed_csv_sha256']:
            raise ValueError(f'Archive member fingerprint mismatch for replicate {index}')
        result[str(index)]=select_one(loss,table,args.output_dir/f'{index}.tsv')
        if result[str(index)]['epoch']!=expected['selected_epoch']:
            raise ValueError(f'Selected epoch mismatch for replicate {index}')
    report={'criterion':'minimum mean training loss across folds, as in the authors’ Figures.ipynb',
            'n_replicates':len(result),'manifest_sha256':sha(args.manifest),'replicates':result,
            'scope':'Reconstruction from caller-supplied author archive members; source archive identity and fold assignments must be checked separately.'}
    with (args.output_dir/'selection_audit.json').open('x') as stream:
        json.dump(report,stream,indent=2,allow_nan=False);stream.write('\n')
    print(json.dumps({'n_replicates':len(result),'output_dir':str(args.output_dir)},indent=2))


if __name__=='__main__':
    main()
