"""Small fixtures for reconstructing the archived CV comparison."""
import csv
import gzip
import hashlib
import io
import struct
import sys
import zipfile
from pathlib import Path

import numpy as np
import pytest
from unittest.mock import patch

BENCHMARK = Path(__file__).resolve().parents[2] / 'benchmark'
sys.path.insert(0,str(BENCHMARK))
from build_dgbyg_pairs import build_rows, read_training, selected_ensemble
from fetch_dgbyg_cv import fetch_member
from select_dgbyg_cv import select_one


def test_selects_minimum_mean_training_loss_from_gzip(tmp_path):
    loss=tmp_path/'loss.npy';np.save(loss,np.array([[3.,1.,0.],[1.,3.,5.]]))
    table=tmp_path/'run.csv.gz'
    with gzip.open(table,'wt',newline='') as f:
        w=csv.writer(f);w.writerow(['','standard_dg_prime','0','1','2'])
        w.writerow([0,10,100,101,102]);w.writerow([1,20,200,201,202])
    out=tmp_path/'selected.tsv';info=select_one(loss,table,out,expected_rows=2)
    assert info['epoch']==0  # epochs 0 and 1 tie; NumPy selects the first
    rows=list(csv.DictReader(out.open(),delimiter='\t'))
    assert [float(r['prediction']) for r in rows]==[100,200]


def test_selection_rejects_duplicate_columns_and_missing_targets(tmp_path):
    loss=tmp_path/'loss.npy';np.save(loss,np.array([[1.,0.]]))
    table=tmp_path/'run.csv';table.write_text(',standard_dg_prime,1,1\n0,10,5,6\n')
    with pytest.raises(ValueError,match='duplicate'):
        select_one(loss,table,tmp_path/'selected.tsv',expected_rows=1)


def test_training_table_must_be_the_pinned_historical_snapshot(tmp_path):
    training=tmp_path/'training.csv';training.write_text('standard_dg_prime\n1\n')
    with pytest.raises(ValueError,match='hash mismatch'):
        read_training(training)


def test_range_fetch_checks_member_identity_and_preserves_existing_output(tmp_path):
    archive=tmp_path/'source.zip'; data=b'author archive member fixture'
    name='dGbyG_results/cross_validation_results/10_fold_cross_validation/0.npy'
    with zipfile.ZipFile(archive,'w',compression=zipfile.ZIP_STORED) as z:
        z.writestr(name,data)
    with zipfile.ZipFile(archive) as z:
        item=z.getinfo(name)
    info={'name':name,'header_offset':item.header_offset,'size':item.file_size,'crc32':item.CRC}
    raw=archive.read_bytes();dest=tmp_path/'0.npy';digest=hashlib.sha256(data).hexdigest()
    with patch('fetch_dgbyg_cv.read_range',side_effect=lambda start,end:io.BytesIO(raw[start:end+1])):
        assert fetch_member(info,digest,dest,compressed=False)=='downloaded'
        assert dest.read_bytes()==data
        assert fetch_member(info,digest,dest,compressed=False)=='reused'
        dest.write_bytes(b'wrong')
        with pytest.raises(ValueError,match='Existing output'):
            fetch_member(info,digest,dest,compressed=False)
        assert dest.read_bytes()==b'wrong'
        wrong=dict(info,name=name+'-wrong')
        with pytest.raises(ValueError,match='member name'):
            fetch_member(wrong,digest,tmp_path/'other.npy',compressed=False)
        invalid=tmp_path/'invalid.npy'
        with pytest.raises(ValueError,match='SHA-256 mismatch'):
            fetch_member(info,'0'*64,invalid,compressed=False)
        assert not invalid.exists() and not list(tmp_path.glob('invalid.npy.*.part'))


def test_ensemble_and_pair_builder_check_identity_and_missing_bridge(tmp_path):
    for i in range(100):
        (tmp_path/f'{i}.tsv').write_text('\tstandard_dg_prime\tprediction\n0\t1\t2\n1\t3\t4\n')
    pred, hashes=selected_ensemble(tmp_path,np.array([1.,3.]))
    assert len(hashes)==100 and pred.tolist()==[2.,4.]
    digest=hashlib.sha256(struct.pack('>dd',2.,4.)).hexdigest()
    manifest={'replicates':{str(i):{'prediction_vector_sha256':digest} for i in range(100)}}
    assert selected_ensemble(tmp_path,np.array([1.,3.]),manifest=manifest)[0].tolist()==[2.,4.]
    (tmp_path/'0.tsv').write_text('\tstandard_dg_prime\tprediction\n0\t1\t2\n1\t3\t5\n')
    with pytest.raises(ValueError,match='fingerprint mismatch'):
        selected_ensemble(tmp_path,np.array([1.,3.]),manifest=manifest)
    identity={'r1':{'author_paper_row':0,'signed_author_scale':-1.,'bridge_I0_noMg_minus_native_kj':2.},
              'r2':{'author_paper_row':1}}
    rows,missing=build_rows(identity,[{'type':'reaction dg'},{'type':'reaction dg'}],pred)
    assert rows==[{'reaction':'r1','paper_row':0,'scale':-1.,'bridge_delta_I0_minus_I025':2.,'dgbyg':2.}]
    assert missing==['r2']
    with pytest.raises(ValueError,match='not a reaction target'):
        build_rows(identity,[{'type':'formation dg'},{'type':'reaction dg'}],pred)
