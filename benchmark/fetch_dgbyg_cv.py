"""Fetch only pinned dGbyG 10-fold archive members using HTTP byte ranges.

Downloads 100 loss arrays and 100 prediction tables (about 6.3 GB uncompressed)
from the authors' 35.5 GB ZIP. Each member is checked against its pinned ZIP
length, CRC32 and SHA-256 before an atomic rename. Existing valid files can be
reused on rerun. The source archive and its licensing are external to MetaG.
"""
import argparse
from datetime import datetime,timezone
import gzip
import hashlib
import json
import math
import os
from pathlib import Path
import struct
import tempfile
import time
import urllib.request
import zlib

HERE=Path(__file__).resolve().parent
URL='https://zenodo.org/records/13885478/files/results_data.zip?download=1'
_MIN_INTERVAL=1.1
_LAST_REQUEST=0.0


def read_range(start,end,*,timeout=90):
    global _LAST_REQUEST
    delay=_MIN_INTERVAL-(time.monotonic()-_LAST_REQUEST)
    if delay>0:time.sleep(delay)
    _LAST_REQUEST=time.monotonic()
    request=urllib.request.Request(URL,headers={'Range':f'bytes={start}-{end}',
        'Accept-Encoding':'identity','User-Agent':'MetaG-dGbyG-reproduction/1.0'})
    response=urllib.request.urlopen(request,timeout=timeout)
    content_range=response.headers.get('Content-Range','')
    if response.status!=206 or not content_range.startswith(f'bytes {start}-{end}/'):
        response.close();raise ValueError('Server did not supply the exact requested range')
    return response


def member_offset(info):
    start=info['header_offset']
    if isinstance(start,bool) or not isinstance(start,int) or start<0:
        raise ValueError('Invalid ZIP member offset')
    with read_range(start,start+29) as response:
        header=response.read()
    if len(header)!=30:raise ValueError('Truncated ZIP local header')
    parts=struct.unpack('<4s5H3I2H',header)
    if parts[0]!=b'PK\x03\x04' or parts[3]!=0 or parts[2]&1:
        raise ValueError('Expected an unencrypted stored ZIP member')
    name_length,extra_length=parts[-2:]
    if name_length<=0 or name_length>4096 or extra_length>65535:
        raise ValueError('Invalid ZIP local header lengths')
    with read_range(start+30,start+29+name_length) as response:
        name=response.read().decode('utf-8')
    if name!=info['name']:
        raise ValueError('ZIP local member name differs from pinned identity')
    return start+30+name_length+extra_length


def inspect_existing(path,*,compressed):
    digest=hashlib.sha256();crc=0;size=0
    with (gzip.open(path,'rb') if compressed else path.open('rb')) as stream:
        for block in iter(lambda:stream.read(1024*1024),b''):
            digest.update(block);crc=zlib.crc32(block,crc);size+=len(block)
    return size,crc&0xffffffff,digest.hexdigest()


def matches(values,info,expected_sha):
    return values==(info['size'],info['crc32'],expected_sha)


def fetch_member(info,expected_sha,dest,*,compressed):
    if dest.exists():
        if not matches(inspect_existing(dest,compressed=compressed),info,expected_sha):
            raise ValueError(f'Existing output does not match pinned member: {dest}')
        return 'reused'
    offset=member_offset(info);size=info['size']
    if isinstance(size,bool) or not isinstance(size,int) or size<=0:raise ValueError('Invalid member size')
    fd,tmpname=tempfile.mkstemp(prefix=dest.name+'.',suffix='.part',dir=dest.parent)
    os.close(fd);tmp=Path(tmpname);digest=hashlib.sha256();crc=0;seen=0
    try:
        with read_range(offset,offset+size-1) as response, (gzip.open(tmp,'wb',compresslevel=1) if compressed else tmp.open('wb')) as stream:
            for block in iter(lambda:response.read(1024*1024),b''):
                seen+=len(block)
                if seen>size:raise ValueError('Member response exceeds pinned length')
                digest.update(block);crc=zlib.crc32(block,crc);stream.write(block)
        if not matches((seen,crc&0xffffffff,digest.hexdigest()),info,expected_sha):
            raise ValueError('Member length, CRC32 or SHA-256 mismatch')
        # A concurrent writer must not be silently replaced.
        os.link(tmp,dest);tmp.unlink()
    finally:
        tmp.unlink(missing_ok=True)
    return 'downloaded'


def main(argv=None):
    global _MIN_INTERVAL
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir',type=Path,required=True)
    parser.add_argument('--manifest',type=Path,default=HERE/'dgbyg_cv_manifest.json')
    parser.add_argument('--replicates',type=int,nargs='+',default=list(range(100)))
    parser.add_argument('--kinds',choices=('npy','csv'),nargs='+',default=['npy','csv'])
    parser.add_argument('--interval',type=float,default=1.1,help='minimum seconds between HTTP requests (1 to 60)')
    args=parser.parse_args(argv)
    if (not math.isfinite(args.interval) or not 1<=args.interval<=60 or not args.replicates or len(set(args.replicates))!=len(args.replicates)
            or any(i<0 or i>=100 for i in args.replicates) or len(set(args.kinds))!=len(args.kinds)):
        raise ValueError('Invalid replicate selection, member kind or request interval')
    manifest=json.loads(args.manifest.read_text())
    if manifest.get('n_replicates')!=100 or set(manifest.get('replicates',{}))!={str(i) for i in range(100)}:
        raise ValueError('Invalid archive manifest')
    _MIN_INTERVAL=args.interval
    args.output_dir.mkdir(parents=True,exist_ok=True);out=[]
    for i in args.replicates:
        record=manifest['replicates'][str(i)]
        for kind in args.kinds:
            info=record[kind+'_archive'];expected=record['loss_sha256' if kind=='npy' else 'uncompressed_csv_sha256']
            expected_name=f'dGbyG_results/cross_validation_results/10_fold_cross_validation/{i}.{kind}'
            if info['name']!=expected_name:raise ValueError('Unexpected member identity')
            dest=args.output_dir/(f'{i}.npy' if kind=='npy' else f'{i}.csv.gz')
            state=fetch_member(info,expected,dest,compressed=(kind=='csv'))
            out.append({'replicate':i,'kind':kind,'state':state,'file':dest.name,'source_sha256':expected})
            print(f'{i} {kind}: {state}',flush=True)
    audit={'created_utc':datetime.now(timezone.utc).isoformat(),'source':URL,
           'manifest_sha256':hashlib.sha256(args.manifest.read_bytes()).hexdigest(),'members':out,
           'scope':'Pinned archive member retrieval and integrity only; no verification of model folds, measurement accuracy or source-data licensing.'}
    name='fetch_audit_'+datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')+'.json'
    with (args.output_dir/name).open('x') as f:json.dump(audit,f,indent=2);f.write('\n')


if __name__=='__main__':main()
