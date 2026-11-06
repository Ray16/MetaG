"""Optionally refit the frozen source curves on one CPU; no quantum inference."""
import argparse
import gzip
import importlib.util
import json
import math
import os
from pathlib import Path
import sys

sys.dont_write_bytecode=True
for key in ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[key]='1'
if hasattr(os,'sched_getaffinity'):
    allowed=os.sched_getaffinity(0)
    os.sched_setaffinity(0,{min(allowed)})


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--packet',type=Path,default=Path(__file__).with_name('aconitate_raw_source_review'))
    parser.add_argument('--profiles',action='store_true',help='Also reoptimize all1539fixedKprofile points.')
    args=parser.parse_args();packet=args.packet.resolve()
    def no_network(event,_args):
        if event in {'socket.connect','socket.getaddrinfo'}:
            raise RuntimeError('No network in source refitting')
    sys.addaudithook(no_network)
    from audit_aconitate_raw import audit
    audit(packet)
    sys.path.insert(0,str(packet))
    spec=importlib.util.spec_from_file_location('frozen_aconitate_fitting',packet/'fit_timecourses.py')
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    raw=json.loads((packet/'calibrated_runs.json').read_text())
    frozen=json.loads(gzip.decompress((packet/'fit_diagnostics.json.gz').read_bytes()))
    n,profiles=0,0;maximum_K_difference=0.;maximum_profile_SSE_difference=0.
    for source,result in zip(raw['runs'],frozen['runs'],strict=True):
        if source['status']!='measured':continue
        t,y=module.arrays(source);sse,fit=module.closed(t,y)
        K=float(fit.x[2]/(1-fit.x[2]));difference=abs(K-result['parameters']['K'])
        assert abs(sse-result['SSE_mM2'])<1e-9 and difference<1e-5
        maximum_K_difference=max(maximum_K_difference,difference);n+=1
        if args.profiles:
            for point in result['K_profile']:
                actual,_=module.profile_minimum(t,y,point['K'])
                difference=abs(actual-point['SSE_mM2']);assert difference<1e-9
                maximum_profile_SSE_difference=max(maximum_profile_SSE_difference,difference);profiles+=1
    print(json.dumps(dict(status='passed',refitted_native_runs=n,reoptimized_profile_points=profiles,
        maximum_K_difference=maximum_K_difference,maximum_profile_SSE_difference_mM2=maximum_profile_SSE_difference,
        CPU_affinity=sorted(os.sched_getaffinity(0)) if hasattr(os,'sched_getaffinity') else None,
        native_threads=1,quantum_inference=False,primary_reference_admission=False,
        scope='Reexecution of the frozen source-fit algorithm; independent experimental validity is not established.'),indent=2))


if __name__=='__main__':main()
