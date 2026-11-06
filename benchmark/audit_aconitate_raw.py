"""Replay the frozen aconitate source arithmetic without inference or optimization.

This checks transcription/calibration bookkeeping and stored fit calculations.
It cannot authenticate the primary experiment or turn conditional fits into
standard Gibbs-energy references.
"""
from __future__ import annotations

import argparse
from collections import Counter
import gzip
import hashlib
import json
import math
from pathlib import Path
import sys


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read(path):
    path = Path(path)
    return json.loads(gzip.decompress(path.read_bytes()) if path.suffix == '.gz' else path.read_bytes())


def close(a, b, tolerance=1e-10):
    assert math.isfinite(a) and math.isfinite(b) and abs(a-b) <= tolerance * max(1., abs(a), abs(b)), (a, b)


def calibrated(run, changes=None):
    changes = changes or {}
    t, values = [], []
    for obs in run['observations']:
        index = obs['source_row']
        t.append(changes.get(('time', index), obs['time_min']))
        row = []
        for name in ('cis', 'trans'):
            cal = run['calibration'][name]
            area = obs[name+'_area']
            if area is None:
                row.append(None)
            else:
                row.append(changes.get((name, index), area)
                           / changes.get(('stdarea', name), cal['area'])
                           * changes.get(('stdconc', name), cal['concentration_mM']))
        values.append(row)
    return t, values


def predicted(parameters, t, fixed_fraction=None):
    if fixed_fraction is None:
        logtotal, initial, fraction, lograte = parameters
    else:
        logtotal, initial, lograte = parameters
        fraction = fixed_fraction
    total = math.exp(logtotal)
    cis = total * (fraction + (initial-fraction) * math.exp(-math.exp(lograte)*t))
    return cis, total-cis


def squared_error(times, values, parameters, fixed_fraction=None):
    return math.fsum((p-v)**2 for t, row in zip(times, values, strict=True)
                     for p, v in zip(predicted(parameters, t, fixed_fraction), row, strict=True)
                     if v is not None)


def profile_kkt(times, values, parameters, K):
    """Independent linear-nuisance stationarity check at a stored profile point."""
    logtotal, initial, lograte = parameters
    total = math.exp(logtotal)
    coefficients = (total*(1-initial), total*initial)
    gradients = [0., 0.]
    fraction = K/(1+K)
    for t, row in zip(times, values, strict=True):
        z = math.exp(-math.exp(lograte)*t)
        w = fraction*(1-z)
        a, b = (w, 1-w), (w+z, 1-w-z)
        for i, observed in enumerate(row):
            if observed is None:
                continue
            residual = a[i]*coefficients[0] + b[i]*coefficients[1] - observed
            gradients[0] += a[i]*residual
            gradients[1] += b[i]*residual
    for coefficient, gradient in zip(coefficients, gradients, strict=True):
        assert gradient >= -1e-8 if coefficient < 1e-8 else abs(gradient) < 1e-8


def audit(packet):
    packet = Path(packet)
    manifest = read(packet/'manifest.json')
    assert {p.relative_to(packet).as_posix() for p in packet.rglob('*') if p.is_file()} == set(manifest['files']) | {'manifest.json'}
    for name, entry in manifest['files'].items():
        path = packet/name
        assert not Path(name).is_absolute() and '..' not in Path(name).parts and not path.is_symlink()
        assert path.stat().st_size == entry['bytes'] and sha(path) == entry['sha256'], name
    data, fits = read(packet/'calibrated_runs.json'), read(packet/'fit_diagnostics.json.gz')
    assert fits['data_sha256'] == sha(packet/'calibrated_runs.json')
    assert fits['protocol_sha256'] == sha(packet/'protocol.json')
    assert fits['analysis_source_sha256'] == sha(packet/'fit_timecourses.py')
    assert data['source_PDF_sha256'] == read(packet/'protocol.json')['source_PDF_sha256']
    assert data['new_primary_references_admitted'] == fits['experimental_labels_admitted'] == 0
    assert not fits['temperature_bridge_constructed'] and fits['source_reference_uncertainty_intervals'] is None
    expected = {(t, b) for t in (70,80,90,97,110) for b in (25,50,75,150)}
    assert {(r['T_C'], r['buffer_mM']) for r in data['runs']} == expected
    assert len(data['runs']) == len(fits['runs']) == 20
    absent = [r for r in data['runs'] if r['status'] != 'measured']
    assert len(absent) == 1 and (absent[0]['T_C'], absent[0]['buffer_mM']) == (110,25)
    assert absent[0]['status'] == 'source_declared_no_experiment'
    counts = Counter()
    profile_entries = 0
    mass_changes, model_K, ratio_K, transcription_widths = [], [], [], []
    for source, result in zip(data['runs'], fits['runs'], strict=True):
        for field in ('table','T_C','buffer_mM','status'):
            assert source[field] == result[field]
        if source['status'] != 'measured':
            continue
        counts['measured_runs'] += 1
        times, values = calibrated(source)
        assert all(b>a for a,b in zip(times,times[1:]))
        counts['timestamps'] += len(times)
        for obs, calculated in zip(source['observations'], values, strict=True):
            for j, name in enumerate(('cis','trans')):
                if calculated[j] is None:
                    assert obs[name+'_mM'] is None and not obs[name+'_area_candidates']
                    counts['blank_'+name] += 1
                else:
                    close(calculated[j], obs[name+'_mM'])
                    assert obs[name+'_area'] in obs[name+'_area_candidates']
                    counts['measured_'+name] += 1
        p = result['parameters']
        parameters = [math.log(p['total_mM']), p['initial_cis_fraction'],
                      p['equilibrium_cis_fraction'], math.log(p['lambda_per_minute'])]
        close(p['K'], p['equilibrium_cis_fraction']/(1-p['equilibrium_cis_fraction']))
        close(squared_error(times,values,parameters), result['SSE_mM2'])
        for t, row, stored, residual in zip(times,values,result['predicted_concentrations_mM'],result['residuals_mM'],strict=True):
            for actual, observed, saved, error in zip(predicted(parameters,t),row,stored,residual,strict=True):
                close(actual,saved)
                if observed is None:
                    assert error is None
                else:
                    close(actual-observed,error)
        paired_totals = [sum(row) for row in values if None not in row]
        change = (paired_totals[-1]-paired_totals[0])/paired_totals[0]
        close(change,result['observed_total']['relative_first_to_last_change'])
        if result['T_C'] == 110:
            mass_changes.append(change)
            model_K.append(p['K'])
        for j,name in enumerate(('cis','trans')):
            single=result['separate_channel_fits'][name]
            error=math.fsum((single['asymptote_mM']+(single['initial_mM']-single['asymptote_mM'])*math.exp(-single['lambda_per_minute']*t)-row[j])**2
                            for t,row in zip(times,values,strict=True) if row[j] is not None)
            close(error,single['SSE_mM2'])
        zero=result['fixed_initial_zero_sensitivity'];lp,fp,lr=zero['parameters']
        close(squared_error(times,values,[lp,0.,fp,lr]),zero['SSE_mM2'])
        variant_Ks=[p['K']]
        for variant in result['transcription_sensitivity']['all_local_combinations']:
            changes={tuple(x['field']):x['value'] for x in variant['choices']}
            tt,yy=calibrated(source,changes)
            close(squared_error(tt,yy,variant['parameters']),variant['SSE_mM2'])
            f=variant['parameters'][2];close(f/(1-f),variant['K']);variant_Ks.append(variant['K'])
        close(min(variant_Ks),result['transcription_sensitivity']['K_range'][0])
        close(max(variant_Ks),result['transcription_sensitivity']['K_range'][1])
        transcription_widths.append(max(variant_Ks)-min(variant_Ks))
        assert len(result['K_profile']) == 81
        n=sum(x is not None for row in values for x in row)
        for i,point in enumerate(result['K_profile']):
            close(point['K'],10**(-3+i/20))
            fraction=point['K']/(1+point['K'])
            error=squared_error(times,values,point['nuisance_parameters'],fraction)
            close(error,point['SSE_mM2'])
            close(n*math.log(max(error/result['SSE_mM2'],1)),point['working_iid_LR'],1e-8)
            profile_kkt(times,values,point['nuisance_parameters'],point['K'])
            profile_entries += 1
        ratio=result['common_recovery_ratio_sensitivity']
        paired=[(t,row[0]/sum(row)) for t,row in zip(times,values,strict=True) if None not in row]
        assert ratio['n_paired_timestamps'] == len(paired)
        counts['paired_ratio_timestamps'] += len(paired)
        close(ratio['K'],ratio['asymptotic_cis_fraction']/(1-ratio['asymptotic_cis_fraction']))
        error=math.fsum((ratio['asymptotic_cis_fraction']+(ratio['initial_cis_fraction']-ratio['asymptotic_cis_fraction'])*math.exp(-ratio['lambda_per_minute']*t)-value)**2 for t,value in paired)
        close(error,ratio['SSE_fraction2'])
        if result['T_C']==110:
            ratio_K.append(ratio['K'])
    assert dict(counts) == dict(measured_runs=19,timestamps=314,measured_cis=309,measured_trans=314,blank_cis=5,paired_ratio_timestamps=309)
    assert profile_entries == 1539
    assert len(data['correction_ledger']) == 105
    for entry in data['correction_ledger']:
        assert sha(packet/entry['panel']) == entry['panel_sha256']
    controls=read(packet/'synthetic_controls.json')
    assert controls['primary_labels_admitted']==0
    close(controls['true_K'],controls['closed_recovery_with_missing_first_cis']['K'])
    close(controls['true_K'],controls['common_multiplicative_loss']['K'])
    assert abs(controls['selective_loss_counterexample']['fitted_ratio_K']-controls['true_K'])>.02
    refinements=read(packet/'profile_refinements.json')
    assert refinements['fit_sha256'] == hashlib.sha256(gzip.decompress((packet/'fit_diagnostics.json.gz').read_bytes())).hexdigest()
    measured_sources=[r for r in data['runs'] if r['status']=='measured']
    measured_results=[r for r in fits['runs'] if r['status']=='measured']
    assert len(refinements['runs']) == 19
    crossings=0
    for source,result,refinement in zip(measured_sources,measured_results,refinements['runs'],strict=True):
        assert (result['T_C'],result['buffer_mM']) == (refinement['T_C'],refinement['buffer_mM'])
        times,values=calibrated(source)
        n=sum(v is not None for row in values for v in row)
        for endpoint in refinement['working_iid_support_endpoints']:
            if endpoint is None:
                continue
            fraction=endpoint['K']/(1+endpoint['K'])
            error=squared_error(times,values,endpoint['nuisance_parameters'],fraction)
            close(error,endpoint['SSE_mM2'])
            close(n*math.log(error/result['SSE_mM2']),endpoint['working_iid_LR'],1e-8)
            close(endpoint['working_iid_LR'],refinements['LR_threshold'],1e-5)
            crossings += 1
    return dict(status='passed',authenticated_files=len(manifest['files']),counts=dict(counts),
                profile_entries_reconstructed=profile_entries,linear_nuisance_KKT_checks=profile_entries,
                adaptive_working_LR_crossings_reconstructed=crossings,
                high_temperature_total_change_range=[min(mass_changes),max(mass_changes)],
                high_temperature_closed_K_range=[min(model_K),max(model_K)],
                high_temperature_ratio_K_range=[min(ratio_K),max(ratio_K)],
                maximum_transcription_K_width=max(transcription_widths),
                primary_labels_admitted=0,source_experiment_independently_authenticated=False,
                source_reference_uncertainty_interval_validated=False,MetaG_accuracy_validated=False)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--packet',type=Path,default=Path(__file__).with_name('aconitate_raw_source_review'))
    args=parser.parse_args()
    def no_network(event,_args):
        if event in {'socket.connect','socket.getaddrinfo'}:
            raise RuntimeError('No network in frozen source arithmetic replay')
    sys.addaudithook(no_network)
    print(json.dumps(audit(args.packet),indent=2))


if __name__=='__main__':main()
