"""Reconstruct a GPDH reference candidate without admitting a thermodynamic label."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import re
import xml.etree.ElementTree as ET


def nominal_exchange(rates, ga_scale, dhap_scale, temperature=298.15):
    """Undo reported concentration-axis scales, not a fitted physical correction.

    The ratio is thermodynamic only if microscopic reversibility, common assay
    environment and the declared chemical/state conventions are established.
    """
    required=['forward_G3P','reverse_DHAP_carbonyl','forward_EtG','reverse_GA_carbonyl']
    values=[rates[k] for k in required]+[ga_scale,dhap_scale,temperature]
    if any(isinstance(v,bool) or not isinstance(v,(int,float)) or not math.isfinite(v) or v<=0
           for v in values):
        raise ValueError('Rates, scales and temperature must be positive finite numbers')
    if ga_scale>1 or dhap_scale>1:
        raise ValueError('Concentration-axis fractions must not exceed one')
    whole=rates['forward_G3P']/rates['reverse_DHAP_carbonyl']
    truncated=rates['forward_EtG']/rates['reverse_GA_carbonyl']
    carbonyl=whole/truncated
    total=carbonyl*ga_scale/dhap_scale
    rt=0.00831446261815324*temperature
    return dict(whole_assay_ratio=whole,truncated_assay_ratio=truncated,
        carbonyl_exchange_ratio=carbonyl,total_substrate_exchange_ratio=total,
        carbonyl_minus_RTln_ratio_kj=-rt*math.log(carbonyl),
        total_minus_RTln_ratio_kj=-rt*math.log(total),
        normalization_change_kj=-rt*math.log(ga_scale/dhap_scale),
        reverse_total_substrate_rates_M_inv_s_inv=dict(
            DHAP=rates['reverse_DHAP_carbonyl']*dhap_scale,
            GA=rates['reverse_GA_carbonyl']*ga_scale))


def verify_primary_sources(facts, folder):
    """Optional authenticated extraction from unchanged original primary XML."""
    roots={}
    for source in facts['source_records']:
        path=Path(folder)/source['filename']
        if hashlib.sha256(path.read_bytes()).hexdigest()!=source['sha256']:
            raise ValueError('Changed primary source: '+source['filename'])
        roots[source['name']]=ET.parse(path)
    rows=roots['Cristobal2024'].find('.//table-wrap[@id="tbl1"]').findall('.//tbody/tr')
    names=['forward_G3P','reverse_DHAP_carbonyl','forward_EtG','reverse_GA_carbonyl']
    for tr,name in zip(rows[:4],names):
        cell=list(tr)[5]
        text=' '.join(''.join(cell.itertext()).split())
        number=float(re.match(r'\(?([\d.]+)',text)[1])
        exponents=list(cell.iter('sup'))
        if exponents:
            assert len(exponents)==1
            number*=10**int(''.join(exponents[0].itertext()).replace('–','-').replace('−','-'))
        if not math.isclose(number,facts['rates_M_inv_s_inv'][name],rel_tol=1e-12):
            raise ValueError('Rate transcription mismatch: '+name)
    paragraphs=[' '.join(''.join(p.itertext()).split()) for p in roots['Reyes2015'].iter('p')]
    assert any('fcar = 0.45' in p and 'DHAP' in p for p in paragraphs)
    assert any('fcar = 0.06' in p and 'GA' in p for p in paragraphs)
    assert any('10 mM TEA' in p and 'I = 0.12' in p for p in paragraphs)
    return True


def audit(packet, primary_dir=None):
    packet=Path(packet)
    facts=json.loads((packet/'source_facts.json').read_text())
    main=nominal_exchange(facts['rates_M_inv_s_inv'],0.06,0.45)
    alternative=nominal_exchange(facts['rates_M_inv_s_inv'],0.06,0.55)
    printed=facts['reported_equilibrium_constants']
    printed_ratio=printed['K1']/printed['K2']
    rate_ratio=main['whole_assay_ratio']
    # Table and the inspected original figure agree, and arithmetic separates
    # the tenfold prose discrepancy from harmless rounding of printed K1.
    assert abs(rate_ratio/printed['K1']-1)<0.02
    assert abs(rate_ratio/facts['conflicting_prose_K1']-1)>0.8
    assert printed['K1']==facts['figure_K1']
    return dict(schema='metag-gpdh-exchange-candidate-audit-v1',
        status='nominal_assay_ratios_reconstructed; thermodynamic_reference_not_admitted',
        equation='(R)-glycerol3-phosphate + glycolaldehyde -> DHAP + ethylene glycol',
        conditions=facts['conditions'],reported_second_order_rate_basis=main,
        printed_equilibrium_basis=dict(carbonyl_exchange_ratio=printed_ratio,
            total_substrate_exchange_ratio=printed_ratio*0.06/0.45),
        alternative_DHAP_normalization_scenario=alternative,
        scenario_difference_kj=alternative['total_minus_RTln_ratio_kj']-main['total_minus_RTln_ratio_kj'],
        scenario_scope='Source decoding alternatives, not a confidence interval or a measured hydration-temperature correction.',
        primary_sources_verified=(verify_primary_sources(facts,primary_dir) if primary_dir is not None else False),
        source_scope='Default replay checks curated source facts and arithmetic. Use --primary-dir for hash-authenticated XML extraction. Original articles are linked, not redistributed.',
        cancellation_scope='Common enzyme-bound cofactor symbols cancel algebraically; a single thermodynamic environment is an additional unverified prerequisite.',
        uncertainty_scope='No numerical confidence interval: full rate covariance, concentration-scale history, differential buffer effects and condition-transfer uncertainties are unresolved.',
        limitations=facts['limitations'],primary_references_admitted=0,
        secondary_thermodynamic_labels_admitted=0,predictions_examined=False)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--packet',type=Path,default=Path(__file__).with_name('gpdh_exchange_reference_review'))
    parser.add_argument('--primary-dir',type=Path)
    parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    result=audit(args.packet,args.primary_dir)
    if args.output:
        with args.output.open('x') as stream:json.dump(result,stream,indent=2,allow_nan=False);stream.write('\n')
    print(json.dumps(result,indent=2,allow_nan=False))


if __name__=='__main__':main()
