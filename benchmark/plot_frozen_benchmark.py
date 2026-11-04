"""Plot the frozen-source benchmark against the Mg-bounded reference.

Requires matplotlib; writes a vector PDF and a 300-dpi PNG. The plot is
strictly descriptive of development-exposed reactions, not an independent test.
"""
import argparse
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

HERE=Path(__file__).resolve().parent


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--predictions',type=Path,default=HERE/'metag_frozen_source_predictions.json')
    parser.add_argument('--reference',type=Path,default=HERE/'reference_mg_bounded.json')
    parser.add_argument('--output-prefix',type=Path,required=True)
    args=parser.parse_args(argv)
    predictions=json.loads(args.predictions.read_text())['predictions']
    labels=json.loads(args.reference.read_text())['2.0']['labels']
    ids=sorted(r for r in labels if predictions.get(r,{}).get('status')=='estimated')
    if not ids:raise ValueError('No common predictions and references')
    reference=np.asarray([labels[r]['exp'] for r in ids],dtype=float)
    predicted=np.asarray([predictions[r]['dG'] for r in ids],dtype=float)
    if not np.isfinite(reference).all() or not np.isfinite(predicted).all():raise ValueError('Nonfinite values')
    plt.rcParams.update({'font.size':9,'axes.labelsize':9,'xtick.labelsize':8,'ytick.labelsize':8,
                         'legend.fontsize':8,'pdf.fonttype':42,'ps.fonttype':42})
    fig,(ax1,ax2)=plt.subplots(1,2,figsize=(7.0,3.15),layout='constrained')
    bound=max(10.,float(np.max(np.abs(np.concatenate((reference,predicted))))));bound=10*np.ceil((bound+5)/10)
    diagonal=np.array([-bound,bound]);ax1.fill_between(diagonal,diagonal-10,diagonal+10,color='0.9',zorder=0,label=r'$\pm 10$ kJ mol$^{-1}$')
    ax1.plot(diagonal,diagonal,color='0.3',lw=1,zorder=1)
    ax1.scatter(reference,predicted,s=13,color='#0072B2',alpha=.72,linewidths=0,zorder=2)
    ax1.set(xlim=(-bound,bound),ylim=(-bound,bound),xlabel=r'Reference $\Delta_r G^{\prime\circ}$ (kJ mol$^{-1}$)',ylabel=r'MetaG $\Delta_r G^{\prime\circ}$ (kJ mol$^{-1}$)')
    ax1.set_aspect('equal',adjustable='box');ax1.legend(loc='upper left',frameon=False)
    for values,name,color,style in ((np.abs(predicted-reference),'MetaG','#0072B2','-'),(np.abs(reference),'Zero predictor','#D55E00','--')):
        ordered=np.sort(values);ax2.step(ordered,np.arange(1,len(ordered)+1)/len(ordered),where='post',label=name,color=color,lw=1.6,linestyle=style)
    ax2.set(xlim=(0,max(float(np.max(np.abs(predicted-reference))),float(np.max(np.abs(reference))))*1.05),ylim=(0,1.02),
            xlabel=r'Absolute error (kJ mol$^{-1}$)',ylabel='Fraction of reactions')
    ax2.legend(loc='lower right',frameon=False)
    for label,axis in zip(('A','B'),(ax1,ax2)):
        axis.text(-.16,1.04,label,transform=axis.transAxes,fontweight='bold',fontsize=10)
        axis.spines[['top','right']].set_visible(False)
    args.output_prefix.parent.mkdir(parents=True,exist_ok=True)
    fig.savefig(args.output_prefix.with_suffix('.pdf'),bbox_inches='tight')
    fig.savefig(args.output_prefix.with_suffix('.png'),dpi=300,bbox_inches='tight')
    plt.close(fig)
    print(f'Plotted {len(ids)} identical-target reaction pairs')

if __name__=='__main__':main()
