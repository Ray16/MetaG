"""Plot a retained-cohort reliability report; no model inference or curve smoothing."""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('report', type=Path)
    parser.add_argument('output_stem', type=Path)
    args = parser.parse_args()
    report = json.loads(args.report.read_text())
    if report['schema'] != 'metag-equal-coverage-reliability-v1':
        raise ValueError('Unsupported reliability report')
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 8,
                         'axes.labelsize': 8, 'axes.titlesize': 9,
                         'pdf.fonttype': 42, 'svg.fonttype': 'none'})
    styles = [('spread', 'CV repeat spread', '#666666', '-'),
              ('disagreement', 'MetaG–comparator disagreement', '#0072B2', '--'),
              ('rank_sum', 'Equal-rank combination', '#D55E00', '-.'),
              ('prediction_magnitude', 'Comparator magnitude control', '#CC79A7', ':')]
    fig, axes = plt.subplots(1, 2, figsize=(7.1, 3.4))
    handles = []
    for name, label, color, line in styles:
        values = report['curves'][name]
        x = [r['coverage'] for r in values]
        a, = axes[0].plot(x, [r['MAE']['comparator'] for r in values],
                         label=label, color=color, linestyle=line, linewidth=1.5)
        handles.append(a)
        axes[1].plot(x, [r['expected_severe_count']['comparator'] for r in values],
                     color=color, linestyle=line, linewidth=1.5)
    axes[0].set_ylabel('Comparator MAE (kJ/mol)')
    axes[1].set_ylabel('Comparator errors >20 kJ/mol (count)')
    axes[1].set_yticks(range(report['full_cohort_comparator_severe_count']+1))
    for axis, letter in zip(axes, ['A', 'B']):
        axis.set_xlim(0, 1.01)
        axis.set_ylim(bottom=0)
        axis.set_xlabel('Retained fraction')
        axis.set_xticks([0, .25, .5, .75, 1])
        axis.grid(alpha=.18, linewidth=.5)
        axis.spines[['top', 'right']].set_visible(False)
        axis.text(.02, .98, letter, transform=axis.transAxes, va='top', fontweight='bold', fontsize=10)
    fig.suptitle(f'Development cohort: {report["n"]} matched reactions', y=.99, fontsize=10)
    fig.legend(handles=handles, loc='lower center', ncol=2, frameon=False, fontsize=8,
               bbox_to_anchor=(.5, -.005))
    fig.tight_layout(rect=(0, .14, 1, .96), w_pad=2.)
    args.output_stem.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output_stem.with_suffix('.pdf'), bbox_inches='tight',
                metadata={'CreationDate': None, 'ModDate': None})
    fig.savefig(args.output_stem.with_suffix('.png'), dpi=250, bbox_inches='tight')
    plt.close(fig)


if __name__ == '__main__':
    main()
