import os
import os.path
import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import subprocess as sp
from operator import itemgetter
from functools import wraps

# import matplotlib
# print(matplotlib.__version__)


os.makedirs(snakemake.output[0], exist_ok=True)

df = pd.read_csv(snakemake.input[0], sep='\t', index_col=[0, 1])

sample_ids = [re.sub(r'\(.*\)$', '', sample_id) for sample_id in df.columns]
# print(sample_ids)

# scoring_features = df.index.get_level_values(0).unique()

scoring_features = [
    '#BSJ_reads(totalRNA)',
    'MCS-detail(species)',
    'MCS-detail(tissue)',
    'MCS-detail(sample)',
    'Multiple Conservation Score(MCS)',
    'Tissue Specifity Index',
    'phyloP(acceptor_in)',
    'phyloP(acceptor_out)',
    'phyloP(donor_in)',
    'phyloP(donor_out)',
    'phastCons(acceptor_in)',
    'phastCons(acceptor_out)',
    'phastCons(donor_in)',
    'phastCons(donor_out)',
    'MAXENT(acceptor)',
    'MAXENT(donor)',
    'MM(acceptor)',
    'MM(donor)',
    'WMM(acceptor)',
    'WMM(donor)',
    'evidence_num_plus',
    'min_dist_of_flanking_RBPs',
    'evidences_score',
    'num_db'
]



def plot_both_simplified_and_original(func):
    '''
    Plot both the simplified and the original figures
    '''

    @wraps(func)
    def wrapper(*args, **kwargs):
        if 'out_prefix' in kwargs:
            out_prefix = kwargs.pop('out_prefix')
        else:
            out_prefix = kwargs.pop('feature_name', args[1])

        simplified_out_prefix = out_prefix
        out_prefix = f'{out_prefix}.bk'

        _ = kwargs.pop('simplify', None)

        func(*args, **kwargs, out_prefix=out_prefix)
        func(*args, **kwargs, out_prefix=simplified_out_prefix, simplify=True)

    return wrapper


def get_effsize_and_conf_interval(df, feature_name):
    all_effsize = df.loc[(feature_name, 'effect_size(two-sided)(R_ranksums)'), :].values
    all_conf_int_high = df.loc[(feature_name, 'conf_int_high(effect_size)(two-sided)(R_ranksums)'), :].values
    all_conf_int_low = df.loc[(feature_name, 'conf_int_low(effect_size)(two-sided)(R_ranksums)'), :].values

    return all_effsize, all_conf_int_high, all_conf_int_low


@plot_both_simplified_and_original
def plot_effsize_and_conf_interval(df, feature_name, simplify=False, width=0.3, y_ticks=None, y_ticks_labels=None, out_prefix=None):
    all_effsize, all_conf_int_high, all_conf_int_low = get_effsize_and_conf_interval(df, feature_name)

    stats = [
        {
            'med': e,
            'q1': e,
            'q3': e,
            'whislo': l,
            'whishi': h
        }
        for e, h, l in zip(all_effsize, all_conf_int_high, all_conf_int_low)
    ]


    fig, ax = plt.subplots(figsize=(13, 1.5))

    ax.bxp(stats, showfliers=False, medianprops={'color': 'black'}, showbox=False, widths=width, capwidths=width)

    if simplify:
        ax.set(xticklabels=[])
        ax.tick_params(bottom=False)
    else:
        ax.set(xticklabels=sample_ids)
        plt.xticks(rotation=90)
        ax.set(ylabel='Effect Size')
        ax.set_title(feature_name)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    if y_ticks:
        if y_ticks_labels:
            plt.yticks(y_ticks, y_ticks_labels)
        else:
            plt.yticks(y_ticks)

    ax.tick_params(axis='y', labelsize=12)

    for tick in ax.get_yticklabels():
        tick.set_fontname('Times New Roman')


    # save figure
    if out_prefix is None:
        out_prefix = feature_name

    svg_out = os.path.join(snakemake.output[0], f'{out_prefix}.svg')
    svg_pdf_out = f'{svg_out}.pdf'

    fig.savefig(svg_out, format='svg', transparent=True, bbox_inches='tight')

    sp.run(['cairosvg', svg_out, '-o', svg_pdf_out])

    plt.close(fig)





plot_effsize_and_conf_interval(df, '#BSJ_reads(totalRNA)', out_prefix='BSJ_reads')
plot_effsize_and_conf_interval(df, 'MCS-detail(species)')
plot_effsize_and_conf_interval(df, 'MCS-detail(tissue)')
plot_effsize_and_conf_interval(df, 'MCS-detail(sample)')
plot_effsize_and_conf_interval(df, 'Multiple Conservation Score(MCS)')
plot_effsize_and_conf_interval(df, 'Tissue Specifity Index')
plot_effsize_and_conf_interval(df, 'phyloP(acceptor_in)')
plot_effsize_and_conf_interval(df, 'phyloP(acceptor_out)')
plot_effsize_and_conf_interval(df, 'phyloP(donor_in)')
plot_effsize_and_conf_interval(df, 'phyloP(donor_out)')
plot_effsize_and_conf_interval(df, 'phastCons(acceptor_in)')
plot_effsize_and_conf_interval(df, 'phastCons(acceptor_out)')
plot_effsize_and_conf_interval(df, 'phastCons(donor_in)')
plot_effsize_and_conf_interval(df, 'phastCons(donor_out)')
plot_effsize_and_conf_interval(df, 'MAXENT(acceptor)')
plot_effsize_and_conf_interval(df, 'MAXENT(donor)')
plot_effsize_and_conf_interval(df, 'MM(acceptor)')
plot_effsize_and_conf_interval(df, 'MM(donor)')
plot_effsize_and_conf_interval(df, 'WMM(acceptor)')
plot_effsize_and_conf_interval(df, 'WMM(donor)')
plot_effsize_and_conf_interval(df, 'evidence_num_plus')
plot_effsize_and_conf_interval(df, 'min_dist_of_flanking_RBPs')
plot_effsize_and_conf_interval(df, 'evidences_score')
plot_effsize_and_conf_interval(df, 'num_db')

