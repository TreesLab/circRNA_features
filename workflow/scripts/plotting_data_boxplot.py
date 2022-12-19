"""
Features:
    - #BSJ_reads(totalRNA)
    - MCS-detail(species)
    - MCS-detail(tissue)
    - MCS-detail(sample)
    - Multiple Conservation Score(MCS)
    - Tissue Specifity Index
    - phyloP(acceptor_in)
    - phyloP(acceptor_out)
    - phyloP(donor_in)
    - phyloP(donor_out)
    - phastCons(acceptor_in)
    - phastCons(acceptor_out)
    - phastCons(donor_in)
    - phastCons(donor_out)
    - MAXENT(acceptor)
    - MAXENT(donor)
    - MM(acceptor)
    - MM(donor)
    - WMM(acceptor)
    - WMM(donor)
    - evidence_num_TransCirc+miRNA/RBP sites across the BSJ
    - evidences_num (translational potential)
    - evidences_score (translational potential)
"""


import os
import os.path
import re
import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
import subprocess as sp
from operator import itemgetter
from functools import wraps



os.makedirs(snakemake.output[0], exist_ok=True)


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
    'evidences_num',
    'evidences_score',
    'num_db'
]


def rename_not_depleted_column(df):
    columns = df.columns
    not_depleted_column = list(filter(lambda col: col.startswith('Not depleted'), columns))[0]
    return df.rename({not_depleted_column: "not_depleted"}, axis=1)


dfs = [
    pd.read_csv(file_, sep='\t').set_index('event_id').iloc[:, 1:].pipe(rename_not_depleted_column)
    for file_ in snakemake.input
]



sample_ids = [
    f'{tissue}({dataset_id})'
    for dataset_id, species, tissue in snakemake.params.groups
]

sample_id_mapper = snakemake.params.sample_id_mapper

new_sample_infos = [sample_id_mapper[id_] for id_ in sample_ids]
new_sample_ids, new_sample_order = zip(*new_sample_infos)
new_sample_ids = [re.sub(r'\(.+\)$', '', id_) for id_ in new_sample_ids]


merged_df = pd.concat([df.assign(sample=id_) for df, id_ in zip(dfs, new_sample_ids)])


sorted_new_sample_ids, _ = zip(*sorted(zip(new_sample_ids, new_sample_order), key=itemgetter(1)))
sorted_new_sample_ids = list(sorted_new_sample_ids)

merged_df.to_csv(os.path.join(snakemake.output[0], 'merged_df.tsv'), sep='\t')

# test

plt.clf()

fig, ax = plt.subplots(figsize=(13, 1.5))


sns.boxplot(
    x='sample',
    y='#BSJ_reads(totalRNA)',
    hue='not_depleted',
    palette={
        0: (1, 1, 1),
        1: sns.color_palette()[7]
    },
    order=sorted_new_sample_ids,
    showfliers=False,
    data=merged_df,
    ax=ax
)


ax.set(xlabel=None)
ax.set(xticklabels=[])
ax.tick_params(bottom=False)
# plt.xticks(ha='right', rotation=30)


ax.set(ylabel=None)
plt.yticks([2] + list(range(5, 30, 5)), [f'tick_{i}' for i in range(6)])

ax.tick_params(axis='y', labelsize=12)

for tick in ax.get_yticklabels():
    tick.set_fontname('Times New Roman')


legend = ax.legend_
# legend.set_title(None)
# legend.texts[0].set_text('Depleted')
# legend.texts[1].set_text('Not depleted')
legend.remove()


sns.despine(top=True, right=True)
# sns.set_style({'font.family': 'Times New Roman'})


# fig.savefig(os.path.join(snakemake.output[0], 'BSJ_reads.svg'), format='svg', transparent=True)
# fig.savefig(os.path.join(snakemake.output[0], 'BSJ_reads.pdf'), format='pdf', transparent=True)
# fig.savefig(os.path.join(snakemake.output[0], 'BSJ_reads.png'), format='png', transparent=True, dpi=300)


svg_out = os.path.join(snakemake.output[0], 'BSJ_reads_test.svg')
svg_pdf_out = os.path.join(snakemake.output[0], 'BSJ_reads_test.svg.pdf')

fig.savefig(svg_out, format='svg', transparent=True, bbox_inches='tight')

sp.run(['cairosvg', svg_out, '-o', svg_pdf_out])



# modularize 

def plot_both_simplified_and_original(func):
    '''
    Plot both the simplified and the original figures
    '''

    @wraps(func)
    def wrapper(*args, **kwargs):
        if 'out_prefix' in kwargs:
            out_prefix = kwargs.pop('out_prefix')
        else:
            out_prefix = kwargs.pop('feature_name', args[0])

        simplified_out_prefix = out_prefix
        out_prefix = f'{out_prefix}.bk'

        _ = kwargs.pop('simplify', None)

        func(*args, **kwargs, out_prefix=out_prefix)
        func(*args, **kwargs, out_prefix=simplified_out_prefix, simplify=True)

    return wrapper


@plot_both_simplified_and_original
def boxplot(feature_name, y_ticks=None, y_ticks_labels=None, out_prefix=None, simplify=False):
    # create figure
    fig, ax = plt.subplots(figsize=(13, 1.5))

    # box plot
    sns.boxplot(
        x='sample',
        y=feature_name,
        hue='not_depleted',
        palette={
            0: (1, 1, 1),
            1: sns.color_palette()[7]
        },
        order=sorted_new_sample_ids,
        showfliers=False,
        data=merged_df,
        ax=ax
    )

    # styling

    ## x-axis
    ax.set(xlabel=None)

    if simplify:
        ax.set(xticklabels=[])
        ax.tick_params(bottom=False)    
    else:
        plt.xticks(rotation=90)

    ## y-axis
    if simplify:
        ax.set(ylabel=None)

    if y_ticks:
        if y_ticks_labels:
            plt.yticks(y_ticks, y_ticks_labels)
        else:
            plt.yticks(y_ticks)

    ax.tick_params(axis='y', labelsize=12)

    for tick in ax.get_yticklabels():
        tick.set_fontname('Times New Roman')

    ## legend
    legend = ax.legend_
    if simplify:
        legend.remove()
    else:
        legend.set_title(None)
        legend.texts[0].set_text('Depleted')
        legend.texts[1].set_text('Not depleted')

    ## figure border
    sns.despine(top=True, right=True)

    # save figure
    if out_prefix is None:
        out_prefix = feature_name

    svg_out = os.path.join(snakemake.output[0], f'{out_prefix}.svg')
    svg_pdf_out = f'{svg_out}.pdf'

    fig.savefig(svg_out, format='svg', transparent=True, bbox_inches='tight')

    sp.run(['cairosvg', svg_out, '-o', svg_pdf_out])

    plt.close(fig)





# test boxplot()
# boxplot('#BSJ_reads(totalRNA)', y_ticks=[2] + list(range(5, 30, 5)), out_prefix='BSJ_reads')
# boxplot('MCS-detail(species)', y_ticks=[1] + list(range(3, 9, 2)))

unit_ticks = [0, 0.5, 1]
unit_ticks_label = ['0', '0.5', '1']


boxplot('#BSJ_reads(totalRNA)', y_ticks=[2] + list(range(5, 30, 5)), out_prefix='BSJ_reads')
boxplot('MCS-detail(species)', y_ticks=[1] + list(range(3, 11, 2)))
boxplot('MCS-detail(tissue)', y_ticks=[1] + list(range(5, 25, 5)))
boxplot('MCS-detail(sample)', y_ticks=[1] + list(range(50, 300, 50)))
boxplot('Multiple Conservation Score(MCS)', y_ticks=[1] + list(range(3, 10, 2)))
boxplot('Tissue Specifity Index', y_ticks=unit_ticks, y_ticks_labels=unit_ticks_label)
boxplot('phyloP(acceptor_in)', y_ticks=list(range(-3, 13, 3)))
boxplot('phyloP(acceptor_out)', y_ticks=list(range(-3, 13, 3)))
boxplot('phyloP(donor_in)', y_ticks=list(range(-3, 13, 3)))
boxplot('phyloP(donor_out)', y_ticks=list(range(-3, 13, 3)))
boxplot('phastCons(acceptor_in)', y_ticks=unit_ticks, y_ticks_labels=unit_ticks_label)
boxplot('phastCons(acceptor_out)', y_ticks=unit_ticks, y_ticks_labels=unit_ticks_label)
boxplot('phastCons(donor_in)', y_ticks=unit_ticks, y_ticks_labels=unit_ticks_label)
boxplot('phastCons(donor_out)', y_ticks=unit_ticks, y_ticks_labels=unit_ticks_label)
boxplot('MAXENT(acceptor)', y_ticks=np.linspace(-13, 17, 4).tolist())
boxplot('MAXENT(donor)', y_ticks=np.linspace(-3, 15, 4).tolist())
boxplot('MM(acceptor)', y_ticks=np.linspace(-11, 19, 4).tolist())
boxplot('MM(donor)', y_ticks=np.linspace(0, 15, 4).tolist())
boxplot('WMM(acceptor)', y_ticks=np.linspace(-18, 21, 4).tolist())
boxplot('WMM(donor)', y_ticks=np.linspace(-3, 15, 4).tolist())
boxplot('evidence_num_plus', y_ticks=list(range(0, 9, 2)))
boxplot('evidences_num', y_ticks=list(range(0, 7, 2)))
boxplot('evidences_score', y_ticks=np.linspace(0, 5, 3).tolist(), y_ticks_labels=['0', '2.5', '5'])
boxplot('num_db', y_ticks=list(range(0, 9, 2)))





