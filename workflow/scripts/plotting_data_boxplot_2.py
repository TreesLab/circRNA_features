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


def get_three_type_df(df):
    not_detected_column = list(filter(lambda col: col.startswith('Not detected'), df.columns))[0]
    not_depleted_column = list(filter(lambda col: col.startswith('Not depleted'), df.columns))[0]
    three_type_df = df.assign(data_type=df[not_detected_column] + (-1) * df[not_depleted_column] + 1)

    return three_type_df


dfs = [
    pd.read_csv(file_, sep='\t').set_index('event_id').pipe(get_three_type_df)
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

    sns.set_palette("bright")
    # box plot
    sns.boxplot(
        x='sample',
        y=feature_name,
        hue='data_type',
        palette={
            0: (0.50, 0.67, 0.33),
            1: (0.86, 0.51, 0.27),
            2: (0.33, 0.68, 0.92)
        },
        saturation=1,
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
        legend.texts[0].set_text('Not depleted')
        legend.texts[1].set_text('Depleted')
        legend.texts[2].set_text('Not detected')

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





boxplot('num_db', y_ticks=list(range(0, 9, 2)))

