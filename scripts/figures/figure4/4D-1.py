#!/usr/bin/env python3

# --- Aim
# Plot Fig. 4d (part 1): ranked correlation plot of genes with STING expression in HEVs
# Shows genes that correlate positively or negatively with STING pathway activation

# --- Define environment
# env: shapely

# --- Load packages

import argparse
import os

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

import pandas as pd
import numpy as np

import seaborn as sns
from adjustText import adjust_text


mpl.rcParams['pdf.fonttype'] = 42


# --- Argparse arguments

parser = argparse.ArgumentParser(description="Fig. 4d (part 1): genes ranked by correlation with STING expression in HEVs",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument("-s", "--samples",
    nargs="*",  # 0 or more values expected => creates a list
    type=str,
    help="List of samples")
parser.add_argument("-r", "--res_dir", help="Path to results directory")
parser.add_argument("--h5ad_dir", help="Path to directory for .h5ad output files")
parser.add_argument("--project_dir", help="Project directory (root)")
parser.add_argument("--sample_sheet", help="Path to sample sheet CSV file", default=None)

args = vars(parser.parse_args())

# Arguments
project_dir = args["project_dir"]
h5ad_dir = args["h5ad_dir"]
res_dir = args["res_dir"]
sample_list = args["samples"]


# --- Load the correlation results

gene_of_i = 'TMEM173'  # STING gene symbol

correlation_dir = os.path.join(res_dir, '05_xenium_mechanisms', '05_gene_gene_correlations', '01_investigate_gene_correlations', f'{gene_of_i}_correlations')
most_recent = correlation_dir

correlation_file = f'All_samples_HEV_{gene_of_i}_gene_correlations.csv'
correlation_df = pd.read_csv(os.path.join(most_recent, 'csv_files', correlation_file), index_col=0)

# Filter for significance
pval_threshold = 0.05
significant_corr_df = correlation_df[correlation_df['p_value'] < pval_threshold]


# --- Plot 1: Plot the correlations

# Add grid behind
plt.rcParams['axes.axisbelow'] = True

fig, ax = plt.subplots(figsize=(5, 2.5))

sns.scatterplot(data=correlation_df, x='rank', y='correlation', color='lightgrey', edgecolor=None, s=2, ax=ax)
sns.scatterplot(data=significant_corr_df, x='rank', y='correlation', color='red', edgecolor=None, s=2, ax=ax)

# Annotate some genes
annotate_genes = {
    'CAMs': ['MCAM', 'ICAM1', 'ICAM2', 'VCAM1', 'BCAM', 'PECAM1'],
    'Chemokines': ['CXCL10', 'CXCL9', 'CCL5', 'CCL2', 'CCL19', 'CCL21', 'CXCL13'],
    'Selectins': ['SELE', 'SELP'],
}

# Define colors for each category
category_colors = {
    'CAMs': 'darkblue',
    'Chemokines': 'darkgreen',
    'Selectins': 'purple',
}

texts = []

for category in annotate_genes:
    color = category_colors[category]
    for gene in annotate_genes[category]:
        if gene not in correlation_df.index:
            continue
        x = correlation_df.loc[gene, 'rank']
        y = correlation_df.loc[gene, 'correlation']
        
        # Add a circle around the annotated point
        ax.scatter(x, y, s=4, facecolors='none', edgecolors='black', linewidths=0.5, zorder=5)
        
        text = ax.annotate(
            gene,
            xy=(x, y),
            fontsize=6,
            color=color,
            ha='center',
            va='center',
            bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor='none', alpha=0.7),
            zorder=10
        )
        texts.append(text)

# Automatically adjust text positions to avoid overlaps
adjust_text(
    texts,
    arrowprops=dict(arrowstyle='-', color='black', lw=0.7, alpha=0.8, shrinkA=5, shrinkB=1),
    ax=ax,
    expand=(1.5, 2.0),
    force_points=(0.3, 0.5),
    force_text=(0.3, 0.5),
    only_move={'text': 'xy'},
    lim=500
)

# Add padding to axis limits to prevent text from being too close to edges
xlim = ax.get_xlim()
x_range = xlim[1] - xlim[0]
ax.set_xlim(xlim[0] - 0.02 * x_range, xlim[1] + 0.02 * x_range)

# Add grid behind plot
ax.grid(which='major', linestyle='-', color='lightgrey')

# Also add a line at the middle rank
mid_rank = len(correlation_df) / 2
ax.axvline(mid_rank, color='lightgray', lw=2, ls='-', zorder=-1)

# Adjust labels
ax.set_xlabel(f'Gene rank',)
ax.set_ylabel('Pearson r', )

# Add a legend for pvalue
legend1 = ax.legend(title=f'p<{pval_threshold}', labels=['False', 'True'], loc='upper left',  frameon=False, markerscale=3, bbox_to_anchor=(1, 1),
                    alignment='left', fontsize=8, title_fontsize=9)

# Add a legend for gene categories (text-only, colored labels)
category_handles = [Line2D([0], [0], marker='None', linestyle='None', label=cat) 
                   for cat in annotate_genes.keys()]

legend2 = ax.legend(handles=category_handles, title='Gene category', loc='lower left', 
         frameon=False, bbox_to_anchor=(1, 0), alignment='left', fontsize=8, title_fontsize=9,
         handlelength=0, handletextpad=0)

# Color the legend text labels
for text, cat in zip(legend2.get_texts(), annotate_genes.keys()):
    text.set_color(category_colors[cat])

# Add back the first legend
ax.add_artist(legend1)

# Add a title
ax.set_title(f'Correlation with {gene_of_i} expression', fontweight='bold')

plt.suptitle(f'All samples - Correlation between mean {gene_of_i} expression and other genes per HEV niche', fontweight='bold', y=1.1)


# --- Save the output

out_dir = os.path.join(res_dir, 'figures')
os.makedirs(out_dir, exist_ok=True)

filename = 'Fig_4D-1.pdf'
plt.savefig(os.path.join(out_dir, filename), bbox_inches='tight')


# --- Finish

print("Script successfully completed")
