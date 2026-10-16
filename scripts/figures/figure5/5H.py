#!/usr/bin/env python3

# --- Aim
# Plot Figure 5H: Pearson correlations of HEV SPP1 expression with cell types across samples
# Accounts for variation in cell composition and HEV SPP1 expression levels

# --- Define environment
# env: xen-hev

# --- Load packages
import argparse
import os
import sys

import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
from scipy import stats
import numpy as np

import seaborn as sns

mpl.rcParams['pdf.fonttype'] = 42


# --- Argparse arguments

parser = argparse.ArgumentParser(description="Fig. 5h: correlation of HEV SPP1 expression with cell type abundance",
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

# Import function script as a module
sys.path.append(os.path.join(project_dir, 'scripts', '05_xenium_mechanisms',))
import plotting_functions


# --- Load plotting colors

color_dir = os.path.join(project_dir, 'resources', 'plotting_colors')

# Use the major level for colors
color_level = 'major'
color_file = f'miniatlas_{color_level}_colors.csv'

color_dict = pd.read_csv(os.path.join(color_dir, color_file), index_col=0).to_dict()['color']


# --- Load the miniatlas classifications

# Decide which annotation levels to use
level = 'subset'

# These are needed to check the alignment between major and minor levels
miniatlas_dir = os.path.join(project_dir, 'resources', 'cell_types', 'miniatlas')
miniatlas_file = 'miniatlas_subset_classifications.csv'

miniatlas_df = pd.read_csv(os.path.join(miniatlas_dir, miniatlas_file),)

# For each cell type in the miniatlas at this level, check if it is in the color_dict
# If not, get the corresponding major cell type and use its color
for cell_type in miniatlas_df[f'miniatlas_{level}'].unique():
    if cell_type not in color_dict.keys():
        # Check that the value is not nan
        if pd.isna(cell_type):
            continue
        print(f'Adding color for {cell_type}')
    # Get the corresponding major cell type
    major_type = miniatlas_df[miniatlas_df[f'miniatlas_{level}'] == cell_type][f'miniatlas_major'].values
    # Find which value is in the color_dict
    major_type_in_dict = [m for m in major_type if m in color_dict.keys()]
    if len(major_type_in_dict) > 0:
        color_dict[cell_type] = color_dict[major_type_in_dict[0]]
    else:
        color_dict[cell_type] = 'lightgrey'

# Also define a dict with shorter names for the subset level cell types
miniatlas_sub_dict = miniatlas_df.set_index('miniatlas_subset')[['miniatlas_subset_short']].to_dict()


# --- Sort by lineage, then by major

lineage_order = ['Epithelial', 'Immune', 'Stromal']

major_cell_types = miniatlas_df[f'miniatlas_major'].dropna().unique().tolist()
minor_cell_types = miniatlas_df[f'miniatlas_minor'].dropna().unique().tolist()
subset_cell_types = miniatlas_df[f'miniatlas_subset'].dropna().unique().tolist()

# Sort based on lineage first (recorded in the miniatlas df)
def get_lineage(cell_type, level='major'):
    if cell_type in miniatlas_df[f'miniatlas_{level}'].values:
        lineage = miniatlas_df[miniatlas_df[f'miniatlas_{level}'] == cell_type]['lineage'].values[0]
    else:
        lineage = 'other'
    return lineage

major_cell_types_sorted = sorted(major_cell_types, key=lambda x: (lineage_order.index(get_lineage(x)) if get_lineage(x) in lineage_order else len(lineage_order), x))


# --- Load SPP1 correlations calculated per sample

# Position legend function
def position_legend(axs):
    """
    Given an array (or list) of axes, find the rightmost right_edge and the highest top_edge
    in figure coordinates, for robust legend placement.
    Returns (max_right_edge, max_top_edge)
    """
    # Accept a single axis as well
    if not isinstance(axs, (list, np.ndarray)):
        axs = [axs]
    fig = axs[0].figure
    fig.canvas.draw()
    right_edges = []
    top_edges = []
    left_edges = []
    bottom_edges = []
    for ax in axs:
        bbox = ax.get_tightbbox(fig.canvas.get_renderer())
        inv_trans = fig.transFigure.inverted()
        bbox_fig = inv_trans.transform(bbox)
        # bbox_fig[0] = (x0, y0) = bottom left; bbox_fig[1] = (x1, y1) = top right
        left_edges.append(bbox_fig[0][0])
        bottom_edges.append(bbox_fig[0][1])
        right_edges.append(bbox_fig[1][0])
        top_edges.append(bbox_fig[1][1])  # y1 is the top in figure coordinates
    max_right_edge = max(right_edges)
    max_top_edge = max(top_edges)
    max_left_edge = max(left_edges)
    max_bottom_edge = max(bottom_edges)
    return max_right_edge, max_top_edge, max_left_edge, max_bottom_edge

# Significance level (threshold)
alpha = 0.05

def analyze_cell_type_correlation(correlations, cell_name):
    """
    Performs a one-sample t-test on Fisher-transformed correlation values.
    """
    # 1. Handle potential r-values of 1.0 or -1.0
    # np.arctanh is infinite at 1 or -1, which breaks the t-test.
    # We clip them to a very large/small value if they exist.
    n = len(correlations)
    if n < 3:
        print("Cannot perform t-test: need at least 3 samples.")
        return
    # Clip values to avoid infinity
    correlations_clipped = np.clip(correlations, -0.99999, 0.99999)
    # 2. Apply Fisher r-to-z transformation
    # This normalizes the correlation values
    z_values = np.arctanh(correlations_clipped)
    # 3. Perform the one-sample t-test
    # H0: The mean of the z-transformed values is 0.
    # Ha: The mean is not 0.
    t_statistic, p_value = stats.ttest_1samp(z_values, popmean=0)
    # 4. Interpret the results
    mean_correlation = np.mean(correlations)
    return t_statistic, p_value, mean_correlation

# Define max distance
max_dist = 200

# Pick a gene of interest
gene_of_interest = 'SPP1'

# Set result directory
hev_dir = os.path.join(res_dir, '05_xenium_mechanisms', '06_cellular_enrichment', f'{gene_of_interest}_vs_cell_abundance')
hev_dir = os.path.join(hev_dir, f'maxdist_{max_dist}um')

most_recent = hev_dir

# Load all samples & merge
sample_res = {}
for sample in sample_list:
    hev_file = f'{sample}_HEV_{gene_of_interest}_and_cell_proportions_corr.csv'
    tmp_df = pd.read_csv(os.path.join(most_recent, 'csv_files', hev_file))
    tmp_df['sample_id'] = sample
    sample_res[sample] = tmp_df

spp1_corr = pd.concat(sample_res.values(), ignore_index=True)

# Sort by median correlation
sorted_cell_types = spp1_corr.groupby('cell_type')['pearson_r'].median().sort_values(ascending=False).index.tolist()
spp1_corr['cell_type'] = pd.Categorical(spp1_corr['cell_type'], categories=sorted_cell_types, ordered=True)


# --- Statistical analysis

stat_df = pd.DataFrame(columns=['cell_type', 't_statistic', 'p_value', 'mean_correlation', 'median_correlation'])

for cell_type in sorted_cell_types:
    spp1_correlations = spp1_corr[spp1_corr['cell_type'] == cell_type]['pearson_r'].values
    t_statistic, p_value, mean_correlation = analyze_cell_type_correlation(spp1_correlations, cell_type)
    median_correlation = np.median(spp1_correlations)
    stat_df = stat_df._append({'cell_type': cell_type, 't_statistic': t_statistic, 'p_value': p_value, 'mean_correlation': mean_correlation, 'median_correlation': median_correlation}, ignore_index=True)

stat_df['significant'] = stat_df['p_value'] < alpha


# --- Plot 1a: Plot correlations as boxplot

whisker_props = {'color': 'black'}
median_props = {'color': 'black'}

fig, ax = plt.subplots(figsize=(3, 2))

sns.boxplot(data=spp1_corr, x='cell_type', y='pearson_r', showfliers=False, linewidth=1, whiskerprops=whisker_props, 
                    medianprops=median_props, width=0.8, palette=color_dict, ax=ax, saturation=1)

sns.stripplot(data=spp1_corr, x='cell_type', y='pearson_r', color='black', jitter=True, size=3,
                    marker='o', ax=ax)

# Automatically resize plot
n_ticks = len(ax.get_xticks())
ax_pos = ax.get_position()
ax.set_position([ax_pos.x0, ax_pos.y0, 0.05*n_ticks, ax_pos.height])

# Modify x tick labels to use short names
ax.set_xticklabels(ax.get_xticklabels(), rotation=90,)
new_labels = [miniatlas_sub_dict['miniatlas_subset_short'].get(label.get_text(), label.get_text()) for label in ax.get_xticklabels()]
ax.set_xticklabels(new_labels)

# Adjust boxplot colors
plotting_functions.beautify_plot(ax)

# Add thicker line at 0
ax.axhline(0, color='lightgray', linewidth=2, zorder=0)

# Add asterisks for significant cell types based on p-value
# Just above the top error bar of each boxplot
ylim = ax.get_ylim()
range_y = ylim[1] - ylim[0]

offset = range_y * 0.003  # 0.5% of the y-axis range

for i, cell_type in enumerate(sorted_cell_types):
    cell_stats = stat_df[stat_df['cell_type'] == cell_type]
    if not cell_stats.empty:
        p_val = cell_stats['p_value'].values[0]
        # Determine significance level
        if p_val < 0.001:
            sig_marker = '***'
        elif p_val < 0.01:
            sig_marker = '**'
        elif p_val < 0.05:
            sig_marker = '*'
        else:
            sig_marker = None
        
        if sig_marker:
            # Get the y position of the top whisker
            cell_data = spp1_corr[spp1_corr['cell_type'] == cell_type]['pearson_r']
            whisker_top = cell_data.max()
            # Place the asterisk slightly above the whisker
            ax.text(i, whisker_top + offset, sig_marker, ha='center', va='bottom', color='black', fontsize=9, fontweight='bold')

# Increase ymax a bit to make space for asterisks
ax.set_ylim(top=ax.get_ylim()[1]*1.1)

# Add legend
labels = major_cell_types_sorted
handles = [
    mpl.patches.Patch(facecolor=color_dict[label], label=label, edgecolor='black', linewidth=0.8)
    for label in labels if label in color_dict
]
legend = ax.legend(handles=handles, bbox_to_anchor=(1, 0.5), loc='center left', frameon=False, handlelength=0.8, fontsize=9)

for legpatch in legend.get_patches():
    fill_color = legpatch.get_facecolor()
    fill_color = list(fill_color)
    darkened_color = fill_color.copy()
    for i in range(len(darkened_color)-1):
        darkened_color[i] *= 0.8
    fill_color[-1] = 0.8 # make the fill color a bit transparent
    legpatch.set_color(tuple(fill_color))
    legpatch.set_edgecolor(tuple(darkened_color))
    legpatch.set_linewidth(1)

ax.set_ylabel('Pearson correlation', fontweight='bold')
ax.set_xlabel('Cell type', fontweight='bold')
ax.set_title(f'HEV {gene_of_interest} expression\nvs\ncell type proportion (<{max_dist} µm from HEV)', fontweight='bold')
ax.yaxis.grid(True, linestyle='-', which='major', color='lightgrey')
ax.xaxis.grid(True, linestyle='-', which='major', color='lightgrey')
ax.set_axisbelow(True)
title_x = ax.get_position().x0 + ax.get_position().width / 2
title_y = ax.get_position().y1 + 0.4
plt.suptitle(f'All samples - Correlation between HEV {gene_of_interest} expression (ACKR1+ endothelial) and cell type abundance (within {max_dist} µm of HEV)', x=title_x, y=title_y, fontweight='bold')


# --- Save output

out_dir = os.path.join(res_dir, 'figures')
os.makedirs(out_dir, exist_ok=True)

filename = 'Fig_5H.pdf'
plt.savefig(os.path.join(out_dir, filename), bbox_inches='tight')


print("Script successfully completed")
