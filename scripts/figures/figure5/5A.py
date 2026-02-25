#!/usr/bin/env python3

# --- Aim
# Plot Figure 5A: CCL19 and CCL21 chemokine expression in HEVs vs distance to TLS
# Shows spatial gradients of T-cell recruiting chemokines relative to TLS structures

# --- Define environment
# env: xen-hev

# --- Load packages
import argparse
import os
from datetime import datetime
import sys
import anndata as ad

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

import pandas as pd
import numpy as np

import seaborn as sns
from scipy import stats
from statsmodels.stats.multitest import multipletests

mpl.rcParams['pdf.fonttype'] = 42


# --- Argparse arguments

parser = argparse.ArgumentParser(description="Import Visium unfiltered data & plot",
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


# --- Load the anndata objects

adata_dict = {}

for sample in sample_list:
    adata = ad.read_h5ad(os.path.join(h5ad_dir, 'xenium', f'{sample}.h5ad'))
    adata_dict[sample] = adata


# --- Load the HEV distance to TLS calculations

# Find most recent results
def find_most_recent(dir):
    # Find most recent results
    subdirs = [d for d in os.listdir(dir) if os.path.isdir(os.path.join(dir, d))]
    # Filter subdirectories to keep only those with date format YYYYmmdd
    date_format = '%Y%m%d'
    valid_subdirs = []
    for subdir in subdirs:
        try:
            datetime.strptime(subdir, date_format)
            valid_subdirs.append(subdir)
        except ValueError:
            pass
    # Sort the valid subdirectories by date and get the most recent one
    if valid_subdirs:
        most_recent_dir = max(valid_subdirs, key=lambda x: datetime.strptime(x, date_format))
        most_recent_dir_path = os.path.join(dir, most_recent_dir)
        print("Most recent results:", most_recent_dir_path)
        return most_recent_dir_path
    else:
        print("No valid directories found.")

distance_dir = os.path.join(res_dir, '05_xenium_mechanisms', '02_identify_HEV_regions')
most_recent_distance = find_most_recent(distance_dir)

distance_dict = {}

for sample in sample_list:
    distance_file = f'{sample}_HEV_niches_info.csv'
    distance_df = pd.read_csv(os.path.join(most_recent_distance, 'csv_files', distance_file))
    # Prepend sample name to niche_id
    distance_df['niche_id'] = sample + '_' + distance_df['niche_id'].astype(str)
    # Make a mapping dict of distance to TLS
    distance_dict[sample] = distance_df.set_index('niche_id')['min_distance_to_TLS'].to_dict()
    # Also load HEV niche barcode and add to adata
    filename = f'{sample}_HEV_niches_cell_to_niche.csv'
    niche_barcode_df = pd.read_csv(os.path.join(most_recent_distance, 'csv_files', filename), index_col=0)
    niche_barcode_df.dropna(subset=['niche_id'], inplace=True)
    # Prepend sample name to niche_id
    niche_barcode_df['niche_id'] = sample + '_' + niche_barcode_df['niche_id'].astype(str)
    adata_dict[sample].obs['niche_id'] = niche_barcode_df['niche_id']

# Extract an unnested distance mapping dict
unnested_distance_dict = {}
for sample in sample_list:
    unnested_distance_dict.update(distance_dict[sample])


# --- Also load cell annotations

# Decide which annotation levels to use
level = 'subset'

for sample in sample_list:
    annotation_file = f'{sample}_{level}_pred_cell_types.csv'
    annotation_dir = os.path.join(project_dir, 'results', '02_cell_type_mapping', 'label_transfer', 'xenium')
    annotations = pd.read_csv(os.path.join(annotation_dir, f'{level}_level', annotation_file), index_col=0)
    # Add to adata
    adata_dict[sample].obs[level] = annotations['pred_cell_type']


# --- For the genes of interest, calculate average expression per HEV niche across all samples

genes_of_interest = ['CCL19', 'CCL21']

# Subset adatas to HEV niches only
for sample in sample_list:
    adata = adata_dict[sample]
    adata_hev = adata[adata.obs['niche_id'].notna()]
    adata_dict[sample] = adata_hev

# Extract a single adata object with all samples
adata_list = list(adata_dict.values())
all_adata = ad.concat(adata_list, merge='same')

average_expression_df = pd.DataFrame()
average_expression_df['niche_id'] = all_adata.obs['niche_id'].unique()
average_expression_df.set_index('niche_id', inplace=True)
average_expression_df['min_distance_to_TLS'] = average_expression_df.index.map(unnested_distance_dict)

for gene_of_interest in genes_of_interest:
    print(f'Processing gene: {gene_of_interest}')
    # Calculate average expression per HEV niche
    avg_expr = pd.DataFrame()
    avg_expr[f'avg_{gene_of_interest}_expression_all_cells'] = all_adata.to_df()[[gene_of_interest]].groupby(all_adata.obs['niche_id']).mean()
    # Add to df
    average_expression_df = average_expression_df.merge(avg_expr, left_index=True, right_index=True)


# --- Bin HEVs into distance bins

# Define bins and labels
bins = [-0.1, 0, 100, 200, np.inf]

labels = []
for i in range(len(bins) - 1):
    left = bins[i]
    right = bins[i + 1]
    if right == np.inf:
        labels.append(f'>{int(left)}')
    elif left == -0.1 and right == 0:
        labels.append('0')
    else:
        labels.append(f'{int(left)}-{int(right)}')

# Now use these labels in pd.cut
average_expression_df['TLS_distance_bin'] = pd.cut(
    average_expression_df['min_distance_to_TLS'],
    bins=bins,
    labels=labels,
    include_lowest=True,
    right=True
)

# Count how many hevs per bin
hev_counts = average_expression_df['TLS_distance_bin'].value_counts().sort_index()


# --- Load the HEV colors

color_dir = os.path.join(project_dir, 'resources', 'plotting_colors')
filename = 'TLS_HEV_colors.csv'

hev_colors_df = pd.read_csv(os.path.join(color_dir, filename), index_col=0)
color_dict = hev_colors_df['color'].to_dict()

# Take the 2 most extreme colors for the distance bins
n_colors = len(labels)
color_list = [color_dict['near_TLS_HEV'], color_dict['non_TLS_HEV']]

cmap = mpl.colors.LinearSegmentedColormap.from_list('HEV_distance_cmap', colors=color_list, N=n_colors)

# Extract palette
color_dict_bins = {label: cmap(i / (n_colors - 1)) for i, label in enumerate(labels)}


# --- Calculate stats for all genes and cell types

stats_results = {}

for gene_of_interest in genes_of_interest:
    stats_results[gene_of_interest] = {}
    
    # Define the columns to analyze for this gene
    columns_to_analyze = [
        f'avg_{gene_of_interest}_expression_all_cells',
    ]
    
    for col in columns_to_analyze:
        # 1. Spearman Correlation (Trend Test)
        corr_df = average_expression_df.dropna(subset=['min_distance_to_TLS', col])
        r_val, p_val_corr = stats.spearmanr(corr_df['min_distance_to_TLS'], corr_df[col])
        
        # 2. Kruskal-Wallis (Group Difference Test)
        groups = []
        for label in labels:
            values = average_expression_df[average_expression_df['TLS_distance_bin'] == label][col].values
            values = values[~np.isnan(values)]
            groups.append(values)
        
        stat, p_val_kw = stats.kruskal(*groups)
        
        # Store results
        stats_results[gene_of_interest][col] = {
            'r': r_val,
            'p_corr': p_val_corr,
            'kw_stat': stat,
            'p_kw': p_val_kw
        }
        
        print(f"\n{gene_of_interest} - {col}:")
        print(f"  Spearman Correlation: r={r_val:.3f}, p={p_val_corr:.2e}")
        print(f"  Kruskal-Wallis: H={stat:.3f}, p={p_val_kw:.2e}")


# --- Post-hoc tests: Mann-Whitney U between groups

posthoc_mode = 'neighbor'  # Options: 'baseline' or 'neighbor'

# Define comparison pairs based on mode
if posthoc_mode == 'baseline':
    # All bins vs first bin
    comparison_pairs = [
        ('0', '0-100'),
        ('0', '100-200'),
        ('0', '>200')
    ]
elif posthoc_mode == 'neighbor':
    # Neighbor vs neighbor
    comparison_pairs = [
        ('0', '0-100'),
        ('0-100', '100-200'),
        ('100-200', '>200')
    ]

# Store results for each gene and column
posthoc_results = {}

for gene_of_interest in genes_of_interest:
    posthoc_results[gene_of_interest] = {}
    
    columns_to_analyze = [
        f'avg_{gene_of_interest}_expression_all_cells',
    ]
    
    for col in columns_to_analyze:
        p_values = []
        
        # 1. Run the pairwise tests
        for group1, group2 in comparison_pairs:
            # Get the data for the two groups
            data1 = average_expression_df[average_expression_df['TLS_distance_bin'] == group1][col].dropna()
            data2 = average_expression_df[average_expression_df['TLS_distance_bin'] == group2][col].dropna()
            
            # Run Mann-Whitney U
            stat, p = stats.mannwhitneyu(data1, data2)
            p_values.append(p)
        
        # 2. Correct for multiple testing (Bonferroni)
        reject, p_adjusted, _, _ = multipletests(p_values, alpha=0.05, method='bonferroni')
        
        # Store results
        posthoc_results[gene_of_interest][col] = {
            'comparison_pairs': comparison_pairs,
            'p_adjusted': p_adjusted,
            'reject': reject
        }
        
        # 3. Print results
        print(f"\n{gene_of_interest} - {col}:")
        for (group1, group2), p_adj, is_sig in zip(comparison_pairs, p_adjusted, reject):
            sig_str = "*" if is_sig else "ns"
            print(f"  {group1} vs {group2}: p_adj = {p_adj:.2e} ({sig_str})")


# --- Functions for plotting

def add_stat_annotation(ax, x1, x2, h, text, fontsize=10):
    """
    Draws a simple horizontal line between x1 and x2 at height h, with text on top.
    """
    line_h = h
    text_h = h
    # Draw the line (just x1 to x2 at height line_h)
    ax.plot([x1, x2], [line_h, line_h], lw=1.0, c='black')
    # Add the text
    ax.text((x1 + x2) * 0.5, text_h, text, ha='center', va='bottom', 
            color='black', fontsize=fontsize)


# --- Plot 1: Make a boxplot of genes of interest expression vs distance to TLS

# Decide on stat display mode
stat_display_mode = 'brackets'  # Options: 'text', 'brackets', 'hybrid'

# Add info whether HEV is close to TLS or not
nrows = len(genes_of_interest)
ncols = 1

# Properties for boxplot
whisker_props = {'color': 'black'}
median_props = {'color': 'black'}

# Set fontsize
plt.rcParams.update({'font.size': 10})

# Set grid behind
plt.rcParams['axes.axisbelow'] = True

fig, axs = plt.subplots(nrows=nrows, ncols=ncols, figsize=(3, 3))

if ncols == 1:
    axs = axs[:, np.newaxis]  # Add an extra dimension for consistent indexing

for n, gene_of_interest in enumerate(genes_of_interest):
    print(f'Plotting gene: {gene_of_interest}')
    # Extract the average expression df for this gene
    spp1_avg = average_expression_df[[f'avg_{gene_of_interest}_expression_all_cells']].copy()
    spp1_avg['TLS_distance_bin'] = average_expression_df['TLS_distance_bin']
    spp1_avg['niche_id'] = average_expression_df.index
    # Plot all cells
    sns.stripplot(data=spp1_avg, x='TLS_distance_bin', y=f'avg_{gene_of_interest}_expression_all_cells', color='black', jitter=True, size=2,
                    marker='o', ax=axs[n][0],)
    sns.boxplot(data=spp1_avg, x='TLS_distance_bin', y=f'avg_{gene_of_interest}_expression_all_cells', ax=axs[n][0], showfliers=False, linewidth=1, linecolor='black',
                    whiskerprops=whisker_props, medianprops=median_props, saturation=1, zorder=10, dodge=False,
                    legend=False, hue='TLS_distance_bin', palette=color_dict_bins)
    axs[n][0].set_title(f'HEV expression\n(all cells)', fontweight='bold', fontsize=10)

    # Add stats to each subplot for this gene based on display mode
    columns_to_analyze = [
        f'avg_{gene_of_interest}_expression_all_cells',
    ]
    
    if stat_display_mode in ['text', 'hybrid']:
        for i, col in enumerate(columns_to_analyze):
            r_val = stats_results[gene_of_interest][col]['r']
            p_val_kw = stats_results[gene_of_interest][col]['p_kw']
            
            # Format KW pvalue
            if p_val_kw < 0.001:
                pval_text = 'p < 0.001'
            else:
                pval_text = f'p = {p_val_kw:.3f}'
            
            axs[n][i].text(0.05, 0.95, 
                    f'r = {r_val:.2f}\nKW {pval_text}', 
                    transform=axs[n][i].transAxes, 
                    ha='left', va='top', fontsize=8,
                    bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor='none', alpha=0.7))
    
    # Define index mapping for bins
    bin_indices = {'0': 0, '0-100': 1, '100-200': 2, '>200': 3}
    
    # Draw brackets if in bracket or hybrid mode
    if stat_display_mode in ['brackets', 'hybrid']:
        for i, col in enumerate(columns_to_analyze):
            ax = axs[n][i]
            
            # Get the post-hoc results for this gene and column
            comparison_pairs_data = posthoc_results[gene_of_interest][col]['comparison_pairs']
            p_adjusted = posthoc_results[gene_of_interest][col]['p_adjusted']
            reject = posthoc_results[gene_of_interest][col]['reject']
            
            # Determine max y for bracket placement
            if stat_display_mode == 'hybrid':
                # In hybrid mode, only consider the first comparison pair
                group1, group2 = comparison_pairs_data[0]
                data1 = average_expression_df[average_expression_df['TLS_distance_bin'] == group1][col].dropna()
                data2 = average_expression_df[average_expression_df['TLS_distance_bin'] == group2][col].dropna()
                y_max = max(data1.max(), data2.max())
            else:
                # In brackets mode, consider all data
                y_max = average_expression_df[col].max()
            h_start = y_max * 1.05
            h_step = y_max * 0.15
            
            # For hybrid mode, only draw the first comparison
            comparisons_to_draw = 1 if stat_display_mode == 'hybrid' else len(comparison_pairs_data)
            
            # Draw brackets
            bracket_count = 0
            for idx, ((group1, group2), p_adj, is_sig) in enumerate(zip(comparison_pairs_data, p_adjusted, reject)):
                if idx >= comparisons_to_draw:
                    break
                # Get x-coordinates
                x1 = bin_indices[group1]
                x2 = bin_indices[group2]
                # Calculate height for this bracket
                h = h_start + (bracket_count * h_step)
                # Define text (*, **, ***, or ns)
                if p_adj < 0.001: sig_text = "***"
                elif p_adj < 0.01: sig_text = "**"
                elif p_adj < 0.05: sig_text = "*"
                else: sig_text = "ns"
                # Draw it on the specific subplot
                add_stat_annotation(ax, x1, x2, h, sig_text, fontsize=8)
                bracket_count += 1
            
            # Adjust ylim to fit brackets
            if bracket_count > 0:
                # In hybrid mode, need extra space for text annotation at top
                ylim_multiplier = 1.6 if stat_display_mode == 'hybrid' else 1.1
                ax.set_ylim(top=h_start + (bracket_count * h_step) * ylim_multiplier)
            
            # Prevent top ytick from being too close to ylim (prevents overlap with count_ax)
            ax.yaxis.set_major_locator(MaxNLocator(nbins='auto', prune='upper'))
    
    for ax in axs[n]:
        plotting_functions.beautify_plot(ax)
        ax.set_ylabel(f'{gene_of_interest}\n(Average)')
        # Automatically resize plots based on number of x-ticks
        ax_pos = ax.get_position()
        ax.set_position([ax_pos.x0, ax_pos.y0, 0.1 * len(ax.get_xticks()), ax_pos.height])
        # Rotate xlabels
        ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
        # Set xlabels
        ax.set_xlabel('Distance to TLS (µm)')
    # Remove xticks if not the last row
    if n < nrows - 1:
        for ax in axs[n]:
            ax.set_xticks([])

for ax in axs.flat:
    # Resize the plots
    ax_pos = ax.get_position()
    ax.set_position([ax_pos.x0, ax_pos.y0, 0.4, 0.4])
    # Add grid
    ax.grid(which='major', axis='y', linestyle='-', color='lightgrey')

plotting_functions.move_plots_x_y(axs, x_offset=0.4, y_offset=0.02, remove_ylabels=False, remove_xlabels=True, remove_title=True)

# Add ax with count bars on top of each ax
for ax in axs[0]:
    count_ax = plt.gcf().add_axes([ax.get_position().x0, ax.get_position().y1 + 0.02, ax.get_position().width, ax.get_position().height * 0.3])
    sns.barplot(x=hev_counts.index, y=hev_counts.values, color='lightgrey', edgecolor='black', ax=count_ax, hue=hev_counts.index, palette=color_dict_bins,
                saturation=1)
    # Annotate counts
    for n, p in enumerate(count_ax.patches):
        count_ax.annotate(f'{int(p.get_height())}', (p.get_x() + p.get_width() / 2., p.get_height()), ha='center', va='bottom', fontsize=8)
    # Adjust ylim to fit counts
    count_ax.set_ylim(0, hev_counts.values.max() * 1.4)
    # Add title from main plot to count_ax instead
    count_ax.set_title(ax.get_title(), fontweight='bold')
    ax.set_title('')  # Remove title from main ax
    # Remove xticks
    count_ax.set_xticks([])
    # Remove top and right spines
    count_ax.spines['top'].set_visible(False)
    count_ax.spines['right'].set_visible(False)
    # Adjust labels
    count_ax.set_xlabel('')
    count_ax.set_ylabel('Count')
    # Beautify
    plotting_functions.beautify_plot(count_ax)

title_x = (axs[0][0].get_position().x0 + axs[0][-1].get_position().x1) / 2
title_y = count_ax.get_position().y1 + 0.4

plt.suptitle(f'All samples - Average {gene_of_interest} Expression (average per HEV niche) by HEV Proximity to TLS', fontweight='bold', y=title_y, x=title_x)


# --- Save the output

out_dir = os.path.join(res_dir, 'figures')
os.makedirs(out_dir, exist_ok=True)

filename = 'Fig_5A.pdf'
plt.savefig(os.path.join(out_dir, filename), bbox_inches='tight')


# --- Finish

print("Script successfully completed")