#!/usr/bin/env python3

# --- Aim
# Plot Fig. 5f: gene expression in HEVs correlated with distance to TLS
# Shows spatial expression gradients of SPP1 as function of proximity to TLS

# --- Define environment
# env: xen-hev

# --- Load packages
import argparse
import os
import sys
import anndata as ad

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from matplotlib.colors import TwoSlopeNorm
from matplotlib.backends.backend_pdf import PdfPages 

import pandas as pd
import numpy as np

import seaborn as sns
from scipy import stats
from statsmodels.stats.multitest import multipletests

mpl.rcParams['pdf.fonttype'] = 42


# --- Argparse arguments

parser = argparse.ArgumentParser(description="Fig. 5f: HEV gene expression vs distance to TLS",
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


# --- Extract some info from the sample sheet

# Use sample_sheet argument if provided, else default
if args.get("sample_sheet"):
    sample_sheet_path = args["sample_sheet"]
else:
    sample_dir = os.path.join(project_dir, 'config')
    sample_sheet_file = 'sample_sheet_xenium.csv'
    sample_sheet_path = os.path.join(sample_dir, sample_sheet_file)

sample_sheet = pd.read_csv(sample_sheet_path)

# Also load tissue metadata
metadata_dir = os.path.join(project_dir, 'resources',)
tissue_metadata_file = 'tissue_metadata.csv'
tissue_metadata = pd.read_csv(os.path.join(metadata_dir, tissue_metadata_file))

# Make a section_id to subtype mapping
section_to_subtype = tissue_metadata.set_index('tissue_id')['clinical_subtype'].to_dict()
# Map section_id to subtype
sample_sheet['subtype'] = sample_sheet['tissue_id'].map(section_to_subtype)

subtype_dict = sample_sheet.set_index('section_id')['subtype'].to_dict()


# --- Load the anndata objects

adata_dict = {}

for sample in sample_list:
    adata = ad.read_h5ad(os.path.join(h5ad_dir, 'xenium', f'{sample}.h5ad'))
    adata_dict[sample] = adata


# --- Load the HEV distance to TLS calculations

distance_dir = os.path.join(res_dir, '05_xenium_mechanisms', '02_identify_HEV_regions')
most_recent_distance = distance_dir

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
    annotation_dir = os.path.join(res_dir, '02_cell_type_mapping', 'label_transfer', 'xenium')
    annotations = pd.read_csv(os.path.join(annotation_dir, f'{level}_level', annotation_file), index_col=0)
    # Add to adata
    adata_dict[sample].obs[level] = annotations['pred_cell_type']


# --- For the gene of interest, calculate average expression per HEV niche across all samples

gene_of_interest = 'SPP1'

# Decide on stat display mode and posthoc mode
stat_display_mode = 'hybrid'  # 'text', 'brackets', 'hybrid'
posthoc_mode = 'baseline'  # 'baseline' (all vs first bin) or 'neighbor' (neighbor vs neighbor)

# Subset adatas to HEV niches only
for sample in sample_list:
    adata = adata_dict[sample]
    adata_hev = adata[adata.obs['niche_id'].notna()]
    adata_dict[sample] = adata_hev

# Extract a single adata object with all samples
adata_list = list(adata_dict.values())
all_adata = ad.concat(adata_list, merge='same')

# Calculate average expression per HEV niche
spp1_avg = pd.DataFrame()

# Calculate average expression in ACKR1+ endothelial cells only
adata_sub = all_adata[all_adata.obs[level] == 'Endothelial ACKR1'].copy()
spp1_avg[f'avg_{gene_of_interest}_expression_ACKR1'] = adata_sub.to_df()[[gene_of_interest]].groupby(adata_sub.obs['niche_id']).mean()

# Add distance to TLS info
spp1_avg['min_distance_to_TLS'] = spp1_avg.index.map(unnested_distance_dict)

# Add niche_id as a column
spp1_avg['niche_id'] = spp1_avg.index


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
spp1_avg['TLS_distance_bin'] = pd.cut(
    spp1_avg['min_distance_to_TLS'],
    bins=bins,
    labels=labels,
    include_lowest=True,
    right=True
)

# Count how many hevs per bin
hev_counts = spp1_avg['TLS_distance_bin'].value_counts().sort_index()

# Count bars: only HEVs with ACKR1+ endothelial cells have a value in the box plot and enter the heatmap means
# (matches the number of points)
hev_counts_plot = (spp1_avg.loc[spp1_avg[f'avg_{gene_of_interest}_expression_ACKR1'].notna(), 'TLS_distance_bin']
                   .value_counts().reindex(hev_counts.index, fill_value=0))


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


# --- Calculate stats for all three groups

# Define the columns to analyze
columns_to_analyze = [
    f'avg_{gene_of_interest}_expression_ACKR1',
]

stats_results = {}

for col in columns_to_analyze:
    # 1. Spearman Correlation (Trend Test)
    corr_df = spp1_avg.dropna(subset=['min_distance_to_TLS', col])
    r_val, p_val_corr = stats.spearmanr(corr_df['min_distance_to_TLS'], corr_df[col])
    
    # 2. Kruskal-Wallis (Group Difference Test)
    groups = []
    for label in labels:
        values = spp1_avg[spp1_avg['TLS_distance_bin'] == label][col].values
        values = values[~np.isnan(values)]
        groups.append(values)
    
    stat, p_val_kw = stats.kruskal(*groups)
    
    # Store results
    stats_results[col] = {
        'r': r_val,
        'p_corr': p_val_corr,
        'kw_stat': stat,
        'p_kw': p_val_kw
    }
    
    print(f"\n{col}:")
    print(f"  Spearman Correlation: r={r_val:.3f}, p={p_val_corr:.2e}")
    print(f"  Kruskal-Wallis: H={stat:.3f}, p={p_val_kw:.2e}")


# --- Post-hoc tests: Mann-Whitney U between groups

# Define comparison pairs based on mode
if posthoc_mode == 'baseline':
    # All bins vs first bin
    baseline = '0'
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

# Store results for each column
posthoc_results = {}

for col in columns_to_analyze:
    p_values = []
    
    # 1. Run the pairwise tests
    for group1, group2 in comparison_pairs:
        # Get the data for the two groups
        data1 = spp1_avg[spp1_avg['TLS_distance_bin'] == group1][col].dropna()
        data2 = spp1_avg[spp1_avg['TLS_distance_bin'] == group2][col].dropna()
        
        # Run Mann-Whitney U
        stat, p = stats.mannwhitneyu(data1, data2)
        p_values.append(p)
    
    # 2. Correct for multiple testing (Bonferroni)
    # 'reject' tells you if it's significant (True/False)
    reject, p_adjusted, _, _ = multipletests(p_values, alpha=0.05, method='bonferroni')
    
    # Store results
    posthoc_results[col] = {
        'comparison_pairs': comparison_pairs,
        'p_adjusted': p_adjusted,
        'reject': reject
    }
    
    # 3. Print results
    print(f"\n{col}:")
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


# --- Plot 1: Make a boxplot of SPP1 expression vs distance to TLS

# Add info whether HEV is close to TLS or not
ncols = 1

# Properties for boxplot
whisker_props = {'color': 'black'}
median_props = {'color': 'black'}

# Set fontsize
plt.rcParams.update({'font.size': 10})

# Set grid behind
plt.rcParams['axes.axisbelow'] = True

fig, axs = plt.subplots(nrows=1, ncols=ncols, figsize=(3, 1.8))

if ncols == 1:
    axs = [axs]

# First plot ACKR1 only
sns.stripplot(data=spp1_avg, x='TLS_distance_bin', y=f'avg_{gene_of_interest}_expression_ACKR1', color='black', jitter=True, size=2,
                marker='o', ax=axs[0],)

sns.boxplot(data=spp1_avg, x='TLS_distance_bin', y=f'avg_{gene_of_interest}_expression_ACKR1', ax=axs[0], showfliers=False, linewidth=1, linecolor='black',
            whiskerprops=whisker_props, medianprops=median_props, saturation=1, zorder=10, dodge=False,
            legend=False, hue='TLS_distance_bin', palette=color_dict_bins)

axs[0].set_title(f'HEV {gene_of_interest}\n(ACKR1+ endothelial)', fontweight='bold', fontsize=10)

ymax = 0

for ax in axs:
    plotting_functions.beautify_plot(ax)
    ax.set_ylabel(f'{gene_of_interest} expression\n(Average)')
    ymax = max(ymax, ax.get_ylim()[1])
    # Automatically resize plots based on number of x-ticks
    ax_pos = ax.get_position()
    ax.set_position([ax_pos.x0, ax_pos.y0, 0.1 * len(ax.get_xticks()), ax_pos.height])
    # Rotate xlabels
    ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
    # Set xlabels
    ax.set_xlabel('Distance to TLS (µm)')

plotting_functions.move_plots_x_y([axs], x_offset=0.4, y_offset=0.0, remove_ylabels=False, remove_xlabels=False, remove_title=False)

# Set same ylims
for ax in axs:
    ax.set_ylim(0, ymax)
    # Add grid
    ax.grid(which='major', axis='y', linestyle='-', color='lightgrey')

# Add ax with count bars on top of each ax
for ax in axs:
    count_ax = plt.gcf().add_axes([ax.get_position().x0, ax.get_position().y1 + 0.02, ax.get_position().width, ax.get_position().height * 0.3])
    sns.barplot(x=hev_counts_plot.index, y=hev_counts_plot.values, color='lightgrey', edgecolor='black', ax=count_ax, hue=hev_counts_plot.index, palette=color_dict_bins,
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

# Add stats to each subplot based on display mode
if stat_display_mode in ['text', 'hybrid']:
    for i, (ax, col) in enumerate(zip(axs, columns_to_analyze)):
        r_val = stats_results[col]['r']
        p_val_kw = stats_results[col]['p_kw']
        
        # Format KW pvalue
        if p_val_kw < 0.001:
            pval_text = 'p < 0.001'
        elif p_val_kw > 0.01:
            pval_text = f'p = {p_val_kw:.2f}'
        else:
            pval_text = f'p = {p_val_kw:.3f}'
        
        ax.text(0.05, 0.95, 
                f'Sp. r = {r_val:.2f}\nKW {pval_text}', 
                transform=ax.transAxes, 
                ha='left', va='top', fontsize=8,
                bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor='none', alpha=0.7))

# Define index mapping for bins
bin_indices = {'0': 0, '0-100': 1, '100-200': 2, '>200': 3}

# Draw brackets if in bracket or hybrid mode
if stat_display_mode in ['brackets', 'hybrid']:
    # Draw brackets for all three subplots
    for i, (ax, col) in enumerate(zip(axs, columns_to_analyze)):
        # Get the post-hoc results for this column
        comparison_pairs = posthoc_results[col]['comparison_pairs']
        p_adjusted = posthoc_results[col]['p_adjusted']
        reject = posthoc_results[col]['reject']
        
        # Get global max for calculating h_step (relative spacing between brackets)
        global_y_max = spp1_avg[col].max()
        h_step = global_y_max * 0.15    # Step up for each subsequent bracket
        
        # For hybrid mode, only draw the first comparison
        comparisons_to_draw = 1 if stat_display_mode == 'hybrid' else len(comparison_pairs)
        
        # Track the current maximum height (for stacking brackets)
        current_max_h = 0
        bracket_count = 0
        
        # Draw brackets ONLY for significant comparisons
        for idx, ((group1, group2), p_adj, is_sig) in enumerate(zip(comparison_pairs, p_adjusted, reject)):
            if idx >= comparisons_to_draw:
                break
            
            # Get x-coordinates
            x1 = bin_indices[group1]
            x2 = bin_indices[group2]
            
            # Calculate the local max for just the groups being compared
            # This allows lower brackets when the max value is elsewhere
            data1 = spp1_avg[spp1_avg['TLS_distance_bin'] == group1][col].dropna()
            data2 = spp1_avg[spp1_avg['TLS_distance_bin'] == group2][col].dropna()
            local_y_max = max(data1.max(), data2.max())
            
            # Start this bracket just above the local max
            h_local_start = local_y_max * 1.05
            
            # But also make sure it's above any previous brackets (to avoid overlap)
            h = max(h_local_start, current_max_h + h_step)
            
            # Update the current max height for next bracket
            current_max_h = h
            
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
            ylim_multiplier = 1.5 if stat_display_mode == 'hybrid' else 1.1
            ax.set_ylim(top=current_max_h * ylim_multiplier)
        
        # Prevent top ytick from being too close to ylim (prevents overlap with count_ax)
        ax.yaxis.set_major_locator(MaxNLocator(nbins='auto', prune='upper'))

title_x = (axs[0].get_position().x0 + axs[-1].get_position().x1) / 2
title_y = count_ax.get_position().y1 + 0.4

plt.suptitle(f'All samples - Average {gene_of_interest} Expression (average per HEV niche) by HEV Proximity to TLS', fontweight='bold', y=title_y, x=title_x)


# --- Now also look at some other genes that we don't expect to change with distance to TLS

# Calculate mean & pct expressing per bin
summarize_df = spp1_avg.groupby('TLS_distance_bin').agg(
    mean_SPP1_ACKR1 = (f'avg_{gene_of_interest}_expression_ACKR1', 'mean'),
    count = ('niche_id', 'count'),
    pct_expressing_ACKR1 = (f'avg_{gene_of_interest}_expression_ACKR1', lambda x: (x > 0).sum() / x.count() * 100),
).reset_index()


# Calculate average gene expression per niche
# Split the genes into LRI genes and HEV markers
gene_dict = {
    'Markers': ['TSPAN7', 'RAMP3', 'MEOX2',],
    'Receptors': ['ACKR1', 'SELE', 'SELP', 'LTBR',],
}

LRI_genes = gene_dict['Receptors'].copy()
HEV_markers = gene_dict['Markers']

genes = [gene for sublist in gene_dict.values() for gene in sublist]

# If gene of interest not in genes, add it
if gene_of_interest not in genes:
    genes.append(gene_of_interest)
    gene_dict['Other'] = [gene_of_interest]
    LRI_genes.append(gene_of_interest) # Add it here too

# Also add SPP1 to other if not already there
if 'SPP1' not in genes:
    genes.append('SPP1')
    if 'Other' in gene_dict:
        gene_dict['Other'].append('SPP1')
    else:
        gene_dict['Other'] = ['SPP1']
    LRI_genes.append('SPP1')

# Reverse the order
genes = genes[::-1]

# Calculate average expression per HEV group
avg_exp_per_HEV = adata_sub.to_df()[genes].groupby(adata_sub.obs['niche_id']).mean()

# Add TLS distance bin info
avg_exp_per_HEV['TLS_distance_bin'] = spp1_avg.set_index('niche_id').loc[avg_exp_per_HEV.index, 'TLS_distance_bin']

# Check that there aren't any niches with NaN
if avg_exp_per_HEV['TLS_distance_bin'].isna().sum() > 0:
    print('Warning: There are niches with NaN TLS distance bin!')

# Now calculate mean & pct expressing per bin for each gene
summary_df = pd.DataFrame()

for gene in genes:
    gene_summary = avg_exp_per_HEV.groupby('TLS_distance_bin').agg(
        mean_expression = (gene, 'mean'),
        pct_expressing = (gene, lambda x: (x > 0).sum() / x.count() * 100),
    ).reset_index()
    # Add a gene column
    gene_summary['gene'] = gene
    summary_df = summary_df._append(gene_summary, ignore_index=True)

# Now also make a df with mean expression as a fraction of the mean expression in the 0 bin
fraction_df = pd.DataFrame()

for gene in genes:
    gene_data = summary_df[summary_df['gene'] == gene].copy()
    # Get mean expression in 0 bin
    mean_0 = gene_data[gene_data['TLS_distance_bin'] == '0']['mean_expression'].values[0]
    # Calculate log2 fold change vs 0
    gene_data['mean_expression_LFC_vs_0'] = np.log2(gene_data['mean_expression'] / mean_0 + 1e-10)  # add small value to avoid log(0)
    # Also calculate the fraction
    gene_data['mean_expression_fraction_vs_0'] = gene_data['mean_expression'] / mean_0
    fraction_df = fraction_df._append(gene_data, ignore_index=True)


# --- Plot 2: Heatmap of gene expression vs distance to TLS for the different genes

hue_key = 'mean_expression_LFC_vs_0'

vmin = fraction_df[hue_key].min()
vmax = fraction_df[hue_key].max()

# Calculate vmin & vmax centered around midpoint
midpoint = 0

# Find the largest distance from the midpoint (in either direction)
abs_max = max(abs(vmax - midpoint), abs(vmin - midpoint))

vmin = midpoint - abs_max
vmax = midpoint + abs_max

norm = TwoSlopeNorm(vmin=vmin, vcenter=midpoint, vmax=vmax)

fig, top_ax = plt.subplots(figsize=(2, 2))

axs = []

for n, (gene_type, gene_list) in enumerate(gene_dict.items()):
    print(f'{gene_type}: {gene_list}')

    if n == 0:
        ax = top_ax
    else:
        ax_width = ax_prev.get_position().width
        ax_height = ax_prev.get_position().height / len(ax_prev.get_yticks()) * len(gene_list)
        ax = plt.gcf().add_axes([ax_prev.get_position().x0, ax_prev.get_position().y0 - ax_height - 0.02, ax_width, ax_height])

    plot_df = fraction_df[fraction_df['gene'].isin(gene_list)].copy()

    g = sns.heatmap(data=plot_df.pivot(index='gene', columns='TLS_distance_bin', values=hue_key),
                    ax=ax, cmap='RdBu_r', norm=norm,
                    linewidths=1, linecolor='white', cbar=False)

    if n == 0:
        # Resize ax based on xticks & yticks
        ax_pos = ax.get_position()
        ax.set_position([ax_pos.x0, ax_pos.y0, 0.15 * len(ax.get_xticks()), 0.12 * len(ax.get_yticks())])

    # Set spines visible
    for _, spine in g.spines.items():
        spine.set_visible(True)

    # Add y label
    ax.set_ylabel(gene_type, fontweight='bold')

    # Rotate yticklabels
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0, ha='right')
    ax_prev = ax
    axs.append(ax)

# Remove xticks from all but the bottom ax
for ax in axs[:-1]:
    ax.set_xticks([])

# Add xlabel
axs[-1].set_xlabel('Distance to TLS (µm)', fontweight='bold')

# Add an ax for colorbar
cbar_ax = plt.gcf().add_axes([axs[-1].get_position().x1 + 0.08, axs[-1].get_position().y0, 0.12*4, 0.12*0.8])

# Add colorbar
cbar = plt.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap='RdBu_r'), cax=cbar_ax, orientation='horizontal')

# Add title
cbar_ax.set_title('Mean expression\n(LFC 0 µm bin)', fontsize=10)

# Set xlim to actual vmin & vmax
vmin = fraction_df[hue_key].min()
vmax = fraction_df[hue_key].max()
cbar.set_ticks([vmin, midpoint])
cbar.set_ticklabels([f'{vmin:.2f}', f'{midpoint:.2f}'])
cbar_ax.set_xlim(vmin, vmax)

# Add countax at the top
count_ax = plt.gcf().add_axes([top_ax.get_position().x0, top_ax.get_position().y1 + 0.02, top_ax.get_position().width, 0.12 * 2])
sns.barplot(x=hev_counts_plot.index, y=hev_counts_plot.values, color='lightgrey', edgecolor='black', ax=count_ax, hue=hev_counts_plot.index, palette=color_dict_bins,
               saturation=1)

# Beautify
plotting_functions.beautify_plot(count_ax)

# Annotate counts
for n, p in enumerate(count_ax.patches):
    count_ax.annotate(f'{int(p.get_height())}', (p.get_x() + p.get_width() / 2., p.get_height()), ha='center', va='bottom', fontsize=8)

# Adjust ylim
count_ax.set_ylim(top=hev_counts.values.max() * 1.4)

# Set top and right spine invisible
count_ax.spines['top'].set_visible(False)
count_ax.spines['right'].set_visible(False)

# Remove xticks
count_ax.set_xticks([])
count_ax.set_xlabel('')
count_ax.set_ylabel('Count')

# Add title
count_ax.set_title('HEV endothelial expression', fontweight='bold', fontsize=11,)

# Adjust main title position
title_x = (top_ax.get_position().x0 + top_ax.get_position().x1) / 2
title_y = count_ax.get_position().y1 + 0.4

plt.suptitle(f'All samples - HEV ACKR1+ Endothelial gene expression (average per HEV niche) by HEV proximity to TLS', fontweight='bold', y=title_y, x=title_x)


# --- Save the output

out_dir = os.path.join(res_dir, 'figures')
os.makedirs(out_dir, exist_ok=True)

filename = 'Fig_5F.pdf'

# Save all plots in single pdf
filename_dir=os.path.join(out_dir, filename)

pp = PdfPages(filename_dir)
fig_nums = plt.get_fignums()

figs = [plt.figure(n) for n in fig_nums]

for fig in figs:
    fig.savefig(pp, format='pdf', bbox_inches='tight')

pp.close()

plt.close('all')


# --- Finish

print("Script successfully completed")