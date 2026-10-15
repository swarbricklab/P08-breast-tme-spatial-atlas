#!/usr/bin/env python3

# --- Aim
# Examine the cellular composition of spots sending LTB signals to HEVs.
# Identify which cell types are involved in LTB-LTBR signaling in the vascular microenvironment.

# --- Define environment
# env: vis-hd

# --- Load packages
import argparse
import os
import re
import sys
from collections.abc import Iterable

from scipy import stats

import seaborn as sns

import anndata as ad
import matplotlib as mpl
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

sc.set_figure_params(facecolor="white", figsize=(8, 8))
sc.settings.verbosity = 3

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
sys.path.append(os.path.join(project_dir, 'scripts'))
from helper_functions import save_multi_image

print(f'Samples in python: {sample_list}')


# --- Extract some info from the sample sheet
# Use sample_sheet argument if provided, else default
if args.get("sample_sheet"):
    sample_sheet_path = args["sample_sheet"]
else:
    sample_sheet_path = os.path.join(project_dir, 'config', 'sample_sheet_visium_hd.csv')

sample_sheet = pd.read_csv(sample_sheet_path)
sample_lut = sample_sheet.set_index('section_id')


# --- Also get some tissue metadata

resource_dir = os.path.join(project_dir, 'resources')
metadata_file = 'tissue_metadata.csv'

tissue_metadata = pd.read_csv(os.path.join(resource_dir, metadata_file))

# Merge the two DataFrames on the common column 'tissue_id'
merged_df = pd.merge(sample_sheet, tissue_metadata, on='tissue_id')
# Create a sample - subtype dictionary
subtype_dict = dict(zip(merged_df['section_id'], merged_df['clinical_subtype']))


# --- Load the stlearn adata objects

resolution = 'square_016um'
threshold = 0.5
threshold_adj = str(threshold).split('.')[-1]
bin_size = 25  # um

lr_dir = os.path.join(res_dir, '04_vascular_microenvironment', '04_ligand_receptor_analysis', '01_run_LR_analysis_stlearn', f'{bin_size}um_bin', 'csv_files')

most_recent_dir_path = lr_dir

adata_dict_stlearn = {}

for sample in sample_list:
    # Load the anndata object
    filename = f'{sample}_HEV_LR_analysis_stlearn.h5ad'
    adata_dict_stlearn[sample] = ad.read_h5ad(os.path.join(most_recent_dir_path, filename))


# --- Load and add the decon results

decon_resolution = resolution

resolution_numbers = re.findall(r'\d+', decon_resolution)[0]
resolution_numbers = int(resolution_numbers)

decon_level = 'subset'
decon_dir = os.path.join(res_dir, '02_cell_type_mapping', 'deconvolution', 'visium_hd', f'{decon_level}_level', decon_resolution,)

for sample in sample_list:
    decon_file = f'{sample}_{decon_level}_{resolution}_q05_cell_abundance_w_sf.csv'
    decon_res = pd.read_csv(os.path.join(decon_dir, decon_file), index_col=0)
    # Convert to proportions within each spot
    decon_res = decon_res.div(decon_res.sum(axis=1), axis=0)
    cell_types = list(decon_res.columns)
    adata_dict_stlearn[sample].obs[cell_types] = decon_res


# --- Also load the revisited spot labels & add to adata
# These distinguish whether a nearby spot is nearby to TLS_HEV or non_TLS_HEV

lr_dir = os.path.join(res_dir, '04_vascular_microenvironment', '04_ligand_receptor_analysis', '03_plot_LR_per_sample', f'{bin_size}um_bin', 'csv_files')

most_recent_dir_path = lr_dir

for sample in sample_list:
    spot_label_file = f'{sample}_HEV_spot_grouping.csv'
    label_df = pd.read_csv(os.path.join(most_recent_dir_path, spot_label_file), index_col=0)
    adata_dict_stlearn[sample].obs['nearby_spots_revisited'] = label_df['nearby_spots_revisited']


# --- Load the miniatlas classifications

# These are needed to check the alignment between major and minor levels
miniatlas_dir = os.path.join(project_dir, 'resources', 'cell_types', 'miniatlas')
miniatlas_file = 'miniatlas_subset_classifications.csv'

miniatlas_df = pd.read_csv(os.path.join(miniatlas_dir, miniatlas_file),)

# Also define a dict with shorter names for the subset level cell types
miniatlas_sub_dict = miniatlas_df.set_index('miniatlas_subset')[['miniatlas_subset_short']].to_dict()

# Order miniatlas_major by lineage
lineage_order = ['Epithelial', 'Immune', 'Stromal']

major_cell_types = miniatlas_df[f'miniatlas_major'].dropna().unique().tolist()

# Sort based on lineage first (recorded in the miniatlas df)
def get_lineage(cell_type, level='major'):
    if cell_type in miniatlas_df[f'miniatlas_{level}'].values:
        lineage = miniatlas_df[miniatlas_df[f'miniatlas_{level}'] == cell_type]['lineage'].values[0]
    else:
        lineage = 'other'
    return lineage

major_cell_types_sorted = sorted(major_cell_types, key=lambda x: (lineage_order.index(get_lineage(x)) if get_lineage(x) in lineage_order else len(lineage_order), x))


# --- Load plotting colors

color_dir = os.path.join(project_dir, 'resources', 'plotting_colors')

# Use the major level for colors
color_level = 'major'
color_file = f'miniatlas_{color_level}_colors.csv'

color_dict = pd.read_csv(os.path.join(color_dir, color_file), index_col=0).to_dict()['color']

# Add the other levels to the color_dict
for cell_type in miniatlas_df[f'miniatlas_{decon_level}'].unique():
    if cell_type not in color_dict.keys():
        # Check that the value is not nan
        if pd.isna(cell_type):
            continue
        print(f'Adding color for {cell_type}')
    # Get the corresponding major cell type
    major_type = miniatlas_df[miniatlas_df[f'miniatlas_{decon_level}'] == cell_type][f'miniatlas_major'].values
    # Find which value is in the color_dict
    major_type_in_dict = [m for m in major_type if m in color_dict.keys()]
    if len(major_type_in_dict) > 0:
        color_dict[cell_type] = color_dict[major_type_in_dict[0]]
    else:
        color_dict[cell_type] = 'lightgrey'


# --- Find the spots that send LTB signals

adata = adata_dict_stlearn[sample_list[0]]

# Extract all significant lrs as list
lr = 'LTB_LTBR'

# Initialize a list to collect results
spot_count_results = []
group_count_results = []

print(f'Analyzing LR: {lr}')

avg_comp_dict = {}

for sample in sample_list:
    adata = adata_dict_stlearn[sample]
    # Check if the LR is in the adata object
    if 'lr_summary' not in adata.uns or lr not in adata.uns['lr_summary'].index:
        print(f'LR {lr} not found in {sample}. Skipping...')
        continue
    # Get the list of all LR pair names
    lr_pairs = adata.uns['lr_summary']
    # Find the index of the LR pair of interest
    if isinstance(lr_pairs, pd.DataFrame):
        lr_index = lr_pairs.index.get_loc(lr)
    else:
        if isinstance(lr_pairs, Iterable):
            lr_index = lr_pairs.index(lr)
        else:
            raise ValueError("lr_pairs should be either a DataFrame or a list of LR pairs.")
    # Get adjusted p-values and scores
    p_adj_matrix = adata.obsm['p_adjs'][:, lr_index]
    sig_mask = p_adj_matrix < 0.05
    barcodes_sig = adata.obs_names[sig_mask]
    if len(barcodes_sig) == 0:
        continue  # No significant spots found
    # Subset to significant barcodes
    adata_sub = adata[adata.obs_names.isin(barcodes_sig)]
    # Subset to nearby spots near TLS-HEV as these appear to be main senders
    # Going to exclude the non-TLS HEV nearby spots bc these will have very different composition
    adata_sub = adata_sub[adata_sub.obs['nearby_spots_revisited'] == 'nearby_near_TLS_HEV']
    # Calculate average composition of these spots
    avg_comp = adata_sub.obs[cell_types].mean()
    # Also calculate the average comp of all nearby spots for comparison
    adata_nearby = adata[adata.obs['nearby_spots_revisited'] == 'nearby_near_TLS_HEV']
    # Remove the spots sending LTB signals
    adata_nearby = adata_nearby[~adata_nearby.obs_names.isin(barcodes_sig)]
    avg_comp_all = adata_nearby.obs[cell_types].mean()
    # Calculate fold change
    fold_change = avg_comp / avg_comp_all
    # Make log2 fold change
    fold_change = np.log2(fold_change.replace(0, np.nan))  # Replace 0 with NaN to avoid -inf
    # Also record the absolute change
    abs_change = avg_comp - avg_comp_all
    # Make a df with all info
    avg_comp_df = pd.DataFrame({
        'cell_type': cell_types,
        'avg_composition_sending_nearby': avg_comp.values,
        'avg_composition_all_nearby': avg_comp_all.values,
        'LFC': fold_change.values,
        'abs_change': abs_change.values,
        'sample_id': sample,
    })
    avg_comp_dict[sample] = avg_comp_df

# Extract a single df for all samples
all_avg_comp_df = pd.concat(avg_comp_dict.values(), ignore_index=True)

# Sort by median LFC
all_avg_comp_df['cell_type'] = pd.Categorical(all_avg_comp_df['cell_type'], 
                                              categories=all_avg_comp_df.groupby('cell_type')['LFC'].median().sort_values(ascending=False).index,
                                              ordered=True)

sorted_cell_types = all_avg_comp_df['cell_type'].cat.categories.tolist()


# --- Statistical testing for enrichment/depletion per cell type

def analyze_cell_type_enrichment(lfc_values, cell_name):
    """
    Performs a one-sample t-test on LFC values to check for significant
    enrichment (mean > 0) or depletion (mean < 0).
    """
    print(f"--- Analyzing {cell_name} ---")
    n = len(lfc_values)
    # A t-test with n=2 is technically possible but has 0 degrees of freedom 
    # for the t-distribution, so it's not meaningful.
    if n < 3: 
        print(f"Cannot perform t-test: only {n} samples. Need at least 3.")
        return None, None, None
    # 1. Perform the one-sample t-test directly on LFC values
    #    Test against a population mean (popmean) of 0.
    t_statistic, p_value = stats.ttest_1samp(lfc_values, popmean=0)
    # 2. Get the mean LFC
    mean_lfc = np.mean(lfc_values)
    print("\n")
    return t_statistic, p_value, mean_lfc

stat_df = pd.DataFrame(columns=['cell_type', 't_statistic', 'p_value', 'mean_correlation'])

for cell_type in sorted_cell_types:
    LFC_df = all_avg_comp_df[all_avg_comp_df['cell_type'] == cell_type]['LFC'].values
    t_statistic, p_value, mean_correlation = analyze_cell_type_enrichment(LFC_df, cell_type)
    stat_df = stat_df._append({'cell_type': cell_type, 't_statistic': t_statistic, 'p_value': p_value, 'mean_correlation': mean_correlation}, ignore_index=True)

alpha = 0.05
stat_df['significant'] = stat_df['p_value'] < alpha


# --- Plot 1: Barplot of LFC as boxplot for all cell types

sc.set_figure_params(fontsize=10)

# Properties for boxplot
whisker_props = {'color': 'black'}
median_props = {'color': 'black'}

fig, ax = plt.subplots(figsize=(5, 3))

sns.boxplot(data=all_avg_comp_df, x='cell_type', y='LFC', showfliers=False, linewidth=1, whiskerprops=whisker_props, 
                        medianprops=median_props, width=0.8, palette=color_dict, ax=ax, saturation=1)

sns.stripplot(data=all_avg_comp_df, x='cell_type', y='LFC', color='black', jitter=True, size=3,
                        marker='o', ax=ax)

# Beautify the plot
for bar in ax.patches: # Adjust the colors of the boxplot
    # Get the RGBA values of the fill color
    fill_color = bar.get_facecolor()
    fill_color = list(fill_color)
    darkened_color = fill_color.copy()
    for i in range(len(darkened_color)-1):
        darkened_color[i] *= 0.8
    fill_color[-1] = 0.8 # make the fill color a bit transparent
    bar.set_color(tuple(fill_color))
    bar.set_edgecolor(tuple(darkened_color))

# Automatically resize based on number of xticks
n_ticks = len(ax.get_xticks())
ax_pos = ax.get_position()
ax.set_position([ax_pos.x0, ax_pos.y0, 0.05*n_ticks, ax_pos.height])

# Rotate xlabels
ax.set_xticklabels(ax.get_xticklabels(), rotation=90,)

# Adjust xticklabels using miniatlas_short notations
new_labels = [miniatlas_sub_dict['miniatlas_subset_short'].get(label.get_text(), label.get_text()) for label in ax.get_xticklabels()]
ax.set_xticklabels(new_labels)

# Add thick gray line at y=0
ax.axhline(0, color='lightgray', linewidth=2, zorder=0)

# Add an asterisk for significant cell types
# Just above the top error bar of each boxplot
ylim = ax.get_ylim()
range_y = ylim[1] - ylim[0]

offset = range_y * 0.003  # 0.5% of the y-axis range

for i, cell_type in enumerate(sorted_cell_types):
    cell_stats = stat_df[stat_df['cell_type'] == cell_type]
    if not cell_stats.empty and cell_stats['significant'].values[0]:
        # Get the y position of the top whisker
        cell_data = all_avg_comp_df[all_avg_comp_df['cell_type'] == cell_type]['LFC']
        whisker_top = cell_data.max()
        # Place the asterisk slightly above the whisker
        ax.text(i, whisker_top + offset, '*', ha='center', va='bottom', color='darkred', fontsize=9, fontweight='bold')

# Increase ymax a bit to make space for asterisks
ax.set_ylim(top=ax.get_ylim()[1]*1.1)

# Add some text at top right of ax explaining * p<0.05 (one-sample t-test on Fisher-transformed r-values)
# Only if there is at least one significant cell type
if stat_df['significant'].any():
    ax.text(0.993, 0.97, '* p<0.05', ha='right', va='top', transform=ax.transAxes, fontsize=8, color='darkred',
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.5))

# Create legend handles (colored boxes with labels)
labels = major_cell_types_sorted
handles = [
    mpl.patches.Patch(facecolor=color_dict[label], label=label, edgecolor='black', linewidth=0.8)
    for label in labels if label in color_dict
]

legend = ax.legend(handles=handles, bbox_to_anchor=(1, 0.5), loc='center left', frameon=False, handlelength=0.8, fontsize=10)

# Make the legend boxes look the same as the barplot
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


lr_adj = lr.replace('_', '→')   # gene names as recorded in the object
ax.set_title(f'Cells enriched in spots sending {lr_adj} signals', fontweight='bold')

# Position title
title_x = ax.get_position().x0 + ax.get_position().width / 2
title_y = ax.get_position().y1 + 0.2

plt.suptitle(f'All samples - Cell type composition of spots sending {lr_adj} signals in nearby_TLS_HEV vs all_TLS_HEV (bin_size={bin_size}µm)', fontsize=10, x=title_x, y=title_y, fontweight='bold')


# --- Plot 2: Same as above but only cell types with positive median LFC
# To make it a bit easier to interpret (we're looking for potential senders)

median_LFC = all_avg_comp_df.groupby('cell_type')['LFC'].median()
pos_cell_types = median_LFC[median_LFC > 0].index.tolist()

# Filter df
all_avg_comp_df_pos = all_avg_comp_df[all_avg_comp_df['cell_type'].isin(pos_cell_types)]

# Reorder the cell types based on median LFC
all_avg_comp_df_pos['cell_type'] = pd.Categorical(all_avg_comp_df_pos['cell_type'],
                                                categories=pos_cell_types,
                                                ordered=True)

fig, ax = plt.subplots(figsize=(5, 3))

sns.boxplot(data=all_avg_comp_df_pos, x='cell_type', y='LFC', showfliers=False, linewidth=1, whiskerprops=whisker_props, 
                        medianprops=median_props, width=0.8, palette=color_dict, ax=ax, saturation=1)

sns.stripplot(data=all_avg_comp_df_pos, x='cell_type', y='LFC', color='black', jitter=True, size=3,
                        marker='o', ax=ax)

# Beautify the plot
for bar in ax.patches: # Adjust the colors of the boxplot
    # Get the RGBA values of the fill color
    fill_color = bar.get_facecolor()
    fill_color = list(fill_color)
    darkened_color = fill_color.copy()
    for i in range(len(darkened_color)-1):
        darkened_color[i] *= 0.8
    fill_color[-1] = 0.8 # make the fill color a bit transparent
    bar.set_color(tuple(fill_color))
    bar.set_edgecolor(tuple(darkened_color))

# Automatically resize based on number of xticks
n_ticks = len(ax.get_xticks())
ax_pos = ax.get_position()
ax.set_position([ax_pos.x0, ax_pos.y0, 0.05*n_ticks, ax_pos.height*0.5])

# Rotate xlabels
ax.set_xticklabels(ax.get_xticklabels(), rotation=90,)

# Adjust xticklabels using miniatlas_short notations
new_labels = [miniatlas_sub_dict['miniatlas_subset_short'].get(label.get_text(), label.get_text()) for label in ax.get_xticklabels()]
ax.set_xticklabels(new_labels)

# Add thick gray line at y=0
ax.axhline(0, color='lightgray', linewidth=2, zorder=0)

# Add an asterisk for significant cell types
# Just above the top error bar of each boxplot
ylim = ax.get_ylim()
range_y = ylim[1] - ylim[0]

offset = range_y * 0.003  # 0.5% of the y-axis range

for i, cell_type in enumerate(pos_cell_types):
    cell_stats = stat_df[stat_df['cell_type'] == cell_type]
    if not cell_stats.empty and cell_stats['significant'].values[0]:
        # Get the y position of the top whisker
        cell_data = all_avg_comp_df[all_avg_comp_df['cell_type'] == cell_type]['LFC']
        whisker_top = cell_data.max()
        # Place the asterisk slightly above the whisker
        ax.text(i, whisker_top + offset, '*', ha='center', va='bottom', color='darkred', fontsize=9, fontweight='bold')

# Calculate equal positioning x & y
offset = 0.04
aspect = ax.get_position().width / ax.get_position().height
x_pos = 1 - offset / aspect
y_pos = 1 - offset

# Add some text at top right of ax explaining * p<0.05 (one-sample t-test on Fisher-transformed r-values)
# Check if there is at least one significant cell type among the positive LFC ones
stat_df_sub = stat_df[stat_df['cell_type'].isin(pos_cell_types)]
if stat_df_sub['significant'].any():
    # Increase ymax a bit to make space for asterisks
    ax.set_ylim(top=ax.get_ylim()[1]*1.1)
    ax.text(x_pos, y_pos, '* p<0.05', ha='right', va='top', transform=ax.transAxes, fontsize=8, color='darkred',
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.5))

# Create legend handles (colored boxes with labels)
# Filter major_cell_types_sorted for major cell types plotted
cell_types_plotted = all_avg_comp_df_pos['cell_type'].unique().tolist()

# Find matching major labels using miniatlas_df
miniatlas_subset_to_major = miniatlas_df.set_index('miniatlas_subset')['miniatlas_major'].to_dict()
major_cell_types_plotted = []
for ct in cell_types_plotted:
    if ct in miniatlas_subset_to_major:
        major_ct = miniatlas_subset_to_major[ct]
        if major_ct not in major_cell_types_plotted:
            major_cell_types_plotted.append(major_ct)

labels = [label for label in major_cell_types_sorted if label in major_cell_types_plotted]

handles = [
    mpl.patches.Patch(facecolor=color_dict[label], label=label, edgecolor='black', linewidth=0.8)
    for label in labels if label in color_dict
]

legend = ax.legend(handles=handles, bbox_to_anchor=(1, 0.5), loc='center left', frameon=False, handlelength=0.8, fontsize=10)

# Make the legend boxes look the same as the barplot
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

lr_adj = lr.replace('_', '→')   # gene names as recorded in the object
ax.set_title(f'Cells enriched in spots sending {lr_adj} signals to HEVs', fontweight='bold')

# Set labels
ax.set_ylabel('LFC', fontweight='bold')
ax.set_xlabel('Cell type', fontweight='bold')

# Position title
title_x = ax.get_position().x0 + ax.get_position().width / 2
title_y = ax.get_position().y1 + 0.2

plt.suptitle(f'All samples - Cell type composition of spots sending {lr_adj} signals in nearby_TLS_HEV vs all_TLS_HEV (bin_size={bin_size}µm)', fontsize=10, x=title_x, y=title_y, fontweight='bold')


# --- Plot 3: Same as Plot 2 but only top 10 cell types for brevity

# Get top 10 cell types by median LFC
top10_cell_types = median_LFC[median_LFC > 0].sort_values(ascending=False).head(10).index.tolist()

# Filter df
all_avg_comp_df_top10 = all_avg_comp_df[all_avg_comp_df['cell_type'].isin(top10_cell_types)]

# Reorder the cell types based on median LFC
all_avg_comp_df_top10['cell_type'] = pd.Categorical(all_avg_comp_df_top10['cell_type'],
                                                categories=top10_cell_types,
                                                ordered=True)

fig, ax = plt.subplots(figsize=(5, 3))

sns.boxplot(data=all_avg_comp_df_top10, x='cell_type', y='LFC', showfliers=False, linewidth=1, whiskerprops=whisker_props, 
                        medianprops=median_props, width=0.8, palette=color_dict, ax=ax, saturation=1)

sns.stripplot(data=all_avg_comp_df_top10, x='cell_type', y='LFC', color='black', jitter=True, size=3,
                        marker='o', ax=ax)

# Beautify the plot
for bar in ax.patches: # Adjust the colors of the boxplot
    # Get the RGBA values of the fill color
    fill_color = bar.get_facecolor()
    fill_color = list(fill_color)
    darkened_color = fill_color.copy()
    for i in range(len(darkened_color)-1):
        darkened_color[i] *= 0.8
    fill_color[-1] = 0.8 # make the fill color a bit transparent
    bar.set_color(tuple(fill_color))
    bar.set_edgecolor(tuple(darkened_color))

# Automatically resize based on number of xticks
n_ticks = len(ax.get_xticks())
ax_pos = ax.get_position()
ax.set_position([ax_pos.x0, ax_pos.y0, 0.05*n_ticks, ax_pos.height*0.5])

# Rotate xlabels
ax.set_xticklabels(ax.get_xticklabels(), rotation=90,)

# Adjust xticklabels using miniatlas_short notations
new_labels = [miniatlas_sub_dict['miniatlas_subset_short'].get(label.get_text(), label.get_text()) for label in ax.get_xticklabels()]
ax.set_xticklabels(new_labels)

# Add thick gray line at y=0
ax.axhline(0, color='lightgray', linewidth=2, zorder=0)

# Add an asterisk for significant cell types
# Just above the top error bar of each boxplot
ylim = ax.get_ylim()
range_y = ylim[1] - ylim[0]

offset = range_y * 0.02  # 2% of the y-axis range

for i, cell_type in enumerate(top10_cell_types):
    cell_stats = stat_df[stat_df['cell_type'] == cell_type]
    if not cell_stats.empty:
        # Get the y position of the top whisker
        cell_data = all_avg_comp_df[all_avg_comp_df['cell_type'] == cell_type]['LFC']
        whisker_top = cell_data.max()
        if cell_stats['significant'].values[0]:
            # Place the asterisk slightly above the whisker
            ax.text(i, whisker_top + offset, '*', ha='center', va='bottom', color='darkred', fontsize=9, fontweight='bold')
        else:
            # Place the p-value above the whisker, rounded to 2 digits
            pval = cell_stats['p_value'].values[0]
            ax.text(i, whisker_top + offset, f"{pval:.3f}".lstrip('0'), ha='center', va='bottom', color='black', fontsize=7)

# Increase ymax a bit to make space for asterisks/p-values
ax.set_ylim(top=ax.get_ylim()[1]*1.2)

# Calculate equal positioning x & y
offset = 0.04
aspect = ax.get_position().width / ax.get_position().height
x_pos = 1 - offset / aspect
y_pos = 1 - offset

# Add some text at top right of ax explaining * p<0.05 (one-sample t-test on Fisher-transformed r-values)
# Check if there is at least one significant cell type among the top 10
stat_df_top10 = stat_df[stat_df['cell_type'].isin(top10_cell_types)]
if stat_df_top10['significant'].any():
    # Increase ymax a bit to make space for asterisks
    ax.set_ylim(top=ax.get_ylim()[1]*1.1)
    ax.text(x_pos, y_pos, '* p<0.05', ha='right', va='top', transform=ax.transAxes, fontsize=8, color='darkred',
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.5))

# Create legend handles (colored boxes with labels)
# Filter major_cell_types_sorted for major cell types plotted
cell_types_plotted_top10 = all_avg_comp_df_top10['cell_type'].unique().tolist()

# Find matching major labels using miniatlas_df
major_cell_types_plotted_top10 = []
for ct in cell_types_plotted_top10:
    if ct in miniatlas_subset_to_major:
        major_ct = miniatlas_subset_to_major[ct]
        if major_ct not in major_cell_types_plotted_top10:
            major_cell_types_plotted_top10.append(major_ct)

labels_top10 = [label for label in major_cell_types_sorted if label in major_cell_types_plotted_top10]

handles_top10 = [
    mpl.patches.Patch(facecolor=color_dict[label], label=label, edgecolor='black', linewidth=0.8)
    for label in labels_top10 if label in color_dict
]

legend = ax.legend(handles=handles_top10, bbox_to_anchor=(1, 0.5), loc='center left', frameon=False, handlelength=0.8, fontsize=10)

# Make the legend boxes look the same as the barplot
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

lr_adj = lr.replace('_', '→')   # gene names as recorded in the object
ax.set_title(f'Cellular enrichment\nBins sending {lr_adj} signals to HEVs', fontweight='bold')

# Set labels
ax.set_ylabel('LFC', fontweight='bold')
ax.set_xlabel('Cell type', fontweight='bold')

# Position title
title_x = ax.get_position().x0 + ax.get_position().width / 2
title_y = ax.get_position().y1 + 0.3

plt.suptitle(f'All samples - Cell type composition of spots sending {lr_adj} signals in nearby_TLS_HEV vs all_TLS_HEV (bin_size={bin_size}µm)', fontsize=10, x=title_x, y=title_y, fontweight='bold')


# --- Save the output

qc_dir = os.path.join(res_dir, '04_vascular_microenvironment', '04_ligand_receptor_analysis', '04_check_composition_LTB_signaling_bins')

os.makedirs(qc_dir, exist_ok=True)

# Save plots
out_dir_plots = os.path.join(qc_dir, 'plots')
os.makedirs(out_dir_plots, exist_ok=True)

# Save plots
filename = f'All_samples_HEV_{lr}_cellular_enrichment.pdf'
save_multi_image(os.path.join(out_dir_plots, filename))

# Also save csv for future reference
out_dir_csv = os.path.join(qc_dir, 'csv_files')
os.makedirs(out_dir_csv, exist_ok=True)

filename = f'All_samples_HEV_{lr}_cellular_enrichment.csv'
all_avg_comp_df.to_csv(os.path.join(out_dir_csv, filename), index=False)

# And save the stat df
filename = f'All_samples_HEV_{lr}_cellular_enrichment_statistics.csv'
stat_df.to_csv(os.path.join(out_dir_csv, filename), index=False)


# --- Finish

print("Script successfully completed")
