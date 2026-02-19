#!/usr/bin/env python3

# --- Aim
# Evaluate and visualize germinal center (GC) marker gene expression in epithelial subsets across samples.
# Filter GC signature genes based on expression and prevalence, and generate comparative dotplots by subtype.


# --- Define environment
# env: spatial-pilot

# --- Load packages
import scanpy as sc
import anndata as ad
import pandas as pd
import json

import numpy as np
import os

from datetime import date

import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib as mpl
import argparse

# Custom functions
from sophiesfunctions.save_multi_image import save_multi_image

mpl.rcParams['pdf.fonttype'] = 42

# --- Arguments

# --- Argparse arguments
parser = argparse.ArgumentParser(description="Import Visium unfiltered data & plot",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)


parser.add_argument("-s", "--samples",
    nargs="*",  # 0 or more values expected => creates a list
    type=str,
    help="List of samples")
parser.add_argument("-r", "--res_dir",
                    help="Path to results directory (optional)")
parser.add_argument("--h5ad_dir", help="Path to directory for .h5ad output files") # Not necessary for this script
parser.add_argument("--project_dir", help="Project directory (root)")
parser.add_argument("--sample_sheet", help="Path to sample sheet CSV file", default=None)

args = vars(parser.parse_args())

# Arguments
h5ad_dir = args["h5ad_dir"]
res_dir = args["res_dir"]
sample_list = args["samples"]
project_dir = args["project_dir"]

print(f'Samples in python: {sample_list}')


# --- Specify layer to use

use_layer = 'log_norm'

# --- Extract some info from the sample sheet

# Use sample_sheet argument if provided, else default
if args.get("sample_sheet"):
    sample_sheet_path = args["sample_sheet"]
else:
    sample_sheet_path = os.path.join(project_dir, 'config', 'sample_sheet.csv')

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


# --- Read in the epithelial path annotations

epithelial_annotations = {}

annotation_dir_base = os.path.join(project_dir, 'resources', 'path_annotations')

for sample in sample_list:
    annotation_dir = os.path.join(annotation_dir_base, sample)
    if os.path.exists(annotation_dir):
        annotation_files = os.listdir(annotation_dir)
        epithelial_file = [file for file in annotation_files if 'epithelial' in file.lower()][0]
        if os.path.exists(os.path.join(annotation_dir, epithelial_file)):
            epithelial_annotations[sample] = pd.read_csv(os.path.join(annotation_dir, epithelial_file), index_col=0)


# --- Load data
# Multiple samples as a dictionary seems like an okay idea
adata_dict = {}

subtype_key = 'subtype'
sample_id_key = 'sample_id'

# Load only samples with epithelial annotations

for sample in epithelial_annotations:
    adata_dict[sample] = ad.read_h5ad(os.path.join(h5ad_dir, 'visium', f'{sample}.h5ad'))
    # Add subtype to obs for easier splitting
    adata_dict[sample].obs[subtype_key] = subtype_dict[sample]
    # Also add sample id for easier splitting
    adata_dict[sample].obs[sample_id_key] = sample
    # Add epithelial annotations
    adata_dict[sample].obs['epithelial'] = epithelial_annotations[sample]


# --- Order based on unique sample ID

ordered_keys = adata_dict.keys()
ordered_keys = sorted(ordered_keys)


# --- Load the reworked annotations
# Summarized by major category

path_layer = 'epithelial'

reworked_anno_dir = os.path.join(project_dir, 'results', 'path_annotations', 'reworked_annotations',)

for sample in epithelial_annotations:
    sample_anno_dir = os.path.join(reworked_anno_dir, sample, 'csv_files', 'reworked_annotations_per_layer')
    if os.path.exists(sample_anno_dir):
        filename = f'{sample}_{path_layer}_reworked.csv'
        path_to_file = os.path.join(sample_anno_dir, filename)
        if os.path.exists(path_to_file):
            print(f'Found reworked annotations for {sample}')
            reworked_anno = pd.read_csv(path_to_file, index_col=0)
            adata_dict[sample].obs[f'{path_layer}_subset'] = reworked_anno
        else:
            adata_dict[sample].obs[f'{path_layer}_subset'] = None


# Also load the way the raw annotations were categorized into major categories
# Function to merge dictionaries of lists without duplicates and sort each list
def merge_and_sort_dictionaries(dicts):
    merged_dict = {}
    for d in dicts:
        for key, value in d.items():
            if key not in merged_dict:
                merged_dict[key] = []
            # Append unique values to the list
            merged_dict[key].extend([v for v in value if v not in merged_dict[key]])
    # Sort each list alphabetically ignoring capitalization
    for key in merged_dict:
        merged_dict[key] = sorted(merged_dict[key], key=str.lower)
    return merged_dict

# Load all JSON files
dicts = []

for sample in epithelial_annotations:
    sample_anno_dir = os.path.join(reworked_anno_dir, sample, 'csv_files', 'reworked_annotations_per_layer')
    if os.path.exists(sample_anno_dir):
        filename = f'{sample}_reworked_annotations_all_layers.json'
        path_to_file = os.path.join(sample_anno_dir, filename)
        if os.path.exists(path_to_file):
            with open(path_to_file, 'r') as json_file:
                data = json.load(json_file)
                if path_layer in data:
                    dicts.append(data[path_layer])

# Merge extracted dictionaries
merged_dict = merge_and_sort_dictionaries(dicts)

# Sort the keys based on the master dictionary for reworking annotations
# Load the desired categories per layer (to get the desired order)
filename = 'desired_categories_per_layer.json'
path_to_file = os.path.join(annotation_dir_base, filename)

with open(path_to_file, 'r') as json_file:
                annotations_of_interest_dict = json.load(json_file)

sorted_keys = annotations_of_interest_dict[path_layer]

# Add _mixed and _other categories
sorted_keys += [f'{path_layer}_mixed', f'{path_layer}_other']

# Sort the dictionary
ordered_major_cat_dict = {key: merged_dict[key] for key in sorted_keys if key in merged_dict}

# Add any remaining keys that were not in key_order
# There shouldn't be any but just to capture mistakes
remaining_keys = {key: merged_dict[key] for key in merged_dict if key not in sorted_keys}
ordered_major_cat_dict.update(remaining_keys)


# --- Convert to adata list & concatenate (easier for plotting)

adata_list = list(adata_dict.values())
all_adata = ad.concat(adata_list, merge='same')
all_adata.obs_names_make_unique()

del adata_list


# --- Prepare for plotting

# Split the adatas by subtype for plotting

all_subtypes = set(subtype_dict.values())

# Define a desired order for plotting the subtypes
subtype_order = ['Luminal A', 'Luminal B (ER+)', 'Luminal B (HER2+)', 'HER2+', 'TNBC']
sorted_subtypes = sorted(all_subtypes, key=lambda x: subtype_order.index(x) if x in subtype_order else float('inf'))

adata_per_subtype = {}

for subtype in sorted_subtypes:
    adata_per_subtype[subtype] = all_adata[all_adata.obs[subtype_key] == subtype]

ncols = len(sorted_subtypes)


# --- Load the GC gene list

gene_dir = os.path.join(project_dir, 'resources', 'gene_lists', 'germinal_center_scoring')
gene_file = 'GC_score_v2.csv'

gene_table = pd.read_csv(os.path.join(gene_dir, gene_file))

if 'genes' in gene_table.columns:
    gene_list = list(gene_table['genes']) # Assuming the genes are in a column named 'genes'
else:
    gene_list = list(gene_table.iloc[:, 0]) # Assuming the genes are in the first column of this table

# Remove any potential nan values
gene_list = [x for x in gene_list if isinstance(x, str) or not np.isnan(x)]

# Limit to 50 genes
top = 100
if len(gene_list) > top:
    list_for_plot = gene_list[:top]


# --- Prepare for plotting split by subtype

# --- Some default parameters for dotplot

# Set some default parameters for all plots
sc.pl.DotPlot.DEFAULT_COLORMAP='Reds'
sc.pl.DotPlot.DEFAULT_DOT_MIN=0
sc.pl.DotPlot.DEFAULT_DOT_MAX=1

# Alpha for subtype color markings
alpha = 0.7

# --- Read in colors per subtype

color_dir = os.path.join(project_dir, 'resources', 'plotting_colors')
color_file = 'SCSubtype_colors.csv'

color_lut = pd.read_csv(os.path.join(color_dir, color_file))
color_lut = color_lut.set_index('subtype')


# --- Package the plot into a function

def plot_gc_dotplot(adata_per_subtype_sub, sorted_subtypes, sample_id_key, gene_list, color_lut, vmax=3, space_per_dot=0.04, yspace=0.02, alpha=0.7, title=None, subtitle=None):
    ncols = 1
    nrows = len(sorted_subtypes)
    sc.set_figure_params(fontsize=10)
    fig, axs = plt.subplots(ncols=ncols, nrows=nrows, figsize=(5, 5))
    main_axs = []
    count_bar_axs = []
    for n, subtype in enumerate(sorted_subtypes):
        plot = sc.pl.DotPlot(
            adata_per_subtype_sub[subtype],
            var_names=gene_list,
            groupby=[sample_id_key],
            layer=use_layer,
            ax=axs[n],
            vmin=0,
            vmax=vmax
        )
        plot.add_totals(size=0.6)
        plot_axes = plot.get_axes()
        main_ax = plot_axes['mainplot_ax']
        xticks = main_ax.get_xticks()
        yticks = main_ax.get_yticks()
        ax_pos = main_ax.get_position()
        main_ax.set_position([ax_pos.x0, ax_pos.y0, len(xticks)*space_per_dot, len(yticks)*space_per_dot])
        main_axs.append(main_ax)
        if n != (nrows-1):
            main_ax.set_xticklabels([])
            plot_axes['size_legend_ax'].set_visible(False)
            plot_axes['color_legend_ax'].set_visible(False)
        count_bar_ax = plot_axes['group_extra_ax']
        count_bar_axs.append(count_bar_ax)
    # Stop axs from overlapping
    for n, ax in enumerate(main_axs):
        if n == 0:
            ax_pos_prev = ax.get_position()
        else:
            ax_pos = ax.get_position()
            new_pos = ax_pos_prev.y0 - ax_pos.height - yspace
            ax.set_position([ax_pos.x0, new_pos, ax_pos.width, ax_pos.height])
            ax_pos_prev = ax.get_position()
    # Position the count bars and set xmax
    xmax = 0
    for n, count_bar_ax in enumerate(count_bar_axs):
        main_ax = main_axs[n]
        main_ax_pos = main_ax.get_position()
        count_bar_ax.set_position([main_ax_pos.x1, main_ax_pos.y0, space_per_dot*3, main_ax_pos.height])
        if count_bar_ax.get_xlim()[1] > xmax:
            xmax = count_bar_ax.get_xlim()[1]
    xmax = xmax*1.1
    for ax in count_bar_axs:
        ax.set_xlim([0, xmax])
    # Reposition legends
    plot_axes = plot.get_axes()
    size_legend_ax = plot_axes['size_legend_ax']
    color_legend_ax = plot_axes['color_legend_ax']
    new_pos = count_bar_ax.get_position().x1 + 0.03
    color_legend_ax.set_position([new_pos, ax_pos_prev.y0, space_per_dot*5, space_per_dot*0.7])
    new_color_pos = color_legend_ax.get_position()
    size_legend_ax.set_position([new_pos, new_color_pos.y1 + 0.12, space_per_dot*5, space_per_dot*2])
    # Add colors for subtypes
    for n, subtype in enumerate(sorted_subtypes):
        color = color_lut.loc[subtype, 'color']
        height = len(main_axs[n].get_yticks()) + 1
        width = 0.4
        xlim = main_axs[n].get_xlim()
        rect = patches.Rectangle((xlim[0], -0.5), width=width, height=height, facecolor=color, alpha=alpha, zorder=-1)
        main_axs[n].add_patch(rect)
        rect = patches.Rectangle((xlim[1]-width, -0.5), width=width, height=height, facecolor=color, alpha=alpha, zorder=-1)
        main_axs[n].add_patch(rect)
    # Add a title to the first plot
    if title:
        main_axs[0].set_title(title, fontweight='bold')
    # Add a title to the entire plot
    title_x = (main_axs[0].get_position().x0 + main_axs[0].get_position().x1)/2
    title_y = main_axs[0].get_position().y1 + 0.1
    if subtitle:
        plt.suptitle(subtitle, y=title_y, x=title_x, fontweight='bold')
    return fig, main_axs


# --- Plot 1a: Expression of GC genes in invasive spots

annotation = 'invasive'

# Subset adata to be just the spots of interest
all_adata_sub = all_adata[all_adata.obs['epithelial_subset'] == annotation]

# Also subset the adata per subtype
adata_per_subtype_sub = {}
for subtype in adata_per_subtype:
    adata_per_subtype_sub[subtype] = adata_per_subtype[subtype][adata_per_subtype[subtype].obs['epithelial_subset'] == annotation]

# Set parameters for plot
space_per_dot = 0.04
yspace = 0.02
vmax=3

max_genes = 400
list_for_plot = gene_list[:max_genes]

# For all genes
plot_gc_dotplot(
    adata_per_subtype_sub, sorted_subtypes, sample_id_key, list_for_plot, color_lut,
    vmax=vmax, space_per_dot=space_per_dot, yspace=yspace, alpha=alpha,
    title=f'GC score top {top} genes',
    subtitle=f'All samples - gene expression GC score in spots annotated as {annotation} (use_layer={use_layer}, vmax={vmax})',
)


# --- Plot 1b: Same as above but for the top 50 genes

top = 50
list_for_plot = gene_list[:top]

plot_gc_dotplot(
    adata_per_subtype_sub, sorted_subtypes, sample_id_key, list_for_plot, color_lut,
    vmax=vmax, space_per_dot=space_per_dot, yspace=yspace, alpha=alpha,
    title=f'GC score top {top} genes',
    subtitle=f'All samples - gene expression GC score in spots annotated as {annotation} (use_layer={use_layer}, vmax={vmax})',
)


# --- For each gene, calculate pct expressing and mean expression per sample

dfs = []

for sample in epithelial_annotations:
    print(f'Processing sample {sample}')
    adata = adata_dict[sample]
    # Subset to invasive spots
    adata_sub = adata[adata.obs['epithelial_subset'] == annotation]
    # Get the expression DataFrame for the selected genes
    expr_df = adata_sub.to_df(layer=use_layer)[gene_list]
    # Calculate percent expressing (fraction > 0) for each gene
    pct_expressing = (expr_df > 0).sum(axis=0) / expr_df.shape[0]
    # Calculate mean expression for each gene
    mean_expression = expr_df.mean(axis=0)
    # Create a DataFrame with the results
    tmp_df = pd.DataFrame({
        'genes': gene_list,
        'pct_expressing': pct_expressing.values,
        'mean_expression': mean_expression.values,
        'sample_id': sample,
    })
    dfs.append(tmp_df)

# Concatenate the DataFrames from all samples into one DataFrame
final_df = pd.concat(dfs, ignore_index=True)

mean_exp_threshold = 1
pct_expressing_threshold = 0.75
min_samples = 3

final_df['keep'] = True
final_df['reason'] = ''

for gene in gene_list:
    gene_rows = final_df[final_df['genes'] == gene]
    mean_exp_exceeds = sum(gene_rows['mean_expression'] > mean_exp_threshold)
    pct_expressing_exceeds = sum(gene_rows['pct_expressing'] > pct_expressing_threshold)
    reasons = []
    if mean_exp_exceeds >= min_samples:
        final_df.loc[final_df['genes'] == gene, 'keep'] = False
        reasons.append('mean_expression')
    if pct_expressing_exceeds >= min_samples:
        final_df.loc[final_df['genes'] == gene, 'keep'] = False
        reasons.append('pct_expressing')
    if reasons:
        final_df.loc[final_df['genes'] == gene, 'reason'] = ';'.join(reasons)

top_genes = final_df[final_df['keep']]

filtered_genes = list(top_genes['genes'].unique())

# Apply maximum cap to filtered list
if len(filtered_genes) > max_genes:
    filtered_genes_plot = filtered_genes[:max_genes]
else:
    filtered_genes_plot = filtered_genes

print(f'Filtered genes: {len(filtered_genes)} out of {len(gene_list)}')


# --- Plot 2a: Plot the genes after filtering

# For filtered genes
plot_gc_dotplot(
    adata_per_subtype_sub, sorted_subtypes, sample_id_key, filtered_genes_plot, color_lut,
    vmax=vmax, space_per_dot=space_per_dot, yspace=yspace, alpha=alpha,
    title=f'GC score - genes after filtering (n={len(filtered_genes)}/{len(gene_list)})',
    subtitle=f'All samples - gene expression GC score in spots annotated as {annotation} (use_layer={use_layer}, vmax={vmax}, mean_exp_threshold={mean_exp_threshold}, pct_expressing_threshold={pct_expressing_threshold}, min_samples={min_samples})',
)


# --- Plot 2b: Plot only the top genes after filtering

top = 50
list_for_plot = filtered_genes[:top]

# For filtered genes
plot_gc_dotplot(
    adata_per_subtype_sub, sorted_subtypes, sample_id_key, list_for_plot, color_lut,
    vmax=vmax, space_per_dot=space_per_dot, yspace=yspace, alpha=alpha,
    title=f'GC score - genes after filtering (top {top})',
    subtitle=f'All samples - gene expression GC score in spots annotated as {annotation} (use_layer={use_layer}, vmax={vmax}, mean_exp_threshold={mean_exp_threshold}, pct_expressing_threshold={pct_expressing_threshold}, min_samples={min_samples})',
)


# --- Plot 2c: Plot only the top genes after filtering

top = 100
list_for_plot = filtered_genes[:top]

# For filtered genes
plot_gc_dotplot(
    adata_per_subtype_sub, sorted_subtypes, sample_id_key, list_for_plot, color_lut,
    vmax=vmax, space_per_dot=space_per_dot, yspace=yspace, alpha=alpha,
    title=f'GC score - genes after filtering (top {top})',
    subtitle=f'All samples - gene expression GC score in spots annotated as {annotation} (use_layer={use_layer}, vmax={vmax}, mean_exp_threshold={mean_exp_threshold}, pct_expressing_threshold={pct_expressing_threshold}, min_samples={min_samples})',
)


# --- Prepare for saving a new gene_table with the top genes

# Filter gene table for genes to keep
top_genes_table = gene_table[gene_table['genes'].isin(filtered_genes)]

# To the original table, also add a column with reason_discarded
gene_table['reason_discarded_invasive'] = ''

# Add the reason for discarding genes
for gene in gene_list:
    gene_rows = final_df[final_df['genes'] == gene]
    reasons = gene_rows['reason'].values[0] if not gene_rows.empty else ''
    gene_table.loc[gene_table['genes'] == gene, 'reason_discarded_invasive'] = reasons

# Optionally, set empty strings to np.nan if you want missing reasons to be NaN
gene_table['reason_discarded_invasive'].replace('', np.nan, inplace=True)

# Also extract a df with filtered genes only
filtered_genes_table = gene_table[gene_table['reason_discarded_invasive'].notna()]

# --- Plot 3a: Plot all genes showing which ones were filtered out

fig, axs = plot_gc_dotplot(
    adata_per_subtype_sub, sorted_subtypes, sample_id_key, gene_list, color_lut,
    vmax=vmax, space_per_dot=space_per_dot, yspace=yspace, alpha=alpha,
    title=f'GC score - genes after filtering (n={len(filtered_genes)}/{len(gene_list)})',
    subtitle=f'All samples - gene expression GC score in spots annotated as {annotation} (use_layer={use_layer}, vmax={vmax}, mean_exp_threshold={mean_exp_threshold}, pct_expressing_threshold={pct_expressing_threshold}, min_samples={min_samples})',
)

# Figure out indeces of filtered genes from the last x
xtick_labels = axs[-1].get_xticklabels()
# Get the text from the xtick labels
xtick_text = [label.get_text() for label in xtick_labels]

# Get the filtered genes from gene_table
removed_genes = gene_table[gene_table['reason_discarded_invasive'].notna()]['genes'].tolist()

# Get the indices of the filtered genes in the xtick text
filtered_indices = [i for i, gene in enumerate(xtick_text) if gene in removed_genes]

# Draw a red rectangle at the filtered indices on each ax
for ax in axs:
    for idx in filtered_indices:
        x = ax.get_xticks()[idx]
        y_min, y_max = ax.get_ylim()
        width = 0.9  # Width of the rectangle
        rect = patches.Rectangle((x - width/2, y_min), width, y_max - y_min, linewidth=0, facecolor='red', zorder=-1, alpha=0.3)
        ax.add_patch(rect)

# Also give the xtick text a red highlight
for idx in filtered_indices:
    xtick_labels[idx].set_color('red')
    xtick_labels[idx].set_fontweight('bold')

axs[-1].set_xticklabels(xtick_labels)


# --- Plot 3b: Same as above but for the top 50 genes

top = 50
list_for_plot = gene_list[:top]

# Calculate what percentage was filtered out in the top 50 genes
n_filtered = sum(filtered_genes_table['genes'].isin(list_for_plot))

fig, axs = plot_gc_dotplot(
    adata_per_subtype_sub, sorted_subtypes, sample_id_key, list_for_plot, color_lut,
    vmax=vmax, space_per_dot=space_per_dot, yspace=yspace, alpha=alpha,
    title=f'GC score - top {top} genes filtered ({n_filtered}/{len(list_for_plot)})',
    subtitle=f'All samples - gene expression GC score in spots annotated as {annotation} (use_layer={use_layer}, vmax={vmax}, mean_exp_threshold={mean_exp_threshold}, pct_expressing_threshold={pct_expressing_threshold}, min_samples={min_samples})',
)

# Figure out indeces of filtered genes from the last x
xtick_labels = axs[-1].get_xticklabels()
# Get the text from the xtick labels
xtick_text = [label.get_text() for label in xtick_labels]

# Get the filtered genes from gene_table
removed_genes = gene_table[gene_table['reason_discarded_invasive'].notna()]['genes'].tolist()

# Get the indices of the filtered genes in the xtick text
filtered_indices = [i for i, gene in enumerate(xtick_text) if gene in removed_genes]

# Draw a red rectangle at the filtered indices on each ax
for ax in axs:
    for idx in filtered_indices:
        x = ax.get_xticks()[idx]
        y_min, y_max = ax.get_ylim()
        width = 0.9  # Width of the rectangle
        rect = patches.Rectangle((x - width/2, y_min), width, y_max - y_min, linewidth=0, facecolor='red', zorder=-1, alpha=0.3)
        ax.add_patch(rect)

# Also give the xtick text a red highlight
for idx in filtered_indices:
    xtick_labels[idx].set_color('red')
    xtick_labels[idx].set_fontweight('bold')

axs[-1].set_xticklabels(xtick_labels)

plt.savefig('test.pdf', bbox_inches='tight')

# --- Plot 3c: Same as above but for the top 100 genes

top = 100
list_for_plot = gene_list[:top]

# Calculate what percentage was filtered out in the top 50 genes
n_filtered = sum(filtered_genes_table['genes'].isin(list_for_plot))

fig, axs = plot_gc_dotplot(
    adata_per_subtype_sub, sorted_subtypes, sample_id_key, list_for_plot, color_lut,
    vmax=vmax, space_per_dot=space_per_dot, yspace=yspace, alpha=alpha,
    title=f'GC score - top {top} genes filtered ({n_filtered}/{len(list_for_plot)})',
    subtitle=f'All samples - gene expression GC score in spots annotated as {annotation} (use_layer={use_layer}, vmax={vmax}, mean_exp_threshold={mean_exp_threshold}, pct_expressing_threshold={pct_expressing_threshold}, min_samples={min_samples})',
)

# Figure out indeces of filtered genes from the last x
xtick_labels = axs[-1].get_xticklabels()
# Get the text from the xtick labels
xtick_text = [label.get_text() for label in xtick_labels]

# Get the filtered genes from gene_table
removed_genes = gene_table[gene_table['reason_discarded_invasive'].notna()]['genes'].tolist()

# Get the indices of the filtered genes in the xtick text
filtered_indices = [i for i, gene in enumerate(xtick_text) if gene in removed_genes]

# Draw a red rectangle at the filtered indices on each ax
for ax in axs:
    for idx in filtered_indices:
        x = ax.get_xticks()[idx]
        y_min, y_max = ax.get_ylim()
        width = 0.9  # Width of the rectangle
        rect = patches.Rectangle((x - width/2, y_min), width, y_max - y_min, linewidth=0, facecolor='red', zorder=-1, alpha=0.3)
        ax.add_patch(rect)

# Also give the xtick text a red highlight
for idx in filtered_indices:
    xtick_labels[idx].set_color('red')
    xtick_labels[idx].set_fontweight('bold')

axs[-1].set_xticklabels(xtick_labels)


out_dir = os.path.join(gene_dir, 'GC_score_filtered_genes')
os.makedirs(out_dir, exist_ok=True)
out_file = os.path.join(out_dir, f'GC_score_filtered_genes.csv')
top_genes_table.to_csv(out_file, index=False)
out_file = os.path.join(out_dir, f'GC_score_filtered_genes_all.csv')
gene_table.to_csv(out_file, index=False)
os.makedirs(qc_dir, exist_ok=True)
today = date.today()
today = today.strftime("%Y%m%d")
out_dir = os.path.join(qc_dir, 'filter_GC_score', 'plots', today)
os.makedirs(out_dir, exist_ok=True)
save_multi_image(os.path.join(out_dir, filename))

# --- Save the output to results/03_TLS_niche_discovery/03_extract_niche_signatures
base_results = os.path.join(project_dir, 'results', '03_TLS_niche_discovery', '03_extract_niche_signatures')
os.makedirs(base_results, exist_ok=True)

out_file = os.path.join(base_results, 'GC_score_filtered_genes.csv')
top_genes_table.to_csv(out_file, index=False)

out_file = os.path.join(base_results, 'GC_score_filtered_genes_all.csv')
gene_table.to_csv(out_file, index=False)

# Also save the plots
today = date.today().strftime("%Y%m%d")
plot_dir = os.path.join(base_results, 'filter_GC_score', 'plots', today)
os.makedirs(plot_dir, exist_ok=True)
filename = 'All_samples_GC_score_filtering.pdf'
save_multi_image(os.path.join(plot_dir, filename))

# Also save the values per sample
out_dir_csv = os.path.join(qc_dir, 'filter_GC_score', 'csv_files', today)
os.makedirs(out_dir_csv, exist_ok=True)

filename = 'All_samples_mean_pct_expressing_per_gene.csv'
final_df.to_csv(os.path.join(out_dir_csv, filename), index=False)


# --- Finish

print("Script successfully completed")
