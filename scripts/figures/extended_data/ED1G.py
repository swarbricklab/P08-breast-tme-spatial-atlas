#!/usr/bin/env python3

# --- Aim
# Molecular validation of the pathology annotations.
# For each pathologist-defined region (across all four annotation layers), find the
# differentially expressed genes (one region vs. the rest), run GO Biological Process 2025
# term enrichment on the up-regulated genes, and draw the enrichment dot plot
# (top terms per region; dot size = combined enrichment score, color = -log10 adjusted P)
# together with a bar of the total spot count per annotation.
# This creates Extended Data Fig. 1g

# --- Define environment
# env: gseapy
# Note: concatenates all 33 Visium sections for the cohort-wide differential expression
# (peak ~33 GB).

# --- Load packages

import os
import sys
import json
import argparse

import numpy as np
import pandas as pd

import scanpy as sc
import anndata as ad
import gseapy as gp

import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.ticker import FuncFormatter

mpl.rcParams['pdf.fonttype'] = 42

# Repo-local helper functions (scripts/helper_functions.py + scripts/figures/dotplot_functions.py)
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))
from helper_functions import save_multi_image


# --- Argparse arguments

parser = argparse.ArgumentParser(description="Molecular validation of pathology annotations (DEG + GO enrichment dot plot)",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument("-s", "--samples",
    nargs="*",  # 0 or more values expected => creates a list
    type=str,
    help="List of samples")
parser.add_argument("-r", "--res_dir", help="Path to results directory")
parser.add_argument("--h5ad_dir", help="Path to directory for .h5ad output files")
parser.add_argument("--project_dir", help="Project directory (root)")
parser.add_argument("--sample_sheet", help="Path to sample sheet CSV file", default=None)
parser.add_argument("--top_n", help="Number of top GO terms to plot per annotation category", type=int, default=5)
parser.add_argument("--gene_set", help="Enrichr gene set library (a .gmt in resources/gene_lists/enrichment_analysis)",
                    default='GO_Biological_Process_2025')

args = vars(parser.parse_args())

# Arguments
project_dir = args["project_dir"]
h5ad_dir = args["h5ad_dir"]
res_dir = args["res_dir"]
sample_list = args["samples"]
top_n = args["top_n"]
gene_set = args["gene_set"]

# Default h5ad directory
if h5ad_dir is None:
    h5ad_dir = os.path.join(project_dir, 'data', 'anndata_objects')

print(f'Samples in python: {sample_list}')
print(f'Top N terms: {top_n}')

# Import the shared dot plot helpers
sys.path.append(os.path.join(project_dir, 'scripts', 'figures'))
import dotplot_functions


# --- Specify layer to use

use_layer = 'log_norm'

annotation_layers = ['epithelial', 'immune', 'stroma', 'additional_features']


# --- Load the pathology annotations

pathology_dir = os.path.join(project_dir, 'data', 'pathology_annotations', 'visium')

annotation_dict = {}

for sample in sample_list:
    filename = f'{sample}_pathology_annotations.csv'
    path_to_file = os.path.join(pathology_dir, filename)
    if os.path.exists(path_to_file):
        pathology_anno = pd.read_csv(path_to_file, index_col=0)
        annotation_dict[sample] = pathology_anno[annotation_layers].copy()
        print(f"Loaded pathology annotations for {sample}")
    else:
        print(f"Warning: No pathology annotations found for {sample}")

# Keep only samples that also have an AnnData object
annotated_samples = [
    sample for sample in annotation_dict
    if os.path.exists(os.path.join(h5ad_dir, 'visium', f'{sample}.h5ad'))
]
annotated_samples = sorted(annotated_samples)

print(f'Samples with annotations and AnnData: {annotated_samples}')


# --- Load the AnnData objects and attach the annotations
# Only the log-normalized expression and the annotations are needed, so each object is reduced to
# those (log_norm as .X) as it is read; keeping the 33 full objects needed ~190 GB.

adata_list = []

for sample in annotated_samples:
    adata = ad.read_h5ad(os.path.join(h5ad_dir, 'visium', f'{sample}.h5ad'))
    obs = pd.DataFrame(index=adata.obs_names)
    # Record the sample of origin
    obs['sample_id'] = sample
    # Attach the per-layer annotations (index aligned on spot barcodes)
    for layer in annotation_layers:
        obs[layer] = annotation_dict[sample][layer]
    adata_list.append(ad.AnnData(X=adata.layers[use_layer], obs=obs, var=pd.DataFrame(index=adata.var_names)))
    del adata


# --- Concatenate into a single object

all_adata = ad.concat(adata_list, merge='same')
all_adata.obs_names_make_unique()

del adata_list


# --- Differential expression per annotation layer (one region vs. rest)

for layer in annotation_layers:
    print(f'Find top genes {layer}')
    rank_genes_key = f'rank_genes_{layer}'
    # Drop annotations that occur in a single spot (cannot test one vs. rest reliably)
    single_spot_clusters = all_adata.obs[layer].value_counts()[all_adata.obs[layer].value_counts() == 1].index
    single_spot_clusters = list(single_spot_clusters.unique())
    all_adata.obs[f'{layer}_2'] = all_adata.obs[f'{layer}']
    all_adata.obs[f'{layer}_2'] = np.where(all_adata.obs[f'{layer}_2'].isin(single_spot_clusters), None, all_adata.obs[f'{layer}_2'])
    # Spots without an annotation for this layer form the reference group
    all_adata.obs[f'{layer}_2'] = np.where(all_adata.obs[f'{layer}_2'].isna(), f'non_{layer}', all_adata.obs[f'{layer}_2'])
    # Rank genes (Wilcoxon, group vs. rest)
    sc.tl.rank_genes_groups(all_adata, f'{layer}_2', use_raw=False, method='wilcoxon',   # .X = use_layer
                            key_added=rank_genes_key, pts=True)


# --- Extract the differential genes per layer

deg_dict = {}

for layer in annotation_layers:
    rank_genes_key = f'rank_genes_{layer}'
    result = all_adata.uns[rank_genes_key]
    groups = result['names'].dtype.names
    reference = result['params']['reference']
    dfs = []
    for group in groups:
        df = pd.DataFrame({
            'group': group,
            'reference': reference,
            'genes': result['names'][group],
            'pvals': result['pvals'][group],
            'logfoldchanges': result['logfoldchanges'][group],
            'pts': result['pts'][group],
            'pvals_adj': result['pvals_adj'][group],
            'z-score': result['scores'][group],
        })
        dfs.append(df)
    deg_dict[layer] = pd.concat(dfs, ignore_index=True)


# --- Prepare for saving

out_dir = os.path.join(res_dir, 'figures', 'extended_data')
os.makedirs(out_dir, exist_ok=True)

csv_dir = os.path.join(out_dir, 'csv_files')
os.makedirs(csv_dir, exist_ok=True)

# Save the marker genes of all layers in one table
marker_genes_df = pd.concat(deg_dict, names=['annotation_layer']).reset_index(level=0)
marker_genes_df.to_csv(os.path.join(csv_dir, 'ED1G_pathology_marker_genes.csv'), index=False)

# Keep a record of which samples were used
pd.DataFrame({'sample_id': annotated_samples}).to_csv(os.path.join(csv_dir, 'ED1G_samples_included.csv'), index=False)


# --- Load the desired category order per layer

annotation_dir_base = os.path.join(project_dir, 'resources', 'pathology_annotations')
with open(os.path.join(annotation_dir_base, 'pathology_categories_per_layer.json'), 'r') as json_file:
    desired_categories = json.load(json_file)


# --- Count spots per category across all samples

spot_counts_per_category = {}

for layer in annotation_layers:
    spot_counts_per_category[layer] = {}
    for sample in annotation_dict:
        value_counts = annotation_dict[sample][layer].dropna().value_counts()
        for category, count in value_counts.items():
            spot_counts_per_category[layer][category] = spot_counts_per_category[layer].get(category, 0) + count

print(f'Spot counts per category: {spot_counts_per_category}')


# --- GO term enrichment on the up-regulated genes per annotation category

# Load the gene set library from file (a downloaded Enrichr .gmt)
gene_set_dir = os.path.join(project_dir, 'resources', 'gene_lists', 'enrichment_analysis')
go_bio = gp.read_gmt(os.path.join(gene_set_dir, f'{gene_set}.gmt'))

enr_dict_full = {}

for layer in annotation_layers:
    enr_dict_full[layer] = {}
    deg_df = deg_dict[layer]
    if layer not in desired_categories:
        continue
    for group in desired_categories[layer]:
        deg_df_sub = deg_df[deg_df['group'] == group]
        # Significant up-regulated genes (adjusted P < 0.05, log fold change > 0.4)
        degs_sig = deg_df_sub[deg_df_sub['pvals_adj'] < 0.05]
        degs_up = degs_sig[degs_sig['logfoldchanges'] > 0.4]
        if len(degs_up) > 1:
            enr_up = gp.enrichr(degs_up['genes'], gene_sets=go_bio, outdir=None)
            # Strip the trailing GO identifier from the term names
            enr_up.res2d.Term = enr_up.res2d.Term.str.split(r" \(GO").str[0]
            enr_dict_full[layer][group] = enr_up.res2d

# Save the enrichment results of all annotation categories in one table
enr_df = pd.concat(
    {(layer, group): enr_dict_full[layer][group] for layer in enr_dict_full for group in enr_dict_full[layer]},
    names=['annotation_layer', 'group'],
).reset_index(level=[0, 1])
enr_df.to_csv(os.path.join(csv_dir, 'ED1G_pathology_GO_enrichment.csv'), index=False)


# --- Identify the top terms per category and filter

# First pass: collect the union of top N terms across all categories (order preserved)
all_top_terms_ordered = []
all_top_terms_set = set()

for layer in annotation_layers:
    for group in enr_dict_full[layer]:
        df_top = enr_dict_full[layer][group].nsmallest(top_n, 'Adjusted P-value')
        for term in df_top['Term'].tolist():
            if term not in all_top_terms_set:
                all_top_terms_ordered.append(term)
                all_top_terms_set.add(term)

# Second pass: keep only the top terms in each category, in a consistent order
enr_dict = {}

for layer in annotation_layers:
    enr_dict[layer] = {}
    for group in enr_dict_full[layer]:
        df = enr_dict_full[layer][group]
        df_filtered = df[df['Term'].isin(all_top_terms_set)].copy()
        df_filtered['term_order'] = df_filtered['Term'].map({term: i for i, term in enumerate(all_top_terms_ordered)})
        df_filtered = df_filtered.sort_values('term_order').drop('term_order', axis=1)
        enr_dict[layer][group] = df_filtered


# --- Order the categories per layer (only those with enrichment results)

ordered_major_cat_dict = {}

for layer in annotation_layers:
    if layer in desired_categories:
        categories_with_data = [cat for cat in desired_categories[layer] if cat in enr_dict[layer]]
        if categories_with_data:
            ordered_major_cat_dict[layer] = categories_with_data


# --- Colors per annotation layer

color_dir = os.path.join(project_dir, 'resources', 'plotting_colors')
path_color_lut = pd.read_csv(os.path.join(color_dir, 'path_colors.csv'), index_col=0)

top_cat_colors = {layer: path_color_lut.loc[layer, 'color_adj'] for layer in annotation_layers}


# --- Build the plotting dataframes

dot_plot_df_dict = {}
count_df_x_dict = {}
count_df_y_dict = {}

# Order of terms used for the shared axis
all_terms_unique = all_top_terms_ordered

for layer in ordered_major_cat_dict:
    dot_plot_df = pd.DataFrame()
    for group in ordered_major_cat_dict[layer]:
        if group in enr_dict[layer]:
            df = enr_dict[layer][group]
            for idx, row in df.iterrows():
                adj_pval = row['Adjusted P-value']
                plot_row = pd.DataFrame({
                    'category': [group],
                    'term': [row['Term']],
                    'adj_pval': [adj_pval],
                    'neg_log10_pval': [-np.log10(adj_pval)],
                    'combined_score': [row['Combined Score']],
                })
                dot_plot_df = pd.concat([dot_plot_df, plot_row], ignore_index=True)
    # Set term as an ordered categorical to keep the axes aligned across layers
    dot_plot_df['term'] = pd.Categorical(dot_plot_df['term'], categories=all_terms_unique, ordered=True)
    dot_plot_df_dict[layer] = dot_plot_df
    # X-axis: number of terms per category
    count_df_x = dot_plot_df.groupby('category').size().to_frame()
    count_df_x.columns = ['count']
    count_df_x = count_df_x.reset_index()
    count_df_x.index = count_df_x['category']
    count_df_x_dict[layer] = count_df_x
    # Y-axis: how many categories each term appears in
    count_df_y = dot_plot_df.groupby('term').size().to_frame()
    count_df_y.columns = ['count']
    count_df_y_dict[layer] = count_df_y


# --- Plotting helpers

def add_size_legend(dotplot_ax, target_ax, size='proportion', size_title='Proportion (per layer)', limit_circles=False, labelspacing=0.5,
                    handletextpad=0.5, fontsize=8, loc='center', bbox_to_anchor=(0.5, 0.5)):
    # Build a size legend from the seaborn scatter legend handles
    handles, labels = dotplot_ax.get_legend_handles_labels()
    index = labels.index(size)
    size_labels = labels[index + 1:]
    size_handles = handles[index + 1:]
    # Optionally limit to 3 circles (smallest, middle, largest)
    if limit_circles:
        if len(size_handles) > 3:
            size_indices = [0, len(size_handles) // 2, -1]
            size_handles = [size_handles[i] for i in size_indices]
            size_labels = [size_labels[i] for i in size_indices]
    # Make the marker circles gray
    for handle in size_handles:
        handle.set_markerfacecolor('lightgray')
        handle.set_linewidth(0.8)
    # Remove the original legend
    dotplot_ax.legend().remove()
    # Add the new size legend
    legend = target_ax.legend(size_handles, size_labels, loc=loc, labelspacing=labelspacing, bbox_to_anchor=bbox_to_anchor,
                              ncol=1, frameon=False, title=f'{size_title}', handletextpad=handletextpad, fontsize=fontsize,
                              alignment='left')
    legend.get_title().set_fontsize(fontsize)
    legend.get_title().set_fontweight('bold')
    for txt in legend.get_texts():
        txt.set_ha("left")
    return legend

def format_func(value, tick_number):
    if value >= 1000:
        return f'{value/1000:.0f}k'
    return f'{value:.0f}'


# --- Shared plotting parameters

space_per_dot = 0.04
hue = 'neg_log10_pval'
size = 'combined_score'

# Keep the color scale and size scale consistent across all layers
vmin = 0
vmax = max([df['neg_log10_pval'].max() for df in dot_plot_df_dict.values()])
size_min = 0

available_layers_with_data = list(ordered_major_cat_dict.keys())

# Cap the combined score for the dot size
all_combined_scores = []
for layer in dot_plot_df_dict:
    all_combined_scores.extend(dot_plot_df_dict[layer]['combined_score'].tolist())

max_combined_score = max(all_combined_scores) if all_combined_scores else 100

comb_score_limit = 800
if max_combined_score > comb_score_limit:
    max_combined_score = comb_score_limit

dot_plot_df_dict_limited = {}
for layer in dot_plot_df_dict:
    df_limited = dot_plot_df_dict[layer].copy()
    df_limited['combined_score'] = df_limited['combined_score'].apply(lambda x: min(x, comb_score_limit))
    dot_plot_df_dict_limited[layer] = df_limited

caps = ['tls', 'lvi']


# --- Plot 1: enrichment dot plot (categories on x-axis, terms on y-axis, spot-count bars above)

axs, grid_spec = dotplot_functions.prepare_dotplot(ordered_major_cat_dict, height=6, width=6)

for n, layer in enumerate(available_layers_with_data):
    dot_plot_df = dot_plot_df_dict_limited[layer]
    if len(dot_plot_df) == 0:
        continue
    ax = axs[n]
    sns.scatterplot(
        data=dot_plot_df,
        x='category',
        y='term',
        hue='neg_log10_pval',
        size='combined_score',
        sizes=(0, 400),
        palette='Reds',
        edgecolor='black',
        linewidth=0.8,
        ax=ax,
        hue_norm=(vmin, vmax),
        size_norm=(size_min, max_combined_score),
    )
    ax.set_title(layer.capitalize().replace('_', ' '), fontweight='bold', fontsize=10)
    ax.set_xlabel('')
    ax.set_ylabel('')
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')
    # Show all terms on the y-axis (even those absent from this layer)
    ax.set_yticks(range(len(all_terms_unique)))
    ax.set_yticklabels(all_terms_unique)
    # Size the axis based on the number of categories and terms
    dim_x = len(dot_plot_df['category'].unique())
    dim_y = len(all_terms_unique)
    ax_pos = ax.get_position()
    ax.set_position([ax_pos.x0, ax_pos.y0, space_per_dot*dim_x, space_per_dot*dim_y*0.75])
    num_cats_y = len(all_terms_unique)
    num_cats_x = len(dot_plot_df['category'].unique())
    ax.set_xlim(-0.75, num_cats_x-1+0.75)
    ax.set_ylim(num_cats_y-1+0.75, -0.75)
    if ax.get_legend():
        ax.get_legend().remove()

# Move plots so they don't overlap
dotplot_functions.move_plots(axs, offset=0.01, remove_yticks=False)

# Spot-count bars above each dot plot
bar_axs = []

all_counts = []
for layer in available_layers_with_data:
    counts = [spot_counts_per_category[layer].get(cat, 0) for cat in ordered_major_cat_dict[layer]]
    all_counts.extend(counts)
max_count_all = max(all_counts) if all_counts else 0

for n, layer in enumerate(available_layers_with_data):
    ax = axs[n]
    ax_pos = ax.get_position()
    categories_in_plot = list(ordered_major_cat_dict[layer])
    counts = [spot_counts_per_category[layer].get(cat, 0) for cat in categories_in_plot]
    bar_height = 0.08
    bar_ax = plt.gcf().add_axes([ax_pos.x0, ax_pos.y1 + 0.01, ax_pos.width, bar_height])
    bar_ax.set_facecolor('none')
    x_positions = range(len(categories_in_plot))
    color = top_cat_colors.get(layer, 'lightgray')
    bar_ax.bar(x_positions, counts, color=color, edgecolor='black', linewidth=0.8, alpha=0.9)
    for i, count in enumerate(counts):
        label = f'{count/1000:.1f}k' if count >= 1000 else str(count)
        bar_ax.text(i, count, label, ha='center', va='bottom', fontsize=6)
    bar_ax.set_xlim(-0.75, len(categories_in_plot) - 0.25)
    bar_ax.set_ylim(0, max_count_all * 1.1)
    bar_ax.set_xticks([])
    bar_ax.spines['top'].set_visible(False)
    bar_ax.spines['right'].set_visible(False)
    if n != 0:
        bar_ax.set_yticks([bar_ax.get_yticks()[1]])
    if n == 0:
        bar_ax.set_ylabel('# spots', fontsize=8, rotation=90)
        bar_ax.tick_params(axis='y', labelsize=8)
        if max_count_all > 1000:
            bar_ax.yaxis.set_major_formatter(FuncFormatter(format_func))
    else:
        bar_ax.set_yticklabels([])
        bar_ax.tick_params(axis='y', labelsize=8)
    # Move the layer title from the dot plot up onto the bar plot
    bar_ax.set_title(ax.get_title(), fontweight='bold', fontsize=10)
    ax.set_title('')
    bar_axs.append(bar_ax)

for ax in axs:
    ax.grid(True, which='both', axis='both', linestyle='-', linewidth=0.5, alpha=0.5)
    ax.set_axisbelow(True)

# Color bar and size legend on the last plot
if len(available_layers_with_data) > 0:
    last_ax = axs[-1]
    norm = plt.Normalize(vmin, vmax)
    sm = plt.cm.ScalarMappable(cmap='Reds', norm=norm)
    sm.set_array([])
    fig = plt.gcf()
    cbar_height = space_per_dot * 5
    cax = fig.add_axes([last_ax.get_position().x1 + 0.05, last_ax.get_position().y0, 0.02, cbar_height])
    fig.colorbar(sm, cax=cax, orientation='vertical')
    cax.set_ylabel('-log10(Adj. P-value)', fontsize=10, rotation=90, labelpad=5)
    cax.yaxis.set_label_position("left")
    legend = add_size_legend(dotplot_ax=last_ax, target_ax=cax, size='combined_score',
                             size_title='Combined\nscore', limit_circles=True,
                             labelspacing=1, fontsize=10, handletextpad=0.5, loc='lower left', bbox_to_anchor=(-2, 1.1))
    legend.get_title().set_fontsize(10)
    legend.get_title().set_fontweight('normal')

# Clean the x tick labels
for ax in axs:
    new_labels = []
    for label in ax.get_xticklabels():
        text = label.get_text()
        text = text.upper() if text in caps else text.capitalize()
        text = text.replace('_', ' ')
        new_labels.append(text)
    ax.set_xticklabels(new_labels, rotation=90, ha='center')

title_x = (axs[0].get_position().x0 + axs[-1].get_position().x1) / 2
title_y = axs[0].get_position().y1 + 0.2
plt.suptitle(f'Gene enrichment analysis - top {top_n} terms per annotation category (combined score capped at {comb_score_limit})',
             x=title_x, y=title_y, fontweight='bold', fontsize=12)


# --- Save the plots

filename = 'Fig_ED1G.pdf'
save_multi_image(os.path.join(out_dir, filename))
plt.close('all')

print(f"Figure saved to {os.path.join(out_dir, filename)}")
print("Script successfully completed")
