#!/usr/bin/env python3

# --- Aim
# Investigate gene expression correlations with SPP1 or TMEM173 (STING) expression in HEV niches.
# Identify gene-gene relationships.

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

from scipy.stats import pearsonr
from scipy.spatial import cKDTree

import scanpy as sc

import pandas as pd
import numpy as np

import seaborn as sns

# Import custom functions
from sophiesfunctions.save_multi_image import save_multi_image

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
sys.path.append(os.path.join(project_dir, 'scripts', '05_xenium_mechanisms'))
import plotting_functions


# --- Extract some info from the sample sheet

# Use sample_sheet argument if provided, else default
if args.get("sample_sheet"):
    sample_sheet_path = args["sample_sheet"]
else:
    sample_dir = os.path.join(project_dir, 'config')
    sample_sheet_file = 'sample_sheet.csv'
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


# --- Load the HEV annotations

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

distance_dir = os.path.join(res_dir, '05_xenium_mechanisms', '02_identify_HEV_regions',)

most_recent_distance = find_most_recent(distance_dir)

for sample in sample_list:
    # Also load HEV niche barcode and add to adata
    filename = f'{sample}_HEV_niches_cell_to_niche.csv'
    niche_barcode_df = pd.read_csv(os.path.join(most_recent_distance, 'csv_files', filename), index_col=0)
    niche_barcode_df.dropna(subset=['niche_id'], inplace=True)
    # Prepend sample name to niche_id
    niche_barcode_df['niche_id'] = sample + '_' + niche_barcode_df['niche_id'].astype(str)
    adata_dict[sample].obs['niche_id'] = niche_barcode_df['niche_id']
    # Also load the distance calculations
    distance_file = f'{sample}_HEV_niches_info.csv'
    distance_df = pd.read_csv(os.path.join(most_recent_distance, 'csv_files', distance_file))
    # Prepend sample name to niche_id
    distance_df['niche_id'] = sample + '_' + distance_df['niche_id'].astype(str)
    # Make a mapping dict of distance to TLS
    distance_dict = distance_df.set_index('niche_id')['close_to_TLS'].to_dict()
    # Map to adata
    adata_dict[sample].obs['close_to_TLS'] = adata_dict[sample].obs['niche_id'].map(distance_dict)


# --- Also load cell annotations

# Decide which annotation levels to use
level = 'subset'

for sample in sample_list:
    annotation_dir = os.path.join(project_dir, 'results', '02_cell_type_mapping', 'label_transfer', 'xenium', f'{level}_level')
    annotation_file = f'{sample}_{level}_pred_cell_types.csv'
    annotations = pd.read_csv(os.path.join(annotation_dir, annotation_file), index_col=0)
    # Add to adata
    adata_dict[sample].obs[level] = annotations['pred_cell_type']


# --- Look at gene expression in HEV niches

adata_dict_HEV = {}

# Subset adatas to HEV only
for sample in sample_list:
    adata = adata_dict[sample]
    adata_HEV = adata[adata.obs['niche_id'].notna()].copy()
    adata_dict_HEV[sample] = adata_HEV

# Extract a single adata object with all samples
adata_list = list(adata_dict_HEV.values())
all_adata = ad.concat(adata_list, merge='same')


# --- Create an adata with mean expression per niche

mean_adata = all_adata.to_df().groupby(all_adata.obs['niche_id']).mean()

mean_adata = ad.AnnData(mean_adata)

mean_adata.obs['niche_id'] = mean_adata.obs.index


# --- Calculate correlations between STING expression and all other genes

# Pick a gene of interest
gene_of_i = 'TMEM173'  # STING gene symbol

correlation_df = pd.DataFrame(columns=['gene', 'correlation', 'p_value'])

for gene in mean_adata.var_names:
    if gene == gene_of_i:
        continue
    # Calculate Pearson correlation
    corr, p_val = pearsonr(mean_adata[:, gene_of_i].X.flatten(), mean_adata[:, gene].X.flatten())
    correlation_df = correlation_df._append({'gene': gene, 'correlation': corr, 'p_value': p_val}, ignore_index=True)

# Sort by correlation
correlation_df = correlation_df.sort_values(by='correlation', ascending=False).reset_index(drop=True)

# Drop any NaN values
correlation_df = correlation_df.dropna(subset=['correlation', 'p_value'])
correlation_df = correlation_df.set_index('gene')

# Add rank
correlation_df['rank'] = range(1, len(correlation_df) + 1)

# Filter significant correlations
pval_threshold = 0.05
significant_corr_df = correlation_df[correlation_df['p_value'] < pval_threshold].copy()
correlation_df['significant'] = correlation_df['p_value'] < pval_threshold


# --- Plot 1: Plot the correlations

# Add grid behind
plt.rcParams['axes.axisbelow'] = True

fig, ax = plt.subplots(figsize=(8, 4))

sns.scatterplot(data=correlation_df, x='rank', y='correlation', color='lightgrey', edgecolor=None, s=2, ax=ax)
sns.scatterplot(data=significant_corr_df, x='rank', y='correlation', color='red', edgecolor=None, s=2, ax=ax)

# Annotate some genes
annotate_genes = ['SELE', 'SELP', 'CXCL10', 'CXCL9', 'CCL5', 'CCL2', 'CCL19', 'CCL21', 'CXCL13',
                    'MCAM', 'ICAM1', 'ICAM2', 'VCAM1', 'BCAM', 'PECAM1']

annotated_positions = []
min_dist = 15  # minimum pixel distance between annotations

for gene in annotate_genes:
    if gene not in correlation_df.index:
        continue
    x = correlation_df.loc[gene, 'rank']
    y = correlation_df.loc[gene, 'correlation']
    xytext = np.array([5, 5])
    # Check for overlap with previous annotations
    for prev_x, prev_y, prev_xytext in annotated_positions:
        # Calculate screen position of previous annotation
        prev_screen = ax.transData.transform((prev_x, prev_y)) + prev_xytext
        curr_screen = ax.transData.transform((x, y)) + xytext
        dist = np.linalg.norm(curr_screen - prev_screen)
        if dist < min_dist:
            # If too close, shift the annotation further
            xytext += np.array([min_dist, min_dist])
    ax.annotate(
        gene,
        xy=(x, y),
        xytext=tuple(xytext),
        textcoords='offset points',
        fontsize=8,
        arrowprops=dict(arrowstyle='->', color='black', lw=0.7)
    )
    annotated_positions.append((x, y, xytext))

# Add grid behind plot
ax.grid(which='major', linestyle='-', color='lightgrey')

# Also add a line at the middle rank
mid_rank = len(correlation_df) / 2
ax.axvline(mid_rank, color='lightgray', lw=3, ls='-', zorder=-1)

# Adjust labels
ax.set_xlabel(f'Gene rank (by correlation with {gene_of_i} expression)',)
ax.set_ylabel('Pearson correlation coefficient', )

# Add a legend for pvalue
ax.legend(title=f'p<{pval_threshold}', labels=['False', 'True'], loc='center left',  frameon=False, markerscale=3, bbox_to_anchor=(1, 0.5))

# Add a title
ax.set_title(f'Correlation mean {gene_of_i} expression & other genes per HEV', fontweight='bold')

plt.suptitle(f'All samples - Correlation between mean {gene_of_i} expression and other genes per HEV niche', fontweight='bold', y=1.02)


# --- Find cells within 30 microns of HEV niches

def find_cells_within_distance_min(df, group_key="HEV_group", x_col="x", y_col="y", max_distance=200):
    """
    For each HEV group, find all cell barcodes within max_distance (µm) of any cell assigned to that group.
    The distance recorded is the minimum distance from any cell in the HEV niche.
    Returns a DataFrame with columns: cell_barcode, HEV_group, distance
    """
    results = []
    # Get all cells assigned to a HEV group (the niche members)
    group_df = df.dropna(subset=[group_key])
    # For each HEV group
    for group in group_df[group_key].unique():
        niche_cells = group_df[group_df[group_key] == group]
        niche_coords = niche_cells[[x_col, y_col]].values
        # Build KDTree for the niche
        tree = cKDTree(niche_coords)
        # For all cells, query the minimum distance to any cell in the niche
        all_coords = df[[x_col, y_col]].values
        min_distances, _ = tree.query(all_coords)
        # Find all cells within max_distance
        within_mask = min_distances <= max_distance
        for idx in np.where(within_mask)[0]:
            cell_barcode = df.index[idx]
            results.append({
                "cell_barcode": cell_barcode,
                "HEV_group": group,
                "distance": min_distances[idx]
            })
    return pd.DataFrame(results)

# Calculate for a variety of thresholds
mean_sting_per_niche = mean_adata.to_df()[[gene_of_i]]

dists = [0, 10, 20, 30, 40, 50, 100, 200]

results_dict = {}

for max_dist in dists:
    print(f'Calculating for max distance {max_dist} µm')
    group_df_list = []
    for sample in sample_list:
        adata = adata_dict[sample]
        coords = adata.obs[['x_centroid', 'y_centroid']].copy()
        coords.columns = ['x', 'y']
        groups = adata.obs['niche_id']
        df = pd.DataFrame(coords, columns=["x", "y"])
        df["HEV_group"] = groups.values
        df_results = find_cells_within_distance_min(df, group_key="HEV_group", x_col="x", y_col="y", max_distance=max_dist)
        # Precompute cell type annotations
        cell_type_ann = adata.obs[level]
        niche_barcodes_dict = adata.obs.groupby('niche_id').groups
        for HEV_group in adata.obs['niche_id'].dropna().unique():
            group_mask = df_results['HEV_group'] == HEV_group
            cell_barcodes_in_group = set(df_results.loc[group_mask, 'cell_barcode'])
            # Add niche barcodes
            cell_barcodes_in_group.update(niche_barcodes_dict.get(HEV_group, []))
            # Extract cell types
            cell_types = cell_type_ann.loc[list(cell_barcodes_in_group)]
            cell_type_counts = cell_types.value_counts()
            cell_type_proportions = cell_type_counts / cell_type_counts.sum()
            group_df_list.append(pd.DataFrame({
                'niche_id': HEV_group,
                'cell_type': cell_type_proportions.index,
                'proportion': cell_type_proportions.values
            }))
    cell_proportions_df = pd.concat(group_df_list, ignore_index=True)
    # Pivot the df to have cell types as columns and niche id as rows
    cell_type_pivot_df = cell_proportions_df.pivot(index='niche_id', columns='cell_type', values='proportion').fillna(0)
    # --- Calculate correlation with cell type composition per niche
    # Make sure indexes match
    cell_type_pivot_df = cell_type_pivot_df.loc[mean_sting_per_niche.index]
    # Calculate correlations between STING expression and cell type proportions
    correlation_ct_df = pd.DataFrame(columns=['cell_type', 'correlation', 'p_value'])
    for cell_type in cell_type_pivot_df.columns:
        # Calculate Pearson correlation
        corr, p_val = pearsonr(mean_sting_per_niche[gene_of_i], cell_type_pivot_df[cell_type].values)
        correlation_ct_df = correlation_ct_df._append({'cell_type': cell_type, 'correlation': corr, 'p_value': p_val}, ignore_index=True)
    # Sort by correlation
    correlation_ct_df = correlation_ct_df.sort_values(by='correlation', ascending=False).reset_index(drop=True)
    # Add column for significance
    correlation_ct_df['significant'] = correlation_ct_df['p_value'] <= pval_threshold
    results_dict[max_dist] = correlation_ct_df

# Extract a single df with all results
all_results = []

for max_dist, df in results_dict.items():
    df_copy = df.copy()
    df_copy['max_dist'] = max_dist
    all_results.append(df_copy)

distance_corr = pd.concat(all_results, ignore_index=True)

# Make a df for plotting with distance as rows and cell types as columns, values are correlation coefficients
plot_df = distance_corr.pivot(index='max_dist', columns='cell_type', values='correlation')


# --- Prepare for plotting heatmap

color_dir = os.path.join(project_dir, 'resources', 'plotting_colors')

# Use the major level for colors
color_level = 'major'
color_file = f'miniatlas_{color_level}_colors.csv'

color_dict_cells = pd.read_csv(os.path.join(color_dir, color_file), index_col=0).to_dict()['color']

# Load the miniatlas classifications
miniatlas_dir = os.path.join(project_dir, 'resources', 'cell_types', 'miniatlas')
miniatlas_file = 'miniatlas_subset_classifications.csv'
miniatlas_df = pd.read_csv(os.path.join(miniatlas_dir, miniatlas_file),)

# For each cell type in the miniatlas at this level, check if it is in the color_dict
# If not, get the corresponding major cell type and use its color
for cell_type in miniatlas_df[f'miniatlas_{level}'].unique():
    if cell_type not in color_dict_cells.keys():
        if pd.isna(cell_type):
            continue
        print(f'Adding color for {cell_type}')
    major_type = miniatlas_df[miniatlas_df[f'miniatlas_{level}'] == cell_type][f'miniatlas_major'].values
    major_type_in_dict = [m for m in major_type if m in color_dict_cells.keys()]
    if len(major_type_in_dict) > 0:
        color_dict_cells[cell_type] = color_dict_cells[major_type_in_dict[0]]
    else:
        color_dict_cells[cell_type] = 'lightgrey'

# Also define a dict with shorter names for the subset level cell types
miniatlas_sub_dict = miniatlas_df.set_index('miniatlas_subset')[['miniatlas_subset_short']].to_dict()

# Make a row colors df
cell_types_df = pd.DataFrame({'cell_type': distance_corr['cell_type'].unique()})
cell_types_df['major'] = cell_types_df['cell_type'].apply(lambda x: miniatlas_df[miniatlas_df[f'miniatlas_{level}'] == x]['miniatlas_major'].values[0] if x in miniatlas_df[f'miniatlas_{level}'].values else 'other')
row_colors = cell_types_df['major'].map(color_dict_cells).tolist()
cell_types_df = cell_types_df.set_index('cell_type')
row_colors = pd.Series(row_colors, index=cell_types_df.index)
row_colors.name = 'Major type'

# Sort by lineage, then by major
lineage_order = ['Epithelial', 'Immune', 'Stromal']
major_cell_types = miniatlas_df[f'miniatlas_major'].dropna().unique().tolist()
minor_cell_types = miniatlas_df[f'miniatlas_minor'].dropna().unique().tolist()
subset_cell_types = miniatlas_df[f'miniatlas_subset'].dropna().unique().tolist()

def get_lineage(cell_type, level='major'):
    if cell_type in miniatlas_df[f'miniatlas_{level}'].values:
        lineage = miniatlas_df[miniatlas_df[f'miniatlas_{level}'] == cell_type]['lineage'].values[0]
    else:
        lineage = 'other'
    return lineage

major_cell_types_sorted = sorted(major_cell_types, key=lambda x: (lineage_order.index(get_lineage(x)) if get_lineage(x) in lineage_order else len(lineage_order), x))

def position_legend(axs):
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
        left_edges.append(bbox_fig[0][0])
        bottom_edges.append(bbox_fig[0][1])
        right_edges.append(bbox_fig[1][0])
        top_edges.append(bbox_fig[1][1])
    max_right_edge = max(right_edges)
    max_top_edge = max(top_edges)
    max_left_edge = max(left_edges)
    max_bottom_edge = max(bottom_edges)
    return max_right_edge, max_top_edge, max_left_edge, max_bottom_edge


# --- Plot 2: Heatmap of correlations between STING expression and cell type proportions at different distances

# Drop columns with NaN values
plot_df = plot_df.dropna(axis=1, how='any')

# Make the plot
g = sns.clustermap(plot_df, row_cluster=False, col_cluster=True, cmap='RdBu_r', linewidths=0.5, figsize=(4, 4), annot=False,
            colors_ratio=(0.01, 0.1), linecolor='white', xticklabels=True, yticklabels=True, fmt='g', col_colors=row_colors,
            center=0, cbar_kws={'orientation': 'vertical'})

axs = plotting_functions.beautify_heatmap(g, row_colors=False, col_colors=True, cbar_title=f'Pearson r', label_size=10, fontsize_title=11, title='',
                cbar_pad=0.1, cbar_height_square=4, cbar_width_square=0.8, space_per_square=0.05, width_to_height_ratio=1,
                colors_div=2,)

g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0)
ytick_texts = [tick.get_text() for tick in g.ax_heatmap.get_yticklabels()]
xtick_texts = [tick.get_text() for tick in g.ax_heatmap.get_xticklabels()]

# Order plot df according to plot order
plot_df.index = plot_df.index.astype(str)
plot_df.columns = plot_df.columns.astype(str)
plot_df_reordered = plot_df.loc[ytick_texts, xtick_texts]

# Adjust labels
heatmap_ax = g.ax_heatmap
heatmap_ax.set_xlabel('Cell type', fontweight='bold')
heatmap_ax.set_ylabel('HEV distance\n(µm)', fontweight='bold')
new_labels = [miniatlas_sub_dict['miniatlas_subset_short'].get(label.get_text(), label.get_text()) for label in g.ax_heatmap.get_xticklabels()]
g.ax_heatmap.set_xticklabels(new_labels, rotation=90)

# Adjust colorbar position
plotting_functions.position_cbar(g.ax_heatmap, g.ax_cbar, offset=0.1)

# Adjust dendrograms
plotting_functions.adjust_dendrograms(g, num_square=2, title=f'HEV {gene_of_i} expression vs cell type abundance')

g.ax_col_dendrogram.set_title(g.ax_col_dendrogram.get_title(), fontsize=12, fontweight='bold')
cbar_lim = g.ax_cbar.get_ylim()
max_val = max(abs(cbar_lim[0]), abs(cbar_lim[1]))
center = 0

for i in range(plot_df_reordered.shape[0]):
    for j in range(plot_df_reordered.shape[1]):
        p_val = distance_corr[
            (distance_corr['cell_type'] == plot_df_reordered.columns[j]) &
            (distance_corr['max_dist'].astype(str) == plot_df_reordered.index[i])
        ]['p_value'].values[0]
        if p_val < 0.05:
            arr = heatmap_ax.collections[0].get_array()
            square_color = arr[i, j] if arr.ndim == 2 else arr[i * plot_df_reordered.shape[1] + j]
            if hasattr(square_color, 'item'):
                square_color = square_color.item()
            extreme_distance = abs(square_color - center) / max_val
            if extreme_distance > 0.6:
                text_color = 'white'
            else:
                text_color = 'black'
            heatmap_ax.text(j + 0.5, i + 0.7, '*', color=text_color, ha='center', va='center', fontsize=12, fontweight='bold')

# Add a legend
fig = plt.gcf()
cbar_ax = g.ax_cbar
right_edge, top_edge, left_edge, bottom_edge = position_legend([cbar_ax])
plotting_functions.add_color_legend(fig, labels=major_cell_types_sorted, color_dict=color_dict_cells, title='Major cell type',
                                x=right_edge+0.05, y=bottom_edge, offset=0, spacing_per_label=0.03, loc='lower left')
col_colors_ax = g.ax_col_colors
ax_pos = col_colors_ax.get_position()
title_x = (ax_pos.x0 + ax_pos.x1) / 2
title_y = ax_pos.y1 + 0.3
plt.suptitle(f'All samples - Correlations between HEV {gene_of_i} expression and cell type abundance (within x µm of HEV, asterisk = p<0.05)', x=title_x, y=title_y, fontweight='bold')


# --- Save the output

out_dir = os.path.join(res_dir, '05_xenium_mechanisms', '05_gene_gene_correlations', f'{gene_of_i}_correlations')

today = datetime.now().strftime('%Y%m%d')
out_dir_plots = os.path.join(out_dir, today, 'plots')
os.makedirs(out_dir_plots, exist_ok=True)

# Save the figures
filename = f'All_samples_HEV_{gene_of_i}_gene_correlations.pdf'
save_multi_image(os.path.join(out_dir_plots, filename))

# Also save the correlation data
out_dir_csv = os.path.join(out_dir, today, 'csv_files')
os.makedirs(out_dir_csv, exist_ok=True)

correlation_df.to_csv(os.path.join(out_dir_csv, f'All_samples_HEV_{gene_of_i}_gene_correlations.csv'))
distance_corr.to_csv(os.path.join(out_dir_csv, f'All_samples_HEV_{gene_of_i}_cell_type_correlations_at_distances.csv'), index=False)


# --- Finish

print("Script successfully completed")