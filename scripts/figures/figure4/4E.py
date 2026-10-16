#!/usr/bin/env python3

# --- Aim
# Plot Figure 4E: adhesion-receptor expression across cell types near HEVs
# Shows the scaled expression of the integrin / SPN / CD38 receptors that bind the
# HEV-expressed adhesion ligands, across HEV-enriched cell types within <50um of HEVs

# --- Define environment
# env: xen-hev

# --- Load packages

import argparse
import os

import anndata as ad
import scanpy as sc

import pandas as pd
import numpy as np
from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import linkage, dendrogram, leaves_list

import seaborn as sns

import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.patches as mpatches

mpl.rcParams['pdf.fonttype'] = 42


# --- Argparse arguments

parser = argparse.ArgumentParser(description="Fig. 4e: adhesion receptor expression across cell types near HEVs",
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


# --- Load the anndata objects

adata_dict = {}

for sample in sample_list:
    adata = ad.read_h5ad(os.path.join(h5ad_dir, 'xenium', f'{sample}.h5ad'))
    adata_dict[sample] = adata


# --- Load the tangram annotations

# Decide which annotation level to use
level = 'subset'

# Initialize a df to store all annotations
anno_df = pd.DataFrame()

for sample in sample_list:
    tmp_df = pd.DataFrame()
    annotation_file = f'{sample}_{level}_pred_cell_types.csv'
    annotation_dir = os.path.join(res_dir, '02_cell_type_mapping', 'label_transfer', 'xenium')
    annotations = pd.read_csv(os.path.join(annotation_dir, f'{level}_level', annotation_file), index_col=0)
    # Add to adata
    adata_dict[sample].obs[level] = annotations['pred_cell_type']
    # Prepend sample id to index
    annotations.index = [f'{sample}_{idx}' for idx in annotations.index]
    tmp_df[level] = annotations['pred_cell_type']
    anno_df = pd.concat([anno_df, tmp_df])

# Add sample_id column
# Split on the last _ and take the first part to get the sample id
anno_df['sample_id'] = [idx.rsplit('_', 1)[0] for idx in anno_df.index]


# --- Load plotting colors

color_dir = os.path.join(project_dir, 'resources', 'plotting_colors')

# Use the major level for colors
color_level = 'major'
color_file = f'miniatlas_{color_level}_colors.csv'

color_dict = pd.read_csv(os.path.join(color_dir, color_file), index_col=0).to_dict()['color']


# --- Load the miniatlas classifications

# These are needed to check the alignment between major and minor levels
miniatlas_dir = os.path.join(project_dir, 'resources', 'cell_types', 'miniatlas')
miniatlas_file = 'miniatlas_subset_classifications.csv'

miniatlas_df = pd.read_csv(os.path.join(miniatlas_dir, miniatlas_file),)

# Add the other levels to the color_dict
for level in ['minor', 'subset']:
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


# --- Load the HEV distances and bin

def bin_distances(df, column_name, bin_size, max_distance):
    # Define the bin edges
    bins = np.arange(0, max_distance + bin_size, bin_size)
    # Use pd.cut to assign each value in the column to a bin
    df['binned'] = pd.cut(df[column_name], bins=bins, right=False, labels=[f'{bins[i]}-{bins[i+1]}' for i in range(len(bins)-1)], include_lowest=True)
    return df['binned']

# Bin distances
bin_size = 25
column_name = 'nearest_distance'
max_distance = 200  # in microns

# Decide whether to look at TLS-HEV only or all HEV
file_extension = 'closest_HEV_distances' # for all HEVs

distance_dir = os.path.join(res_dir, '05_xenium_mechanisms', '03_HEV_microenvironment', '01_calculate_distance')
most_recent = distance_dir

# Load closest distance to any HEV (TLS-associated or not)
for sample in sample_list:
    filename = f'{sample}_{file_extension}.csv'
    tmp_df = pd.read_csv(os.path.join(most_recent, 'csv_files', filename), index_col=0)
    adata_dict[sample].obs['HEV_distance'] = tmp_df['nearest_distance']
    # Also bin the distances
    tmp_df['binned'] = bin_distances(tmp_df, column_name=column_name, bin_size=bin_size, max_distance=max_distance)
    adata_dict[sample].obs['distance_bin'] = tmp_df['binned']


# --- Calculate mean proportions per distance bin per cell type

mean_prop_dict = {}

# Initialize a df to store all distances
dist_df = pd.DataFrame()
for sample in sample_list:
    tmp_df = pd.DataFrame()
    tmp_df['HEV_distance'] = adata_dict[sample].obs['HEV_distance']
    tmp_df['distance_bin'] = adata_dict[sample].obs['distance_bin']
    # Add sample_id column
    tmp_df['sample_id'] = sample
    # And add the cell type annotation
    tmp_df[level] = adata_dict[sample].obs[level]
    dist_df = pd.concat([dist_df, tmp_df])

# Remove any NaN distance bins
dist_df = dist_df[~dist_df['distance_bin'].isna()]

# Convert distance_bin to string
dist_df['distance_bin'] = dist_df['distance_bin'].astype(str)

# Convert distance_bin to ordered categorical
dist_df['distance_bin'] = pd.Categorical(dist_df['distance_bin'], categories=[f'{i}-{i+bin_size}' for i in range(0, max_distance, bin_size)], ordered=True)

# Sort by distance_bin
dist_df = dist_df.sort_values('distance_bin')

# Count number of cells per distance bin per sample
cell_counts = dist_df.groupby(['sample_id', 'distance_bin', level]).size().unstack(fill_value=0)

# Reset multiindex
cell_counts = cell_counts.reset_index()

# Make a cell proportion df
all_cell_types = dist_df[level].unique().tolist()
prop_df = cell_counts.copy()

# Calculate proportions per bin per cell type 
prop_df[all_cell_types] = prop_df[all_cell_types].div(prop_df[all_cell_types].sum(axis=1), axis=0)

# Calculate mean per cell type across samples
mean_prop = prop_df.groupby('distance_bin')[all_cell_types].mean().reset_index()
mean_prop_dict[level] = mean_prop


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

# Now sort minor level by corresponding major level
def get_major(cell_type, level, upper_level='major'):
    level_clean = level.replace('_stepwise', '')
    if cell_type in miniatlas_df[f'miniatlas_{level_clean}'].values:
        major = miniatlas_df[miniatlas_df[f'miniatlas_{level_clean}'] == cell_type][f'miniatlas_{upper_level}'].values[0]
    else:
        major = 'other'
    return major

minor_cell_types_sorted = sorted(minor_cell_types, key=lambda x: (major_cell_types_sorted.index(get_major(x, level='minor')) if get_major(x, level='minor') in major_cell_types_sorted else len(major_cell_types_sorted), x))
subset_cell_types_sorted = sorted(subset_cell_types, key=lambda x: (major_cell_types_sorted.index(get_major(x, level='subset')) if get_major(x, level='subset') in major_cell_types_sorted else len(major_cell_types_sorted), x))

cell_type_orders = {
    'major': major_cell_types_sorted,
    'minor': minor_cell_types_sorted,
    'subset': subset_cell_types_sorted
}


# --- Extract a single adata object with all samples

# Extract a single adata object with all samples
adata_list = list(adata_dict.values())
all_adata = ad.concat(adata_list, merge='same')

all_adata.obs_names_make_unique()


# --- Some functions for plotting

# Add legends
def position_legend(ax):
    """
    Given a single axis, return the bottom, top, left, and right edges
    in figure coordinates, for robust legend placement.
    Returns (bottom_edge, top_edge, left_edge, right_edge)
    """
    fig = ax.figure
    fig.canvas.draw()
    bbox = ax.get_tightbbox(fig.canvas.get_renderer())
    inv_trans = fig.transFigure.inverted()
    bbox_fig = inv_trans.transform(bbox)
    # bbox_fig: [[x0, y0], [x1, y1]]
    left_edge = bbox_fig[0][0]
    right_edge = bbox_fig[1][0]
    bottom_edge = bbox_fig[0][1]
    top_edge = bbox_fig[1][1]
    return bottom_edge, top_edge, left_edge, right_edge

# Function to add legend
def add_color_legend(fig, labels, color_dict, title='', x=0, y=0, offset=0.01, spacing_per_label=0.05, loc='center left', fontsize=10):
    # Create legend handles (colored boxes with labels)
    handles = [
        mpatches.Patch(facecolor=color_dict[label], label=label, edgecolor='black', linewidth=0.8)
        for label in labels if label in color_dict
    ]
    # Estimate height needed
    est_height = spacing_per_label * len(handles)
    bbox_anchor = (x + offset, y)
    # Add the legend
    legend = fig.legend(
        handles=handles,
        title=title,
        loc=loc,
        bbox_to_anchor=bbox_anchor,
        borderaxespad=0,
        frameon=False,
        handlelength=0.8,
        bbox_transform=fig.transFigure,  # <-- This is the key fix!
        fontsize=fontsize,
    )
    legend.get_title().set_fontweight('bold')
    legend.get_title().set_fontsize(fontsize)  # Set the title fontsize
    legend._legend_box.align = "left"
    # Return updated vertical position for next stacked legend
    return legend



# --- Calculate slope based on mean proportions
# To sort the cell types by increasing or decreasing proportions with distance

def calculate_slope(df, cell_type, bin_size):
    # Extract x and y values
    x = np.array([int(bin.split('-')[0]) + bin_size/2 for bin in df['distance_bin']])
    y = df[cell_type].values
    # Fit a linear regression model (1st degree polynomial)
    coeffs = np.polyfit(x, y, deg=1)
    slope = coeffs[0]
    return slope

slope_dict = {}
slope_sorted_cell_types = {}

mean_prop = mean_prop_dict[level]
all_cell_types = [col for col in mean_prop.columns if col != 'distance_bin']

# Calculate slopes for all cell types
slope_dict[level] = {}
for cell_type in all_cell_types:
    slope = calculate_slope(mean_prop, cell_type, bin_size)
    slope_dict[level][cell_type] = slope

# Sort cell types by slope
sorted_cell_types = sorted(slope_dict[level], key=slope_dict[level].get)

# Update the cell type order in adata
all_adata.obs[level] = pd.Categorical(all_adata.obs[level], categories=sorted_cell_types, ordered=True)
slope_sorted_cell_types[level] = sorted_cell_types


# --- Prepare for plotting

# Define a max distance - consider only cells within this distance
max_hev_dist = 50

# Filter the adata to only include cells within this distance
all_adata_filtered = all_adata[all_adata.obs['HEV_distance'] <= max_hev_dist].copy()

all_adata_filtered.obs[level].value_counts()

# Set fontsize
sc.set_figure_params(fontsize=9)

# Disable the grid
plt.rcParams['axes.grid'] = False

# Define shared color scale
vmin = 0
vmax = 1

# Adhesion receptors expressed by immune cells that bind the HEV-expressed adhesion
# ligands, grouped by integrin family. Genes absent from the Xenium panel are dropped.
receptor_families_dict = {
    'β2-integrin\nfamily': ['ITGAL', 'ITGAX', 'ITGB2'],
    'α4-integrin': ['ITGA4', 'ITGA9', 'ITGB1'],
    'Non-integrin': ['SPN', 'CD38'],
}

for family, genes in receptor_families_dict.items():
    genes_present = [g for g in genes if g in all_adata_filtered.var_names]
    if len(genes_present) < len(genes):
        print(f'Warning: Missing receptor genes in {family}: {set(genes) - set(genes_present)}')
    receptor_families_dict[family] = genes_present

receptor_genes = [g for genes in receptor_families_dict.values() for g in genes]

# Genes scored in the pseudobulk matrixplot below
unique_genes = receptor_genes


# --- Define HEV cell types per annotation level

HEV_cell_types_dict = {
    'major': ['Endothelial', 'PVL'],
    'minor': ['Endothelial ACKR1', 'PVL Immature'],
    'subset': ['Endothelial ACKR1', 'PVL Immature s1'],
}


# --- Plot per annotation level

# Remove any stepwise part
level_clean = level.replace('_stepwise', '')

HEV_cell_types = HEV_cell_types_dict[level_clean]

print(f'Plotting for level: {level}')
sorted_cell_types = slope_sorted_cell_types[level]
all_adata_filtered.obs[level] = pd.Categorical(all_adata_filtered.obs[level], categories=sorted_cell_types, ordered=True)


# --- Restrict to the HEV-enriched (negative-slope) cell types

negative_slope_cell_types = [cell for cell in sorted_cell_types if slope_dict[level][cell] < 0]
all_adata_sub = all_adata_filtered[all_adata_filtered.obs[level].isin(negative_slope_cell_types)].copy()

# Pseudobulk scaled-expression matrix (receptors x cell types) via matrixplot
dp = sc.pl.matrixplot(
    all_adata_sub,
    var_names=unique_genes,
    groupby=level,
    cmap='viridis',
    vmin=vmin,
    vmax=vmax,
    show=False,
    swap_axes=True,
    standard_scale='var',
)

heatmap_ax = dp['mainplot_ax']

# Extract values df from plot
values_df = heatmap_ax.collections[0].get_array().reshape(-1, 1)

# Create a df with cell types and genes
genes = [tick.get_text() for tick in heatmap_ax.get_yticklabels()]
cell_types = [tick.get_text() for tick in heatmap_ax.get_xticklabels()]

# values_df is your array, genes and cell_types are your labels
heatmap_values = values_df.flatten()
df = pd.DataFrame(
    heatmap_values.reshape(len(genes), len(cell_types)),
    index=genes,
    columns=cell_types
)

# Order cell types by cumulative scaled receptor expression
sum_df_rec = df.loc[receptor_genes, :].sum(axis=0).sort_values(ascending=False)
sorted_cell_types_by_sum_rec = sum_df_rec.index.tolist()
sum_df_rec = pd.DataFrame(sum_df_rec, columns=['sum'])

# Close the plot
plt.close()


# --- Receptors-only heatmap (Figure 4E)
# The published panel shows only the receptors (integrins / SPN / CD38) that bind the
# HEV-expressed adhesion ligands, not the ligands themselves. Receptor rows are
# hierarchically clustered; cell-type columns are kept in cumulative-receptor-expression
# order, with the HEV cell types fixed on the left.

sum_cmap = 'Purples'
slope_cmap = 'bwr'

# Cluster the receptor rows across all HEV-enriched (negative-slope) cell types
other_cell_types = [ct for ct in sorted_cell_types_by_sum_rec if ct not in HEV_cell_types]
all_neg_slope = HEV_cell_types + other_cell_types

if len(receptor_genes) >= 3 and len(all_neg_slope) >= 2:
    row_matrix = df.loc[receptor_genes, all_neg_slope].values
    row_linkage = linkage(pdist(row_matrix), method='average')
    clustered_genes = [receptor_genes[i] for i in leaves_list(row_linkage)]
else:
    row_linkage = None
    clustered_genes = receptor_genes

# Layout parameters
space_per_dot = 0.04
row_dendro_size = space_per_dot          # row dendrogram width
n_hev = len(HEV_cell_types)
n_genes = len(clustered_genes)
hev_width = n_hev * space_per_dot * 0.7
heatmap_height = n_genes * space_per_dot
hev_x = 0.32
hev_y = 0.15

# Major-lineage colours restricted to the lineages present among the plotted cell types
major_cell_types_sorted_neg = [ct for ct in major_cell_types_sorted
                               if ct in [get_major(cell, level=level) for cell in negative_slope_cell_types]]

fig = plt.figure(figsize=(6, 5))

n_other = len(other_cell_types)
other_width = n_other * space_per_dot * 0.7

# HEV heatmap (fixed on the left)
hev_ax = fig.add_axes([hev_x, hev_y, hev_width, heatmap_height])
sns.heatmap(
    df.loc[clustered_genes, HEV_cell_types],
    ax=hev_ax, cmap='viridis', vmin=vmin, vmax=vmax, cbar=False,
    yticklabels=True, xticklabels=True, linewidths=1, linecolor='gray',
)
hev_ax.set_yticklabels(hev_ax.get_yticklabels(), rotation=0)
hev_ax.set_xticklabels(hev_ax.get_xticklabels(), rotation=90)
for spine in hev_ax.spines.values():
    spine.set_linewidth(1); spine.set_visible(True)

# Row dendrogram — placed to the left of the widest yticklabel of the HEV heatmap
hev_left_edge = position_legend(hev_ax)[2]
row_dendro_ax = fig.add_axes([hev_left_edge - row_dendro_size - 0.005, hev_y, row_dendro_size, heatmap_height])
if row_linkage is not None:
    dendrogram(row_linkage, ax=row_dendro_ax, orientation='left', no_labels=True,
               color_threshold=-1, link_color_func=lambda k: 'black')
    for coll in row_dendro_ax.collections:
        coll.set_linewidth(0.6)
row_dendro_ax.set_xticks([]); row_dendro_ax.set_yticks([])
row_dendro_ax.invert_yaxis()
for spine in row_dendro_ax.spines.values():
    spine.set_visible(False)

# Non-HEV heatmap (columns in gene-sum order, no clustering)
other_x = hev_x + hev_width + 0.01
other_ax = fig.add_axes([other_x, hev_y, other_width, heatmap_height])
sns.heatmap(
    df.loc[clustered_genes, other_cell_types],
    ax=other_ax, cmap='viridis', vmin=vmin, vmax=vmax, cbar=False,
    yticklabels=False, xticklabels=True, linewidths=1, linecolor='gray',
)
other_ax.set_yticks([])
other_ax.set_xticklabels(other_ax.get_xticklabels(), rotation=90)
for spine in other_ax.spines.values():
    spine.set_linewidth(1); spine.set_visible(True)

# Scaled expression colorbar
cbar_ax = fig.add_axes([other_x + other_width + 0.03, hev_y,
                        space_per_dot * 5 * 0.75, space_per_dot * 0.8 * 0.9])
mpl.colorbar.ColorbarBase(cbar_ax, cmap='viridis',
                          norm=mpl.colors.Normalize(vmin=vmin, vmax=vmax), orientation='horizontal')
cbar_ax.set_xlim(vmin, vmax)
cbar_ax.set_xticks([vmin, vmax])
cbar_ax.set_xticklabels([f'{vmin:.1f}', f'{vmax:.1f}'])
cbar_ax.set_title('Scaled expr.\nin group')

# Slope / gene-sum / cell-count annotation rows above the heatmaps
metadata_y0 = hev_y + heatmap_height + 0.01
heatmap_positions = [(hev_x, hev_width, HEV_cell_types, True),
                     (other_x, other_width, other_cell_types, False)]

all_slope_values = list(slope_dict[level].values())
slope_vmax = max(abs(min(all_slope_values)), abs(max(all_slope_values)))

slope_axs = []
for x, width, cts, is_hev in heatmap_positions:
    slope_ax = fig.add_axes([x, metadata_y0, width, 0.04])
    slope_values = [slope_dict[level][ct] for ct in cts]
    slope_plot_df = pd.DataFrame(slope_values, index=cts, columns=['slope'])
    sns.heatmap(slope_plot_df.T, ax=slope_ax, cmap=slope_cmap, cbar=False, annot=False,
                linewidths=1, linecolor='gray', yticklabels=False,
                vmin=-slope_vmax, vmax=slope_vmax)
    slope_ax.set_xticks([])
    if is_hev:
        slope_ax.set_yticks([0.5])
        slope_ax.set_yticklabels(['Slope'], rotation=0, va='center')
    for spine in slope_ax.spines.values():
        spine.set_linewidth(1); spine.set_visible(True)
    slope_axs.append(slope_ax)

cbar_ax_slope = fig.add_axes([slope_axs[-1].get_position().x1 + 0.05, slope_axs[-1].get_position().y0,
                              cbar_ax.get_position().width, cbar_ax.get_position().height])
mpl.colorbar.ColorbarBase(cbar_ax_slope, cmap='bwr',
                          norm=mpl.colors.Normalize(vmin=-slope_vmax, vmax=slope_vmax), orientation='horizontal')
cbar_ax_slope.set_xlim(min(all_slope_values), 0)
cbar_ax_slope.set_xticks([min(all_slope_values), 0])
cbar_ax_slope.set_xticklabels([f'{min(all_slope_values):.1e}', '0'], fontsize=8)

sum_axs = []
for i, (x, width, cts, is_hev) in enumerate(heatmap_positions):
    slope_ax = slope_axs[i]
    sum_ax = fig.add_axes([x, slope_ax.get_position().y1 + 0.01, width, 0.04])
    sum_axs.append(sum_ax)
    plot_sum_df = sum_df_rec.reindex(cts)
    max_sum = sum_df_rec['sum'].max()
    sns.heatmap(plot_sum_df.T, ax=sum_ax, cmap=sum_cmap, vmin=0, vmax=max_sum, cbar=False,
                annot=False, linewidths=1.3, linecolor='gray', yticklabels=False)
    sum_ax.set_xticks([])
    if is_hev:
        sum_ax.set_yticks([0.5])
        sum_ax.set_yticklabels(['Gene sum'], rotation=0, va='center')
    for spine in sum_ax.spines.values():
        spine.set_linewidth(1); spine.set_visible(True)

max_sum_val = sum_df_rec['sum'].max()
cbar_ax_sum = fig.add_axes([sum_axs[-1].get_position().x1 + 0.05, sum_axs[-1].get_position().y0,
                            cbar_ax.get_position().width, cbar_ax.get_position().height])
mpl.colorbar.ColorbarBase(cbar_ax_sum, cmap=sum_cmap,
                          norm=mpl.colors.Normalize(vmin=0, vmax=max_sum_val), orientation='horizontal')
cbar_ax_sum.set_xlim(0, max_sum_val)
cbar_ax_sum.set_xticks([0, max_sum_val])
cbar_ax_sum.set_xticklabels(['0.0', f'{max_sum_val:.1f}'], fontsize=8)
cbar_ax_sum.xaxis.set_ticks_position('top')

bar_axs = []
for i, (x, width, cts, is_hev) in enumerate(heatmap_positions):
    sum_ax = sum_axs[i]
    sum_ax_pos = sum_ax.get_position()
    bar_ax = fig.add_axes([sum_ax_pos.x0, sum_ax_pos.y1 + 0.01, sum_ax_pos.width, sum_ax_pos.height * 3])
    bar_axs.append(bar_ax)
    cell_counts = all_adata_filtered.obs[level].value_counts().reindex(cts).fillna(0).astype(int)
    x_pos = np.arange(len(cts))
    color_list = [color_dict[c] if c in color_dict else 'lightgrey' for c in cts]
    bar_ax.bar(x_pos, cell_counts[cts], color=color_list, linewidth=0.5, edgecolor='black')
    bar_ax.set_xlim(sum_ax.get_xlim()[0] - 0.5, sum_ax.get_xlim()[1] - 0.5)
    bar_ax.spines['top'].set_visible(False)
    bar_ax.spines['right'].set_visible(False)
    bar_ax.spines['left'].set_visible(False)
    bar_ax.set_xticks([]); bar_ax.set_yticks([])
    for rect, count in zip(bar_ax.patches, cell_counts[cts]):
        h = rect.get_height()
        if h > 0:
            bar_ax.text(rect.get_x() + rect.get_width() / 2,
                        h + bar_ax.get_ylim()[1] * 0.02,
                        f'{count/1000:.1f}',
                        ha='center', va='bottom', fontsize=7, rotation=90)

max_ylim = max([b.get_ylim()[1] for b in bar_axs])
for b in bar_axs:
    b.set_ylim(0, max_ylim * 1.3)

bar_axs[-1].set_title(f'Receptor expression <{max_hev_dist}µm from HEVs', fontweight='bold', fontsize=10)
if level_clean in ['subset', 'minor']:
    bar_axs[0].set_title('HEV')
bar_axs[0].spines['top'].set_visible(True)
bar_axs[0].set_ylabel('# cells\n(x10³)', fontsize=8)

# Major-type colour legend
fig.canvas.draw()
right_edges = [position_legend(cb_ax)[3] for cb_ax in [cbar_ax, cbar_ax_slope, cbar_ax_sum]]
x_legend = max(right_edges)
y_legend = sum_axs[-1].get_position().y1

add_color_legend(
    fig, labels=major_cell_types_sorted_neg, color_dict=color_dict, title='',
    x=x_legend, y=y_legend, offset=0.01, spacing_per_label=0.025,
    loc='upper left', fontsize=8,
)

# Replace subset-level cell-type labels with their short names
if 'subset' in level:
    for xtick_ax in [hev_ax, other_ax]:
        new_labels = [miniatlas_sub_dict['miniatlas_subset_short'].get(label.get_text(), label.get_text())
                      for label in xtick_ax.get_xticklabels()]
        xtick_ax.set_xticklabels(new_labels)


# --- Save the output

out_dir = os.path.join(res_dir, 'figures')
os.makedirs(out_dir, exist_ok=True)

filename = 'Fig_4E.pdf'
plt.savefig(os.path.join(out_dir, filename), bbox_inches='tight')


# --- Finish

print("Script successfully completed")
