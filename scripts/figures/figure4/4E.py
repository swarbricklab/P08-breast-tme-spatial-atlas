#!/usr/bin/env python3

# --- Aim
# Plot gene expression of all cell types within x µm of HEVs
# Specifically for LRi analysis - find likely interacting cell types for given LR pairs active in HEV niches (identified in Visium data)
# This script specifically looks at interactions involving cell adhesion molecules

# --- Define environment
# env: xen-hev

# --- Load packages

import argparse
import os
from datetime import datetime

import anndata as ad
import scanpy as sc

import pandas as pd
import numpy as np

import seaborn as sns

import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.patches as mpatches

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
    annotation_dir = os.path.join(project_dir, 'results', '02_cell_type_mapping', 'label_transfer', 'xenium')
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
most_recent = find_most_recent(distance_dir)

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

# Automatic resizing function
def auto_resize_dotplot(ax, space_per_dot=0.1, width_ratio=1):
    """
    Automatically resize the dotplot axes based on the number of genes.
    """
    # Get the current position of the axes
    pos = ax.get_position()
    # Resize based on xticks & yticks
    xticks = len(ax.get_xticks())
    yticks = len(ax.get_yticks())
    # Calculate new width based on number of genes
    new_width = xticks * space_per_dot * width_ratio
    new_height = yticks * space_per_dot
    # Set the new position
    ax.set_position([pos.x0, pos.y0, new_width, new_height])

# Function for repositioning the dotplot axes
def reposition_dotplot_axs(axs, y_offset=0.1):
    for n, ax in enumerate(axs):
        # Get the current position of the axes
        pos = ax.get_position()
        if n == 0:
            ypos = pos.y0  # Keep the first plot at the original position
        else:
            ypos = ax_pos_prev.y0 - pos.height - y_offset  # Offset from the previous plot
        # Set the new position with an offset
        ax.set_position([pos.x0, ypos, pos.width, pos.height])
        ax_pos_prev = ax.get_position()  # Store the position for the next iteration
        # Remove the xticks for all but the last plot
        if n < len(axs) - 1:
            ax.set_xticks([])

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

def join_lr_multi(axs, lr_pairs, offset=0.3, offset_start=0.2):
    """
    Draws brackets connecting ligands and receptors, extending the bracket
    further left if a ligand has multiple receptors or a receptor has
    multiple ligands.

    This version correctly handles labels that appear on multiple axes,
    counting and offsetting them independently per axis.

    Args:
        axs (list | Axes): The matplotlib Axes object(s).
        lr_pairs (list[tuple]): A list of (ligand, receptor) string pairs.
        offset (float): Horizontal distance from the label's left edge to the end of the bracket.
        offset_start (float): Horizontal distance from the longest label to the main vertical line.
    """
    if not isinstance(axs, (list, tuple, np.ndarray)):
        axs = [axs]
    # 1. Build a lookup that can handle labels on MULTIPLE axes.
    # label_to_info['L1'] -> [{'ax': ax1, 'idx': 2, ...}, {'ax': ax2, 'idx': 4, ...}]
    label_to_info = {}
    all_labels = set(lr for pair in lr_pairs for lr in pair)
    if not all_labels:
        return 0
    renderer = axs[0].get_figure().canvas.get_renderer()
    for ax in axs:
        inv = ax.transData.inverted()
        for idx, label in enumerate(ax.get_yticklabels()):
            text = label.get_text()
            if text in all_labels:
                if text not in label_to_info:
                    label_to_info[text] = []
                info = {
                    'ax': ax,
                    'idx': idx,
                    'bbox': label.get_window_extent(renderer),
                    'inv': inv
                }
                label_to_info[text].append(info)
    # 2. Resolve pairs by finding common axes.
    # A pair (L1, R1) is only valid if L1 and R1 are found on the same axis.
    resolved_pairs = []
    for lr1, lr2 in lr_pairs:
        info_list1 = label_to_info.get(lr1, [])
        info_list2 = label_to_info.get(lr2, [])
        # Find all valid (L1, R1) instances on a common axis
        for info1 in info_list1:
            for info2 in info_list2:
                if info1['ax'] == info2['ax']:
                    resolved_pairs.append({
                        'lr1': lr1,
                        'lr2': lr2,
                        'ax': info1['ax'],
                        'info1': info1, # All info for lr1 (ax, idx, bbox, inv)
                        'info2': info2  # All info for lr2 (ax, idx, bbox, inv)
                    })
    if not resolved_pairs:
        return 0 # No valid pairs found
    # 3. Count occurrences PER AXIS, using only resolved pairs.
    ligand_counts_per_ax = {ax: {} for ax in axs}
    receptor_counts_per_ax = {ax: {} for ax in axs}
    for pair in resolved_pairs:
        ax, lr1, lr2 = pair['ax'], pair['lr1'], pair['lr2']
        ligand_counts_per_ax[ax].setdefault(lr1, 0)
        receptor_counts_per_ax[ax].setdefault(lr2, 0)
        ligand_counts_per_ax[ax][lr1] += 1
        receptor_counts_per_ax[ax][lr2] += 1
    # 4. Determine the main vertical alignment line (Global)
    all_x0_coords = []
    for pair in resolved_pairs:
        all_x0_coords.append(pair['info1']['inv'].transform((pair['info1']['bbox'].x0, 0))[0])
        all_x0_coords.append(pair['info2']['inv'].transform((pair['info2']['bbox'].x0, 0))[0])
    if not all_x0_coords:
        return 0
    min_x0 = min(all_x0_coords)
    x_left_align = min_x0 - offset_start - offset
    # 5. Prepare for drawing with PER-AXIS seen counters.
    ligand_seen_per_ax = {ax: {lig: 0 for lig in ax_counts} for ax, ax_counts in ligand_counts_per_ax.items()}
    receptor_seen_per_ax = {ax: {rec: 0 for rec in ax_counts} for ax, ax_counts in receptor_counts_per_ax.items()}
    drawn_verticals_per_ax = {ax: [] for ax in axs}
    indent_step = offset * 2
    # 6. Main drawing loop (using resolved pairs).
    for pair in resolved_pairs:
        # Get all resolved data for this specific pair instance
        ax = pair['ax']
        lr1, lr2 = pair['lr1'], pair['lr2']
        info1, info2 = pair['info1'], pair['info2']
        idx_lig, idx_rec = info1['idx'], info2['idx']
        bbox_lig, bbox_rec = info1['bbox'], info2['bbox']
        inv_lig, inv_rec = info1['inv'], info2['inv'] # inv_lig will equal inv_rec
        # Calculate horizontal end-points
        x_lig_end = inv_lig.transform((bbox_lig.x0, 0))[0] - offset
        x_rec_end = inv_rec.transform((bbox_rec.x0, 0))[0] - offset
        # Get PER-AXIS counts for vertical spread
        n_total_lig = ligand_counts_per_ax[ax].get(lr1, 0)
        n_total_rec = receptor_counts_per_ax[ax].get(lr2, 0)
        n_seen_lig = ligand_seen_per_ax[ax].get(lr1, 0)
        n_seen_rec = receptor_seen_per_ax[ax].get(lr2, 0)
        # Calculate vertical positions (using per-axis counts)
        y_spread_lig, y_spread_rec = 0.3, 0.3
        y_lig_data = idx_lig + 0.5 + ((n_seen_lig - (n_total_lig - 1) / 2) * y_spread_lig if n_total_lig > 1 else 0)
        y_rec_data = idx_rec + 0.5
        if n_total_rec > 1:
            order = (n_total_rec - 1) - n_seen_rec
            y_rec_data += (order - (n_total_rec - 1) / 2) * y_spread_rec
        # Iteratively find a non-overlapping x-position (per-axis check)
        y_span = tuple(sorted((y_lig_data, y_rec_data)))
        current_x = x_left_align
        if n_total_lig > 1 or n_total_rec > 1:
            current_x -= indent_step
        drawn_verticals = drawn_verticals_per_ax[ax] # Check on common axis
        while True:
            has_overlap = any(
                max(y_span[0], s['span'][0]) < min(y_span[1], s['span'][1])
                for s in drawn_verticals if s['x'] == current_x
            )
            if not has_overlap:
                x_left_draw = current_x
                break
            else:
                current_x -= indent_step
        # Draw the bracket on the common axis
        x_points = [x_lig_end, x_left_draw, x_left_draw, x_rec_end]
        y_points = [y_lig_data, y_lig_data, y_rec_data, y_rec_data]
        ax.plot(x_points, y_points, color='black', linewidth=2, zorder=10, clip_on=False, solid_joinstyle='miter')
        # Record the segment that was just drawn for this axis
        drawn_verticals.append({'x': x_left_draw, 'span': y_span})
        # Update seen counters PER AXIS
        ligand_seen_per_ax[ax][lr1] += 1
        receptor_seen_per_ax[ax][lr2] += 1
    # 7. Collapse all drawn verticals to find the furthest point
    all_drawn_verticals = [d for verticals in drawn_verticals_per_ax.values() for d in verticals]
    if all_drawn_verticals:
        furthest_x = min(d['x'] for d in all_drawn_verticals)
    else:
        furthest_x = x_left_align 
    return furthest_x


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

adhesion_genes = ['ICAM1', 'ICAM2', 'VCAM1', 'PECAM1', 'MCAM']

lr_pairs = [('ICAM1', 'ITGAL'), ('ICAM1', 'ITGAX'), ('ICAM1', 'ITGB2'), ('ICAM1', 'SPN'), ('ICAM2', 'ITGAL'), ('ICAM2', 'ITGB2'), 
            ('VCAM1', 'ITGA4'), ('VCAM1', 'ITGA9'), ('VCAM1', 'ITGB1'), ('VCAM1', 'ITGB2'), ('PECAM1', 'CD38'), ('MCAM', 'MCAM')]

# Create a dict with all receptors per adhesion gene
top_genes_dict = {}

for gene in adhesion_genes:
    receptors = [rec for lig, rec in lr_pairs if lig == gene]
    # Sort receptors alphabetically
    receptors = sorted(receptors)
    # If gene is also in receptors, filter
    if gene in receptors:
        receptors = [rec for rec in receptors if rec != gene]
    top_genes_dict[gene] = [gene] + receptors

# Now also sort the lr_pairs based on the order in top_genes_dict
sorted_lr_pairs = []
for gene in adhesion_genes:
    receptors = [rec for lig, rec in lr_pairs if lig == gene]
    # Sort receptors alphabetically
    receptors = sorted(receptors)
    for rec in receptors:
        sorted_lr_pairs.append((gene, rec))

lr_pairs = sorted_lr_pairs

if len(lr_pairs) == 0:
    top_genes_dict = {'Adhesion genes': adhesion_genes}

# Check if genes are in the object
for zone, genes in top_genes_dict.items():
    genes_present = [gene for gene in genes if gene in all_adata_filtered.var_names]
    if len(genes_present) < len(genes):
        missing_genes = set(genes) - set(genes_present)
        print(f'Warning: The following genes are missing from the data and will be skipped in the {zone} plot: {missing_genes}')
    top_genes_dict[zone] = genes_present

# Also identify missing genes in lr_pairs
for lig, rec in lr_pairs:
    if lig not in all_adata_filtered.var_names:
        print(f'Warning: Ligand {lig} is missing from the data and will be skipped in the ligand-receptor brackets.')
        lr_pairs.remove((lig, rec))
    if rec not in all_adata_filtered.var_names:
        print(f'Warning: Receptor {rec} is missing from the data and will be skipped in the ligand-receptor brackets.')
        lr_pairs.remove((lig, rec))

lr_present = []
# Filter LR pairs to only include those present in any list in top_genes_dict
for zone, genes in top_genes_dict.items():
    for lig, rec in lr_pairs:
        if lig in genes and rec in genes:
            lr_present.append((lig, rec))

lr_pairs = lr_present

n_keys = len(top_genes_dict)
total_genes = sum([len(genes) for genes in top_genes_dict.values()])


# --- Function to get heatmap values as a df

def plot_heatmap(adata, ax, gene, annotation_key, vmin=0, vmax=1, last_ax=False, scale='var'):
    dp = sc.pl.matrixplot(
        adata,
        var_names=gene,
        groupby=annotation_key,
        ax=ax,
        cmap='viridis',
        vmin=vmin,
        vmax=vmax,
        show=False,
        swap_axes=True,
        standard_scale=scale,
    )
    # Adjust the size of the axes
    auto_resize_dotplot(dp['mainplot_ax'], space_per_dot=0.04, width_ratio=0.7)
    heatmap_ax = dp['mainplot_ax']
    # Adjust the ylims
    buffer = 0
    dp['mainplot_ax'].set_ylim(len(dp['mainplot_ax'].get_yticks()) + buffer, -buffer)
    # If not the last plot, make all other axes invisible
    if not last_ax:
        for name, dp_ax in dp.items():
            if not name == 'mainplot_ax':
                dp_ax.set_visible(False)
    cbar_ax = dp['color_legend_ax']
    return heatmap_ax, cbar_ax

def get_heatmap_values(all_adata, gene, level, all_bins, sorted_cell_types, mean_prop, nrows, vmin=0, vmax=None, scale='var'):
    heatmap_values = pd.DataFrame(index=sorted_cell_types, columns=all_bins)
    fig, axs = plt.subplots(ncols=1, nrows=nrows, figsize=(5, 5))
    heatmap_axs = []
    for i, n_bin in enumerate(all_bins):
        ax = axs[i]
        # Subset adata to bin 
        adata_sub = all_adata[all_adata.obs['distance_bin'] == n_bin].copy()
        heatmap_ax, cbar_ax = plot_heatmap(adata_sub, ax, gene, annotation_key=level, vmin=vmin, vmax=vmax, 
                                            scale=scale, last_ax=False)
        # Set bin as yticklabel
        heatmap_ax.set_yticklabels([n_bin])
        # Extract zscore df from plot
        zscore_df = heatmap_ax.collections[0].get_array().reshape(-1, 1)
        # Add to heatmap_values df
        heatmap_values[n_bin] = zscore_df.flatten()
        heatmap_ax.set_visible(False)
    # Make long format
    heatmap_values = heatmap_values.reset_index().melt(id_vars='index', var_name='distance_bin', value_name='z_score')
    # Rename the first column
    heatmap_values = heatmap_values.rename(columns={'index': 'cell_type'})
    # Add the proportions to the zscore df
    heatmap_values = heatmap_values.merge(mean_prop.melt(id_vars='distance_bin', var_name='cell_type', value_name='proportion'), on=['distance_bin', 'cell_type'])
    # Close this figure
    plt.close(fig)
    return heatmap_values

# Get the non adhesion molecule genes in top_genes_dict
non_adhesion_genes = []
for genes in top_genes_dict.values():
    for gene in genes:
        if gene not in adhesion_genes:
            non_adhesion_genes.append(gene)

# Get all unique genes in top_genes_dict
unique_genes = set()
for genes in top_genes_dict.values():
    for gene in genes:
        unique_genes.add(gene)

unique_genes = sorted(list(unique_genes))


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


# --- Plot 2: Same as above but only cell types with negative slope

negative_slope_cell_types = [cell for cell in sorted_cell_types if slope_dict[level][cell] < 0]
all_adata_sub = all_adata_filtered[all_adata_filtered.obs[level].isin(negative_slope_cell_types)].copy()

#  Plot 0: Make a heatmap to calculate total scores per cell type
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

# Calculate sum per cell type
sum_df = df.loc[non_adhesion_genes,:].sum(axis=0).sort_values(ascending=False)

# Sort cell types by sum
sorted_cell_types_by_sum = sum_df.index.tolist()

# Reorder the cell types in adata
all_adata_sub.obs[level] = pd.Categorical(all_adata_sub.obs[level], categories=sorted_cell_types_by_sum, ordered=True)

# Close the plot
plt.close()

# -- Make the actual plot

ncols = 2
sum_cmap = 'Purples'

fig, axs = plt.subplots(ncols=ncols, nrows=len(top_genes_dict), figsize=(5, 5))

if n_keys == 1:
    axs = [axs]

heatmap_axs_HEV = []
heatmap_axs = []

# First plot just the HEVs
for i, (zone, genes) in enumerate(top_genes_dict.items()):
    # Use the first column for HEV cell types only
    plot_df = df.loc[genes, HEV_cell_types]
    ax = axs[i][0]
    g = sns.heatmap(
        plot_df,
        ax=ax,
        cmap='viridis',
        vmin=vmin,
        vmax=vmax,
        cbar=False,
        yticklabels=True,
        xticklabels=True,
        linewidths=1,
        linecolor='gray',
    )
    # Adjust the size of the axes
    auto_resize_dotplot(ax, space_per_dot=0.04, width_ratio=0.7)
    # Adjust the ylims
    buffer = 0
    ax.set_ylim(len(ax.get_yticks()) + buffer, -buffer)
    # Set spines visible
    for spine in ax.spines.values():
        spine.set_linewidth(1)  # Change to your preferred thickness
        spine.set_visible(True)  # Change to your preferred thickness
    # Make sure the yticklabels are not rotated
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0)
    # Also make sure the xticklabels are rotated
    ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
    heatmap_axs_HEV.append(ax)

# Then plot the remaining cell types
for i, (zone, genes) in enumerate(top_genes_dict.items()):
    other_cell_types = [ct for ct in sorted_cell_types_by_sum if ct not in HEV_cell_types]
    plot_df = df.loc[genes, other_cell_types]
    ax = axs[i][1]
    g = sns.heatmap(
        plot_df,
        ax=ax,
        cmap='viridis',
        vmin=vmin,
        vmax=vmax,
        cbar=False,
        yticklabels=True,
        xticklabels=True,
        linewidths=1,
        linecolor='gray',
    )
    # Adjust the size of the axes
    auto_resize_dotplot(ax, space_per_dot=0.04, width_ratio=0.7)
    # Adjust the ylims
    buffer = 0
    ax.set_ylim(len(ax.get_yticks()) + buffer, -buffer)
    # Set spines visible
    for spine in ax.spines.values():
        spine.set_linewidth(1)  # Change to your preferred thickness
        spine.set_visible(True)  # Change to your preferred thickness
    # Remove yticks
    ax.set_yticks([])
    # Also make sure the xticklabels are rotated
    ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
    heatmap_axs.append(ax)

# Reposition the dotplot axes
reposition_dotplot_axs(heatmap_axs, y_offset=0.01)
reposition_dotplot_axs(heatmap_axs_HEV, y_offset=0.01)

# Move the heatmap_axs right next to the HEV heatmaps
for n, ax in enumerate(heatmap_axs):
    pos_hev = heatmap_axs_HEV[n].get_position()
    pos = ax.get_position()
    ax.set_position([pos_hev.x1 + 0.01, pos.y0, pos.width, pos.height])

# Add a colorbar ax next to the last heatmap_ax
cbar_ax = fig.add_axes([heatmap_axs[-1].get_position().x1 + 0.03, heatmap_axs[-1].get_position().y0, 0.04*4, 0.03])

# Add colorbar to the new ax
cbar = mpl.colorbar.ColorbarBase(cbar_ax, cmap='viridis', norm=mpl.colors.Normalize(vmin=vmin, vmax=vmax), orientation='horizontal')

# Set xlim to vmin and vmax
cbar_ax.set_xlim(vmin, vmax)
cbar_ax.set_xticks([vmin, vmax])
cbar_ax.set_xticklabels([f'{vmin:.1f}', f'{vmax:.1f}'])
cbar_ax.set_title('Scaled expr.\nin group')

# Position the cbar
ax_pos = heatmap_axs[-1].get_position()
space_per_dot = 0.04
cbar_height = space_per_dot * 0.8
cbar_width = space_per_dot * 5
cbar_ax.set_position([ax_pos.x1 + 0.03, ax_pos.y0, cbar_width*0.75, cbar_height*0.9])

if len(lr_pairs) > 0:
    # Join the ligand-receptor pairs
    x_left_global = join_lr_multi(heatmap_axs_HEV, lr_pairs, offset=0.25)
    # Draw rectangles from top to bottom of each ax
    for n, ax in enumerate(heatmap_axs_HEV):
        # Add vertical text
        label = list(top_genes_dict.keys())[n]
        face_color = 'lightskyblue'
        # Convert x_left_global from data to figure coordinates
        x_fig, _ = fig.transFigure.inverted().transform(ax.transData.transform((x_left_global, 0)))
        pos = ax.get_position()
        y0 = pos.y0
        y1 = pos.y1
        height = abs(y1 - y0)
        width = 0.02  # Width in figure coordinates
        x1 = x_fig - 0.02  # Center the rectangle on x_left_global
        highlight_rect = mpl.patches.Rectangle(
            (x1-width, y0), width, height,
            transform=fig.transFigure,
            facecolor=face_color,
            edgecolor='black',
            linewidth=0,
            alpha=1
        )
        fig.patches.append(highlight_rect)
        # Calculate midpoint y
        y_text = (y0 + y1) / 2
        x_text = x1-width - 0.01
        if 'zone' in label:
            label = label.replace(' zone', '\nzone')
        fig.text(
            x_text, y_text,
            label,
            ha='right', va='center',
            color='black',
        )

# Add axs for slope for each column
slope_cmap = 'bwr'
slope_axs = []

for i in range(ncols):
    ref_ax = heatmap_axs_HEV[0] if i == 0 else heatmap_axs[0]
    # Add an ax for slope
    slope_ax = fig.add_axes([ref_ax.get_position().x0, ref_ax.get_position().y1 + 0.01, ref_ax.get_position().width, 0.04])
    # Get the xticklabels from the last ax to ensure same order
    last_ax = heatmap_axs_HEV[-1] if i == 0 else heatmap_axs[-1]
    cell_type_order = [tick.get_text() for tick in last_ax.get_xticklabels()]
    # Extract all values
    all_values = [v for v in slope_dict[level].values()]
    # Plot slopes as a heatmap
    slope_values = [slope_dict[level][cell_type] for cell_type in cell_type_order]
    slope_plot_df = pd.DataFrame(slope_values, index=cell_type_order, columns=['slope'])
    vmax_cbar = max(abs(min(all_values)), abs(max(all_values)))
    vmin_cbar = -vmax_cbar
    sns.heatmap(slope_plot_df.T, ax=slope_ax, cmap=slope_cmap, cbar=False, annot=False, linewidths=1, linecolor='gray', yticklabels=False,
                vmin=vmin_cbar, vmax=vmax_cbar)
    # Remove xticks
    slope_ax.set_xticks([])
    # Add ytick with slope
    if i == 0:
        slope_ax.set_yticks([0.5])
        slope_ax.set_yticklabels(['Slope'], rotation=0, va='center')
    # Set spines visible
    for spine in slope_ax.spines.values():
        spine.set_linewidth(1)  # Change to your preferred thickness
        spine.set_visible(True)  # Change to your preferred thickness
    if i == ncols - 1:
        # Add a cbar ax next to the slope_ax
        cbar_ax_slope = fig.add_axes([slope_ax.get_position().x1 + 0.05, slope_ax.get_position().y0, cbar_ax.get_position().width, cbar_ax.get_position().height])
        # Add colorbar to the new ax
        cbar_slope = mpl.colorbar.ColorbarBase(cbar_ax_slope, cmap='bwr', norm=mpl.colors.Normalize(vmin=vmin_cbar, vmax=vmax_cbar), orientation='horizontal')
        # Set xlim to actual min and ma
        vmin_cbar = min(all_values)
        vmax_cbar = 0
        cbar_ax_slope.set_xlim(vmin_cbar, vmax_cbar)
        cbar_ax_slope.set_xticks([vmin_cbar, vmax_cbar])
        # Use scientific not
        cbar_ax_slope.set_xticklabels([f'{vmin_cbar:.1e}', f'{vmax_cbar:.1e}'], fontsize=8)
    slope_axs.append(slope_ax)

# Then add an ax showing summed gene values
sum_axs = []

for i in range(ncols):
    slope_ax = slope_axs[i]
    # Add an ax for sum
    sum_ax = fig.add_axes([slope_ax.get_position().x0, slope_ax.get_position().y1 + 0.01, slope_ax.get_position().width, 0.04])
    sum_axs.append(sum_ax)
    sum_df = pd.DataFrame(sum_df, index=sorted_cell_types_by_sum, columns=['sum'])
    # Get the xticklabels from the last ax to ensure same order
    last_ax = heatmap_axs_HEV[-1] if i == 0 else heatmap_axs[-1]
    cell_type_order = [tick.get_text() for tick in last_ax.get_xticklabels()]
    plot_sum_df = sum_df.reindex(cell_type_order)
    max_sum = sum_df['sum'].max()
    # Add heatmap for summed gene values
    sns.heatmap(plot_sum_df.T, ax=sum_ax, cmap=sum_cmap, vmin=0, vmax=max_sum, cbar=False, annot=False, linewidths=1.3, linecolor='gray', yticklabels=False)
    # Remove xticks
    sum_ax.set_xticks([])
    # Add ytick with 'Sum'
    if i == 0:
        sum_ax.set_yticks([0.5])
        sum_ax.set_yticklabels(['Gene sum'], rotation=0, va='center')
    # Set spines visible
    for spine in sum_ax.spines.values():
        spine.set_linewidth(1)  # Change to your preferred thickness
        spine.set_visible(True)  # Change to your preferred thickness
    if i == ncols - 1:
        # Add cbar ax next to the sum_ax
        cbar_ax_sum = fig.add_axes([sum_ax.get_position().x1 + 0.05, sum_ax.get_position().y0, cbar_ax.get_position().width, cbar_ax.get_position().height])
        # Add colorbar to the new ax
        vmax_cbar = sum_df['sum'].max()
        cbar_sum = mpl.colorbar.ColorbarBase(cbar_ax_sum, cmap=sum_cmap, norm=mpl.colors.Normalize(vmin=0, vmax=vmax_cbar), orientation='horizontal')
        # Set xlim to actual min and ma
        vmin_cbar = 0
        cbar_ax_sum.set_xlim(vmin_cbar, vmax_cbar)
        cbar_ax_sum.set_xticks([vmin_cbar, vmax_cbar])
        # Round to first decimal
        cbar_ax_sum.set_xticklabels([f'{vmin_cbar:.1f}', f'{vmax_cbar:.1f}'], fontsize=8)
        # Set the x axis at the top
        cbar_ax_sum.xaxis.set_ticks_position('top')

# Add barplot showing total cell numbers on top of top heatmap ax
bar_axs = []

for i in range(ncols):
    top_ax = sum_axs[i]
    bar_ax = plt.axes([0.1, 0.1, 0.2, 0.2])
    bar_axs.append(bar_ax)
    top_ax_pos = top_ax.get_position()
    bar_ax.set_position([top_ax_pos.x0, top_ax_pos.y1 + 0.01, top_ax_pos.width, top_ax_pos.height*3])
    # Get the xticklabels from the last ax to ensure same order
    last_ax = heatmap_axs_HEV[-1] if i == 0 else heatmap_axs[-1]
    cell_type_order = [tick.get_text() for tick in last_ax.get_xticklabels()]    # Get cell counts
    cell_counts = all_adata_filtered.obs[level].value_counts().reindex(cell_type_order).fillna(0).astype(int)
    # Plot bars
    x = np.arange(len(cell_type_order))
    # Extract colors for the sorted cells
    color_list = [color_dict[cell] if cell in color_dict else 'lightgrey' for cell in cell_type_order]
    bar_ax.bar(x, cell_counts[cell_type_order], color=color_list, linewidth=0.5, edgecolor='black')
    # Set xlim same as heatmap
    bar_ax.set_xlim(top_ax.get_xlim()[0] - 0.5, top_ax.get_xlim()[1] - 0.5)
    # Set spines inivisible
    bar_ax.spines['top'].set_visible(False)
    bar_ax.spines['right'].set_visible(False)
    bar_ax.spines['left'].set_visible(False)
    # Remove xticks & yticks
    bar_ax.set_xticks([])
    bar_ax.set_yticks([])
    # Annotate the bars with cell counts
    for rect, count in zip(bar_ax.patches, cell_counts[cell_type_order]):
        height = rect.get_height()
        if height > 0:
            # Divide by 1000 
            # Round to 1 decimal
            count = f'{count/1000:.1f}'
            bar_ax.text(
                rect.get_x() + rect.get_width() / 2, height + bar_ax.get_ylim()[1]*0.02,
                f'{count}',
                ha='center', va='bottom',
                fontsize=7, rotation=90
            )

# Set equal ylim for both bar_axs
max_ylim = max([bar_ax.get_ylim()[1] for bar_ax in bar_axs])
for bar_ax in bar_axs:
    bar_ax.set_ylim(0, max_ylim)

# Add a title
bar_ax.set_title(f'Gene expression <{max_hev_dist}µm from HEVs', fontweight='bold', fontsize=10)
if level_clean == 'subset':
    bar_axs[0].set_title(f'HEV',)
# Set top ax spien visible
bar_axs[0].spines['top'].set_visible(True)

# Add a color legend
bottom_edge, top_edge, left_edge, right_edge = position_legend(cbar_ax)
x_legend = left_edge

bottom_edge, top_edge, left_edge, right_edge = position_legend(heatmap_axs[0])
y_legend = top_edge - 0.03

# Add legend only if minimum number of genes 
# Filter major_cell_types_sorted to only those present in negative_slope_cell_types
major_cell_types_sorted_neg = [ct for ct in major_cell_types_sorted if ct in [get_major(cell, level=level) for cell in negative_slope_cell_types]]

if total_genes >= 7:
    add_color_legend(
        fig,
        labels=major_cell_types_sorted_neg,
        color_dict=color_dict,
        title='',
        x=x_legend,
        y=y_legend,
        offset=0.01,
        spacing_per_label=0.025,
        loc='upper left',
        fontsize=8
    )

# Add ylabel
bar_axs[0].set_ylabel('# cells\n(x10³)', fontsize=8)

# Adjust the ylim to make space for the annotations
for bar_ax in bar_axs:
    bar_ax.set_ylim(0, bar_ax.get_ylim()[1]*1.3)

# If subset level, replace xticklabels of the last heatmap ax by the short names
if 'subset' in level:
    last_heatmap_ax = heatmap_axs[-1]
    xticklabels = last_heatmap_ax.get_xticklabels()
    new_labels = []
    for label in xticklabels:
        cell_type = label.get_text()
        if cell_type in miniatlas_sub_dict['miniatlas_subset_short']:
            new_label = miniatlas_sub_dict['miniatlas_subset_short'][cell_type]
        else:
            new_label = cell_type
        new_labels.append(new_label)
    last_heatmap_ax.set_xticklabels(new_labels)
    # Also for HEV heatmap
    last_heatmap_ax_hev = heatmap_axs_HEV[-1]
    xticklabels_hev = last_heatmap_ax_hev.get_xticklabels()
    new_labels_hev = []
    for label in xticklabels_hev:
        cell_type = label.get_text()
        if cell_type in miniatlas_sub_dict['miniatlas_subset_short']:
            new_label = miniatlas_sub_dict['miniatlas_subset_short'][cell_type]
        else:
            new_label = cell_type
        new_labels_hev.append(new_label)
    last_heatmap_ax_hev.set_xticklabels(new_labels_hev)

# Position the title
bar_ax_pos = bar_axs[-1].get_position()
title_x = (bar_ax_pos.x0 + bar_ax_pos.x1) / 2
title_y = bar_ax_pos.y1 + 0.1

plt.suptitle(f'All samples - Gene expression near (<{max_hev_dist}µm) HEVs', x=title_x, y=title_y, fontweight='bold')


# --- Save the output

out_dir = os.path.join(res_dir, 'figures')
os.makedirs(out_dir, exist_ok=True)

filename = 'Fig_4E.pdf'
plt.savefig(os.path.join(out_dir, filename), bbox_inches='tight')


# --- Finish

print("Script successfully completed")