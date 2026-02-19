#!/usr/bin/env python3

# --- Aim
# Classify vascular structures based on dominant endothelial cell type from deconvolution results.
# Annotate vessel types for downstream spatial and compositional analyses.

# --- Define environment
# env: vis-hd

# --- Load packages
import argparse
import math
import os
from datetime import date, datetime

import anndata as ad
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
import squidpy as sq
import seaborn as sns
from matplotlib.colors import ListedColormap
from matplotlib.colors import to_rgb
import itertools
import re
import numpy as np
import sys

# Import custom functions
from sophiesfunctions.auto_crop import auto_crop_scanpy
from sophiesfunctions.save_multi_image import save_multi_image

sc.set_figure_params(facecolor="white", figsize=(8, 8))
sc.settings.verbosity = 3

mpl.rcParams['pdf.fonttype'] = 42


# --- Argparse arguments

parser = argparse.ArgumentParser(description="Import Visium data & do basic QC",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

# Required
parser.add_argument("--res_dir", help="Path to results directory")
parser.add_argument("--h5ad_dir", help="Path to directory for .h5ad output files")
parser.add_argument("--project_dir", help="Project directory (root)")
parser.add_argument("--sample_id", help="Sample ID")

parser.add_argument("--min_count", help="Min count for QC", type=int)
parser.add_argument("--min_gene", help="Min gene for QC", type=int)
parser.add_argument("--min_spots", help="Minimal spots per gene", type=int)

args = vars(parser.parse_args())


# --- Arguments

res_dir = args["res_dir"]
h5ad_dir=args["h5ad_dir"]
sample_id=args["sample_id"]
project_dir=args["project_dir"]

# Import function script as a module
sys.path.append(os.path.join(project_dir, 'scripts', 'TLS_model', '05_HEV_analysis', '04_cellular_composition',))
import plotting_functions

# --- Get some info from the sample sheet

# Extract the external_id
sample_dir = os.path.join(project_dir, 'config')

sample_sheet_file = 'sample_sheet.csv'

sample_sheet = pd.read_csv(os.path.join(sample_dir, sample_sheet_file))
external_id = sample_sheet[sample_sheet['section_id'] == sample_id]['external_id'].values[0]
tissue_id = sample_sheet[sample_sheet['section_id'] == sample_id]['tissue_id'].values[0]


# --- Also get some tissue metadata

resource_dir = os.path.join(project_dir, 'resources')
metadata_file = 'tissue_metadata.csv'

tissue_metadata = pd.read_csv(os.path.join(resource_dir, metadata_file))

# Create a sample - subtype dictionary
subtype = tissue_metadata[tissue_metadata['tissue_id'] == str(tissue_id)]['clinical_subtype'].values[0]


# --- Load the adata object

use_data = 'log_norm'
# Need both the 48 & 16 res
resolution = 'square_016um'

adata = ad.read_h5ad(os.path.join(h5ad_dir, resolution, 'raw', f'{sample_id}_raw.h5ad'))


# --- Calculate crop coords for optimal presentation

coords = auto_crop_scanpy(adata, border_size=0.1)


# --- Load decon results

decon_resolution = resolution

resolution_numbers = re.findall(r'\d+', decon_resolution)[0]
resolution_numbers = int(resolution_numbers)

decon_level = 'minor'
ncells = 3

decon_dir = os.path.join(project_dir, 'resources', 'deconvolution', 'cell2location', f'{decon_level}_level', sample_id, decon_resolution,
                         f'n_cells_{ncells}', 'model_output')

decon_file = 'q05_cell_abundance_w_sf.csv'

decon_res = pd.read_csv(os.path.join(decon_dir, decon_file), index_col=0)

cell_types = list(decon_res.columns)
adata.obs[cell_types] = decon_res


# --- Load miniatlas classifications

# To make it easier to visualize trends by lineage
class_dir = os.path.join(project_dir, 'resources', 'cell_types', 'miniatlas')
class_file = 'miniatlas_subset_classifications.csv'

class_df = pd.read_csv(os.path.join(class_dir, class_file))

# Get endothelial classes
endothelial = list(class_df[class_df['miniatlas_major'] == 'Endothelial'][f'miniatlas_{decon_level}'].unique())
endothelial = sorted(endothelial)
 

# --- Plot 1: Spatial plots endothelial cells

ncols = len(endothelial)
nrows = 1

fig, axs = plt.subplots(ncols=ncols, nrows=nrows, figsize=(7*ncols, 5))

for n, cell_type in enumerate(endothelial):
    ax = axs[n]
    sc.pl.spatial(adata, color=cell_type, bw=False, alpha_img=0.7, crop_coord=coords, 
                title=f'{sample_id} ({subtype})\n{cell_type}',
                cmap='magma', ax=ax, vmin=0, vmax='p99.2')
    ax.set_title(ax.get_title(), fontweight='bold')

plt.suptitle(f'{sample_id} ({subtype}) - Endothelial decon results (decon_level={decon_level}, decon_res={decon_resolution}, ncells={ncells})', 
             y=1.1, x=0.5, fontweight='bold')


# --- Normalize the cell numbers & plot in single plot

tmp_df = adata.obs[endothelial]
df_normalized = (tmp_df - tmp_df.min()) / (tmp_df.max() - tmp_df.min())

adata.obs[endothelial] = df_normalized


# --- Plot 2: Endothelial after normalizing

# Set a vmax
vmax=0.1

fig, axs = plt.subplots(ncols=ncols, nrows=nrows, figsize=(7*ncols, 5))

for n, cell_type in enumerate(endothelial):
    ax = axs[n]
    sc.pl.spatial(adata, color=cell_type, bw=False, alpha_img=0.7, crop_coord=coords, 
                title=f'{sample_id} ({subtype})\n{cell_type}',
                cmap='magma', ax=ax, vmax=vmax)
    ax.set_title(ax.get_title(), fontweight='bold')

plt.suptitle(f'{sample_id} ({subtype}) - Endothelial decon results normalized from 0-1 (decon_level={decon_level}, decon_res={decon_resolution}, ncells={ncells}, vmax={vmax})', 
             y=1.1, x=0.5, fontweight='bold')


# --- Plot 3: All endothelial in a single plot

# Use function from cell
white_spacing=20

# Use function from cell2location to make colormap with increasing transparency gradient
# To be able to plot multiple scores in single plot
def create_colormap(R, G, B):
    spacing = int(white_spacing * 2.55)
    N = 255
    M = 3
    alphas = np.concatenate([[0] * spacing * M,
                                np.linspace(0, 1.0, (N - spacing) * M)])
    vals = np.ones((N * M, 4))
    #         vals[:, 0] = np.linspace(1, R / 255, N * M)
    #         vals[:, 1] = np.linspace(1, G / 255, N * M)
    #         vals[:, 2] = np.linspace(1, B / 255, N * M)
    for i, color in enumerate([R, G, B]):
        vals[:, i] = color / 255
    vals[:, 3] = alphas
    return ListedColormap(vals)

# Also load the colors for the endothelial types
colors_file = os.path.join(project_dir, 'resources', 'plotting_colors', 'endothelial_colors.csv')
colors_df = pd.read_csv(colors_file)

# Make a color_dict
color_dict = colors_df.set_index('endothelial_type')['color'].to_dict()

# Also make a color list (ordered by endothelial list)
color_list = [color_dict[endo] for endo in endothelial]

color_list_rgb = []
for value in color_dict.values():
    rgb = tuple(int(c * 255) for c in to_rgb(value))
    color_list_rgb.append(rgb)

cmap_list = []

for rgb in color_list_rgb:
    cmap = create_colormap(rgb[0], rgb[1], rgb[2])
    cmap_list.append(cmap)

# Plotting
ncols = len(endothelial)+1
vmin = 0

fig, axs = plt.subplots(ncols=ncols, nrows=nrows, figsize=(7*ncols, 5))

# First plotting separately
# Then all in 1 plot
for n, cell_type in enumerate(endothelial):
    ax = axs[n]
    sc.pl.spatial(adata, color=cell_type, bw=False, alpha_img=0.5, crop_coord=coords, 
                title=f'{sample_id} ({subtype})\n{cell_type}',
                cmap=cmap_list[n], ax=ax, vmin=vmin, vmax=vmax,)
    ax.set_title(ax.get_title(), fontweight='bold')

ax = axs[n+1]
# Then all in 1 plot
for n, cell_type in enumerate(endothelial):
    if n == 0:
        alpha_img = 0.7
    else:
        alpha_img = 0
    sc.pl.spatial(adata, color=cell_type, bw=False, alpha_img=alpha_img, crop_coord=coords, 
                title=f'{sample_id} ({subtype})\n{cell_type}',
                cmap=cmap_list[n], ax=ax, vmin=vmin, vmax=vmax, colorbar_loc=None)

ax.set_title('Endothelial subtypes', fontweight='bold')

# Add colorbars per score
for i, label in enumerate(endothelial):
    cmap = cmap_list[i]
    norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
    spacing = 0.13
    mid_y = 0.5
    num_cbars = len(endothelial)
    y_pos = mid_y + (i - (num_cbars - 1) / 2) * spacing  # Distribute around the midpoint
    cbar_ax = ax.inset_axes([1.1, y_pos, 0.3, 0.05])
    cbar = fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), 
                        cax=cbar_ax, orientation='horizontal')
    cbar.ax.set_ylabel(f'{label}', rotation=0, ha='left', va='center', fontsize=18)
    cbar.ax.yaxis.set_label_position("right")
    cbar.ax.yaxis.set_label_coords(1.05, 1.05)
    cbar.ax.tick_params(labelsize=15)

plt.suptitle(f'{sample_id} ({subtype}) - Endothelial decon results normalized from 0-1 (decon_level={decon_level}, decon_res={decon_resolution}, ncells={ncells}, vmax={vmax})', 
             y=1.1, x=0.5, fontweight='bold')


# --- Prepare for extracting groups

# Get a scalefactor for converting pixel coords to um
spot_diameter_pixel = adata.uns['spatial'][sample_id]['scalefactors']['spot_diameter_fullres']

resolution_numbers = re.findall(r'\d+', resolution)[0]
resolution_numbers = int(resolution_numbers)

spot_diameter_um = resolution_numbers # Square size is whatever the resolution is 
scalefactor_um = spot_diameter_pixel/spot_diameter_um

# Vis hd spots have no distance between them, so distance from center to center is the square size
spot_dist_um = resolution_numbers
spot_dist_pixel = spot_dist_um * scalefactor_um

# But Vis hd spots are squares, I also want to catch the diagonal neighbors
# So need to adjust the radius for this (will not catch these neighbors if I use the length of the square's edge)

# Lucky I remember how to do pythagoras
c = pow(spot_dist_pixel/2, 2) + pow(spot_dist_pixel/2, 2)
c = math.sqrt(c)

radius_max = c*2

# Round up to whole pixel
radius_max = math.ceil(radius_max)

# In theory this should work, but for some reason it doesn't in some areas??
# Make it a bit larger and then it works on all spots
radius_max = radius_max * 1.1

# Make connectivity matrix
sq.gr.spatial_neighbors(adata, coord_type='generic', radius=(0, radius_max))


# --- Define some functions for grouping

# Function to expand existing groups with lower threshold
def expand_existing_groups(
    adata, 
    group_key='groups', 
    score_key='score', 
    lower_thresh=0.2, 
    obs_key='expanded_groups'
):
    """
    Expands existing groups in `adata.obs[group_key]` by including neighboring spots
    with a score above `lower_thresh`.
    Parameters:
    - adata: AnnData object
    - group_key: Column in `adata.obs` containing initial groups
    - score_key: Column in `adata.obs` containing the score values
    - lower_thresh: Minimum score to include neighbors in a group
    - obs_key: Column name to store the expanded group annotation
    Returns:
    - Updated `adata.obs[obs_key]` with expanded group labels.
    """
    # Get spatial connectivity matrix
    connectivity_matrix = adata.obsp['spatial_connectivities'].tocsc()
    # Get spot names and scores
    spot_names = list(adata.obs_names)
    scores = adata.obs[score_key]
    # Initialize a mapping of spots to their groups
    spot_to_group = adata.obs[group_key].copy()
    # Identify NaN spots (not initially assigned to any group)
    unassigned_spots = set(spot_to_group[spot_to_group.isna()].index)
    # Create a mapping of existing groups
    group_dict = {}
    for spot, group in spot_to_group.dropna().items():
        if group not in group_dict:
            group_dict[group] = set()
        group_dict[group].add(spot)
    # Expand groups
    merging = True
    while merging:
        merging = False
        for group, members in group_dict.items():
            # Collect neighbors of all group members
            neighbors = set()
            for spot in members:
                spot_idx = spot_names.index(spot)
                neighbor_idxs = list(adata.obsp['spatial_connectivities'][spot_idx].indices)
                # Only consider neighbors that are not yet assigned to any group or have a score above threshold
                for neighbor_idx in neighbor_idxs:
                    neighbor_spot = spot_names[neighbor_idx]
                    if spot_to_group[neighbor_spot] is np.nan:  # Unassigned spot
                        if scores[neighbor_spot] >= lower_thresh:
                            neighbors.add(neighbor_spot)
                    elif spot_to_group[neighbor_spot] != group:  # Neighbor belongs to a different group
                        if scores[neighbor_spot] >= lower_thresh:
                            neighbors.add(neighbor_spot)
            # If new valid members found, add them to the group
            new_members = {n for n in neighbors if n in unassigned_spots and scores[n] >= lower_thresh}
            if new_members:
                group_dict[group].update(new_members)
                unassigned_spots.difference_update(new_members)
                merging = True  # Continue expanding
    # Assign expanded groups to adata
    expanded_series = pd.Series(np.nan, index=adata.obs_names, dtype="object")
    for group, members in group_dict.items():
        expanded_series.loc[list(members)] = group
    adata.obs[obs_key] = expanded_series
    return adata.obs[obs_key]


# squidpy's spatial_connectivities just saves which spots are neighbors
# so there will be overlap between the sublists (if a spot is a neighbor of multiple other spots, it appears multiple times)
# therefore merge the sublists until all values are unique

# Define a function to do this
def extract_neighbor_groups(connected_groups, adata, obs_key='group'):
    pooled = [set(subList) for subList in connected_groups]
    merging = True
    while merging:
        merging=False
        for i,group in enumerate(pooled):
            merged = next((g for g in pooled[i+1:] if g.intersection(group)),None)
            if not merged: continue
            group.update(merged)
            pooled.remove(merged)
            merging = True
    # Convert the list of sets to a list of lists
    list_of_lists = [list(s) for s in pooled]
    # Generate group names
    groups = [f'group{i+1}' for i in range(len(list_of_lists))]
    # Flatten the list of lists and pair each item with its corresponding group
    data_list = []
    for group, sublist in zip(groups, list_of_lists):
        for item in sublist:
            data_list.append({'index': item, 'group': group})
    # Create a DataFrame
    df = pd.DataFrame(data_list)
    # Sort DataFrame based on the 'index' column
    df_sorted = df.sort_values(by='index')
    # Transfer the barcode id to the index from the anndata object
    df_sorted.index = adata.obs_names
    # Tack the group annotation onto the anndata object
    adata.obs[obs_key] = df_sorted['group']
    # Return the grouping
    return list_of_lists

# And define a function that incorporates everything
def make_neighbor_groups(adata, cell_type, high_thresh, lower_thresh, min_spots=2):
    # Subset adata
    adata_sub = adata[adata.obs[cell_type] > high_thresh]
    # Make connectivity matrix
    sq.gr.spatial_neighbors(adata_sub, coord_type='generic', radius=(0, radius_max))
    # Extract connected spots as groups
    connected_groups = []
    for i, neighbors in enumerate(adata_sub.obsp['spatial_connectivities'].toarray()):
        group = [j for j, is_neighbor in enumerate(neighbors) if is_neighbor]
        group.append(i)  # Include the spot itself in the group
        connected_groups.append(group)
    # Extract groups
    list_of_lists = extract_neighbor_groups(connected_groups, adata_sub, obs_key='group')
    adata.obs['group'] = adata_sub.obs['group']
    # Expand groups to neighboring spots with lower values
    lower_thresh = 0.05
    expand_existing_groups(adata, group_key='group', score_key=cell_type, lower_thresh=lower_thresh, obs_key='expanded_groups')
    adata_sub.obs['expanded_groups'] = adata.obs['expanded_groups']
    # Find groups that are too small
    group_key = 'expanded_groups'
    too_small = adata_sub.obs[group_key].value_counts()[adata_sub.obs[group_key].value_counts() < min_spots].index
    # Make a temporary df
    temp_df = adata_sub.obs[[group_key]]
    # Set small groups as too_small
    too_small_value = 'too_small'
    temp_df[group_key] = np.where(temp_df[group_key].isin(too_small), too_small_value, temp_df[group_key])
    # Save singletons as separat df
    too_small = temp_df[temp_df[group_key] == too_small_value]
    # Remove the singletons for now
    temp_df = temp_df[temp_df[group_key] != too_small_value]
    # Update group names based on size
    group_counts = temp_df[group_key].value_counts()
    sorted_groups = group_counts.sort_values(ascending=False).index
    group_map = {group: f"group{i+1}" for i, group in enumerate(sorted_groups)}
    temp_df[group_key] = temp_df[group_key].map(group_map)
    # Add results to a new column
    adata_sub.obs['group_too_small_removed'] = temp_df[group_key]
    # Add groups + singletons to another column
    group_too_small = pd.concat([temp_df, too_small])
    adata_sub.obs['group_too_small'] = group_too_small
    return adata_sub.obs['group_too_small_removed']


# --- First make groups with high threshold

high_thresh = 0.15
lower_thresh = 0.05

for cell_type in endothelial:
    group_key = f'{cell_type}_groups'
    adata.obs[group_key] = make_neighbor_groups(adata, cell_type, high_thresh=high_thresh, lower_thresh=lower_thresh)

# Save these initial assignments for later
group_keys = [item + '_groups' for item in endothelial]

pre_conflict = adata.obs[group_keys]

# Investigate how many spots have conflicts
pre_conflict = pre_conflict.dropna(how='all')

# Rename columns for ease
pre_conflict.columns = endothelial

# Step 1: Convert to assignment matrix
assign = pre_conflict.notna().astype(int)

# Step 2: Initialize output matrix
categories = assign.columns
result = pd.DataFrame(0, index=categories, columns=categories)

# Step 3: Iterate over rows and populate the co-assignment counts
for _, row in assign.iterrows():
    assigned = row[row == 1].index.tolist()
    if len(assigned) == 1:
        # Unique assignment → only count (i, i)
        cat = assigned[0]
        result.loc[cat, cat] += 1
    else:
        # Multi-assignment → count all unique unordered pairs (i, j) where i ≠ j
        for i, j in itertools.combinations(assigned, 2):
            result.loc[i, j] += 1
            result.loc[j, i] += 1

conflict_df = result

# --- Plot 4a: Plot groups

# Plotting
ncols = len(endothelial)

fig, axs = plt.subplots(ncols=ncols, nrows=nrows, figsize=(5*ncols, 5))

# First plotting separately
# Then all in 1 plot
for n, cell_type in enumerate(endothelial):
    ax = axs[n]
    group_key = f'{cell_type}_groups'
    adata_sub = adata[adata.obs[group_key].notna()]
    sc.pl.spatial(adata_sub, color=group_key, bw=False, alpha_img=0.5, crop_coord=coords, 
                title=f'{cell_type}\nSpatial groups', legend_loc=None, ax=ax,)
    ax.set_title(ax.get_title(), fontweight='bold')

plt.suptitle(f'{sample_id} ({subtype}) - Endothelial decon spatial groups (high_thresh={high_thresh}, lower_thresh={lower_thresh})', 
             y=1.1, x=0.5, fontweight='bold')


# --- Plot 4b: Plot conflicts in assignment as a confusion matrix

sc.set_figure_params(fontsize=10)

# Disable the grid
plt.rcParams['axes.grid'] = False

# Make a new col_colors_df
col_colors = pd.DataFrame(index=endothelial)
col_colors['endothelial_type'] = col_colors.index.map(color_dict)

# Make heatmap
g = sns.clustermap(conflict_df, row_cluster=False, col_cluster=False, cmap='Blues', linewidths=0.5, figsize=(4, 4), annot=True,
               colors_ratio=(0.01, 0.1), linecolor='black', xticklabels=True, yticklabels=True, fmt='g', col_colors=col_colors,)

axs = plotting_functions.beautify_heatmap(g, row_colors=False, col_colors=True, cbar_title='', label_size=10, fontsize_title=11, title='',
                    cbar_pad=0.1, cbar_height_square=len(endothelial), cbar_width_square=0.3, space_per_square=0.1, width_to_height_ratio=1,
                    colors_div=3,)

heatmap_ax = g.ax_heatmap

# Move y-axis labels to the left
heatmap_ax.yaxis.set_label_position("left")  # Move axis label
heatmap_ax.yaxis.tick_left()  # Move tick labels

g.ax_col_colors.set_title(f'Vessel assignment per spot', fontweight='bold')

# Add labels
heatmap_ax.set_xlabel('')
heatmap_ax.set_ylabel('')

# Calculate the midpoint to position the title
title_x = (heatmap_ax.get_position().x0 + heatmap_ax.get_position().x1)/2
title_y = heatmap_ax.get_position().y1 + 0.2

plt.suptitle(f'{sample_id} ({subtype}) - Number of spots with conflicting assignment', y=title_y, x=title_x, fontweight='bold')


# --- Plot 4c: Same as above but proportions

rowsum = conflict_df.sum(axis=1)
prop_df = conflict_df.div(rowsum, axis=0)

# Make heatmap
g = sns.clustermap(prop_df.round(2), row_cluster=False, col_cluster=False, cmap='Blues', linewidths=0.5, figsize=(4, 4), annot=True,
               colors_ratio=(0.01, 0.1), linecolor='black', xticklabels=True, yticklabels=True, fmt='g', col_colors=col_colors,)

axs = plotting_functions.beautify_heatmap(g, row_colors=False, col_colors=True, cbar_title='', label_size=10, fontsize_title=11, title='',
                    cbar_pad=0.1, cbar_height_square=len(endothelial), cbar_width_square=0.3, space_per_square=0.1, width_to_height_ratio=1,
                    colors_div=3,)

heatmap_ax = g.ax_heatmap

# Move y-axis labels to the left
heatmap_ax.yaxis.set_label_position("left")  # Move axis label
heatmap_ax.yaxis.tick_left()  # Move tick labels

g.ax_col_colors.set_title(f'Vessel assignment per spot\nProportion (row)', fontweight='bold')

# Add labels
heatmap_ax.set_xlabel('')
heatmap_ax.set_ylabel('')

# Calculate the midpoint to position the title
title_x = (heatmap_ax.get_position().x0 + heatmap_ax.get_position().x1)/2
title_y = heatmap_ax.get_position().y1 + 0.2

# Remove yticks col color ax
g.ax_col_colors.set_yticks([])

plt.suptitle(f'{sample_id} ({subtype}) - Proportion of spots (per row) with conflicting assignment', y=title_y, x=title_x, fontweight='bold')


# --- Also look at conflict on the vessel (i.e. spatial group) level

# Make an empty df to add counts  to
result = pd.DataFrame(0, index=endothelial, columns=endothelial)

for col in categories:
    col_groups = pre_conflict[col].dropna().unique()
    for group in col_groups:
        # Find barcodes (rows) where this group appears in this column
        barcodes = pre_conflict[pre_conflict[col] == group].index
        # Subset the original df to just these barcodes
        df_subset = pre_conflict.loc[barcodes]
        # Check which columns have any non-null values in these barcodes
        active_columns = df_subset.notna().any(axis=0)
        assigned_columns = active_columns[active_columns].index.tolist()
        if len(assigned_columns) == 1:
            # Only assigned in this column
            result.loc[col, col] += 1
        else:
            for other_col in assigned_columns:
                if other_col != col:
                    result.loc[col, other_col] += 1

vessel_conflict = result


df_sub = pre_conflict[pre_conflict['Endothelial ACKR1'].notna()]
len(df_sub['Endothelial ACKR1'].unique())

df_sub_sub = df_sub[df_sub['Endothelial RGS5'].notna()]
len(df_sub_sub['Endothelial ACKR1'].unique())


df_sub = pre_conflict[pre_conflict['Endothelial CXCL12'].notna()]
len(df_sub['Endothelial CXCL12'].unique())

df_sub_sub = df_sub[df_sub['Endothelial RGS5'].notna()]
len(df_sub_sub['Endothelial CXCL12'].unique())


df_sub = pre_conflict[pre_conflict['Endothelial RGS5'].notna()]
len(df_sub['Endothelial RGS5'].unique())

df_sub_sub = df_sub[df_sub['Endothelial CXCL12'].notna()]
len(df_sub_sub['Endothelial RGS5'].unique())


# --- Plot 3d: Confusion matrix vessel level

# Make heatmap
g = sns.clustermap(vessel_conflict, row_cluster=False, col_cluster=False, cmap='Blues', linewidths=0.5, figsize=(4, 4), annot=True,
               colors_ratio=(0.01, 0.1), linecolor='black', xticklabels=True, yticklabels=True, fmt='g', col_colors=col_colors,)

axs = plotting_functions.beautify_heatmap(g, row_colors=False, col_colors=True, cbar_title='', label_size=10, fontsize_title=11, title='',
                    cbar_pad=0.1, cbar_height_square=len(endothelial), cbar_width_square=0.3, space_per_square=0.1, width_to_height_ratio=1,
                    colors_div=3,)

heatmap_ax = g.ax_heatmap

# Move y-axis labels to the left
heatmap_ax.yaxis.set_label_position("left")  # Move axis label
heatmap_ax.yaxis.tick_left()  # Move tick labels

g.ax_col_colors.set_title(f'Assignment per vessel', fontweight='bold')

# Add labels
heatmap_ax.set_xlabel('')
heatmap_ax.set_ylabel('')

# Calculate the midpoint to position the title
title_x = (heatmap_ax.get_position().x0 + heatmap_ax.get_position().x1)/2
title_y = heatmap_ax.get_position().y1 + 0.2

plt.suptitle(f'{sample_id} ({subtype}) - Number of vessels with conflicting assignment', y=title_y, x=title_x, fontweight='bold')


# --- Plot 3e: Same as above but with proportions

rowsum = vessel_conflict.sum(axis=1)
prop_df = vessel_conflict.div(rowsum, axis=0)

# Make heatmap
g = sns.clustermap(prop_df.round(2), row_cluster=False, col_cluster=False, cmap='Blues', linewidths=0.5, figsize=(4, 4), annot=True,
               colors_ratio=(0.01, 0.1), linecolor='black', xticklabels=True, yticklabels=True, fmt='g', col_colors=col_colors,)

axs = plotting_functions.beautify_heatmap(g, row_colors=False, col_colors=True, cbar_title='', label_size=10, fontsize_title=11, title='',
                    cbar_pad=0.1, cbar_height_square=len(endothelial), cbar_width_square=0.3, space_per_square=0.1, width_to_height_ratio=1,
                    colors_div=3,)

heatmap_ax = g.ax_heatmap

# Move y-axis labels to the left
heatmap_ax.yaxis.set_label_position("left")  # Move axis label
heatmap_ax.yaxis.tick_left()  # Move tick labels

g.ax_col_colors.set_title(f'Assignment per vessel\nProportions (row)', fontweight='bold')

# Add labels
heatmap_ax.set_xlabel('')
heatmap_ax.set_ylabel('')

# Calculate the midpoint to position the title
title_x = (heatmap_ax.get_position().x0 + heatmap_ax.get_position().x1)/2
title_y = heatmap_ax.get_position().y1 + 0.2

# Remove yticks col color ax
g.ax_col_colors.set_yticks([])

plt.suptitle(f'{sample_id} ({subtype}) - Proportion of vessels with conflicting assignment (proportions per row)', y=title_y, x=title_x, fontweight='bold')


# --- Deal with overlap - single vessels with multiple groups: assign dominant

group_keys = [f"{item}_groups" for item in endothelial]

tmp_df = adata.obs[group_keys]

def resolve_group_conflicts(df):
    # Count how many barcodes are in each group across all columns
    group_sizes = df.apply(pd.Series.value_counts).sum(axis=1).to_dict()
    # Lists to collect conflict status per row
    conflict_flags = []
    def keep_largest(row):
        groups = row.dropna().unique()
        if len(groups) <= 1:
            conflict_flags.append('no_conflict')
            return row  # no conflict
        # There’s a conflict
        conflict_flags.append('conflict')
        # Find largest group
        largest_group = max(groups, key=lambda g: group_sizes.get(g, 0))
        return row.apply(lambda g: g if g == largest_group else None)
    # Apply row-wise conflict resolution
    df_resolved = df.apply(keep_largest, axis=1)
    # Add conflict status column
    df_resolved['conflict_status'] = conflict_flags
    return df_resolved

def clean_and_rename_groups(df_resolved, min_group_size=2):
    # Exclude the conflict_status column
    group_df = df_resolved.drop(columns='conflict_status', errors='ignore')
    # Step 1: Filter out groups that are too small, per column
    df_filtered = group_df.copy()
    for col in df_filtered.columns:
        counts = df_filtered[col].value_counts()
        valid_groups = counts[counts >= min_group_size].index
        df_filtered[col] = df_filtered[col].apply(lambda g: g if g in valid_groups else None)
    # Step 2: Get unique group names (from filtered), sorted by frequency
    updated_counts = df_filtered.apply(pd.Series.value_counts).sum(axis=1)
    sorted_groups = updated_counts.sort_values(ascending=False).index.tolist()
    # Step 3: Generate new names: group1, group2, ...
    group_renaming = {old: f"group{i+1}" for i, old in enumerate(sorted_groups)}
    # Step 4: Apply renaming
    df_renamed = df_filtered.applymap(lambda g: group_renaming.get(g) if g in group_renaming else None)
    # Step 5: Add back conflict_status if it was in the original df
    if 'conflict_status' in df_resolved.columns:
        df_renamed['conflict_status'] = df_resolved['conflict_status']
    return df_renamed

df_resolved = resolve_group_conflicts(tmp_df.dropna(how='all'))
df_renamed = clean_and_rename_groups(df_resolved)

# Drop 'conflict_status' if present (we only care about group assignment columns)
group_df = df_renamed.drop(columns='conflict_status', errors='ignore')

# For each row, find the column where the value is not None
assigned_from = group_df.apply(lambda row: row.first_valid_index(), axis=1)

# Add as new column
df_renamed_with_source = df_renamed.copy()
df_renamed_with_source['assigned_from'] = assigned_from
df_renamed_with_source['assigned_from'] = df_renamed_with_source['assigned_from'].str.removesuffix('_groups')

df_renamed_with_source['assigned_from'] = pd.Categorical(df_renamed_with_source['assigned_from'], ordered=True, categories=endothelial)

# Add to adata - recording vessel type
adata.obs['vessel_type'] = df_renamed_with_source['assigned_from']


# --- Plot 5: Overview of conflicts and how dealt with

adata.obs[df_renamed.columns] = df_renamed

# Plotting
ncols = len(endothelial)+3

fig, axs = plt.subplots(ncols=ncols, nrows=nrows, figsize=(5*ncols, 5))

# First plotting separately
# Then all in 1 plot
for n, cell_type in enumerate(endothelial):
    ax = axs[n]
    group_key = f'{cell_type}_groups'
    adata_sub = adata[adata.obs[group_key].notna()]
    sc.pl.spatial(adata_sub, color=group_key, bw=False, alpha_img=0.5, crop_coord=coords, 
                title=f'{cell_type}\nSpatial groups', legend_loc=None,
                palette=color_list, ax=ax,)
    ax.set_title(ax.get_title(), fontweight='bold')

# And plot whether or not there was conflict
ax = axs[n+1]
group_key = 'conflict_status'
adata_sub = adata[adata.obs[group_key].notna()]
sc.pl.spatial(adata_sub, color=group_key, bw=False, alpha_img=0.5, crop_coord=coords, 
            title=f'{cell_type}\nGrouping conflicts', ax=ax,)
ax.set_title(ax.get_title(), fontweight='bold')

# Make the next ax invisible to accommodate legend
axs[n+2].set_visible(False)

# And plot the vessel type for each group
ax = axs[n+3]
group_key = 'vessel_type'
adata_sub = adata[adata.obs[group_key].notna()]
sc.pl.spatial(adata_sub, color=group_key, bw=False, alpha_img=0.5, crop_coord=coords, 
            title=f'{cell_type}\nVessel_type', ax=ax, palette=color_list)
ax.set_title(ax.get_title(), fontweight='bold')

# Move a bit closer to previous ax
ax_prev_pos = axs[n+1].get_position()
ax_pos = ax.get_position()

ax.set_position([ax_prev_pos.x1+ax_prev_pos.width/2, ax_pos.y0, ax_pos.width, ax_pos.height])

plt.suptitle(f'{sample_id} ({subtype}) - Resolve grouping conflicts (assign largest group) (high_thresh={high_thresh}, lower_thresh={lower_thresh})', 
             y=1.1, x=0.5, fontweight='bold')


# --- Plot 6: Plots of vessel number + vessel size

# Make font a bit smaller
plt.rcParams.update({'font.size': 10})

# Enable gridlines behind the plot for all plots
plt.rcParams['axes.axisbelow'] = True

# Also make a color dict
color_dict = {}

for n, endo_type in enumerate(endothelial):
    color_dict[endo_type] = color_list[n]

# Count number of groups per category
unique_counts = df_renamed[group_keys].apply(lambda col: col.dropna().nunique()).to_frame()
unique_counts.reset_index(inplace=True)
unique_counts.columns = ['group_type', 'counts']

unique_counts['endothelial_type'] = unique_counts['group_type'].str.removesuffix('_groups')

# Calculate group size per category
group_counts_df = df_renamed[group_keys].melt(var_name='group_type', value_name='group') \
                    .dropna() \
                    .value_counts(subset=['group_type', 'group']) \
                    .reset_index(name='count')

# Remove the suffix
group_counts_df['endothelial_type'] = group_counts_df['group_type'].str.removesuffix('_groups')

# Plotting
nrows = 1
ncols = 2
fig, axs = plt.subplots(ncols=ncols, nrows=nrows, figsize=(3*ncols, 2.5))

# -- First barplot showing number of groups
n=0
ax = axs[n]
sns.barplot(unique_counts, x='endothelial_type', y='counts', order=endothelial, palette=color_dict, ax=ax, hue='endothelial_type',
            saturation=1)

ax.set_xticklabels(ax.get_xticklabels(), rotation=90)

# Beautify the plot
plotting_functions.beautify_plot(ax)
# Annotate the bars
plotting_functions.annotate_bars(ax, unique_counts, x='endothelial_type', y='counts', plot_order=endothelial, annotation_size=9)

ax.set_ylabel('# groups')
ax.set_xlabel('')
ax.set_title(f'{sample_id}\nNumber of groups', fontweight='bold')

# -- Then boxplot showing groups size

# Boxplot size
# Properties for boxplot
whisker_props = {'color': 'black'}
median_props = {'color': 'black'}

n+=1
ax = axs[n]

sns.boxplot(group_counts_df, x='endothelial_type', y='count', order=endothelial, palette=color_dict, ax=ax, hue='endothelial_type',
            saturation=1, whiskerprops=whisker_props, medianprops=median_props, fliersize=1)

# Beautify the plot
plotting_functions.beautify_plot(ax)

ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
ax.set_ylabel('Size (# spots)')
ax.set_xlabel('')
ax.set_title(f'{sample_id}\nGroup size', fontweight='bold')

plt.subplots_adjust(wspace=0.3)

plt.suptitle(f'{sample_id} ({subtype}) - Grouping overview (high_thresh={high_thresh}, lower_thresh={lower_thresh})', 
             y=1.2, x=0.5, fontweight='bold')


# --- Save the output

threshold = 0.5
threshold_adj = str(threshold).split('.')[-1]

qc_dir = os.path.join(project_dir, 'results', 'TLS_model')

today = date.today()
today = today.strftime("%Y%m%d")

out_dir = os.path.join(qc_dir, 'other_signatures', 'HEV_analysis', 'per_sample', 'square_016um', f'threshold_{threshold_adj}', 'other_vessels',
                       'classify_vessels')

os.makedirs(out_dir, exist_ok=True)

# Save plots
out_dir_plots = os.path.join(out_dir, 'plots', today)
os.makedirs(out_dir_plots, exist_ok=True)

# Also save the groupings as a .csv
out_dir_csv = os.path.join(out_dir, 'csv_files', today)
os.makedirs(out_dir_csv, exist_ok=True)

# Save plots
filename = f'{sample_id}_vessel_classification.pdf'

save_multi_image(os.path.join(out_dir_plots, filename))

# And save the dfs
filename = f'{sample_id}_vessel_classification.csv'
# Trim column suffix
df_renamed.columns = [col.removesuffix('_groups') for col in df_renamed.columns]
df_renamed.to_csv(os.path.join(out_dir_csv, filename), index=True)

filename = f'{sample_id}_vessel_size.csv'
group_counts_df.to_csv(os.path.join(out_dir_csv, filename), index=False)

filename = f'{sample_id}_vessel_numbers.csv'
unique_counts.to_csv(os.path.join(out_dir_csv, filename), index=False)

# Also save the conflict dfs
filename = f'{sample_id}_vessel_conflict_spot_level.csv'
conflict_df.to_csv(os.path.join(out_dir_csv, filename), index=True)

filename = f'{sample_id}_vessel_conflict_vessel_level.csv'
vessel_conflict.to_csv(os.path.join(out_dir_csv, filename), index=True)

# And save the normalized endothelial numbers
filename = f'{sample_id}_endothelial_numbers_normalized.csv'
print(f"Saving endothelial numbers for columns: {endothelial}")
print(f"Shape of data to save: {adata.obs[endothelial].shape}")
print(f"First few rows:\n{adata.obs[endothelial].head()}")
adata.obs[endothelial].to_csv(os.path.join(out_dir_csv, filename), index=True)


# --- Finish

print("Script successfully completed")
