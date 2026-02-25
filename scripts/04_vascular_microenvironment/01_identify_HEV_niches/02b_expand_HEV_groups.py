#!/usr/bin/env python3

# --- Aim
# Expand initial HEV spot groups by including neighboring spots with lower HEV scores.
# Capture the full vessel structure of HEVs for more accurate spatial representation.

# --- Define environment
# env: vis-hd

# --- Load packages
import argparse
import json
import math
import os
import sys
from datetime import date, datetime

import anndata as ad
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from matplotlib.colors import ListedColormap, to_rgb
import re
import squidpy as sq

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
parser.add_argument("--sample_sheet", help="Path to sample sheet CSV file", default=None)

parser.add_argument("--min_count", help="Min count for QC", type=int)
parser.add_argument("--min_gene", help="Min gene for QC", type=int)
parser.add_argument("--min_spots", help="Minimal spots per gene", type=int)

args = vars(parser.parse_args())


# --- Arguments

res_dir = args["res_dir"]
h5ad_dir=args["h5ad_dir"]
sample_id=args["sample_id"]
project_dir=args["project_dir"]
sys.path.append(os.path.join(project_dir, 'scripts'))
from helper_functions import auto_crop_scanpy, save_multi_image


# --- Get some info from the sample sheet
# Use sample_sheet argument if provided, else default
if args.get("sample_sheet"):
    sample_sheet_path = args["sample_sheet"]
else:
    sample_dir = os.path.join(project_dir, 'config')
    sample_sheet_file = 'sample_sheet.csv'
    sample_sheet_path = os.path.join(sample_dir, sample_sheet_file)

sample_sheet = pd.read_csv(sample_sheet_path)
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
resolutions = ['square_048um', 'square_016um']

adata_dict = {}

for resolution in resolutions:
    adata_dict[resolution] = ad.read_h5ad(os.path.join(h5ad_dir, 'visium_hd', resolution, f'{sample_id}.h5ad'))


# --- Calculate crop coords for optimal presentation
coords_dict  = {}

for resolution in resolutions:
    coords = auto_crop_scanpy(adata_dict[resolution], border_size=0.1)
    coords_dict[resolution] = coords


# --- Load HEV groups

# Function to find the most recent results
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
    else:
        print("No valid directories found.")
    return most_recent_dir_path

# Find most recent results
threshold = 0.5
threshold_adj = str(threshold).split('.')[-1]

TLS_group_dir = os.path.join(res_dir, '04_vascular_microenvironment', '01_identify_HEV_niches', '02a_group_HEV_bins', resolution, f'threshold_{threshold_adj}', 'csv_files')

most_recent_dir_path = find_most_recent(TLS_group_dir)

# Read in the results & add to anndata
score_of_i = 'HEVs'
TLS_file_name = f'{sample_id}_{score_of_i}_groups.csv'
TLS_results = pd.read_csv(os.path.join(most_recent_dir_path, TLS_file_name), index_col=0)

# Load the groups that also include groups considered too small (may be large enough when expanded)
group_key = 'group'

adata_dict['square_016um'].obs[group_key] = TLS_results[group_key]


# --- And load the scoring results

clusters_of_i = {
    'B_cell_zone': 'b_cell_zone_top_genes_filtered.csv',
    'T_cell_zone': 't_cell_zone_top_genes_filtered.csv',
    'HEVs': 'hevs_filtered_endothelial_genes.csv'
}

# Define the TLS zone names
TLS_zone_names = ['T_cell_zone', 'B_cell_zone']

score_keys_dict = {}

for score_name, filename in clusters_of_i.items():
    out_name = f'{score_name}_score'
    if score_name in TLS_zone_names:
        resolution = 'square_048um'
    else:
        resolution = 'square_016um'
    GC_res_dir = os.path.join(res_dir, '04_vascular_microenvironment', '01_identify_HEV_niches', '01_score_signature', out_name, resolution)
    most_recent_dir_path = find_most_recent(GC_res_dir)
    # Read results and add to adata
    GC_file_name = f'{sample_id}_{out_name}.csv'
    GC_results = pd.read_csv(os.path.join(most_recent_dir_path, 'csv_files', GC_file_name), index_col=0)
    # Choose number of genes used for scoring
    top_gene = 20
    GC_key = f'{score_name}_score_{top_gene}'
    score_keys_dict[score_name] = GC_key
    adata_dict[resolution].obs[GC_key] = GC_results[GC_key]


# Also make colormaps per score
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
    for i, color in enumerate([R, G, B]):
        vals[:, i] = color / 255
    vals[:, 3] = alphas
    return ListedColormap(vals)

# Define some colors per score
color_dict = {
    'GC': mpl.colors.cnames["cyan"],
    'B_cell_zone': mpl.colors.cnames["cornflowerblue"],
    'T_cell_zone': mpl.colors.cnames["crimson"],
    'HEVs': mpl.colors.cnames["indigo"],    
}

color_dict_rgb = {}
for key, value in color_dict.items():
    rgb = tuple(int(c * 255) for c in to_rgb(value))
    color_dict_rgb[key] = rgb

cmap_dict = {}

for score_name in clusters_of_i:
    rgb = color_dict_rgb[score_name]
    cmap = create_colormap(rgb[0], rgb[1], rgb[2])
    cmap_dict[score_name] = cmap


# --- Plot 1: Plot HEV groups and TLS scores for sanity

# Plotting
ncols = 2
nrows = 1

sc.set_figure_params(fontsize=15)

fig, axs = plt.subplots(nrows, ncols, figsize=(8*ncols, 6.5*nrows))

# -- First plot TLS groups & HEV classifications
n = 0
ax = axs[n]

sc.pl.spatial(adata_dict['square_016um'], color=group_key, bw=False, alpha_img=0.7, crop_coord=coords_dict['square_016um'], title=f'{sample_id} ({subtype})\nHEV groups > {threshold}',
                palette=['indigo'], ax=ax, legend_loc=None)

ax.set_title(ax.get_title(), fontweight='bold')


# -- Then the score results
n += 1
ax = axs[n]

vmax = 0.8

# Then all scores in single plot
# Use scanpy because squidpy ignores info abt alpha in cmap
score_names = list(clusters_of_i.keys())
score_names_r = score_names[::-1]

for i, score_name in enumerate(score_names_r):
    if i == 0:
        alpha_img = 0.7
    else: 
        alpha_img = 0
    if score_name in TLS_zone_names:
        resolution = 'square_048um'
    else:
        resolution = 'square_016um'
    cmap = cmap_dict[score_name]
    sc.pl.spatial(adata_dict[resolution], color=score_keys_dict[score_name], alpha_img=alpha_img, cmap=cmap,
        crop_coord=coords_dict[resolution], title=f'{sample_id} ({subtype})\nTLS zones', ax=ax, vmin=0, vmax=vmax, colorbar_loc=None,)

ax.set_title(ax.get_title(), fontweight='bold')

# Add colorbars per score
for i, (label, cmap) in enumerate(cmap_dict.items()):
    norm = mpl.colors.Normalize(vmin=0, vmax=vmax)
    spacing = 0.13
    mid_y = 0.5
    num_cbars = len(cmap_dict)
    y_pos = mid_y + (i - (num_cbars - 1) / 2) * spacing  # Distribute around the midpoint
    cbar_ax = ax.inset_axes([1.1, y_pos, 0.3, 0.05])
    cbar = fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), 
                        cax=cbar_ax, orientation='horizontal')
    cbar.ax.set_ylabel(f'{score_keys_dict[label]}', rotation=0, ha='left', va='center', fontsize=18)
    cbar.ax.yaxis.set_label_position("right")
    cbar.ax.yaxis.set_label_coords(1.05, 1.05)
    cbar.ax.tick_params(labelsize=15)

plt.suptitle(f'{sample_id} ({subtype}) - HEV groups & TLS zone scores (top_gene={top_gene}, threshold={threshold}, resolution={resolution})', y=1.1, fontweight='bold', x=0.5)


# --- Make connectivity matrix

# Use squidpy to make the connectivity matrix (figure out which spots are neighbors)
resolution = 'square_016um'
adata = adata_dict[resolution]

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

# Use squidpy to define which TLS-high spots are neighbors
sq.gr.spatial_neighbors(adata, coord_type='generic', radius=(0, radius_max))


# --- Expand existing groups to neighboring spots with lower scores

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


# --- Apply the function

lower_thresh = 0.3

expanded_groups = expand_existing_groups(
    adata,
    group_key=group_key,
    score_key=score_keys_dict[score_of_i],
    lower_thresh=lower_thresh,
    obs_key='expanded_groups',
)

# Ensure the numbering is still correct
group_counts = expanded_groups.value_counts()
sorted_groups = group_counts.sort_values(ascending=False).index
group_map = {group: f"group{i+1}" for i, group in enumerate(sorted_groups)}

temp_df = pd.DataFrame(expanded_groups)
temp_df['group_updated'] = temp_df['expanded_groups'].map(group_map)

# Add to adata
adata.obs['expanded_groups'] = temp_df['group_updated']


# --- Merge neighboring groups
# The function above adds spots to the existing groups but doesn't merge pre-existing groups that become connected through the added spots
# So after expanding the groups need to go back and make sure all connected groups are merged

# Subset adata to just spots in the expanded groups
adata_sub = adata[adata.obs['expanded_groups'].notna()]

# Use squidpy to define which spots are neighbors
sq.gr.spatial_neighbors(adata_sub, coord_type='generic', radius=(0, radius_max))

# Extract connected spots as groups
connected_groups = []
for i, neighbors in enumerate(adata_sub.obsp['spatial_connectivities'].toarray()):
    group = [j for j, is_neighbor in enumerate(neighbors) if is_neighbor]
    group.append(i)  # Include the spot itself in the group
    connected_groups.append(group)

# squidpy's spatial_connectivities just saves which spots are neighbors
# so there will be overlap between the sublists (if a spot is a neighbor of multiple other spots, it appears multiple times)
# therefore merge the sublists until all values are unique

# Define a function to do this
def extract_neighbor_groups(connected_groups, adata=adata_sub, obs_key='group'):
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

# Extract groups
group_key = 'group'
list_of_lists = extract_neighbor_groups(connected_groups, obs_key=group_key)

# Also filter out groups that are too small
min_spots = 2

# Find groups that are too small
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
adata_sub.obs['expanded_groups'] = temp_df[group_key]

# Add groups + singletons to another column
group_too_small = pd.concat([temp_df, too_small])

adata_sub.obs['group_too_small'] = group_too_small

# Overwrite expanded groups in original adata
adata.obs['expanded_groups'] = adata_sub.obs['expanded_groups']


# --- Plot 2: Plot the expanded groups

adata.obs['above_low_thresh'] = np.where(
    adata.obs[score_keys_dict[score_of_i]] > lower_thresh,
    f'>{lower_thresh}',
    None,
)

adata.obs['spot_added'] = np.where(
    adata.obs['expanded_groups'].notna() & adata.obs[group_key].isna(),
    f'>{lower_thresh} neighbor',
    None,
)

# Plotting
ncols = 4
nrows = 1

fig, axs = plt.subplots(nrows, ncols, figsize=(8*ncols, 6.5*nrows))

# -- First plot HEV groups & spots above lower thresh
n = 0
ax = axs[n]

sc.pl.spatial(adata_dict['square_016um'], color='above_low_thresh', bw=False, alpha_img=0.7, crop_coord=coords_dict['square_016um'], title=f'{sample_id} ({subtype})\nHEV groups > {threshold}',
                palette=['darkorange'], ax=ax,)

sc.pl.spatial(adata_dict['square_016um'], color=group_key, bw=False, alpha_img=0, crop_coord=coords_dict['square_016um'], title=f'{sample_id} ({subtype})\nHEV groups > {threshold}\nHEV score > {lower_thresh}',
                palette=['indigo'], ax=ax, legend_loc=None)

ax.set_title(ax.get_title(), fontweight='bold')


# -- Then plot the expanded groups
n += 1
ax = axs[n]

sc.pl.spatial(adata_dict['square_016um'], color='expanded_groups', bw=False, alpha_img=0.7, crop_coord=coords_dict['square_016um'], title=f'{sample_id} ({subtype})\nHEV groups expanded > {lower_thresh}',
                ax=ax, legend_loc=None, palette=['indigo'],)

ax.set_title(ax.get_title(), fontweight='bold')


# -- Then plot the groups
n += 1
ax = axs[n]

# Make a copy so it will generate distinct colors
adata_dict['square_016um'].obs['expanded_groups_copy'] = adata_dict['square_016um'].obs['expanded_groups']

sc.pl.spatial(adata_dict['square_016um'], color='expanded_groups_copy', bw=False, alpha_img=0.7, crop_coord=coords_dict['square_016um'], title=f'{sample_id} ({subtype})\nHEV groups expanded > {threshold}',
                ax=ax, legend_loc=None)

ax.set_title(ax.get_title(), fontweight='bold')

plt.suptitle(f'{sample_id} ({subtype}) - HEV groups expanded to >{lower_thresh} neighbors (top_gene={top_gene}, threshold={threshold}, resolution={resolution})', y=1.1, fontweight='bold', x=0.5)


# -- Then plot the original groups + spots added
n += 1
ax = axs[n]

sc.pl.spatial(adata_dict['square_016um'], color='spot_added', bw=False, alpha_img=0.7, crop_coord=coords_dict['square_016um'], title=f'{sample_id} ({subtype})\nHEV groups > {threshold}',
                palette=['darkorange'], ax=ax,)

sc.pl.spatial(adata_dict['square_016um'], color=group_key, bw=False, alpha_img=0, crop_coord=coords_dict['square_016um'], title=f'{sample_id} ({subtype})\nHEV groups + spots added',
                palette=['indigo'], ax=ax, legend_loc=None)

ax.set_title(ax.get_title(), fontweight='bold')

plt.suptitle(f'{sample_id} ({subtype}) - HEV groups expanded to >{lower_thresh} neighbors (top_gene={top_gene}, threshold={threshold}, resolution={resolution})', y=1.1, fontweight='bold', x=0.5)


# --- Save output

qc_dir = os.path.join(res_dir, '04_vascular_microenvironment', '01_identify_HEV_niches', '02b_expand_HEV_groups')

threshold_adj = str(threshold).split('.')[-1]
lower_thresh_adj = str(lower_thresh).split('.')[-1]

today = date.today()
today = today.strftime("%Y%m%d")

out_dir = os.path.join(qc_dir, resolution, f'threshold_{threshold_adj}')

os.makedirs(out_dir, exist_ok=True)

# Save plots
out_dir_plots = os.path.join(out_dir, 'plots', 'groups', today)
os.makedirs(out_dir_plots, exist_ok=True)

# Also save the groupings as a .csv
out_dir_csv = os.path.join(out_dir, 'csv_files', today)
os.makedirs(out_dir_csv, exist_ok=True)

# Save plots
filename = f'{sample_id}_{score_of_i}_expanded_groups.pdf'

save_multi_image(os.path.join(out_dir_plots, filename))
plt.close('all')

# And save the groupings as a csv
to_save_df = adata.obs[[score_keys_dict[score_of_i], group_key, 'expanded_groups', 'spot_added']].dropna(subset='expanded_groups')

filename = f'{sample_id}_{score_of_i}_groups.csv'
to_save_df.to_csv(os.path.join(out_dir_csv, filename))

# Save the grouping settings as a JSON file
# Save in a dict
settings = {'score_threshold': threshold,
            'lower_threshold': lower_thresh,
            'score_keys': score_keys_dict[score_of_i],
            'resolution': resolution,
            }

with open(os.path.join(out_dir_csv, f'{sample_id}_grouping_settings.json'), 'w') as json_file:
    json.dump(settings, json_file, indent=4)  # indent=4 makes the JSON file more readable

