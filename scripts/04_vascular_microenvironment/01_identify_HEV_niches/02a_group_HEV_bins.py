#!/usr/bin/env python3

# --- Aim
# Group spatial spots with high HEV scores into contiguous regions to delineate individual HEVs.
# Enable quantification of HEV locations and facilitate downstream spatial analyses.

# --- Define environment
# env: vis-hd

# --- Load packages
import argparse
import json
import math
import os
import re
import sys

import anndata as ad
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import squidpy as sq
from matplotlib.colors import LinearSegmentedColormap, ListedColormap, to_rgb
from mycolorpy import colorlist as mcp

sc.set_figure_params(facecolor="white", figsize=(8, 8))
sc.settings.verbosity = 3


# --- Argparse arguments

parser = argparse.ArgumentParser(description="Import Visium data & do basic QC",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

# Required

parser.add_argument("--res_dir", help="Path to results directory")
parser.add_argument("--h5ad_dir", help="Path to directory for .h5ad output files")
parser.add_argument("--project_dir", help="Project directory (root)")
parser.add_argument("--sample_id", help="Sample ID")
parser.add_argument("--sample_sheet", help="Path to sample sheet CSV file", default=None)

args = vars(parser.parse_args())


# --- Arguments

res_dir = args["res_dir"]
h5ad_dir=args["h5ad_dir"]
sample_id=args["sample_id"]
project_dir=args["project_dir"]
sys.path.append(os.path.join(project_dir, 'scripts'))
from helper_functions import auto_crop_scanpy, auto_crop, save_multi_image



# --- Get some info from the sample sheet
# Use sample_sheet argument if provided, else default
if args.get("sample_sheet"):
    sample_sheet_path = args["sample_sheet"]
else:
    sample_dir = os.path.join(project_dir, 'config')
    sample_sheet_file = 'sample_sheet_visium_hd.csv'
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
# Going to make groups based on the 16 µm res
resolution = 'square_016um'

adata = ad.read_h5ad(os.path.join(h5ad_dir, 'visium_hd', resolution, f'{sample_id}.h5ad'))


# --- Calculate crop coords for optimal presentation
coords = auto_crop_scanpy(adata, border_size=0.1)

# Also calculate some coords for use w/ squidpy
coord_squidpy = auto_crop(adata, border_size=0.1)


# --- Load TLS results

clusters_of_i = {
    'GC': 'germinal_center_filtered_epithelial_genes.csv',
    'B_cell_zone': 'b_cell_zone_top_genes_filtered.csv',
    'T_cell_zone': 't_cell_zone_top_genes_filtered.csv',
    'HEVs': 'hevs_filtered_endothelial_genes.csv'
}

# Define HEV score as score of interest
score_of_i = 'HEVs'

# Define the TLS zone names
TLS_zone_names = ['T_cell_zone', 'B_cell_zone']

score_keys_dict = {}

for score_name, filename in clusters_of_i.items():
    out_name = f'{score_name}_score'
    GC_res_dir = os.path.join(res_dir, '04_vascular_microenvironment', '01_identify_HEV_niches', '01_score_signature', out_name, resolution)
    most_recent_dir_path = GC_res_dir
    # Read results and add to adata
    GC_file_name = f'{sample_id}_{out_name}.csv'
    GC_results = pd.read_csv(os.path.join(most_recent_dir_path, 'csv_files', GC_file_name), index_col=0)
    # Choose number of genes used for scoring
    top_gene = 20
    GC_key = f'{score_name}_score_{top_gene}'
    score_keys_dict[score_name] = GC_key
    adata.obs[GC_key] = GC_results[GC_key]


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
    #         vals[:, 0] = np.linspace(1, R / 255, N * M)
    #         vals[:, 1] = np.linspace(1, G / 255, N * M)
    #         vals[:, 2] = np.linspace(1, B / 255, N * M)
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


# --- Prepare for saving the plots

qc_dir = os.path.join(res_dir, '04_vascular_microenvironment', '01_identify_HEV_niches', '02a_group_HEV_bins', resolution)
os.makedirs(qc_dir, exist_ok=True)


# --- Prepare for grouping: figure out which spots above threshold

threshold = 0.5

TLS_keys = [score_keys_dict.get(key, None) for key in TLS_zone_names]

adata.obs['TLS_zone_above_thresh'] = adata.obs[TLS_keys].gt(threshold).any(axis=1)
adata.obs['TLS_zone_above_thresh'] = adata.obs['TLS_zone_above_thresh'].astype('category')

# Also make a filtered object with HEV above threshold
adata.obs[f'{score_of_i}_above_thresh'] = adata.obs[score_keys_dict[score_of_i]] > threshold
adata.obs[f'{score_of_i}_above_thresh'] = adata.obs[f'{score_of_i}_above_thresh'].astype('category')

adata_sub = adata[adata.obs[f'{score_of_i}_above_thresh'] == True]


# --- Plot 1: Plot the TLS model score + spots above a certain threshold

# Define custom colormap for plotting TLS model results
# Set everything under certain threshold as the lowest color of a color map

cmap = "YlOrRd"
color_list = mcp.gen_color(cmap=cmap, n=11)
color_list_extracted = ['lightgray'] + color_list[2:]

n_bins = 100
cmap_name = 'my_list'
cmap_adj = LinearSegmentedColormap.from_list(cmap_name, color_list_extracted, N=n_bins)
cmap_adj.set_under(color='lightgray', alpha=0.4)


# Plotting
sc.set_figure_params(fontsize=20)

nrows = 1
ncols = 3

fig, axs = plt.subplots(nrows, ncols, figsize=(8*ncols, 6.5*nrows))

# Spots above threshold
sc.pl.spatial(adata, color='TLS_zone_above_thresh', groups=[True], bw=False, alpha_img=0.7, crop_coord=coords, title=f'{sample_id} ({subtype})\nTLS zone score > {threshold}', ax=axs[0],
              legend_loc=None, palette=['#1f77b4'])

# Spots above threshold
sc.pl.spatial(adata, color='TLS_zone_above_thresh', groups=[True], bw=False, alpha_img=0.7, crop_coord=coords, title=f'{sample_id} ({subtype})\nTLS zone score > {threshold}', ax=axs[1],
              legend_loc=None, palette=['#1f77b4'])
sc.pl.spatial(adata, color=f'{score_of_i}_above_thresh', groups=[True], bw=False, alpha_img=0, crop_coord=coords, title=f'{sample_id} ({subtype})\n{score_keys_dict[score_of_i]} > {threshold}', ax=axs[1],
              legend_loc=None, palette=['indigo'])


# score results
vmax = 0.8

# Then all scores in single plot
# Use scanpy because squidpy ignores info abt alpha in cmap
score_names = list(clusters_of_i.keys())
score_names_r = score_names[::-1]

for n, score_name in enumerate(score_names_r):
    if n == 0:
        alpha_img = 0.7
    else: 
        alpha_img = 0
    cmap = cmap_dict[score_name]
    sc.pl.spatial(adata, color=f'{score_name}_score_{top_gene}', alpha_img=alpha_img, cmap=cmap,
        crop_coord=coords, title=f'{sample_id} ({subtype})\nTLS zones', ax=axs[2], vmin=0, vmax=vmax, colorbar_loc=None,)


# Add colorbars per score
for i, (label, cmap) in enumerate(cmap_dict.items()):
    norm = mpl.colors.Normalize(vmin=0, vmax=vmax)
    spacing = 0.13
    mid_y = 0.5
    num_cbars = len(cmap_dict)
    y_pos = mid_y + (i - (num_cbars - 1) / 2) * spacing  # Distribute around the midpoint
    cbar_ax = axs[2].inset_axes([1.1, y_pos, 0.3, 0.05])
    cbar = fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), 
                        cax=cbar_ax, orientation='horizontal')
    cbar.ax.set_ylabel(f'{score_keys_dict[label]}', rotation=0, ha='left', va='center', fontsize=18)
    cbar.ax.yaxis.set_label_position("right")
    cbar.ax.yaxis.set_label_coords(1.05, 1.05)
    cbar.ax.tick_params(labelsize=15)

plt.suptitle(f'{sample_id} ({subtype}) - TLS_zone scores (top_gene={top_gene}, threshold={threshold}, resolution={resolution})', y=1.1, fontweight='bold', x=0.5)


# --- Define groups of neighboring spots

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

# Enlarge the neighbour radius slightly so all spots are captured; the exact
# radius occasionally misses spots near region boundaries.
radius_max = radius_max * 1.1

# Use squidpy to define which TLS-high spots are neighbors
sq.gr.spatial_neighbors(adata_sub, coord_type='generic', radius=(0, radius_max))

adata_sub.obs['in_tissue'] = adata_sub.obs['in_tissue'].astype('category')


# --- Plot 2: Plot scanpy's neighborhood graph result

sc.set_figure_params(fontsize=15, figsize=(6,6))

sq.pl.spatial_scatter(
    adata_sub,
    color="in_tissue",
    connectivity_key="spatial_connectivities",
    edges_color="black",
    shape='square',
    edges_width=0.4,
    size=0.5,
    palette='Dark2',
    crop_coord=coord_squidpy,
)

plt.suptitle(f'{sample_id} ({subtype}) - neighborhood graph of spots with {score_keys_dict[score_of_i]} score >{threshold}) (resolution={resolution})', y=1.1, fontweight='bold', x=0.5)


# --- Now annotate the spots as groups

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
list_of_lists = extract_neighbor_groups(connected_groups, obs_key='group')


# --- Plot 2b: plot the neighborhood graph again after grouping

sc.set_figure_params(fontsize=15, figsize=(6,6))

group_key = 'group'

squidpy_size = 4
squidpy_edge = 0.5

sq.pl.spatial_scatter(adata_sub, color=group_key, connectivity_key="spatial_connectivities", edges_color="black", crop_coord=coord_squidpy,
                    palette='viridis', edges_width=squidpy_edge, title='Squidpy spatial connectivities', img_alpha=0.5, outline=True,
                    outline_color=('black', 'black'), outline_width=(0, 0.05), shape='square')

plt.suptitle(f'{sample_id} ({subtype}) - neighborhood graph of spots with {score_keys_dict[score_of_i]} score >{threshold})', y=1.1, fontweight='bold', x=0.5)


# --- Plot 3: Plot the grouping of high GC spots

sc.set_figure_params(fontsize=12, figsize=(6,6))

# Convert to categorical to avoid a downstream plotting error
adata_sub.obs[group_key] = adata_sub.obs[group_key].astype('category')

sq.pl.spatial_scatter(adata_sub, color=group_key, img_alpha=1, crop_coord=coord_squidpy, title=f'{sample_id} - groups of spots >{threshold}', 
              shape='square', legend_loc=None)

plt.suptitle(f'{sample_id} ({subtype}) - grouping of spots with {score_keys_dict[score_of_i]} score >{threshold})', y=1.05, fontweight='bold', x=0.5)


# --- Only keep groups with a minimum number of spots 

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

# Remove the singletons
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


# --- Plot 5: Groups after removing small groups

sc.set_figure_params(fontsize=12)

# Convert to categorical to avoid a downstream plotting error
adata_sub.obs['group_too_small'] = adata_sub.obs['group_too_small'].astype('category')

sq.pl.spatial_scatter(adata_sub, color='group_too_small', img_alpha=0.5, crop_coord=coord_squidpy, title=f'{sample_id} - groups > {threshold}', shape='square',
                      legend_loc=None, figsize=(10,10))

plt.suptitle(f'{sample_id} ({subtype}) - grouping of spots with TLS zone score >{threshold} - groups too small removed (min_spots={min_spots})', y=1.05, fontweight='bold', x=0.5)


# --- Save everything so far

threshold_adj = str(threshold).split('.')[-1]
out_dir = os.path.join(qc_dir, f'threshold_{threshold_adj}',)

os.makedirs(out_dir, exist_ok=True)

# Save plots
out_dir_plots = os.path.join(out_dir, 'plots', 'groups')
os.makedirs(out_dir_plots, exist_ok=True)

# Also save the groupings as a .csv
out_dir_csv = os.path.join(out_dir, 'csv_files')
os.makedirs(out_dir_csv, exist_ok=True)

# Save plots
filename = f'{sample_id}_{score_of_i}_groups.pdf'

save_multi_image(os.path.join(out_dir_plots, filename))
plt.close('all')

# And save the groupings as a csv
to_save_df = adata_sub.obs[[score_keys_dict[score_of_i], 'group', 'group_too_small', 'group_too_small_removed']]

filename = f'{sample_id}_{score_of_i}_groups.csv'
to_save_df.to_csv(os.path.join(out_dir_csv, filename))

# Save the grouping settings as a JSON file
# Save in a dict
settings = {'score_threshold': threshold,
            'min_spots': min_spots,
            'score_keys': score_keys_dict[score_of_i],
            }

with open(os.path.join(out_dir_csv, f'{sample_id}_grouping_settings.json'), 'w') as json_file:
    json.dump(settings, json_file, indent=4)  # indent=4 makes the JSON file more readable


# --- Finish

print("Script successfully completed")
