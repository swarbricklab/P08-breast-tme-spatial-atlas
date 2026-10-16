#!/usr/bin/env python3

# --- Aim
# Plot the expression of different LRs across different labels used for LRi analysis
# To see if lack of interaction is because 1 or both components are not expressed
# Also creates spot grouping CSV file needed by downstream scripts

# --- Define environment
# env: vis-hd

# --- Load packages
import argparse
import os
import sys

import anndata as ad
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import squidpy as sq
import re
import math

sc.set_figure_params(facecolor="white", figsize=(8, 8))
sc.settings.verbosity = 3

mpl.rcParams['pdf.fonttype'] = 42

# --- Argparse arguments

parser = argparse.ArgumentParser(description="Plot LR gene expression and create spot grouping files",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

# Required
parser.add_argument("--res_dir", help="Path to results directory")
parser.add_argument("--h5ad_dir", help="Path to directory for .h5ad output files")
parser.add_argument("--project_dir", help="Project directory (root)")
parser.add_argument("--sample_id", help="Sample ID")
parser.add_argument("--sample_sheet", help="Path to sample sheet CSV file", default=None)
parser.add_argument("--bin_size", type=int, default=25,
                    help="Bin size (µm) of the 01_run_LR_analysis_stlearn.py / 02_analyze_LR_all.py runs to plot")

args = vars(parser.parse_args())


# --- Arguments

res_dir = args["res_dir"]
h5ad_dir = args["h5ad_dir"]
sample_id = args["sample_id"]
project_dir = args["project_dir"]
sys.path.append(os.path.join(project_dir, 'scripts'))
from helper_functions import auto_crop_scanpy, auto_crop, save_multi_image


# --- Get some info from the sample sheet

# Use sample_sheet argument if provided, else default
if args.get("sample_sheet"):
    sample_sheet_path = args["sample_sheet"]
else:
    sample_sheet_path = os.path.join(project_dir, 'config', 'sample_sheet_visium_hd.csv')

sample_sheet = pd.read_csv(sample_sheet_path)
tissue_id = sample_sheet[sample_sheet['section_id'] == sample_id]['tissue_id'].values[0]


# --- Also get some tissue metadata

resource_dir = os.path.join(project_dir, 'resources')
metadata_file = 'tissue_metadata.csv'

tissue_metadata = pd.read_csv(os.path.join(resource_dir, metadata_file))

# Create a sample - subtype dictionary
subtype = tissue_metadata[tissue_metadata['tissue_id'] == str(tissue_id)]['clinical_subtype'].values[0]


# --- Load the adata object

resolution = 'square_016um'

adata = ad.read_h5ad(os.path.join(h5ad_dir, 'visium_hd', resolution, f'{sample_id}.h5ad'))

# The adata files have log-normalized data in the X layer (no separate layer)
use_data = None


# --- Calculate crop coords for optimal presentation

coords = auto_crop_scanpy(adata, border_size=0.1)

# Also for squidpy
coords_sq = auto_crop(adata, border_size=0.1)


# --- Load the stlearn adata

# To make it easier to visualize trends by lineage
class_dir = os.path.join(project_dir, 'resources', 'cell_types', 'miniatlas')
class_file = 'miniatlas_subset_classifications.csv'

class_df = pd.read_csv(os.path.join(class_dir, class_file))


# --- Load HEV groups

bin_size = args["bin_size"]

lr_sample_dir = os.path.join(res_dir, '04_vascular_microenvironment', '04_ligand_receptor_analysis',
                              '01_run_LR_analysis_stlearn', f'{bin_size}um_bin', 'csv_files')

most_recent_dir_path = lr_sample_dir

adata_file = f'{sample_id}_HEV_LR_analysis_stlearn.h5ad'
stlearn_adata = ad.read_h5ad(os.path.join(most_recent_dir_path, adata_file))


# --- Also load the LRs with LFCs

lr_all_dir = os.path.join(res_dir, '04_vascular_microenvironment', '04_ligand_receptor_analysis',
                           '02_analyze_LR_all', f'{bin_size}um_bin', 'csv_files')

most_recent_dir_path = lr_all_dir

# Use the spatial group level results
filename = 'All_samples_group_counts.csv'

lr_summary = pd.read_csv(os.path.join(most_recent_dir_path, filename))

# Filter to LRs that are significant in all samples
n_samples = lr_summary['n_samples_sig'].max()
lr_summary = lr_summary[lr_summary['n_samples_sig'] == n_samples]

# Then sort by LFC
lr_summary = lr_summary.sort_values(by='LFC', ascending=False)


# --- Prepare for saving

out_dir = os.path.join(res_dir, '04_vascular_microenvironment', '04_ligand_receptor_analysis',
                       '03_plot_LR_per_sample', f'{bin_size}um_bin')

os.makedirs(out_dir, exist_ok=True)

# Save plots
out_dir_plots = os.path.join(out_dir, 'plots')
os.makedirs(out_dir_plots, exist_ok=True)

# Also save the groupings as a .csv
out_dir_csv = os.path.join(out_dir, 'csv_files')
os.makedirs(out_dir_csv, exist_ok=True)


# --- Extract spot labels

spot_labels = stlearn_adata.obs['cell_type']

# Add the spot labels to the adata object
adata.obs['spot_labels'] = spot_labels


# --- Define a color_dict for the spot labels

# Get a colordict to make sure all colors align across samples
colors = stlearn_adata.uns['cell_type_colors']
labels = stlearn_adata.obs['cell_type'].cat.categories

color_dict = {}
for label, color in zip(labels, colors):
    color_dict[label] = color


# --- Redefine groups also including the nearby spots
# Then afterwards will asign near TLS vs non TLS based on dominant annotation

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

# Do spatial gropuping per sample
adata_sub = adata[adata.obs['spot_labels'].notna()]

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

# Extract connected spots as groups
connected_groups = []
for i, neighbors in enumerate(adata_sub.obsp['spatial_connectivities'].toarray()):
    group = [j for j, is_neighbor in enumerate(neighbors) if is_neighbor]
    group.append(i)  # Include the spot itself in the group
    connected_groups.append(group)

# Extract groups
list_of_lists = extract_neighbor_groups(connected_groups, adata=adata_sub, obs_key='group')

# Name groups according to size
group_counts = adata_sub.obs['group'].value_counts()
group_counts = group_counts[group_counts > 1]  # Only keep groups with more than one spot

group_names = {group: f'Group {i+1}' for i, group in enumerate(group_counts.index)}

# Assign group names to the adata object
adata_sub.obs['group'] = adata_sub.obs['group'].map(group_names)

# Add to original adata
adata.obs['new_grouping'] = adata_sub.obs['group']


# --- Assign categories to groups based on dominant HEV label

# Remove the nearby_spots label
HEV_labels = [l for l in labels if l != 'nearby_spots']

# Define a function to assign categories based on dominant HEV label
def assign_HEV_category(group, adata=adata_sub, obs_key='cell_type'):
    # Get the HEV labels for the spots in the group
    adata_sub = adata[adata.obs['new_grouping'] == group]
    # Count the occurrences of each HEV label in the group
    hev_counts = adata_sub.obs[obs_key].value_counts()
    # Figure out which of the HEV labels is dominant
    hev_counts = hev_counts[hev_counts.index.isin(HEV_labels)]
    # If there are no HEV labels in the group, return 'unknown'
    if hev_counts.empty:
        dominant_label = 'unknown'
    else:
        dominant_label = hev_counts.idxmax()
    return dominant_label

# Assign dominant HEV label to each group
adata_sub = adata[adata.obs['new_grouping'].notna()]

adata_sub.obs['dominant_HEV_label'] = adata_sub.obs['new_grouping'].apply(assign_HEV_category, adata=adata_sub, obs_key='spot_labels')
adata.obs['dominant_HEV_label'] = adata_sub.obs['dominant_HEV_label']


# --- Adjust nearby_spots labels

# Now add the dominant_HEV_label to the nearby_spots
adata.obs['nearby_spots_revisited'] = np.where(adata.obs['spot_labels'] == 'nearby_spots',
                                            'nearby_' + adata.obs['dominant_HEV_label'],    
                                            adata.obs['spot_labels'])

categories = ['near_TLS_HEV', 'non_TLS_HEV']

# Derive lighter colors for the nearby spots
def derive_lighter_colors(og_color, factor=0.5):
    rgba = mpl.colors.to_rgba(og_color)
    lighter_color = tuple(min(1, c + (1 - c) * factor) for c in rgba[:3]) + (rgba[3],)
    return lighter_color

for cat in categories:
    og_color = color_dict[cat]
    new_color = derive_lighter_colors(og_color, factor=0.5)
    # Convert to hex
    new_color_hex = mpl.colors.to_hex(new_color)
    color_dict[f'nearby_{cat}'] = new_color_hex


# --- Plot 1: Spatial plot of the spot labels

sc.set_figure_params(fontsize=15)

ncols = 3
nrows = 1

fig, axs = plt.subplots(ncols=ncols, nrows=nrows, figsize=(8*ncols, 8*nrows))

# -- First plot the spot labels
n=0
ax = axs[n]

# Extract a listed colormap from color_dict
color_list = [color_dict[label] for label in adata.obs['spot_labels'].astype('category').cat.categories if label in color_dict]
# Create a colormap from the list
cmap = mpl.colors.ListedColormap(color_list)

adata_sub = adata[adata.obs['spot_labels'].notna()]
sq.pl.spatial_scatter(adata_sub, color='spot_labels', img_alpha=0.5, crop_coord=coords_sq,   
                title=f'{sample_id} ({subtype})\nSpot labels', shape='square',
                palette=cmap, ax=ax, legend_loc='right margin')

ax.set_title(ax.get_title(), fontweight='bold')

# -- Then plot the new grouping
n+=1
ax = axs[n]

# Extract a listed colormap from color_dict
color_list = [color_dict[label] for label in adata.obs['dominant_HEV_label'].astype('category').cat.categories if label in color_dict]
# Create a colormap from the list
cmap = mpl.colors.ListedColormap(color_list)

adata_sub = adata[adata.obs['dominant_HEV_label'].notna()]
sq.pl.spatial_scatter(adata_sub, color='dominant_HEV_label', img_alpha=0.5, crop_coord=coords_sq,
                title=f'{sample_id} ({subtype})\nSpot grouping', shape='square',
                palette=cmap, ax=ax, legend_loc='right margin')

ax.set_title(ax.get_title(), fontweight='bold')

# -- Then plot nearby spots adjusted
n+=1
ax = axs[n]

# Extract a listed colormap from color_dict
color_list = [color_dict[label] for label in adata.obs['nearby_spots_revisited'].astype('category').cat.categories if label in color_dict]
# Create a colormap from the list
cmap = mpl.colors.ListedColormap(color_list)

adata_sub = adata[adata.obs['nearby_spots_revisited'].notna()]
sq.pl.spatial_scatter(adata_sub, color='nearby_spots_revisited', img_alpha=0.5, crop_coord=coords_sq,
                title=f'{sample_id} ({subtype})\nSpot grouping', shape='square',
                palette=cmap, ax=ax, legend_loc='right margin')

ax.set_title(ax.get_title(), fontweight='bold')

plt.suptitle(f'{sample_id} ({subtype}) - Spot labels & grouping used for LRi stlearn', y=1, x=0.5, fontweight='bold')

plt.subplots_adjust(wspace=0.6)


# --- Save this first plot with high res

out_file = os.path.join(out_dir_plots, f'{sample_id}_spot_labels_grouping_hires.pdf')
plt.savefig(out_file, dpi=1500, bbox_inches='tight')


# --- Plot the genes of the top & bottom 10 LRs

top = 10

top_lrs = lr_summary.head(top)['lr_pair'].tolist()
bottom_lrs = lr_summary.tail(top)['lr_pair'].tolist()

all_lrs = top_lrs + bottom_lrs

# Make a dict of LR pairs and genes involved
top_genes = {}
for lr in all_lrs:
    gene_a, gene_b = lr.split('_')
    top_genes[lr] = [gene_a, gene_b]

# Extract the LFCs from the lr_summary
lr_lfc = lr_summary.set_index('lr_pair')['LFC'].to_dict()

# Find the highest & lowest LFCs for the LRs of interest
max_lfc = max(lr_lfc[lr] for lr in lr_lfc)
min_lfc = min(lr_lfc[lr] for lr in lr_lfc)

most_extreme = max(abs(max_lfc), abs(min_lfc))

# Make a colormap for the LFCs
cmap = mpl.cm.coolwarm
norm = mpl.colors.Normalize(vmin=-most_extreme, vmax=most_extreme)


# --- Plot 2: Dotplot of genes corresponding to top & bottom LRs

sc.set_figure_params(fontsize=10)

dp = sc.pl.dotplot(adata, var_names=top_genes, groupby=['nearby_spots_revisited'], layer=use_data, return_fig=True)
dp.add_totals().style(dot_edge_color='black', dot_edge_lw=0.5, cmap='Reds').show()

# Get axes
axs = plt.gcf().get_axes()
# Get the legend ax & adjust the legend
dots_ax = axs[4]
dots_ax.set_title('Fraction of spots\nin group (%)', y=0.5, fontsize=8.5)

# Color the bars according to the group label
bar_ax = axs[1]
main_ax = axs[0]

for i, bar in enumerate(bar_ax.patches):
    # Get the group label from the bar's label
    group_label = main_ax.get_yticklabels()[i].get_text()
    # Get the color for the group label
    color = color_dict.get(group_label, 'lightgray')  # Default to light gray if not found
    bar.set_facecolor(color)

# Add a dotted line after the second group to separate TLS from non-TLS
# Get ylabels texts
y_labels = [label.get_text() for label in main_ax.get_yticklabels()]

# Find the index where it goes from near_TLS to non_TLS
# This is the first non_TLS label
for i, label in enumerate(y_labels):
    if 'near_TLS' not in label:
        line_index = i
        break

main_ax.axhline(y=line_index, color='lightgray', linestyle=(0, (3,3)), linewidth=1, zorder=-1)

# Add a dotted line to separate the top & bottom LRs
main_ax.axvline(x=top*2, color='lightgray', linestyle='-', linewidth=2)

# Add small rectangles on the left to indicate the group labels
for i, label in enumerate(main_ax.get_yticklabels()):
    # Make sure the first & the last extend to the edges of the plot
    if i == 0:
        y_start = -0.6
        width = 0.6
    elif i == len(main_ax.get_yticklabels()) - 1:
        y_start = 0
        width = 0.6
    else:
        width = 0
        y_start = 0
    group_label = label.get_text()
    if group_label in color_dict:
        color = color_dict[group_label]
        # Add a rectangle to the left of the y-tick label
        main_ax.add_patch(mpl.patches.Rectangle((-0.5, i+0.1+y_start), 0.5, 0.8 + width, color=color, transform=main_ax.transData))

# Desired rectangle height in pixels
rect_height_px = 5

# Get the axes bounding box in display (pixel) coordinates
bbox = main_ax.get_window_extent()

# Convert pixel height to axes coordinates
height_axes = rect_height_px / bbox.height

# Add rectangles showing LFCs
for i, lr in enumerate(top_genes):
    lfc = lr_lfc[lr]
    x_pos = i * 2
    color = cmap(norm(lfc))
    # Rectangle just above the axes
    rect = mpl.patches.Rectangle(
        (x_pos + 0.25, 1.01),  # x in data coords, y just above axes
        1.4, height_axes,      # width and height in data coords and axes coords
        transform=mpl.transforms.blended_transform_factory(main_ax.transData, main_ax.transAxes),
        color=color, alpha=1, clip_on=False, zorder=-10,
    )
    main_ax.add_patch(rect)

# Add a colorbar for the LFCs
# Duplicate the other cbar ax
cbar_ax = axs[5]
new_cbar_ax = plt.gcf().add_axes([0.92, 0.15, 0.02, 0.7])  # [left, bottom, width, height]
cbar = mpl.colorbar.ColorbarBase(new_cbar_ax, cmap=cmap, norm=norm,
                                 orientation='horizontal')

# Set position based on other cbar
cbar_pos = cbar_ax.get_position()
new_cbar_ax.set_position([cbar_pos.x0, cbar_pos.y1 + cbar_pos.height*11.5, cbar_pos.width, cbar_pos.height])

new_cbar_ax.set_title('LFC\n(TLS / non-TLS)', fontsize=8.5)
cbar.ax.tick_params(labelsize=8.5)

# Limit the x-axis to extreme values
# Find the highest & lowest LFCs for the LRs of interest
max_lfc = max(lr_lfc[lr] for lr in top_genes)
min_lfc = min(lr_lfc[lr] for lr in top_genes)

new_cbar_ax.set_xlim(min_lfc, max_lfc)

plt.suptitle(f'{sample_id} ({subtype}) - Dotplot of top & bottom {top} LRs (found in min {n_samples} samples)', y=1.2, x=0.5, fontweight='bold')


# --- Now plot some LRs of interest

lr_list = lr_summary['lr_pair'].tolist()

# Find all LRs that contain particular substrings
substrings = ['CCL', 'CCR', 'CXCL', 'CXCR','SEL', 'SPP1', 'LTB']
lr_interest = []

for lr in lr_list:
    if any(sub in lr for sub in substrings):
        lr_interest.append(lr)

# Make a dict of LR pairs and genes involved
lr_interest_genes = {}
for lr in lr_interest:
    gene_a, gene_b = lr.split('_')
    lr_interest_genes[lr] = [gene_a, gene_b]


# --- Plot 3: Dotplot LRs of interest

dp = sc.pl.dotplot(adata, var_names=lr_interest_genes, groupby=['nearby_spots_revisited'], layer=use_data, return_fig=True)
dp.add_totals().style(dot_edge_color='black', dot_edge_lw=0.5, cmap='Reds').show()

# Get axes
axs = plt.gcf().get_axes()
# Get the legend ax & adjust the legend
dots_ax = axs[4]
dots_ax.set_title('Fraction of spots\nin group (%)', y=0.5, fontsize=8.5)

# Color the bars according to the group label
bar_ax = axs[1]
main_ax = axs[0]

for i, bar in enumerate(bar_ax.patches):
    # Get the group label from the bar's label
    group_label = main_ax.get_yticklabels()[i].get_text()
    # Get the color for the group label
    color = color_dict.get(group_label, 'lightgray')  # Default to light gray if not found
    bar.set_facecolor(color)

# Add a dotted line after the second group to separate TLS from non-TLS
# Get ylabels texts
y_labels = [label.get_text() for label in main_ax.get_yticklabels()]

# Find the index where it goes from near_TLS to non_TLS
# This is the first non_TLS label
for i, label in enumerate(y_labels):
    if 'near_TLS' not in label:
        line_index = i
        break

main_ax.axhline(y=line_index, color='lightgray', linestyle=(0, (3,3)), linewidth=1, zorder=-1)

# Add small rectangles on the left to indicate the group labels
for i, label in enumerate(main_ax.get_yticklabels()):
    # Make sure the first & the last extend to the edges of the plot
    if i == 0:
        y_start = -0.6
        width = 0.6
    elif i == len(main_ax.get_yticklabels()) - 1:
        y_start = 0
        width = 0.6
    else:
        width = 0
        y_start = 0
    group_label = label.get_text()
    if group_label in color_dict:
        color = color_dict[group_label]
        # Add a rectangle to the left of the y-tick label
        main_ax.add_patch(mpl.patches.Rectangle((-0.5, i+0.1+y_start), 0.5, 0.8 + width, color=color, transform=main_ax.transData))

# Add rectangles showing LFCs
first_neg_LFC = None

for i, lr in enumerate(lr_interest):
    lfc = lr_lfc[lr]
    if lfc < 0 and first_neg_LFC is None:
        first_neg_LFC = i # Store the index of the first negative LFC
    x_pos = i * 2
    color = cmap(norm(lfc))
    # Rectangle just above the axes
    rect = mpl.patches.Rectangle(
        (x_pos + 0.25, 1.01),  # x in data coords, y just above axes
        1.4, height_axes,      # width and height in data coords and axes coords
        transform=mpl.transforms.blended_transform_factory(main_ax.transData, main_ax.transAxes),
        color=color, alpha=1, clip_on=False, zorder=-10,
    )
    main_ax.add_patch(rect)

# Add a colorbar for the LFCs
# Duplicate the other cbar ax
cbar_ax = axs[5]
new_cbar_ax = plt.gcf().add_axes([0.92, 0.15, 0.02, 0.7])  # [left, bottom, width, height]
cbar = mpl.colorbar.ColorbarBase(new_cbar_ax, cmap=cmap, norm=norm,
                                 orientation='horizontal')

# Set position based on other cbar
cbar_pos = cbar_ax.get_position()
new_cbar_ax.set_position([cbar_pos.x0, cbar_pos.y1 + cbar_pos.height*11.5, cbar_pos.width, cbar_pos.height])

new_cbar_ax.set_title('LFC\n(TLS / non-TLS)', fontsize=8.5)
cbar.ax.tick_params(labelsize=8.5)

# Limit the x-axis to extreme values
# Find the highest & lowest LFCs for the LRs of interest
max_lfc = max(lr_lfc[lr] for lr in lr_interest_genes)
min_lfc = min(lr_lfc[lr] for lr in lr_interest_genes)

new_cbar_ax.set_xlim(min_lfc, max_lfc)

# Add a vertical line to indicate where the LFC becomes negative
if first_neg_LFC is not None:
    main_ax.axvline(x=first_neg_LFC*2-0.05, color='lightgray', linestyle='-', linewidth=2)

plt.suptitle(f'{sample_id} ({subtype}) - Dotplot LRs of interest ({substrings}) - sorted by LFC (found in min {n_samples} samples)', y=1.2, x=0.5, fontweight='bold')


# ---  Now also average gene expression by spatial group

# Extract all genes in the lr_list as a single unnested list
# Just to speed things up
gene_list = [gene for lr in lr_list for gene in lr.split('_')]

# Remove duplicates
gene_list = list(set(gene_list))

# Get an unnested list of genes
expr_df = adata_sub.to_df(layer=use_data)[gene_list]

# Add composite group column
adata_sub = adata_sub[adata_sub.obs['new_grouping'].notna()]
expr_df['group'] = adata_sub.obs['new_grouping'].values

# And also add the spot labels
expr_df['nearby_spots_revisited'] = adata_sub.obs['nearby_spots_revisited'].values

# Compute mean expression per composite group
grouped_expr = expr_df.groupby(['group', 'nearby_spots_revisited']).mean()

# Drop any rows with NaN values
grouped_expr = grouped_expr.dropna()

adata_grouped = ad.AnnData(X=grouped_expr.values)

# Collapse the multiindex to a single index
grouped_expr.index = grouped_expr.index.map(lambda x: f"{x[0]}-{x[1]}")
adata_grouped.obs_names = grouped_expr.index
adata_grouped.var_names = grouped_expr.columns
vessel_classes_df = grouped_expr.index.to_series().str.split('-', expand=True)[1]

adata_grouped.obs['nearby_spots_revisited'] = vessel_classes_df


# --- Plot 4: Spatial plot showing grouping

sc.set_figure_params(fontsize=15)

ncols = 2
nrows = 1

fig, axs = plt.subplots(ncols=ncols, nrows=nrows, figsize=(8*ncols, 8*nrows))

# -- First plot the revisited spot labels
n=0
ax = axs[n]

adata_sub = adata[adata.obs['nearby_spots_revisited'].notna()]
sc.pl.spatial(adata_sub, color='nearby_spots_revisited', bw=False, alpha_img=0.5, crop_coord=coords,
                title=f'{sample_id} ({subtype})\nSpot grouping',    
                palette=color_dict, ax=ax, legend_loc='right margin')

ax.set_title(ax.get_title(), fontweight='bold')


# -- Then plot the new grouping
n+=1
ax = axs[n]

adata_sub = adata[adata.obs['new_grouping'].notna()]
sc.pl.spatial(adata_sub, color='new_grouping', bw=False, alpha_img=0.5, crop_coord=coords,
                title=f'{sample_id} ({subtype})\nSpot grouping',    
                palette='Set1', ax=ax, legend_loc=None)

ax.set_title(ax.get_title(), fontweight='bold')

plt.subplots_adjust(wspace=0.6)

plt.suptitle(f'{sample_id} ({subtype}) - Spot labels & grouping for collapsing of gene expression', y=1, x=0.5, fontweight='bold')


# --- Plot 5: Dotplot of top & bottom LRs averaged across spatial groups

sc.set_figure_params(fontsize=10)

dp = sc.pl.dotplot(adata_grouped, var_names=top_genes, groupby=['nearby_spots_revisited'], return_fig=True)
dp.add_totals().style(dot_edge_color='black', dot_edge_lw=0.5, cmap='Reds').show()

# Get axes
axs = plt.gcf().get_axes()

# Get the legend ax & adjust the legend
dots_ax = axs[4]
dots_ax.set_title('Fraction of spatial groups\nexpressing (%)', y=0.5, fontsize=8.5)

# Get the legend ax & adjust the legend
cbar_ax = axs[5]
cbar_ax.set_title('Mean expression\nacross spatial groups', y=1, fontsize=8.5)

# Color the bars according to the group label
bar_ax = axs[1]
main_ax = axs[0]

for i, bar in enumerate(bar_ax.patches):
    # Get the group label from the bar's label
    group_label = main_ax.get_yticklabels()[i].get_text()
    # Get the color for the group label
    color = color_dict.get(group_label, 'lightgray')  # Default to light gray if not found
    bar.set_facecolor(color)

# Add a dotted line after the second group to separate TLS from non-TLS
# Get ylabels texts
y_labels = [label.get_text() for label in main_ax.get_yticklabels()]

# Find the index where it goes from near_TLS to non_TLS
# This is the first non_TLS label
for i, label in enumerate(y_labels):
    if 'near_TLS' not in label:
        line_index = i
        break

main_ax.axhline(y=line_index, color='lightgray', linestyle=(0, (3,3)), linewidth=1, zorder=-1)

# Add a dotted line to separate the top & bottom LRs
main_ax.axvline(x=top*2, color='lightgray', linestyle='-', linewidth=2)

# Add small rectangles on the left to indicate the group labels
for i, label in enumerate(main_ax.get_yticklabels()):
    # Make sure the first & the last extend to the edges of the plot
    if i == 0:
        y_start = -0.6
        width = 0.6
    elif i == len(main_ax.get_yticklabels()) - 1:
        y_start = 0
        width = 0.6
    else:
        width = 0
        y_start = 0
    group_label = label.get_text()
    if group_label in color_dict:
        color = color_dict[group_label]
        # Add a rectangle to the left of the y-tick label
        main_ax.add_patch(mpl.patches.Rectangle((-0.5, i+0.1+y_start), 0.5, 0.8 + width, color=color, transform=main_ax.transData))

# Add rectangles showing LFCs
for i, lr in enumerate(top_genes):
    lfc = lr_lfc[lr]
    x_pos = i * 2
    color = cmap(norm(lfc))
    # Rectangle just above the axes
    rect = mpl.patches.Rectangle(
        (x_pos + 0.25, 1.01),  # x in data coords, y just above axes
        1.4, height_axes,      # width and height in data coords and axes coords
        transform=mpl.transforms.blended_transform_factory(main_ax.transData, main_ax.transAxes),
        color=color, alpha=1, clip_on=False, zorder=-10,
    )
    main_ax.add_patch(rect)

# Add a colorbar for the LFCs
# Duplicate the other cbar ax
cbar_ax = axs[5]
new_cbar_ax = plt.gcf().add_axes([0.92, 0.15, 0.02, 0.7])  # [left, bottom, width, height]
cbar = mpl.colorbar.ColorbarBase(new_cbar_ax, cmap=cmap, norm=norm,
                                 orientation='horizontal')

# Set position based on other cbar
cbar_pos = cbar_ax.get_position()
new_cbar_ax.set_position([cbar_pos.x0, cbar_pos.y1 + cbar_pos.height*11.5, cbar_pos.width, cbar_pos.height])

new_cbar_ax.set_title('LFC\n(TLS / non-TLS)', fontsize=8.5)
cbar.ax.tick_params(labelsize=8.5)

# Limit the x-axis to extreme values
# Find the highest & lowest LFCs for the LRs of interest
max_lfc = max(lr_lfc[lr] for lr in top_genes)
min_lfc = min(lr_lfc[lr] for lr in top_genes)

new_cbar_ax.set_xlim(min_lfc, max_lfc)

plt.suptitle(f'{sample_id} ({subtype}) - Dotplot of top & bottom {top} LRs - expression averaged across spatial groups (found in min {n_samples} samples)', y=1.2, x=0.5, fontweight='bold')


# --- Plot 6: Dotplot LRs of interest

dp = sc.pl.dotplot(adata_grouped, var_names=lr_interest_genes, groupby=['nearby_spots_revisited'], return_fig=True)
dp.add_totals().style(dot_edge_color='black', dot_edge_lw=0.5, cmap='Reds').show()

# Get axes
axs = plt.gcf().get_axes()
# Get the legend ax & adjust the legend
dots_ax = axs[4]
dots_ax.set_title('Fraction of spatial groups\nexpressing (%)', y=0.5, fontsize=8.5)

# Get the legend ax & adjust the legend
cbar_ax = axs[5]
cbar_ax.set_title('Mean expression\nacross spatial groups', y=1, fontsize=8.5)

# Color the bars according to the group label
bar_ax = axs[1]
main_ax = axs[0]

for i, bar in enumerate(bar_ax.patches):
    # Get the group label from the bar's label
    group_label = main_ax.get_yticklabels()[i].get_text()
    # Get the color for the group label
    color = color_dict.get(group_label, 'lightgray')  # Default to light gray if not found
    bar.set_facecolor(color)

# Add a dotted line after the second group to separate TLS from non-TLS
# Get ylabels texts
y_labels = [label.get_text() for label in main_ax.get_yticklabels()]

# Find the index where it goes from near_TLS to non_TLS
# This is the first non_TLS label
for i, label in enumerate(y_labels):
    if 'near_TLS' not in label:
        line_index = i
        break

main_ax.axhline(y=line_index, color='lightgray', linestyle=(0, (3,3)), linewidth=1, zorder=-1)

# Add small rectangles on the left to indicate the group labels
for i, label in enumerate(main_ax.get_yticklabels()):
    # Make sure the first & the last extend to the edges of the plot
    if i == 0:
        y_start = -0.6
        width = 0.6
    elif i == len(main_ax.get_yticklabels()) - 1:
        y_start = 0
        width = 0.6
    else:
        width = 0
        y_start = 0
    group_label = label.get_text()
    if group_label in color_dict:
        color = color_dict[group_label]
        # Add a rectangle to the left of the y-tick label
        main_ax.add_patch(mpl.patches.Rectangle((-0.5, i+0.1+y_start), 0.5, 0.8 + width, color=color, transform=main_ax.transData))

# Add rectangles showing LFC
first_neg_LFC = None

for i, lr in enumerate(lr_interest):
    lfc = lr_lfc[lr]
    if lfc < 0 and first_neg_LFC is None:
        first_neg_LFC = i  # Store the index of the first negative LFC
    x_pos = i * 2
    color = cmap(norm(lfc))
    # Rectangle just above the axes
    rect = mpl.patches.Rectangle(
        (x_pos + 0.25, 1.01),  # x in data coords, y just above axes
        1.4, height_axes,      # width and height in data coords and axes coords
        transform=mpl.transforms.blended_transform_factory(main_ax.transData, main_ax.transAxes),
        color=color, alpha=1, clip_on=False, zorder=-10,
    )
    main_ax.add_patch(rect)

# Add a colorbar for the LFCs
# Duplicate the other cbar ax
new_cbar_ax = plt.gcf().add_axes([0.92, 0.15, 0.02, 0.7])  # [left, bottom, width, height]
cbar = mpl.colorbar.ColorbarBase(new_cbar_ax, cmap=cmap, norm=norm,
                                 orientation='horizontal')

# Set position based on other cbar
cbar_pos = cbar_ax.get_position()
new_cbar_ax.set_position([cbar_pos.x0, cbar_pos.y1 + cbar_pos.height*11, cbar_pos.width, cbar_pos.height])

new_cbar_ax.set_title('LFC\n(TLS / non-TLS)', fontsize=8.5)
cbar.ax.tick_params(labelsize=8.5)

# Limit the x-axis to extreme values
# Find the highest & lowest LFCs for the LRs of interest
max_lfc = max(lr_lfc[lr] for lr in lr_interest_genes)
min_lfc = min(lr_lfc[lr] for lr in lr_interest_genes)
new_cbar_ax.set_xlim(min_lfc, max_lfc)

# Add a vertical line to indicate where the LFC becomes negative
if first_neg_LFC is not None:
    main_ax.axvline(x=first_neg_LFC*2-0.05, color='lightgray', linestyle='-', linewidth=2)

plt.suptitle(f'{sample_id} ({subtype}) - Dotplot LRs of interest ({substrings}) - expression averaged across spatial groups (found in min {n_samples} samples)', y=1.2, x=0.5, fontweight='bold')


# --- Save the output

# Save plots
filename = f'{sample_id}_HEV_LR_gene_expression.pdf'
save_multi_image(os.path.join(out_dir_plots, filename))

# -- Save the grouped anndata object for convenience

filename = f'{sample_id}_HEV_LR_analysis_stlearn_spatial_groups.h5ad'
ad.AnnData.write(adata_grouped, filename=os.path.join(out_dir_csv, filename), compression='gzip')

# Also save the grouping info for convenience
to_save = adata.obs[['spot_labels', 'nearby_spots_revisited', 'new_grouping']]

# Drop any rows with NaN values
to_save = to_save.dropna()
to_save.to_csv(os.path.join(out_dir_csv, f'{sample_id}_HEV_spot_grouping.csv'), index=True)


# --- Finish

print("Script successfully completed")
