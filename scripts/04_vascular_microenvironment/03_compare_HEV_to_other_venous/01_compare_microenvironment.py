#!/usr/bin/env python3

# --- Aim
# Bin spatial distances and compare pathology annotations across distance bins.
# Contrast microenvironmental features of vessels overlapping with HEVs versus non-HEV vessels.

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
import numpy as np
import pandas as pd
import scanpy as sc
import scanpy.plotting.palettes as sc_palettes
import seaborn as sns
import scipy.cluster.hierarchy as sch
import sys
import re
from scipy.spatial import cKDTree
import squidpy as sq
from matplotlib.colors import ListedColormap

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

args = vars(parser.parse_args())


# --- Arguments

res_dir = args["res_dir"]
h5ad_dir=args["h5ad_dir"]
sample_id=args["sample_id"]
project_dir=args["project_dir"]

sys.path.append(os.path.join(project_dir, 'scripts'))
from helper_functions import auto_crop_scanpy, save_multi_image, generate_color_variations, auto_crop

# Import function script as a module
sys.path.append(os.path.join(project_dir, 'scripts', '04_vascular_microenvironment'))
import plotting_functions


# --- Get some info from the sample sheet
# Use --sample_sheet if provided, otherwise default to config/sample_sheet_visium_hd.csv
if args["sample_sheet"] is not None:
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

use_data = 'log_norm'
# Need both the 48 & 16 res
resolutions = ['square_016um', 'square_008um']

adata_dict = {}

for resolution in resolutions:
    adata_dict[resolution] = ad.read_h5ad(os.path.join(h5ad_dir, 'visium_hd', resolution, f'{sample_id}.h5ad'))


# --- Calculate crop coords for optimal presentation
coords_dict  = {}

for resolution in resolutions:
    coords = auto_crop_scanpy(adata_dict[resolution], border_size=0.1)
    coords_dict[resolution] = coords

# Also calculate some coords for squidpy
coords_dict_sq  = {}

for resolution in resolutions:
    coords = auto_crop(adata_dict[resolution], border_size=0.1)
    coords_dict_sq[resolution] = coords


# --- Load miniatlas classifications

decon_level = 'minor'

# To make it easier to visualize trends by lineage
class_dir = os.path.join(project_dir, 'resources', 'cell_types', 'miniatlas')
class_file = 'miniatlas_subset_classifications.csv'

class_df = pd.read_csv(os.path.join(class_dir, class_file))

# Get endothelial classes
endothelial = list(class_df[class_df['miniatlas_major'] == 'Endothelial'][f'miniatlas_{decon_level}'].unique())
endothelial = sorted(endothelial)


# --- Load vessel groups

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

# Specify path
TLS_group_dir = os.path.join(project_dir, 'results', '04_vascular_microenvironment', '02_identify_all_vasculature', '01_classify_vessels', 'csv_files')

most_recent_dir_path = find_most_recent(TLS_group_dir)

base_res = 'square_016um'

filename = f'{sample_id}_vessel_classification.csv'
vessel_classes = pd.read_csv(os.path.join(most_recent_dir_path, filename), index_col=0)
adata_dict[base_res].obs[vessel_classes.columns] = vessel_classes

# Collapse to single vessel class per spot
collapsed_col = vessel_classes[endothelial].apply(lambda row: row.first_valid_index() if row.first_valid_index() is not None else None, axis=1)
adata_dict[base_res].obs['vessel_class'] = collapsed_col


# --- And load the distance calculations

dist_dir = os.path.join(project_dir, 'results', '04_vascular_microenvironment', '02_identify_all_vasculature', '03_calculate_distances', 'square_008um', 'csv_files')
most_recent = find_most_recent(dist_dir)

distance_dict = {}

for cell_type in endothelial:
    dist_file = f'{sample_id}_distance_closest_{cell_type.replace(" ", "_")}_per_spot.csv'
    distance_df = pd.read_csv(os.path.join(most_recent, 'min_dist', dist_file), index_col=0)
    distance_dict[cell_type] = distance_df
    # Also add to adata
    adata_dict[base_res].obs[f'dist_{cell_type}'] = distance_df['nearest_016um_dist_um']


# --- Load the pathology annotations

annotation_dir = os.path.join(project_dir, 'data', 'pathology_annotations', 'visium_hd')
filename = f'{sample_id}_pathology_annotations.csv'

annotation_key = 'path_anno'

annotations = pd.read_csv(os.path.join(annotation_dir, filename), index_col=0)
adata_dict['square_008um'].obs[annotation_key] = annotations


# --- Load HEV overlap information

# Define which endothelial type want to look at
endo_of_interest = 'Endothelial ACKR1'
out_name = endo_of_interest.split(' ')[-1]

TLS_group_dir = os.path.join(res_dir, '04_vascular_microenvironment', '02_identify_all_vasculature', '02_investigate_vessels', 'csv_files')

most_recent_dir_path = find_most_recent(TLS_group_dir)

filename = f'{sample_id}_vessel_HEV_overlap_spot_level.csv'
overlap_df = pd.read_csv(os.path.join(most_recent_dir_path, filename), index_col=0)

adata_dict[base_res].obs['overlap'] = overlap_df['overlap']

# Also load the file with proportion of overlap
filename = f'{sample_id}_vessel_HEV_overlap.csv'
proportion_df = pd.read_csv(os.path.join(most_recent_dir_path, filename), index_col=0)

# Identify groups with at least 50% overlap
min_overlap = 0.5

# Filter for the endothelial type of interest
proportion_df_endo = proportion_df[proportion_df['group_type'] == endo_of_interest].copy()

# Separate groups with overlap above vs below threshold (or NaN)
proportion_df_high = proportion_df_endo[proportion_df_endo['pct_overlap'] > min_overlap]
# Include both low overlap and NaN (vessels with 0 or no overlap data)
proportion_df_low = proportion_df_endo[(proportion_df_endo['pct_overlap'] <= min_overlap) | (proportion_df_endo['pct_overlap'].isna())]

high_overlap_groups = proportion_df_high['group'].to_list()
low_overlap_groups = proportion_df_low['group'].to_list()

# Initialize the column
adata_dict[base_res].obs['majority_overlap'] = np.nan

# Set 'overlap' for vessels with high overlap
barcode_groups_high = adata_dict[base_res].obs[adata_dict[base_res].obs[endo_of_interest].isin(high_overlap_groups)].index
adata_dict[base_res].obs.loc[barcode_groups_high, 'majority_overlap'] = 'overlap'

# Set vessel type name for vessels with low/no overlap
barcode_groups_low = adata_dict[base_res].obs[adata_dict[base_res].obs[endo_of_interest].isin(low_overlap_groups)].index
adata_dict[base_res].obs.loc[barcode_groups_low, 'majority_overlap'] = endo_of_interest


# --- Get some colors for the pathology annotations

adata = adata_dict['square_008um']

annotation_labels = {
    'epithelial': ['invasive', 'in_situ', 'benign', 'normal'],
    'immune': ['lymphocytes'],
    'stroma': ['stroma_fibrous', 'stroma_adipose'],
    'additional_features': ['necrosis', 'blood_vessel', 'LVI']
}

all_colors = sc_palettes.default_20

# Remove red to be more colorblind-friendly
del all_colors[3] 

# Extract all pathology categories
all_cats = list(adata.obs[annotation_key].dropna().unique())

colors = sc_palettes.default_20[:len(annotation_labels)]

top_cat_colors = {}
for n, cat in enumerate(annotation_labels):
    top_cat_colors[cat] = colors[n]

# Get some color variations of the original colors (so variations of the same annotation are a similar color)
color_dict_patho = {}

# And also define a plot order
plot_order = []

for top_cat, labels in annotation_labels.items():
    base_color = top_cat_colors[top_cat]
    # Check how many of the labels are found in the object
    num_cols = sum(item in set(all_cats) for item in labels)
    # Get the matching items
    matching_items = [item for item in labels if item in set(all_cats)]
    plot_order += matching_items
    if num_cols > 4:
        step_size=0.1
    else:
        step_size=0.2
    variations = generate_color_variations(base_color, num_cols, step_size=step_size)
    for n, item in enumerate(matching_items):
        color_dict_patho[item] = variations[n]

# Set the plot_order
adata.obs[annotation_key] = pd.Categorical(adata.obs[annotation_key], categories=plot_order, ordered=True)


# --- Also set colors for endothelial types

colors_file = os.path.join(project_dir, 'resources', 'plotting_colors', 'endothelial_colors.csv')
colors_df = pd.read_csv(colors_file)

# Make a color_dict
color_dict = colors_df.set_index('endothelial_type')['color'].to_dict()

# Add a color for overlap
color_dict['overlap'] = 'indigo'


# --- Plot 1: Spatial overview with path annotations

adata = adata_dict[base_res]

# Plotting
sc.set_figure_params(fontsize=15)

ncols = 2
nrows = 1

fig, axs = plt.subplots(ncols=ncols, nrows=nrows, figsize=(7*ncols, 5*nrows))

# -- First the vessel type of interest
n=0
ax = axs[n]

adata_sub = adata[adata.obs[endo_of_interest].notna()]

sc.pl.spatial(adata_sub, color=endo_of_interest, bw=False, alpha_img=0.5, crop_coord=coords_dict[base_res], 
            title=f'{sample_id} ({subtype})\n{endo_of_interest}',
            palette=[color_dict[endo_of_interest]], ax=ax, legend_loc=None)

ax.set_title(ax.get_title(), fontweight='bold')

# -- Then plot overlap with HEVs
n+=1
ax = axs[n]

adata_sub = adata[adata.obs['majority_overlap'].notna()]

sc.pl.spatial(adata_sub, color='majority_overlap', bw=False, alpha_img=0.5, crop_coord=coords_dict[base_res], 
            title=f'{sample_id} ({subtype})\nHEV overlap (>{min_overlap} spots)',
            palette=color_dict, ax=ax,)

ax.set_title(ax.get_title(), fontweight='bold')

plt.suptitle(f'{sample_id} ({subtype}) - {endo_of_interest} overlap with HEVs (minimum overlap >{min_overlap} spots)', y=1.2, x=0.5, fontweight='bold')


# --- Plot 2: Barplot of overlap

sc.set_figure_params(fontsize=10)

# Disable the grid
plt.rcParams['axes.grid'] = False

# Prepare for plotting
overlap_ratio = proportion_df[proportion_df['group_type'] == endo_of_interest]
overlap_ratio = overlap_ratio.set_index('group')['pct_overlap']
group_classification = overlap_ratio > 0.5
summary_df = group_classification.value_counts().rename_axis('majority_overlap').reset_index(name='num_groups')
summary_df['majority_overlap'] = summary_df['majority_overlap'].replace({True: 'overlap', False: 'no_overlap'})

color_dict['no_overlap'] = 'lightgray'

# Plotting
nrows = 1
ncols = 1
fig, ax = plt.subplots(ncols=ncols, nrows=nrows, figsize=(1.5*ncols, 2.5))

sns.barplot(data=summary_df, x='majority_overlap', y='num_groups', palette=color_dict, hue='majority_overlap', ax=ax,)

plotting_functions.beautify_plot(ax)
plotting_functions.annotate_bars(ax, summary_df, x='majority_overlap', y='num_groups', plot_order=['no_overlap', 'overlap'], margin=1.1, annotation_size=9)

ax.set_xlabel('')
ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
ax.set_ylabel('Number of vessels')

ax.set_title(f'{endo_of_interest}\nHEV overlap', fontweight='bold')

plt.suptitle(f'{sample_id} ({subtype}) - {endo_of_interest} overlap with HEVs (minimum overlap >{min_overlap} spots)', y=1.2, x=0.5, fontweight='bold')


# --- Figure out pixel to um scalefactor

target_res = 'square_008um'
resolution_numbers = re.findall(r'\d+', target_res)[0]
resolution_numbers = int(resolution_numbers)

spot_diameter_um = resolution_numbers # Square size is whatever the resolution is 

# Get a scalefactor for converting pixel coords to um
spot_diameter_pixel = adata_dict[target_res].uns['spatial'][sample_id]['scalefactors']['spot_diameter_fullres']

scalefactor_um = spot_diameter_um/spot_diameter_pixel


# --- Calculate distances on the fly

# Function to calculate distances
def calculate_distance(base_res, target_res, group_key):
    # Extract shortened versions of resolution
    base_res_short = base_res.split('_')[1]
    target_res_short = target_res.split('_')[1]
    adata = adata_dict[base_res]
    coords = adata.obsm['spatial']  # Extract the spot coords
    groups = adata.obs[group_key]  # Group labels
    # Convert to DataFrame for easier manipulation
    df = pd.DataFrame(coords, columns=["x", "y"])
    df["vessel_group"] = groups.values
    # Add barcodes to the 16µm DataFrame
    df[f"barcode_{base_res_short}"] = adata.obs.index
    df = df.dropna(subset='vessel_group')
    # Extract coordinates and barcodes for each group
    grouped_data_16 = {
        group: df[df["vessel_group"] == group][["x", "y", f"barcode_{base_res_short}"]].reset_index(drop=True)
        for group in df["vessel_group"].unique()
    }
    # --- Extract decon_resolution spot coordinates
    adata_48um = adata_dict[target_res]
    df_48um = pd.DataFrame(adata_48um.obsm['spatial'], columns=["x", "y"])
    df_48um["barcode"] = adata_48um.obs.index
    # --- Compute the minimum distance to any 16 µm spot in each group
    results = []
    for group, group_df_16 in grouped_data_16.items():
        if group_df_16.empty:
            continue  # Skip empty groups
        coords_16um = group_df_16[["x", "y"]].to_numpy()
        barcodes_16um = group_df_16[f"barcode_{base_res_short}"].to_numpy()
        # Build KDTree for the current group
        tree = cKDTree(coords_16um)
        # Query the tree to find the closest 16 µm spot for each 48 µm spot
        min_distances, nearest_idx = tree.query(df_48um[['x', 'y']].values)
        # Get coordinates and barcodes of the closest 16 µm spot
        nearest_coords = coords_16um[nearest_idx]
        nearest_barcodes = barcodes_16um[nearest_idx]
        # Store results in a DataFrame, including barcodes
        df_temp = pd.DataFrame({
            f"barcode_{target_res_short}": df_48um["barcode"],  # Save barcode
            "vessel_group": group,
            f"x_{target_res_short}": df_48um["x"],
            f"y_{target_res_short}": df_48um["y"],
            f"nearest_{base_res_short}_x": nearest_coords[:, 0],
            f"nearest_{base_res_short}_y": nearest_coords[:, 1],
            f"nearest_{base_res_short}_barcode": nearest_barcodes,
            f"nearest_{base_res_short}_dist": min_distances
        })
        results.append(df_temp)
    # Combine results for all groups
    df_results = pd.concat(results, ignore_index=True)
    # --- Convert pixel distances to um
    df_results[f'nearest_{base_res_short}_dist_um'] = df_results[f'nearest_{base_res_short}_dist'] * scalefactor_um
    return df_results

# Save results in a dict
group_key = endo_of_interest
print(f'Calculating distances for {group_key}')
dist_df = calculate_distance(base_res, target_res, group_key)


# Also get a df with smallest distance to any vessel within a vessel type
base_res_short = base_res.split('_')[1]
target_res_short = target_res.split('_')[1]

# Now calculate smallest distance for overlap vs no overlap
# Get group names as sets
groups_with_overlap = set(group_classification[group_classification == True].index)
groups_without_overlap = set(group_classification[group_classification == False].index)

df_sub = dist_df[dist_df['vessel_group'].isin(groups_with_overlap)]
min_dist_overlap = df_sub.sort_values(f"nearest_{base_res_short}_dist_um").drop_duplicates(f"barcode_{target_res_short}", keep="first")

df_sub = dist_df[dist_df['vessel_group'].isin(groups_without_overlap)]
min_dist_no_overlap = df_sub.sort_values(f"nearest_{base_res_short}_dist_um").drop_duplicates(f"barcode_{target_res_short}", keep="first")

# Add to adata
adata_dict[target_res].obs[f'min_dist_overlap'] = min_dist_overlap.set_index(f'barcode_{target_res_short}')['nearest_016um_dist_um']
adata_dict[target_res].obs[f'min_dist_no_overlap'] = min_dist_no_overlap.set_index(f'barcode_{target_res_short}')['nearest_016um_dist_um']


# Make some plots to check it worked
# --- Plot 3: Distance calculations per group

sc.set_figure_params(fontsize=15)

ncols = 3
nrows = 1

fig, axs = plt.subplots(ncols=ncols, nrows=nrows, figsize=(8*ncols, 6*nrows),)

adata_sub = adata_dict[base_res][adata_dict[base_res].obs['majority_overlap'].notna()]
adata_sub.obs['majority_overlap'].unique()

ax = axs[0]
sc.pl.spatial(adata_sub, color='majority_overlap', bw=False, alpha_img=0.5, crop_coord=coords_dict[base_res], 
            title=f'{sample_id} ({subtype})\nHEV overlap (>{min_overlap} spots)',
            palette=color_dict, ax=ax,)

ax.set_title(ax.get_title(), fontweight='bold')

for n, overlap in enumerate(['no_overlap', 'overlap']):
    n+=1
    plot_val = overlap
    if overlap == 'no_overlap':
        plot_val = endo_of_interest
    # Set ax
    ax = axs[n]
    # Plot all HEV groups
    sq.pl.spatial_scatter(adata_dict[target_res], color=f'min_dist_{overlap}', img_alpha=0.7,
                    crop_coord=coords_dict_sq[target_res], title=f'{sample_id} ({subtype})\nHEV classification', ax=ax,
                    shape='square', vmax=400)
    sq.pl.spatial_scatter(adata_sub, color='majority_overlap', img_alpha=0, groups=[plot_val],
                    crop_coord=coords_dict_sq[target_res], title=f'{sample_id} ({subtype})\n{plot_val}\nDistance (µm)', ax=ax,
                    shape='square', palette=ListedColormap(['white']), legend_loc=None)
    ax.set_title(ax.get_title(), fontweight='bold')


plt.suptitle(f'{sample_id} ({subtype}) - vessel groups - distance (um) to each {target_res} spot (threshold={threshold})', y=1.15, fontweight='bold', x=0.5)


# --- Now bin the distances

group_types = ['no_overlap', 'overlap']

distance_dict['no_overlap'] = min_dist_no_overlap
distance_dict['overlap'] = min_dist_overlap

def bin_distances(df, column_name, bin_size, max_distance):
    # Define the bin edges
    bins = np.arange(0, max_distance + bin_size, bin_size)
    # Use pd.cut to assign each value in the column to a bin
    df['binned'] = pd.cut(df[column_name], bins=bins, right=False, labels=[f'{bins[i]}-{bins[i+1]}' for i in range(len(bins)-1)], include_lowest=True)
    return df['binned']

# Bin distances
bin_size = 25
column_name = 'nearest_016um_dist_um'
max_distance=400

for cell_type in group_types:
    df_min_dist = distance_dict[cell_type]
    df_min_dist['binned'] = bin_distances(df_min_dist, column_name=column_name, bin_size=bin_size, max_distance=max_distance)
    # Add to adata
    adata_dict[target_res].obs[f'dist_bin_{cell_type}'] = df_min_dist.set_index(f'barcode_{target_res_short}')['binned']


# --- Plot 4: Plot distances after binning

adata_dict[base_res].obs['overlap_adj'] = adata_dict[base_res].obs['majority_overlap']
adata_dict[base_res].obs['overlap_adj'] = adata_dict[base_res].obs['overlap_adj'].replace(endo_of_interest, 'no_overlap')

ncols = len(group_types)
nrows = 1

plt.figure(figsize=(8*ncols, 6*nrows))

for n, cell_type in enumerate(group_types):
    # add a new subplot iteratively
    ax = plt.subplot(nrows, ncols, n + 1)
    # Add legend only to the last one
    if n == (ncols-1):
        legend_loc = 'right margin'
    else:
        legend_loc = None
    adata_sub = adata_dict['square_008um'][adata_dict['square_008um'].obs[f'dist_bin_{cell_type}'].notna()]
    sq.pl.spatial_scatter(adata_sub, color=f'dist_bin_{cell_type}', img_alpha=0.7,
                    crop_coord=coords_dict_sq[target_res], title=f'{sample_id} ({subtype})\nHEV classification', ax=ax,
                    shape='square', palette='viridis', legend_loc=legend_loc)
    # sc.pl.spatial(adata_sub, color=f'dist_bin_{cell_type}', bw=False, alpha_img=0.7, crop_coord=coords_dict['square_008um'], 
    #                 title=f'{sample_id} ({subtype})\n{cell_type}',
    #                 palette='viridis', ax=ax, legend_loc=legend_loc)
    adata_sub = adata_dict[base_res][adata_dict[base_res].obs['overlap_adj'].notna()]
    sq.pl.spatial_scatter(adata_sub, color='overlap_adj', img_alpha=0, groups=[cell_type],
                    crop_coord=coords_dict_sq[base_res], title=f'{sample_id} ({subtype})\n{endo_of_interest} {cell_type}\nDistance (µm)', ax=ax,
                    shape='square', palette=ListedColormap(['white']), legend_loc=None)
    ax.set_title(ax.get_title(), fontweight='bold')

plt.suptitle(f'{sample_id} ({subtype}) - Distances after binning (resolution=square_008um)', y=1.1, fontweight='bold', x=0.5)


# --- Prepare for plotting composition across distance

composition_dict = {}

for cell_type in group_types:
    # Extract a df
    cell_df = adata_dict['square_008um'].obs[[f'dist_bin_{cell_type}'] + [annotation_key]]
    # Drop all spots that are too far away
    cell_df = cell_df.dropna(subset=f'dist_bin_{cell_type}')
    # Count annotations per bin
    counts = cell_df.groupby([f'dist_bin_{cell_type}', annotation_key]).size().unstack(fill_value=0)
    # Convert to proportions
    proportions = counts.div(counts.sum(axis=1), axis=0)
    composition_dict[cell_type] = proportions


# --- Plot 5: Cell type composition across distance

# Disable the grid
plt.rcParams['axes.grid'] = False

nrows = 1
ncols = len(group_types)

fig, axs = plt.subplots(ncols=ncols, nrows=nrows, figsize=(5*ncols, 4*nrows), sharey=True)

for n, cell_type in enumerate(group_types):
    ax = axs[n]
    composition = composition_dict[cell_type][plot_order[::-1]]
    composition.plot.area(ax=ax, color=color_dict_patho, alpha=0.8, linewidth=0)
    if n == (ncols-1):
        # Move legend
        legend = ax.legend()
        legend.set_title('')
        # Automatically calculate how many columns
        ncol_legend = math.ceil(len(plot_order)/10)
        ax.legend(handlelength=0.8, loc='center left', bbox_to_anchor=(1, 0.5), alignment='left', frameon=False, markerscale=2, reverse=True, ncols=ncol_legend)
    else:
        # Remove legend otherwise
        ax.get_legend().remove()
    # Make sure all xticklabels are displayed
    ax.set_xticklabels(ax.get_xticklabels(), rotation=90)  # Rotate x-tick labels if needed
    ax.set_xticks(range(len(composition)))  # Ensure all x-ticks are placed
    ax.set_xticklabels(composition.index)  # Set the x-tick labels to match the DataFrame index
    # Adjust lims
    ax.set_ylim(0, 1)
    ax.set_xlim(0, ax.get_xlim()[1]-0.75)
    # Set xlabel
    ax.set_xlabel('Distance (µm)')
    ax.set_ylabel('Proportion')
    # Set title
    ax.set_title(f'{endo_of_interest}\n{cell_type} HEV', fontweight='bold')

# Adjust positioning & size
for n, ax in enumerate(axs):
    ax_pos = ax.get_position()
    # Shift the subplot to the desired location + set the width so it's equal for all plots
    if n == 0:
        ax.set_position([ax_pos.x0, ax_pos.y0, ax_pos.width, ax_pos.height])
    else:
        new_pos = ax_pos_prev.x1 + 0.02
        ax.set_position([new_pos, ax_pos.y0, ax_pos.width, ax_pos.height])
    ax_pos_prev = ax.get_position()

plt.suptitle(f'{sample_id} ({subtype}) - Pathology annotations across vessel distance bins stratified by overlap HEV (min_pct_spots >{min_overlap}) (threshold={threshold}, bin_size={bin_size}, max_distance={max_distance})', y=1.15, fontweight='bold', x=0.5)


# --- Save output

qc_dir = os.path.join(project_dir, 'results', '04_vascular_microenvironment', '03_compare_HEV_to_other_venous', '01_compare_microenvironment')

today = date.today()
today = today.strftime("%Y%m%d")

out_dir = os.path.join(qc_dir)

os.makedirs(out_dir, exist_ok=True)

# Save plots
out_dir_plots = os.path.join(out_dir, 'plots', today)
os.makedirs(out_dir_plots, exist_ok=True)

# Also save the groupings as a .csv
out_dir_csv = os.path.join(out_dir, 'csv_files', today)
os.makedirs(out_dir_csv, exist_ok=True)

# Save plots
filename = f'{sample_id}_{out_name}_composition_vs_distance.pdf'

save_multi_image(os.path.join(out_dir_plots, filename))
plt.close('all')

# # And save the distance bins
# filename = f'{sample_id}_vessel_bins.csv'
# df_min_dist.to_csv(os.path.join(out_dir_csv, filename), index=False)

# And save the composition results
for cell_type in group_types:
    composition = composition_dict[cell_type]
    filename = f'{sample_id}_{out_name}_{cell_type.replace(" ", "_")}_composition_vs_distance.csv'
    composition.to_csv(os.path.join(out_dir_csv, filename), index=True)


# --- Finish

print("Script successfully completed")
