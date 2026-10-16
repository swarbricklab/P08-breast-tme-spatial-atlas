#!/usr/bin/env python3

# --- Aim
# Plot PVL composition across distance bins to look at tendencies of PVLs association with HEVs vs non-HEV vessels


# --- Define environment
# env: vis-hd

# --- Load packages
import argparse
import os

import anndata as ad
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import squidpy as sq
import seaborn as sns
from matplotlib.colors import ListedColormap
import sys
import re
from scipy.spatial import cKDTree

sc.set_figure_params(facecolor="white", figsize=(8, 8))
sc.settings.verbosity = 3

mpl.rcParams['pdf.fonttype'] = 42


# --- Argparse arguments


parser = argparse.ArgumentParser(description="PVL composition across distance from HEV vs non-HEV venous vessels",
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

# Import function script as a module
sys.path.append(os.path.join(project_dir, 'scripts', '04_vascular_microenvironment'))
import plotting_functions


# --- Extract some info from the sample sheet
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
# Need both the 48 & 16 res
base_res = 'square_016um'
decon_resolution = 'square_016um'
resolutions = [base_res]

if decon_resolution not in resolutions:
    resolutions.append(decon_resolution)

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

# Decide which decon level to use
decon_level = 'subset'

# Also get PVL classes
PVL_classes = list(class_df[class_df['miniatlas_major'] == 'PVL'][f'miniatlas_{decon_level}'].unique())
PVL_classes = sorted(PVL_classes)


# --- Load vessel groups

threshold = 0.5
threshold_adj = str(threshold).split('.')[-1]

# Define which endothelial type want to look at
endo_of_interest = 'Endothelial ACKR1'
out_name = endo_of_interest.split(' ')[-1]

# Specify path
TLS_group_dir = os.path.join(res_dir, '04_vascular_microenvironment', '02_identify_all_vasculature', '01_classify_vessels', 'csv_files')

most_recent_dir_path = TLS_group_dir

filename = f'{sample_id}_vessel_classification.csv'
vessel_classes = pd.read_csv(os.path.join(most_recent_dir_path, filename), index_col=0)
adata_dict[base_res].obs[vessel_classes.columns] = vessel_classes

# Collapse to single vessel class per spot
collapsed_col = vessel_classes[endothelial].apply(lambda row: row.first_valid_index() if row.first_valid_index() is not None else None, axis=1)
adata_dict[base_res].obs['vessel_class'] = collapsed_col


# --- Load the deconvolution results

resolution_numbers = re.findall(r'\d+', decon_resolution)[0]
resolution_numbers = int(resolution_numbers)

decon_dir = os.path.join(res_dir, '02_cell_type_mapping', 'deconvolution', 'visium_hd', f'{decon_level}_level', decon_resolution,)

decon_file = f'{sample_id}_{decon_level}_{resolution}_q05_cell_abundance_w_sf.csv'

decon_res = pd.read_csv(os.path.join(decon_dir, decon_file), index_col=0)

# Only load the PVL classes
adata_dict[decon_resolution].obs[PVL_classes] = decon_res[PVL_classes]


# --- Load HEV overlap information

# Define which endothelial type want to look at
endo_of_interest = 'Endothelial ACKR1'
out_name = endo_of_interest.split(' ')[-1]

TLS_group_dir = os.path.join(res_dir, '04_vascular_microenvironment', '02_identify_all_vasculature', '02_investigate_vessels', 'csv_files')

most_recent_dir_path = TLS_group_dir

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


# --- Also get some colors

# Define some colors per endothelial type
colors_file = os.path.join(project_dir, 'resources', 'plotting_colors', 'endothelial_colors.csv')
colors_df = pd.read_csv(colors_file)

# Make a color_dict
color_dict = colors_df.set_index('endothelial_type')['color'].to_dict()

# Add a color for overlap
color_dict['overlap'] = 'indigo'


# --- Plot 1: Plot endo of i, HEV overlap, distances after binning

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

target_res = decon_resolution
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


# --- Now bin the distances

group_types = ['no_overlap', 'overlap']

distance_dict = {}

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


# --- Plot 3: Plot distances after binning

sc.set_figure_params(fontsize=15)

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
    adata_sub = adata_dict[target_res][adata_dict[target_res].obs[f'dist_bin_{cell_type}'].notna()]
    sq.pl.spatial_scatter(adata_sub, color=f'dist_bin_{cell_type}', img_alpha=0.7,
                    crop_coord=coords_dict_sq[target_res], title=f'{sample_id} ({subtype})\nHEV classification', ax=ax,
                    shape='square', palette='viridis', legend_loc=legend_loc)
    adata_sub = adata_dict[base_res][adata_dict[base_res].obs['overlap_adj'].notna()]
    sq.pl.spatial_scatter(adata_sub, color='overlap_adj', img_alpha=0, groups=[cell_type],
                    crop_coord=coords_dict_sq[base_res], title=f'{sample_id} ({subtype})\n{endo_of_interest} {cell_type}\nDistance (µm)', ax=ax,
                    shape='square', palette=ListedColormap(['white']), legend_loc=None)
    ax.set_title(ax.get_title(), fontweight='bold')

plt.suptitle(f'{sample_id} ({subtype}) - Distances after binning (resolution={target_res})', y=1.1, fontweight='bold', x=0.5)


# --- Plot 4: Plot PVL association between HEV overlap vs no overlap

sc.set_figure_params(fontsize=10)
# Set grid behind plot
plt.rcParams['axes.axisbelow'] = True

# Update color_dict
color_dict['no_overlap'] = color_dict[endo_of_interest]

# Also limit the bins to make the plots smaller
bin_limit = 200

nrows = len(group_types)
ncols = len(PVL_classes)

fig, axs = plt.subplots(ncols=ncols, nrows=nrows, figsize=(2, 2.2), squeeze=False)

for n, PVL_class in enumerate(PVL_classes):
    for i, cell_type in enumerate(group_types):
        # add a new subplot iteratively
        ax = axs[i, n]
        cell_df = adata_dict[target_res].obs[[f'dist_bin_{cell_type}'] + PVL_classes]
        # Melt to long format: one row per (bin, subtype, PVL value)
        temp = cell_df[[f'dist_bin_{cell_type}', PVL_class]].copy()
        temp = temp.rename(columns={f'dist_bin_{cell_type}': 'distance_bin', PVL_class: 'PVL_value'})
        temp['endothelial_subtype'] = cell_type
        df_long = temp
        # Drop rows with missing bins or PVL values
        df_long = df_long.dropna(subset=['distance_bin', 'PVL_value'])
        # Order bins numerically
        bin_order = sorted(df_long['distance_bin'].unique(), key=lambda x: int(x.split('-')[0]))
        df_long['distance_bin'] = pd.Categorical(df_long['distance_bin'], categories=bin_order, ordered=True)
        # Limit the bins
        bin_order = [b for b in bin_order if int(b.split('-')[1]) <= bin_limit]
        df_long = df_long[df_long['distance_bin'].isin(bin_order)]
        # Plot each subtype separately
        ax = axs[i, n]  # Access the correct subplot
        sub_df = df_long[df_long['endothelial_subtype'] == cell_type]
        # Calculate mean per bin
        sub_df = sub_df.groupby('distance_bin')['PVL_value'].mean().reset_index()
        sub_df.set_index('distance_bin').plot.area(ax=ax, color=color_dict[cell_type], legend=False, alpha=0.5)
        # Make sure all xticklabels are displayed
        ax.set_xticklabels(ax.get_xticklabels(), rotation=90)  # Rotate x-tick labels if needed
        ax.set_xticks(range(len(bin_order)))  # Ensure all x-ticks are placed
        ax.set_xticklabels(bin_order)  # Set the x-tick labels to match the DataFrame index
        ax.set_title(f'{PVL_class}', fontweight='bold', fontsize=10)
        # # Rotate xtick labels 90
        # ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
        # Adjust labels
        ax.set_xlabel('Distance (µm)')
        ax.set_ylabel('')
        # Adjust x-limits to fit the bins
        ax.set_xlim(0, len(bin_order) - 1)

ymax = 0

for ax in axs.ravel():
    if ax.get_ylim()[1] > ymax:
        ymax = ax.get_ylim()[1]

# Adjust positioning & size
for n, ax in enumerate(axs.ravel()):
    ax_pos = ax.get_position()
    ax.set_position([ax_pos.x0, ax_pos.y0, len(bin_order)*0.08, ax_pos.height])
    ax.set_ylim(0, ymax)  # Set the same y-limits for all subplots

plotting_functions.move_plots_x_y(axs, x_offset=0.3, y_offset=0.05, remove_ylabels=False)

# Move legend
ax = axs[0,-1]
plotting_functions.add_legend(ax=ax, ordered_list=group_types, color_lut=color_dict, title=f'{endo_of_interest}\nHEV overlap', 
                               loc='center left', y_offset=0, x_offset=1)

# Make the facecolor of the legend handles more transparent
for handle in ax.get_legend().legend_handles:
    facecolor = handle.get_facecolor()
    handle.set_facecolor(handle.get_facecolor()[:3] + (0.5,))  # Add transparency to the facecolor
    # Set the edgecolor to be the same as the facecolor
    handle.set_edgecolor(facecolor)
    # Set the linewidth to 2
    handle.set_linewidth(1.25)

# Get the left-most x-position among all subplots
lefts = [ax.get_position().x0 for ax in axs.flat]
left_min = min(lefts)

# Add shared y-label just to the left of all subplots
fig.text(left_min - 0.3, 0.5, 'Average proportion', va='center', ha='right',
         rotation='vertical')

# Calculate x to center
title_x = (axs[0][0].get_position().x0 + axs[0][-1].get_position().x1)/2

plt.suptitle(f'{sample_id} ({subtype}) - PVL composition across distance bins (threshold={threshold}, decon_res={decon_resolution}, decon_level={decon_level}, (minimum overlap >{min_overlap} spots))',
             y=1.2, fontweight='bold', x=title_x)


# --- Save output

qc_dir = res_dir

out_dir = os.path.join(qc_dir, '04_vascular_microenvironment', '03_compare_HEV_to_other_venous', '03_compare_PVL_association',)

os.makedirs(out_dir, exist_ok=True)

# Save plots
out_dir_plots = os.path.join(out_dir, 'plots')
os.makedirs(out_dir_plots, exist_ok=True)

# Save plots
filename = f'{sample_id}_HEV_vs_non_HEV_PVL_composition_vs_distance.pdf'

save_multi_image(os.path.join(out_dir_plots, filename))
plt.close('all')

# Also save csv file for convenience
out_dir_csv = os.path.join(out_dir, 'csv_files')
os.makedirs(out_dir_csv, exist_ok=True)

group_keys = [f'dist_bin_{cell_type}' for cell_type in group_types]
cell_df = adata_dict[decon_resolution].obs[group_keys + PVL_classes]

# Drop rows with NaN values in the group keys
cell_df = cell_df.dropna(subset=group_keys)

# Save as gzipped csv
cell_df.to_csv(os.path.join(out_dir_csv, f'{sample_id}_PVL_composition_vs_distance.csv'), index=True)


# --- Finish

print("Script successfully completed")
