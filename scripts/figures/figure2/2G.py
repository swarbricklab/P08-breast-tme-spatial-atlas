#!/usr/bin/env python3

# --- Aim
# Plot selected GC examples for Figure 2G
# Shows GCs detected by signature
# This creates Figure 2G

# --- Define environment
# env: spatial-pilot

# --- Load packages
import argparse
import os
from datetime import datetime
import sys

import scanpy as sc
from matplotlib_scalebar.scalebar import ScaleBar
import anndata as ad
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import tifffile

sc.set_figure_params(facecolor="white", figsize=(8, 8))
sc.settings.verbosity = 3


# --- Argparse arguments
parser = argparse.ArgumentParser(description="Plot GC examples for Figure 2G",
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

print(f'Samples in python: {sample_list}')

sys.path.append(os.path.join(project_dir, 'scripts'))
from helper_functions import auto_crop_scanpy


# --- Load sample sheet

if args.get("sample_sheet"):
    sample_sheet_path = args["sample_sheet"]
else:
    sample_sheet_path = os.path.join(project_dir, 'config', 'sample_sheet_visium.csv')

sample_sheet = pd.read_csv(sample_sheet_path)
sample_lut = sample_sheet.set_index('section_id')


# --- Also get some tissue metadata

resource_dir = os.path.join(project_dir, 'resources')
metadata_file = 'tissue_metadata.csv'

tissue_metadata = pd.read_csv(os.path.join(resource_dir, metadata_file))

# Merge the two DataFrames on the common column 'tissue_id'
merged_df = pd.merge(sample_sheet, tissue_metadata, on='tissue_id')

# Create a sample - subtype dictionary
subtype_dict = dict(zip(merged_df['section_id'], merged_df['clinical_subtype']))


# --- Function to find the most recent results
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
        most_recent_dir_path = None
    return most_recent_dir_path


# --- Define thresholds and version
GC_threshold = 0.7
GC_threshold_adj = str(GC_threshold).split('.')[-1]


# --- Load anndata objects 
adata_dict = {}

for sample in sample_list:
    adata_dict[sample] = ad.read_h5ad(os.path.join(h5ad_dir, 'visium', f'{sample}.h5ad'))

# Add some info for plotting
for sample in sample_list:
    adata_dict[sample].obs['sample_id'] = sample
    adata_dict[sample].obs['subtype'] = subtype_dict[sample]


# --- Calculate crop coords for optimal presentation

coords = {}
for sample in adata_dict.keys():
    coords[sample] = auto_crop_scanpy(adata_dict[sample], border_size=0.1)


# --- Load the GC group annotations

GC_group_dir = os.path.join(res_dir, '03_TLS_niche_discovery', '04_map_niche_signatures', '02_group_high_score_areas', 'GC_score')

most_recent_dir_path = find_most_recent(GC_group_dir)

# Read in the results & add to anndata
group_key = 'group_too_small_removed'

for sample in sample_list:
    GC_file_name = f'{sample}_GC_groups.csv'
    GC_file_path = os.path.join(most_recent_dir_path, 'csv_files', GC_file_name)
    if os.path.exists(GC_file_path):
        GC_results = pd.read_csv(GC_file_path, index_col=0)
        adata_dict[sample].obs[group_key] = GC_results[group_key]
    else:
        print(f'Warning: {GC_file_path} not found')


# --- Also load the GC score results

GC_res_dir = os.path.join(res_dir, '03_TLS_niche_discovery', '04_map_niche_signatures',
                          '01_score_signatures', 'GC_score')

most_recent_dir_path = find_most_recent(GC_res_dir)

num_genes = 20
GC_key = f'GC_score_{num_genes}'

# Read in the results & add to anndata
for sample in sample_list:
    GC_file_name = f'{sample}_GC_score.csv'
    GC_file_path = os.path.join(most_recent_dir_path, 'csv_files', GC_file_name)
    if os.path.exists(GC_file_path):
        GC_results = pd.read_csv(GC_file_path, index_col=0)
        adata_dict[sample].obs[GC_key] = GC_results[GC_key]
    else:
        print(f'Warning: {GC_file_path} not found')


# --- Predefine GCs to plot

# These are the selected examples for the figure
selected_non_overlapping_per_sample = {
    'P22_Vis': ['group5', 'group6'],
    'P24_Vis': ['group2', 'group4'],
    'P18_Vis': ['group1'],
    'P29_Vis_Rep1': ['group1'],
}

selected_samples = list(selected_non_overlapping_per_sample.keys())

# Filter sample_list to selected
sample_list = [s for s in sample_list if s in selected_samples]

print(f'Selected samples for plotting: {sample_list}')


# --- Function to move plots

def move_plots_x_y(axs, x_offset_small=0.02, x_offset_large=0.10, y_offset=0.04, remove_ylabels=True, remove_xlabels=True, remove_title=True):
    for i, axs_row in enumerate(axs):
        for j, ax in enumerate(axs_row):
            if j > 0:
                if remove_ylabels:
                    ax.set_ylabel('')
                    ax.set_yticklabels('')
            if i != (len(axs) - 1):
                if remove_xlabels:
                    ax.set_xlabel('')
                    ax.set_xticklabels('')
            if i != 0:
                if remove_title:
                    ax.set_title('')
            ax_pos = ax.get_position()
            if i == 0:
                new_pos_y = ax_pos.y0
            else:
                new_pos_y = ax_pos_row_prev.y0 - ax_pos.height - y_offset
            if j == 0:
                new_pos_x = ax_pos.x0
            else:
                # Use large offset after every 2nd column
                if j % 2 == 0:
                    offset = x_offset_large
                else:
                    offset = x_offset_small
                new_pos_x = ax_pos_prev.x0 + ax_pos_prev.width + offset
            ax.set_position([new_pos_x, new_pos_y, ax_pos.width, ax_pos.height])
            ax_pos_prev = ax.get_position()
        ax_pos_row_prev = axs_row[0].get_position()


# --- Prepare for plotting

# Calculate total GCs to plot
total_non_overlapping_GCs = sum([len(v) for v in selected_non_overlapping_per_sample.values()])

# Define ncols & nrows
ncols = 3 * 2  # Will make 2 plots per GC: H&E & GC score
nrows = int(np.ceil(total_non_overlapping_GCs*2 / ncols))

# Define a vmax
vmax = 0.8

# Define a value for the scalebar
scalebar_val = 100  # in microns
scalebar_units = 'µm'

vis_spot_diameter = 55
microns_per_pixel = 0.4426  # Garvan imaging


# Function to crop the image with safe coordinates
def crop_image(img, coords):
    # Clip the coordinates to ensure they are within the image boundaries
    x_min = max(coords[0], 0)
    x_max = min(coords[1], img.shape[1])
    y_min = max(coords[2], 0)
    y_max = min(coords[3], img.shape[0])
    # Crop the image using the clipped coordinates
    crop_img = img[y_min:y_max, x_min:x_max]
    return crop_img, (x_min, x_max, y_min, y_max)


# --- Load the full resolution images and plot

fig = plt.figure(figsize=(5, 5))

plot_index = 0
for i, sample in enumerate(selected_samples):
    if sample not in sample_list:
        continue
    print(f'Plotting sample: {sample}')
    adata = adata_dict[sample]
    spot_diameter_fulres = adata.uns['spatial'][sample]['scalefactors']['spot_diameter_fullres']
    scale_factor = adata.uns['spatial'][sample]['scalefactors']['tissue_hires_scalef']
    
    # Find and load the fullres image
    image_dir = os.path.join(project_dir, 'data', 'images', 'he', 'aligned_images', 'visium')
    image_file = f'{sample}.ome.tif'
    image_file_path = os.path.join(image_dir, image_file)
    
    if not os.path.exists(image_file_path):
        print(f'No fullres image found for sample {sample}, skipping...')
        continue
    
    print(f'Fullres image: {image_file_path}')
    
    # Load the original high-res scan using tifffile
    img = tifffile.imread(image_file_path)
    print('Fullres H&E image successfully loaded')
    
    # Define a key marking where the fullres will be added in adata
    original_hires_key = 'fulres'
    
    # -- Calculate the scale factor
    # Divide the height of the original hires image by the height of the downsampled hires image
    size_fact = img.shape[0]/adata.uns['spatial'][sample]['images']['hires'].shape[0]
    
    # Extract the downsampled hires scale factor
    hires_fact = adata.uns['spatial'][sample]['scalefactors']['tissue_hires_scalef']
    
    # Calculate the scale factor for the original hires image and add to the anndata object
    adata.uns['spatial'][sample]['scalefactors'][f'tissue_{original_hires_key}_scalef'] = size_fact * hires_fact
    
    for j, gc_group in enumerate(selected_non_overlapping_per_sample[sample]):
        print(f'Plotting GC group: {gc_group}')
        
        # Subset the adata to the GC group to calculate crop coords
        adata_group_sub = adata[adata.obs[group_key] == gc_group]
        
        # Make sure border_size is equal despite differences in area
        desired_border = 8 * (vis_spot_diameter/microns_per_pixel)
        
        # Get original coords to calculate width
        coords_group = auto_crop_scanpy(adata_group_sub, border_size=0, spot_diameter=vis_spot_diameter/microns_per_pixel)
        width = coords_group[1] - coords_group[0]
        border_size = desired_border/width
        
        if len(adata_group_sub) == 1:
            border_size *= 0.3
        
        # Calculate new crop coord based on the group
        coords_group = auto_crop_scanpy(adata_group_sub, border_size=border_size)
        
        # -- Prepare an anndata that contains only the cropped area we want to plot
        # Crop the H&E image
        crop_img, adjusted_coords = crop_image(img, coords_group)
        
        # Calculate the height and width of the cropped image
        height_crop = adjusted_coords[3] - adjusted_coords[2]
        width_crop = adjusted_coords[1] - adjusted_coords[0]
        
        # Also check if either x0 or y0 are below 0
        xmin = 0
        if adjusted_coords[0] < xmin:
            xmin = adjusted_coords[0]
        ymin = 0
        if adjusted_coords[2] < ymin:
            ymin = adjusted_coords[2]
        
        # Make a copy of anndata & add the crop in the hires slot
        adata_crop = adata.copy()
        adata_crop.uns['spatial'][sample]['images'][original_hires_key] = crop_img
        
        # Shift the spot coordinates
        spot_coords = adata_crop.obsm['spatial']
        shifted_coords = spot_coords - np.array([adjusted_coords[0], adjusted_coords[2]])
        adata_crop.obsm['spatial'] = shifted_coords
        
        # Plot 1: H&E with GC group highlighted
        ax1 = plt.subplot(nrows, ncols, plot_index + 1)
        sc.pl.spatial(
            adata_crop,
            img_key=original_hires_key,
            color=None,
            ax=ax1,
            title='',
            crop_coord=(xmin, xmin + width_crop, ymin, ymin + height_crop),
        )
        
        # Plot 2: GC score
        ax2 = plt.subplot(nrows, ncols, plot_index + 2)
        sc.pl.spatial(
            adata_crop,
            img_key=original_hires_key,
            color=GC_key,
            alpha=0.8,
            ax=ax2,
            show=False,
            title='',
            vmin=0,
            vmax=vmax,
            cmap='viridis',
            colorbar_loc=None,
            crop_coord=(xmin, xmin + width_crop, ymin, ymin + height_crop),
            spot_size=vis_spot_diameter/microns_per_pixel,
        )
        
        # Add a scalebar to the first plot
        scalebar = ScaleBar(microns_per_pixel, units=scalebar_units, location='lower right', box_alpha=0.5, fixed_value=scalebar_val,
                frameon=False, scale_loc="none", pad=0.3, height_fraction=0.04)
        ax1.add_artist(scalebar)
        plot_index += 2

# Extract the axs
axs = fig.get_axes()

# Adjust size
for ax in axs:
    pos = ax.get_position()
    new_pos = [pos.x0, pos.y0, 0.4, 0.4]
    ax.set_position(new_pos)
    # Also remove x & ylabels
    ax.set_xlabel('')
    ax.set_ylabel('')

# Move plots so they don't overlap
# Pad the axs list with None to make it evenly divisible by ncols
padded_axs = axs + [None] * (nrows * ncols - len(axs))

# Reshape the padded list
axs_grid = np.array(padded_axs).reshape(nrows, ncols)
axs_grid_cleaned = np.array([[ax for ax in row if ax is not None] for row in axs_grid if any(ax is not None for ax in row)])

move_plots_x_y(axs_grid_cleaned, x_offset_small=0.01, x_offset_large=0.03, y_offset=0.03, remove_title=False)

# Add a colorbar to the last plot on the first row
target_ax = axs_grid_cleaned[0][-1]
cbar_ax = fig.add_axes([target_ax.get_position().x1 + 0.01, target_ax.get_position().y0, 0.04, target_ax.get_position().height*0.7])
norm = plt.Normalize(vmin=0, vmax=vmax)
sm = plt.cm.ScalarMappable(cmap='viridis', norm=norm)
sm.set_array([])
cbar = fig.colorbar(sm, cax=cbar_ax)
# Set title to the left of the ax
cbar_ax.set_title('GC\nsignature', loc='left', pad=10)

# Set the yticklabels as min & max
cbar.set_ticks([0, vmax])
cbar.set_ticklabels(['Low', 'High'])

# Add a suptitle
title_x = axs_grid_cleaned[0][0].get_position().x0 + (axs_grid_cleaned[0][-1].get_position().x1 - axs_grid_cleaned[0][0].get_position().x0)/2
title_y = axs_grid_cleaned[0][0].get_position().y1 + 0.2
plt.suptitle('GCs detected by GC signature not overlapping with TLS annotations', x=title_x, y=title_y, fontweight='bold')


# --- Save the output

out_dir = os.path.join(res_dir, 'figures')
os.makedirs(out_dir, exist_ok=True)

filename = 'Fig_2G.pdf'
# Save with high resolution
plt.savefig(os.path.join(out_dir, filename), bbox_inches='tight', dpi=300)

print(f'Figure saved to {os.path.join(out_dir, filename)}')


# --- Finish

print("Script successfully completed")
