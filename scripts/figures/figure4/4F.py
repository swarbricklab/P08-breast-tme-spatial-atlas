#!/usr/bin/env python3

# --- Aim
# Plot Figure 4F: zoomed view showing TLS-HEVs vs stromal HEV-like vessels
# Shows spatial context for HEV classification used in ligand-receptor analysis
# This script is set up for sample P29_HD

# --- Define environment
# env: shapely

# --- Load packages

import argparse
from datetime import datetime
import os

import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np

import tifffile
import squidpy as sq

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgb
from matplotlib.colors import ListedColormap
import matplotlib.patches as mpatches
from matplotlib_scalebar.scalebar import ScaleBar

import geopandas as gpd

mpl.rcParams['pdf.fonttype'] = 42  # To get editable text in Illustrator

# --- Argparse arguments

parser = argparse.ArgumentParser(description="Plot HEV grouping overview for figure 4F",
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
h5ad_dir = args["h5ad_dir"]
sample_id = args["sample_id"]
project_dir = args["project_dir"]

# Override sample_id
sample_id = 'P29_HD'

# --- Extract some info from the sample sheet

# Extract the external_id
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
subtype = tissue_metadata[tissue_metadata['tissue_id'].astype(str) == str(tissue_id)]['clinical_subtype'].values[0]


# --- Load the LR adata object

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

bin_size = 25  # in um

lr_res_dir = os.path.join(res_dir, '04_vascular_microenvironment', '04_ligand_receptor_analysis',
                          '01_run_LR_analysis_stlearn', f'{bin_size}um_bin', 'csv_files')

most_recent_dir_path = find_most_recent(lr_res_dir)

adata_name = f'{sample_id}_HEV_LR_analysis_stlearn.h5ad'
adata = ad.read_h5ad(os.path.join(most_recent_dir_path, adata_name))


# --- Calculate crop coords for optimal presentation

def auto_crop_scanpy(adata, border_size=0.05, spot_diameter=350):
    # get the coordinates of the spots
    coord = adata.obsm['spatial']
    coord_x = coord[:, 0]
    coord_y = coord[:, 1]
    # do some calculations
    x_min = np.min(coord_x)
    x_max = np.max(coord_x)
    y_min = np.min(coord_y)
    y_max = np.max(coord_y)
    x_range = x_max - x_min
    y_range = y_max - y_min
    # define a mid point
    mid_point = np.array([x_range / 2 + x_min, y_range / 2 + y_min])
    mid_point = mid_point.astype(int)
    if x_range == 0 and y_range == 0:
        # Only one spot, define mid_dist and border based on a fixed size
        mid_dist = spot_diameter
        border = int((border_size / 2) * mid_dist)
    else:
        if x_range > y_range:
            mid_dist = int(x_range / 2)
            border = int((border_size / 2) * x_range)
        else:
            mid_dist = int(y_range / 2)
            border = int((border_size / 2) * y_range)
    # define the crop coordinates for the square
    coord_square = [
        mid_point[0] - mid_dist - border,
        mid_point[0] + mid_dist + border,
        mid_point[1] - mid_dist - border,
        mid_point[1] + mid_dist + border,
    ]
    coord_square = tuple(coord_square)
    return coord_square

crop_coord = auto_crop_scanpy(adata, border_size=0.1)


# --- Also load colors for HEV classes

color_dir = os.path.join(project_dir, 'resources', 'plotting_colors')
color_file = 'TLS_HEV_colors.csv'

# Load the color file
color_df = pd.read_csv(os.path.join(color_dir, color_file), index_col=0)

# Make a color_dict
color_dict = {}
for label, color in zip(color_df.index, color_df['color']):
    color_dict[label] = color


# --- Import the original high-resolution scan

# Find the original high-res image
image_dir = os.path.join(project_dir, 'data', 'images', 'he', 'aligned_images', 'visium_hd')
image_file = f'{sample_id}.ome.tif'
image_file_path = os.path.join(image_dir, image_file)

# Load the original high-res scan using tifffile
img = tifffile.imread(image_file_path)

# Add the scan to the anndata object for easy plotting
original_hires_key = 'original_hires'

# -- Calculate the scale factor
# Divide the height of the original hires image by the height of the downsampled hires image
size_fact = img.shape[0]/adata.uns['spatial'][sample_id]['images']['hires'].shape[0]

# Extract the downsampled hires scale factor
hires_fact = adata.uns['spatial'][sample_id]['scalefactors']['tissue_hires_scalef']

# Calculate the scale factor for the original hires image and add to the anndata object
fulres_fact = size_fact * hires_fact

# Get the microns per pixel for the scalebar
microns_per_pixel = adata.uns['spatial'][sample_id]['scalefactors']['microns_per_pixel']


# --- Load the area to crop to

crop_dir = os.path.join(project_dir, 'resources', 'areas_of_interest', 'visium_hd')
filename = f'{sample_id}_HEV_LR_area.geojson'

crop_df = gpd.read_file(os.path.join(crop_dir, filename))

# Extract xmin, xmax, ymin, ymax
min_x, min_y, max_x, max_y = crop_df.total_bounds

# Make sure all are integers
min_x = int(min_x)
max_x = int(max_x)
min_y = int(min_y)
max_y = int(max_y)


# --- Prepare for saving

out_dir = os.path.join(res_dir, 'figures')
os.makedirs(out_dir, exist_ok=True)


# --- Prepare for plotting

nrows = 1
ncols = 1

# Function to move plots so they don't overlap, moving both in x & y direction
def move_plots_x_y(axs, x_offset=0.02, y_offset=0.02, remove_title=True):
    for i, axs_row in enumerate(axs):
        for j, ax in enumerate(axs_row):
            if i != 0:
                if remove_title:
                    ax.set_title('')
            # Adjust the positions of the subplots
            ax_pos = ax.get_position()
            # Shift the subplot to the desired location + set the width so it's equal for all plots
            if i == 0:
                new_pos_y = ax_pos.y0
            else:
                new_pos_y = ax_pos_row_prev.y0 - ax_pos.height - y_offset
            if j == 0:
                new_pos_x = ax_pos.x0
            else:
                new_pos_x = ax_pos_prev.x0 + ax_pos_prev.width + x_offset
            ax.set_position([new_pos_x, new_pos_y, ax_pos.width, ax_pos.height])
            ax_pos_prev = ax.get_position()
        ax_pos_row_prev = axs_row[0].get_position()


# --- Define an area to crop to

# Function to crop the image with safe coordinates
# Bc. sometimes the calculated crop coordinates might be outside of the image area (especially if the original image has been cropped)
def crop_image(img, coords):
    # Clip the coordinates to ensure they are within the image boundaries
    x_min = max(coords[0], 0)
    x_max = min(coords[1], img.shape[1])
    y_min = max(coords[2], 0)
    y_max = min(coords[3], img.shape[0])
    # Crop the image using the clipped coordinates
    crop_img = img[y_min:y_max, x_min:x_max]
    return crop_img, (x_min, x_max, y_min, y_max)

# Function to crop & subset adata (avoid repetition)
def crop_and_subset_anndata(
    adata,
    adata_crop,
    img,
    sample_id,
    original_hires_key,
    zoom_coord_fulres,
    width_fulres,
    fulres_fact
):
    """
    Crop the H&E image and subset the AnnData object to the cropped region.
    Returns the cropped AnnData object (adata_sub), crop image, and crop coordinates.
    """
    coords_group = [
        zoom_coord_fulres[0],
        zoom_coord_fulres[0] + width_fulres,
        zoom_coord_fulres[1],
        zoom_coord_fulres[1] + width_fulres
    ]
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
    adata_crop.uns['spatial'][sample_id]['images'][original_hires_key] = crop_img
    # Shift the spot coordinates
    spot_coords = adata.obsm['spatial']
    shifted_coords = spot_coords - np.array([adjusted_coords[0], adjusted_coords[2]])
    adata_crop.obsm['spatial'] = shifted_coords
    # Subset to just the spots that are actually in the crop (to speed things up)
    xmax = xmin + width_crop
    ymax = ymin + height_crop
    # Also add scale factor
    adata_crop.uns['spatial'][sample_id]['scalefactors'][f'tissue_{original_hires_key}_scalef'] = fulres_fact
    # Extract spatial coordinates (assuming they are in the form of an (N, 2) array)
    adata_crop.obs["pxl_col_in_fullres"] = adata_crop.obsm["spatial"][:, 0] * fulres_fact  # X-coordinates
    adata_crop.obs["pxl_row_in_fullres"] = adata_crop.obsm["spatial"][:, 1] * fulres_fact # Y-coordinates
    # Subset the AnnData object based on the crop region
    mask = (
        (adata_crop.obs["pxl_col_in_fullres"] >= xmin) & 
        (adata_crop.obs["pxl_col_in_fullres"] < xmax) &
        (adata_crop.obs["pxl_row_in_fullres"] >= ymin) & 
        (adata_crop.obs["pxl_row_in_fullres"] < ymax)
    )
    adata_sub = adata_crop[mask]
    return adata_sub, crop_img, (xmin, ymin, width_crop, height_crop)


# --- Also prepare a dict to rename the labels in the plots

rename_dict = {
    'near_TLS_HEV': 'HEV (TLS)',
    'non_TLS_HEV': 'HEV-like (stroma)',
    'nearby_spots': f'Perivascular bins (<{bin_size}µm)' if bin_size > 0 else 'Non-HEV',
}


# --- Plot 1: Overview of scores at multiple resolutions

# Crop the adata
adata_crop = adata.copy()

adata_crop, crop_img, (xmin, ymin, width_crop, height_crop) = crop_and_subset_anndata(
    adata=adata,
    adata_crop=adata_crop,
    img=img,
    sample_id=sample_id,
    original_hires_key=original_hires_key,
    zoom_coord_fulres=(min_x, min_y),
    width_fulres=max_x - min_x,
    fulres_fact=fulres_fact
)

# Set the text size
sc.set_figure_params(fontsize=20)

# Set the scalebar size
scalebar_val = 300  # in microns
scalebar_units = 'µm'

ncols = 2
nrows = 1

# Make figure
fig, axs = plt.subplots(nrows, ncols, figsize=(5, 5))

# First plot the cropped image
ax = axs[0]
ax.imshow(crop_img,)
ax.set_title('Cropped Image')

# Remove xticks & yticks
ax.set_xticks([])
ax.set_yticks([])

# Add scalebar
scalebar = ScaleBar(microns_per_pixel, units=scalebar_units, location='upper left', box_alpha=0.5, fixed_value=scalebar_val,
                    frameon=True, scale_loc="top", pad=0.3, height_fraction=0.02)
ax.add_artist(scalebar)

# Now plot the HEV groups
ax = axs[1]

# Plot the HEV groups
# Prepare a listed colormap for squidpy
all_cats = adata_crop.obs['cell_type'].cat.categories.tolist()
color_list = [to_rgb(color_dict[cat]) for cat in all_cats if cat in color_dict]

# Convert to listed colormap
listed_cmap = ListedColormap(color_list)

sq.pl.spatial_scatter(adata_crop, color='cell_type', crop_coord=(xmin, ymin, width_crop, height_crop),
    ax=ax, img_res_key=original_hires_key, palette=listed_cmap, legend_loc='none', img_alpha=0.5,
    shape='square', alpha=1)

# Add scalebar
scalebar = ScaleBar(microns_per_pixel, units=scalebar_units, location='upper left', box_alpha=0.5, fixed_value=scalebar_val,
                    frameon=True, scale_loc="top", pad=0.3, height_fraction=0.02)
ax.add_artist(scalebar)

# Remove the x & y labels of all axs
for ax in axs.flat:
    ax.set_xlabel('')
    ax.set_ylabel('')

# Resize the plots
for ax in axs.flat:
    ax_pos = ax.get_position()
    # Set the width to be equal for all plots
    ax.set_position([ax_pos.x0, ax_pos.y0, 0.8, 0.8])

move_plots_x_y([axs], x_offset=0.05, y_offset=0.05, remove_title=True)

# Set titles
titles = ['H&E', 'HEV groups']

for ax, title in zip(axs.flat, titles):
    ax.set_title(title, fontweight='bold')

# Add a legend
legend_labels = all_cats

handles = []
for label in legend_labels:
    if label in color_dict:
        color = color_dict[label]
        handles.append(mpatches.Patch(facecolor=color, label=label, edgecolor='black', linewidth=0))

ax.legend(handles=handles, loc='upper center', bbox_to_anchor=(0.5, 0), ncols=2, columnspacing=-4,
               frameon=False, handlelength=0.7, handletextpad=0.5, borderpad=0, fontsize=15)

# Update the labels of the legend
for text in ax.get_legend().get_texts():
    if text.get_text() in rename_dict:
        text.set_text(rename_dict[text.get_text()])

# Position the title
title_x = (axs[0].get_position().x0 + axs[-1].get_position().x1) / 2
title_y = axs[0].get_position().y1 + 0.2

plt.suptitle(f'{sample_id} ({subtype}) - HEV grouping used for LR analysis', y=title_y, x=title_x, fontweight='bold')


# --- Save plot with high resolution

filename = 'Fig_4F.pdf'

plt.savefig(os.path.join(out_dir, filename), bbox_inches='tight')


# --- Finish

print('Script successfully completed')
print(f'Saved plot to: {os.path.join(out_dir, filename)}')
