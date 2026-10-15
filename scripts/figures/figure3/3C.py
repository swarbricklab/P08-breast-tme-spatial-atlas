#!/usr/bin/env python3

# --- Aim
# Plot vessel morphology showing different endothelial types (Figure 3C)
# Shows high-resolution H&E and deconvolution results for different vessel types

# --- Define environment
# env: vis-hd

# --- Load packages
import argparse
import os

import anndata as ad
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib_scalebar.scalebar import ScaleBar
from matplotlib.colors import ListedColormap, to_rgb

import numpy as np
import pandas as pd
import scanpy as sc
import tifffile

sc.set_figure_params(facecolor="white", figsize=(8, 8))
sc.settings.verbosity = 3

mpl.rcParams['pdf.fonttype'] = 42


# --- Argparse arguments
parser = argparse.ArgumentParser(description="Plot vessel morphology for Figure 3C",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument("-s", "--samples",
    nargs="*",  # 0 or more values expected => creates a list
    type=str,
    help="List of samples")
parser.add_argument("--res_dir", help="Path to results directory")
parser.add_argument("--h5ad_dir", help="Path to directory for .h5ad output files")
parser.add_argument("--project_dir", help="Project directory (root)")
parser.add_argument("--sample_sheet", help="Path to sample sheet CSV file", default=None)

args = vars(parser.parse_args())

# Arguments
project_dir = args["project_dir"]
h5ad_dir = args["h5ad_dir"]
res_dir = args["res_dir"]

print(f'Project directory: {project_dir}')


# --- Extract some info from the sample sheet
# Use sample_sheet argument if provided, else default
if args.get("sample_sheet"):
    sample_sheet_path = args["sample_sheet"]
else:
    sample_sheet_path = os.path.join(project_dir, 'config', 'sample_sheet_visium_hd.csv')

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


# --- Define crop areas per vessel type

crop_areas = {
    'CXCL12': {'top_left': (53807, 51960), 'width': 1496, 'sample_id': 'P27_HD'},
    'RGS5': {'top_left': (44403, 33135), 'width': 970, 'sample_id': 'P26_HD'},
    'LYVE1': {'top_left': (14452, 23804), 'width': 1512, 'sample_id': 'P26_HD'},
    'ACKR1': {'top_left': (47836, 18601), 'width': 1313, 'sample_id': 'P32_HD'}
}

# Figure out which samples are required
sample_ids = [coords['sample_id'] for coords in crop_areas.values()]
sample_list = sorted(list(set(sample_ids)))  # Remove duplicates


# --- Load the adata object

# Define which resolution want to use
resolution = 'square_016um'
adata_dict = {}

for sample in sample_list:
    # Load the adata object
    adata = ad.read_h5ad(os.path.join(h5ad_dir, 'visium_hd', resolution, f'{sample}.h5ad'))
    # Add to dictionary
    adata_dict[sample] = adata


# --- Load miniatlas classifications

decon_level = 'minor'

# To make it easier to visualize trends by lineage
class_dir = os.path.join(project_dir, 'resources', 'cell_types', 'miniatlas')
class_file = 'miniatlas_subset_classifications.csv'

class_df = pd.read_csv(os.path.join(class_dir, class_file))

# Get endothelial classes
endothelial = list(class_df[class_df['miniatlas_major'] == 'Endothelial'][f'miniatlas_{decon_level}'].unique())
endothelial = sorted(endothelial)

# Also load the colors for the endothelial types
colors_file = os.path.join(project_dir, 'resources', 'plotting_colors', 'endothelial_colors.csv')
colors_df = pd.read_csv(colors_file)

# Make a color_dict
color_dict = colors_df.set_index('endothelial_type')['color'].to_dict()


# --- Load normalized endothelial numbers

threshold = 0.5
threshold_adj = str(threshold).split('.')[-1]

# Specify path
vessel_dir = os.path.join(res_dir, '04_vascular_microenvironment', '02_identify_all_vasculature',
                          '01_classify_vessels', 'csv_files')

most_recent_dir_path = vessel_dir

for sample in sample_list:
    # And also load the normalized endothelial cell numbers
    filename = f'{sample}_endothelial_numbers_normalized.csv'
    endo_numbers = pd.read_csv(os.path.join(most_recent_dir_path, filename), index_col=0)[endothelial]
    adata_dict[sample].obs[endothelial] = endo_numbers[endothelial]


# --- Generate color maps for endothelial types

# Also make colormaps per score
white_spacing = 20

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

# Create a colormap for each endothelial type
color_list_rgb = []
for value in endothelial:
    rgb = tuple(int(c * 255) for c in to_rgb(color_dict[value]))
    color_list_rgb.append(rgb)

cmap_dict = {}
for n, rgb in enumerate(color_list_rgb):
    cmap = create_colormap(rgb[0], rgb[1], rgb[2])
    cmap_dict[endothelial[n]] = cmap


# --- Some functions for plotting

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

# Function to crop & subset adata (avoid repetition)
def crop_and_subset_anndata(
    adata,
    adata_crop,
    img,
    sample_id,
    original_hires_key,
    zoom_coord_fulres,
    width_fulres,
    fulres_fact,
    aspect_ratio=1
):
    """
    Crop the H&E image and subset the AnnData object to the cropped region.
    Returns the cropped AnnData object (adata_sub), crop image, and crop coordinates.
    """
    coords_group = [
        zoom_coord_fulres[0],
        zoom_coord_fulres[0] + int(width_fulres*aspect_ratio),
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
    # Have a grace region of 10% around the crop to avoid losing spots that are on the edge
    grace_region = 0.1  # 10% grace region
    xmin_grace = xmin - int(width_crop * grace_region)
    xmax_grace = xmax + int(width_crop * grace_region)
    ymin_grace = ymin - int(height_crop * grace_region)
    ymax_grace = ymax + int(height_crop * grace_region)
    mask = (
        (adata_crop.obs["pxl_col_in_fullres"] >= xmin_grace) & 
        (adata_crop.obs["pxl_col_in_fullres"] < xmax_grace) &
        (adata_crop.obs["pxl_row_in_fullres"] >= ymin_grace) & 
        (adata_crop.obs["pxl_row_in_fullres"] < ymax_grace)
    )
    adata_sub = adata_crop[mask]
    return adata_sub, crop_img, (xmin, ymin, width_crop, height_crop)


# --- Make a dict of cropped adatas

adata_dict_crop = {}

for sample in adata_dict.keys():
    adata_dict_crop[sample] = adata_dict[sample].copy()


# --- Adjust the coords so each area is the same

# Find the max width
max_width = 0
for coords in crop_areas.values():
    if coords['width'] > max_width:
        max_width = coords['width']

# Adjust the crop areas to have the same width
for area, coords in crop_areas.items():
    # Calculate the midpoint
    mid_x = coords['top_left'][0] + coords['width'] / 2
    mid_y = coords['top_left'][1] + coords['width'] / 2
    # Adjust the coordinates
    new_x = int(mid_x - max_width / 2)
    new_y = int(mid_y - max_width / 2)
    coords['top_left'] = (new_x, new_y)
    # Update the coords
    crop_areas[area] = {
        'sample_id': coords['sample_id'],
        'top_left': (new_x, new_y),
        'width': max_width
    }


# --- Load fulres images & make some anndatas cropped

cropped_adata_dict = {}

# Define key for adding to anndata
original_hires_key = 'original_hires'

# Because the resolution of 'hires' isn't great, everything becomes blurry once you zoom in
# Therefore use the original image for plotting the zoomed in plots
aspect_ratio = 1

# Find the original high-res image
image_dir = os.path.join(project_dir, 'data', 'images', 'he', 'aligned_images', 'visium_hd')

for n, endo in enumerate(crop_areas):
    # Get the crop coordinates
    coords = crop_areas[endo]
    sample = coords['sample_id']
    top_left = coords['top_left']
    width = coords['width']
    # Get the adata
    adata = adata_dict[sample]
    # Create a fresh copy for each vessel type to avoid overwriting
    adata_crop = adata_dict[sample].copy()
    
    # Load the original high-res scan using tifffile
    image_file = f'{sample}.ome.tif'
    image_file_path = os.path.join(image_dir, image_file)
    img = tifffile.imread(image_file_path)
    
    # -- Calculate the scale factor
    # Divide the height of the original hires image by the height of the downsampled hires image
    size_fact = img.shape[0]/adata.uns['spatial'][sample]['images']['hires'].shape[0]
    # Extract the downsampled hires scale factor
    hires_fact = adata.uns['spatial'][sample]['scalefactors']['tissue_hires_scalef']
    # Calculate the scale factor for the original hires image and add to the anndata object
    fulres_fact = size_fact * hires_fact
    adata.uns['spatial'][sample]['scalefactors'][f'tissue_{original_hires_key}_scalef'] = fulres_fact
    # Crop the H&E image
    adata_sub, crop_img, (xmin, ymin, width_crop, height_crop) = crop_and_subset_anndata(
        adata=adata,
        adata_crop=adata_crop,
        img=img,
        sample_id=sample,
        original_hires_key=original_hires_key,
        zoom_coord_fulres=top_left,
        width_fulres=width,
        fulres_fact=fulres_fact,
        aspect_ratio=aspect_ratio
    )
    cropped_adata_dict[endo] = adata_sub.copy()

# Set the scalebar size
scalebar_val = 50  # in microns
scalebar_units = 'µm'


# --- Plot 1: Overview of vessels

# Set fontsize
sc.set_figure_params(fontsize=18)

ncols = len(crop_areas)
nrows = 2

vmax = 0.15

# Define a plotting order
plot_order = ['CXCL12', 'RGS5', 'ACKR1', 'LYVE1']

# Create a dict with labels
alt_labels = {
    'CXCL12': 'Arterial',
    'RGS5': 'Capillary',
    'ACKR1': 'Venous',
    'LYVE1': 'Lymphatic',
}

fig, axs = plt.subplots(nrows, ncols, figsize=(5, 5))

for n, endo in enumerate(plot_order):
    adata = cropped_adata_dict[endo]
    # Get the matching endothelial subtype
    matching_endo = [subtype for subtype in endothelial if endo in subtype]
    matching_endo = matching_endo[0]  # Take the first match
    # Get the crop coordinates
    coords = crop_areas[endo]
    sample = coords['sample_id']
    top_left = coords['top_left']
    width = coords['width']
    # Plot the cropped H&E image
    crop_coord = [top_left[0], top_left[0] + width*aspect_ratio, top_left[1], top_left[1] + width]
    ax = axs[0][n]
    hires_fact = adata.uns['spatial'][sample]['scalefactors']['tissue_hires_scalef']
    # Plot all endothelial decon results
    sc.pl.spatial(adata, color=None, ax=ax, bw=False, alpha_img=1, title=f'{endo}',
                    crop_coord=(0, int(width*aspect_ratio), 0, width), img_key=original_hires_key)
    ax.set_title(alt_labels[endo], fontweight='bold')
    # Get the microns per pixel for the scalebar
    microns_per_pixel = adata.uns['spatial'][sample]['scalefactors']['microns_per_pixel']
    frameon = False
    scale_loc = 'none'
    scalebar = ScaleBar(microns_per_pixel, units=scalebar_units, location='lower right', box_alpha=0.5, fixed_value=scalebar_val,
                    frameon=frameon, scale_loc=scale_loc, pad=0.3, height_fraction=0.04)
    ax.add_artist(scalebar)
    # Plot the normalized endothelial values
    ax = axs[1][n]
    sc.pl.spatial(adata, color=matching_endo, ax=ax, bw=False, alpha_img=0.7, cmap=cmap_dict[matching_endo], vmax=vmax, vmin=0,
                    crop_coord=(0, int(width*aspect_ratio), 0, width), title=f'', img_key=original_hires_key, colorbar_loc=None,)
    scalebar = ScaleBar(microns_per_pixel, units=scalebar_units, location='lower right', box_alpha=0.5, fixed_value=scalebar_val,
                    frameon=False, scale_loc='none', pad=0.3, height_fraction=0.04)
    ax.add_artist(scalebar)
    # Add colorbar below the plot
    norm = mpl.colors.Normalize(vmin=0, vmax=vmax)
    cbar_ax = ax.inset_axes([0.2, -0.15, 0.6, 0.1])  # Adjust the position of the colorbar
    cbar = fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap_dict[matching_endo]), 
                        cax=cbar_ax, orientation='horizontal')
    cbar.ax.set_xlabel(f'{endo}+', fontsize=15)
    cbar.ax.tick_params(labelsize=15)
    cbar.ax.set_xticks([0, vmax])

# Remove the x & y labels of all axs
for ax in axs.flat:
    ax.set_xlabel('')
    ax.set_ylabel('')

# Resize the plots
for ax in axs.flat:
    ax_pos = ax.get_position()
    # Set the width to be equal for all plots
    ax.set_position([ax_pos.x0, ax_pos.y0, 0.4*aspect_ratio, 0.4])

# Set ylabels
ylabels = ['H&E', 'Deconvolution',]

for i, ax in enumerate(axs[:, 0]):
    ax.set_ylabel(ylabels[i], fontweight='bold')

move_plots_x_y(axs, x_offset=0.03, y_offset=0.03, remove_title=True)

# Calculate position for title
title_x = axs[0, 0].get_position().x0 + (axs[0, -1].get_position().x1 - axs[0, 0].get_position().x0) / 2
title_y = axs[0, 0].get_position().y1 + 0.2  # Slightly above the first row

plt.suptitle(f'Vessel morphology (scalebar={scalebar_val} {scalebar_units})', 
             fontweight='bold', y=title_y, x=title_x, ha='center')


# --- Save the output

out_dir = os.path.join(res_dir, 'figures')
os.makedirs(out_dir, exist_ok=True)

filename = 'Fig_3C.pdf'
plt.savefig(os.path.join(out_dir, filename), bbox_inches='tight')

print(f"Saved plot to: {os.path.join(out_dir, filename)}")
print("Script successfully completed")
