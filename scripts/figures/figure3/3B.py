#!/usr/bin/env python3

# --- Aim
# Plot endothelial deconvolution and pathology annotations in zoomed areas (Figure 3B)
# Shows H&E, deconvolution results, and pathology annotations for different tissue regions
# This script is set up for P29_HD

# --- Define environment
# env: vis-hd

# --- Load packages
import argparse
import os
from datetime import datetime

import anndata as ad
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patheffects as path_effects
from matplotlib_scalebar.scalebar import ScaleBar
from matplotlib.colors import ListedColormap, to_rgb

import numpy as np
import pandas as pd
import scanpy as sc
import scanpy.plotting.palettes as sc_palettes
import squidpy as sq
import tifffile

from sophiesfunctions.color_functions import generate_color_variations

sc.set_figure_params(facecolor="white", figsize=(8, 8))
sc.settings.verbosity = 3

mpl.rcParams['pdf.fonttype'] = 42


# --- Argparse arguments
parser = argparse.ArgumentParser(description="Plot vessel classification with pathology annotations for Figure 3B",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument("--res_dir", help="Path to results directory")
parser.add_argument("--h5ad_dir", help="Path to directory for .h5ad output files")
parser.add_argument("--project_dir", help="Project directory (root)")
parser.add_argument("--sample_id", help="Sample ID")
parser.add_argument("--sample_sheet", help="Path to sample sheet CSV file", default=None)

args = vars(parser.parse_args())

# Arguments
project_dir = args["project_dir"]
h5ad_dir = args["h5ad_dir"]
res_dir = args["res_dir"]
sample_id = args["sample_id"]

print(f'Sample: {sample_id}')


# --- Extract some info from the sample sheet
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

# Merge the two DataFrames on the common column 'tissue_id'
merged_df = pd.merge(sample_sheet, tissue_metadata, on='tissue_id')
subtype = merged_df[merged_df['section_id'] == sample_id]['clinical_subtype'].values[0]


# --- Load the adata objects

# Define which resolutions to use
resolutions = ['square_016um', 'square_008um']
adata_dict = {}

for resolution in resolutions:
    # Load the adata object
    adata = ad.read_h5ad(os.path.join(h5ad_dir, 'visium_hd', resolution, f'{sample_id}.h5ad'))
    # Add to dictionary
    adata_dict[resolution] = adata


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
        return most_recent_dir_path
    else:
        print("No valid directories found.")
        return None

# Find most recent results
resolution = 'square_016um'
threshold = 0.5
threshold_adj = str(threshold).split('.')[-1]

# Specify path
vessel_dir = os.path.join(res_dir, '04_vascular_microenvironment', '02_identify_all_vasculature',
                          '01_classify_vessels', 'csv_files')

most_recent_dir_path = find_most_recent(vessel_dir)

# And also load the normalized endothelial cell numbers
filename = f'{sample_id}_endothelial_numbers_normalized.csv'
endo_numbers = pd.read_csv(os.path.join(most_recent_dir_path, filename), index_col=0)[endothelial]

adata_dict[resolution].obs[endothelial] = endo_numbers[endothelial]


# --- Also load the pathology annotations

annotation_dir = os.path.join(project_dir, 'data', 'pathology_annotations', 'visium_hd')

filename = f'{sample_id}_pathology_annotations.csv'

annotation_key = 'path_anno'

annotations = pd.read_csv(os.path.join(annotation_dir, filename), index_col=0)
adata_dict['square_008um'].obs[annotation_key] = annotations


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

color_dict_path = {}

# And also define a plot order
plot_order = []

for top_cat, labels in annotation_labels.items():
    top_color = top_cat_colors[top_cat]
    # Get the variations
    colors = generate_color_variations(top_color, len(labels))
    # Add to the color_dict
    for n, label in enumerate(labels):
        color_dict_path[label] = colors[n]
        plot_order.append(label)

# Set the plot_order
adata_dict['square_008um'].obs[annotation_key] = pd.Categorical(adata.obs[annotation_key], categories=plot_order, ordered=True)


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


# --- Import the original high-resolution scan

adata = adata_dict[resolution]

# Define key for adding to anndata
original_hires_key = 'original_hires'

# Find the original high-res image
image_dir = os.path.join(project_dir, 'data', 'images', 'visium_hd')
image_file = f'{sample_id}.ome.tif'
image_file_path = os.path.join(image_dir, image_file)

# Load the original high-res scan using tifffile
img = tifffile.imread(image_file_path)

# -- Calculate the scale factor
# Divide the height of the original hires image by the height of the downsampled hires image
size_fact = img.shape[0]/adata.uns['spatial'][sample_id]['images']['hires'].shape[0]

# Extract the downsampled hires scale factor
hires_fact = adata.uns['spatial'][sample_id]['scalefactors']['tissue_hires_scalef']

# Calculate the scale factor for the original hires image and add to the anndata object
adata.uns['spatial'][sample_id]['scalefactors'][f'tissue_{original_hires_key}_scalef'] = size_fact * hires_fact
print(f'Successfully imported fulres image {sample_id}')


# --- Specify crop areas (for sample P29_HD)

# Specify the top left corner of the crop area
crop_areas_fulres = {
    'invasive': [21833, 34836],
    'lymphocyte-rich': [51779, 42489],
    'mixed': [16678, 24170],
}

# Also specify a width
width_fulres = 6702

# Convert to mpl coords
hires_scalef = adata.uns['spatial'][sample_id]['scalefactors'][f'tissue_hires_scalef']

crop_areas = {}

for area, coords in crop_areas_fulres.items():
    coords_adj = [0, 0]
    # Convert the coordinates to mpl coords
    coords_adj[0] = int(coords[0] * hires_scalef)
    coords_adj[1] = int(coords[1] * hires_scalef)
    crop_areas[area] = coords_adj

width = int(width_fulres * hires_scalef)


# --- Function for moving plots

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


# --- Plot: Zoomed-in areas

sc.set_figure_params(fontsize=30)

# Get the microns per pixel for the scalebar
microns_per_pixel = adata.uns['spatial'][sample_id]['scalefactors']['microns_per_pixel']

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

alpha_img = 0.7

nrows = 3
ncols = len(crop_areas_fulres)

# Set the scalebar size
scalebar_val = 200  # in microns
scalebar_units = 'µm'

# List to store which annotations plotted
annotations_plotted = []

# Plotting parameters
vmin = 0
vmax = 0.15

# Initialize the figure
fig, axs = plt.subplots(nrows, ncols, figsize=(5, 5))

for col_index, (crop_area, coords) in enumerate(crop_areas_fulres.items()):
    # Create fresh copies for each crop area to avoid overwriting
    adata_crop = adata.copy()
    adata_crop_008um = adata_dict['square_008um'].copy()
    
    # -- Prepare an anndata that contains only the cropped area we want to plot
    # Calculate coords based on width
    coords_group = [coords[0], coords[0] + width_fulres, coords[1], coords[1] + width_fulres]
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
    scalef_tmp = adata.uns['spatial'][sample_id]['scalefactors'][f'tissue_{original_hires_key}_scalef']
    adata_crop.uns['spatial'][sample_id]['scalefactors'][f'tissue_{original_hires_key}_scalef'] = scalef_tmp
    # Extract spatial coordinates (assuming they are in the form of an (N, 2) array)
    adata_crop.obs["pxl_col_in_fullres"] = adata_crop.obsm["spatial"][:, 0] * scalef_tmp  # X-coordinates
    adata_crop.obs["pxl_row_in_fullres"] = adata_crop.obsm["spatial"][:, 1] * scalef_tmp  # Y-coordinates
    # Subset the AnnData object based on the crop region
    mask = (
        (adata_crop.obs["pxl_col_in_fullres"] >= xmin) &
        (adata_crop.obs["pxl_col_in_fullres"] < xmax) &
        (adata_crop.obs["pxl_row_in_fullres"] >= ymin) &
        (adata_crop.obs["pxl_row_in_fullres"] < ymax)
    )
    adata_sub = adata_crop[mask]
    
    # -- Plot the H&E
    n = 0
    ax = axs[n][col_index]
    title = crop_area.capitalize()
    if crop_area == 'lymphocyte-rich':
        title = 'Lymphocyte-\nrich'
    sc.pl.spatial(adata_sub, color=None, bw=False, alpha_img=1, crop_coord=(xmin, width_crop, ymin, height_crop), title=title,
                  na_in_legend=False, ax=ax, img_key=original_hires_key)
    ax.set_title(ax.get_title(), fontweight='bold')
    # Add the scalebar
    scalebar = ScaleBar(microns_per_pixel, units=scalebar_units, location='lower right', box_alpha=0.5, fixed_value=scalebar_val,
                        frameon=False, scale_loc="none", pad=0.3, height_fraction=0.02)
    ax.add_artist(scalebar)
    
    # -- Plot the endothelial decon
    n += 1
    ax = axs[n][col_index]
    for i, endo_type in enumerate(endothelial):
        if i == 0:
            alpha_img = 0.7
            img_key = original_hires_key
        else:
            alpha_img = 0
            img_key = None
        sc.pl.spatial(adata_sub, color=endo_type, bw=False, alpha_img=alpha_img, crop_coord=(xmin, width_crop, ymin, height_crop), title=f'',
                      na_in_legend=False, ax=ax, colorbar_loc=None, legend_loc=None, img_key=original_hires_key,
                      vmin=vmin, vmax=vmax, cmap=cmap_dict[endo_type])
    # Add the scalebar
    scalebar = ScaleBar(microns_per_pixel, units=scalebar_units, location='lower right', box_alpha=0.5, fixed_value=scalebar_val,
                        frameon=False, scale_loc="none", pad=0.3, height_fraction=0.02)
    ax.add_artist(scalebar)
    
    # -- Plot the pathology annotations
    n += 1
    ax = axs[n][col_index]
    # Prepare crop object for 008um
    adata_crop_008um.uns['spatial'][sample_id]['images'][original_hires_key] = crop_img
    # Shift the spot coordinates
    spot_coords = adata_dict['square_008um'].obsm['spatial']
    shifted_coords = spot_coords - np.array([adjusted_coords[0], adjusted_coords[2]])
    adata_crop_008um.obsm['spatial'] = shifted_coords
    # Subset to just the spots that are actually in the crop (to speed things up)
    xmax = xmin + width_crop
    ymax = ymin + height_crop
    # Also add scale factor
    scalef_tmp = adata.uns['spatial'][sample_id]['scalefactors'][f'tissue_{original_hires_key}_scalef']
    adata_crop_008um.uns['spatial'][sample_id]['scalefactors'][f'tissue_{original_hires_key}_scalef'] = scalef_tmp
    # Extract spatial coordinates (assuming they are in the form of an (N, 2) array)
    adata_crop_008um.obs["pxl_col_in_fullres"] = adata_crop_008um.obsm["spatial"][:, 0] * scalef_tmp  # X-coordinates
    adata_crop_008um.obs["pxl_row_in_fullres"] = adata_crop_008um.obsm["spatial"][:, 1] * scalef_tmp  # Y-coordinates
    # Subset the AnnData object based on the crop region
    mask = (
        (adata_crop_008um.obs["pxl_col_in_fullres"] >= xmin) &
        (adata_crop_008um.obs["pxl_col_in_fullres"] < xmax) &
        (adata_crop_008um.obs["pxl_row_in_fullres"] >= ymin) &
        (adata_crop_008um.obs["pxl_row_in_fullres"] < ymax)
    )
    adata_sub = adata_crop_008um[mask]
    # Plot the pathology annotations
    # Extract a listed colormap for squidpy
    color_list = [to_rgb(color_dict_path[cat]) for cat in plot_order if cat in adata_sub.obs[annotation_key].cat.categories]
    listed_cmap = ListedColormap(color_list)
    sq.pl.spatial_scatter(adata_sub, color=annotation_key, crop_coord=(xmin, ymin, width_crop, height_crop),
                          title=f'',
                          ax=ax, img_res_key=original_hires_key, palette=listed_cmap, legend_loc=None, img_alpha=1,
                          shape='square', alpha=0.7)
    # Add the scalebar
    scalebar = ScaleBar(microns_per_pixel, units=scalebar_units, location='lower right', box_alpha=0.5, fixed_value=scalebar_val,
                        frameon=False, scale_loc="none", pad=0.3, height_fraction=0.02)
    ax.add_artist(scalebar)
    # Keep track of which annotations have been plotted
    annotations_plotted += list(adata_sub.obs[annotation_key].dropna().unique())

# Remove the x & y labels of all axs
for ax in axs.flat:
    ax.set_xlabel('')
    ax.set_ylabel('')

# Resize the plots
for ax in axs.flat:
    ax_pos = ax.get_position()
    # Set the width to be equal for all plots
    ax.set_position([ax_pos.x0, ax_pos.y0, 0.8, 0.8])

move_plots_x_y(axs, x_offset=0.05, y_offset=0.05, remove_title=True)

# Set ylabel for the first column
labels = ['H&E', 'Deconvolution', 'Pathology\nannotations']

for i in range(nrows):
    axs[i][0].set_ylabel(labels[i], fontsize=30, fontweight='bold')

# Calculate the suptitle position
title_x = (axs[0][0].get_position().x0 + axs[0][-1].get_position().x1) / 2
title_y = axs[0][0].get_position().y1 + 0.3  # Adjust y position slightly above the first row

plt.suptitle(f'{sample_id} ({subtype}) - Endothelial decon results & pathology annotations - zoomed areas (vmax={vmax}, scalebar={scalebar_val}{scalebar_units})',
             y=title_y, x=title_x, fontweight='bold',)

# Add colorbars to the last column second row plot
target_ax = axs[1][-1]

# Adjust the labels
endo_name_dict = {
    'Endothelial CXCL12': 'CXCL12+',
    'Endothelial RGS5': 'RGS5+',
    'Endothelial ACKR1': 'ACKR1+',
    'Endothelial Lymphatic LYVE1': 'LYVE1+',
}

# Add colorbars per score
for i, label in enumerate(list(endo_name_dict.keys())[::-1]):
    cmap = cmap_dict[label]
    norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
    spacing = 0.1
    mid_y = 0.5
    num_cbars = len(endo_name_dict)
    y_pos = mid_y + (i - (num_cbars - 1) / 2) * spacing  # Distribute around the midpoint
    cbar_ax = target_ax.inset_axes([1.08, y_pos, 0.25, 0.05])
    cbar = fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap),
                        cax=cbar_ax, orientation='horizontal')
    cbar.ax.set_ylabel(f'{endo_name_dict[label]}', rotation=0, ha='left', va='center', fontsize=20)
    cbar.ax.yaxis.set_label_position("right")
    cbar.ax.yaxis.set_label_coords(1.1, 1.05)
    cbar.ax.tick_params(labelsize=15)
    # Set xticks as min & max
    cbar.ax.set_xticks([vmin, vmax])
    cbar.ax.set_xticklabels([f'Low', f'High'])
    # Keep xticklabels only on the bottom plot
    if i != 0:
        cbar.ax.set_xticklabels([])
    # Give the first one a title
    if i == num_cbars - 1:
        cbar.ax.set_title('Endothelial type', fontsize=20, loc='left', pad=15)


# Add a pathology legend to the third row last column plot
target_ax = axs[2][-1]
annotations_plotted = list(set(annotations_plotted))  # Remove duplicates

# Sort according to plot_order
annotations_plotted = [cat for cat in plot_order if cat in annotations_plotted]

# Create a legend for the pathology annotations
handles = []
for cat in annotations_plotted:
    if cat in color_dict_path:
        color = color_dict_path[cat]
        # Adjust the label (remove _ and capitalize)
        label = cat.replace('_', ' ').capitalize()
        handle = mpl.patches.Patch(color=color, label=label, alpha=0.7)
        handles.append(handle)

target_ax.legend(handles=handles, loc='center left', bbox_to_anchor=(1.01, 0.5), fontsize=20,
                 handlelength=0.7, frameon=False)


# --- Save the output

out_dir = os.path.join(project_dir, 'results', 'figures')
os.makedirs(out_dir, exist_ok=True)

filename = 'Fig_3B.pdf'
plt.savefig(os.path.join(out_dir, filename), bbox_inches='tight')

print(f"Saved plot to: {os.path.join(out_dir, filename)}")
print("Script successfully completed")
