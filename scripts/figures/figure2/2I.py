#!/usr/bin/env python3

# --- Aim
# Plot the TLS niche signatures on lymph node sample P12_LN_Vis
# Shows H&E, GC annotations (from pathology), and niche signatures
# With overview + 2 zoomed areas (Primary follicle, Secondary follicle)
# This creates Figure 2I

# --- Define environment
# env: shapely

# --- Load packages

import os
from datetime import datetime
import ast
import argparse
import sys

import scanpy as sc
import anndata as ad
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

import matplotlib as mpl
from matplotlib.colors import to_rgb
from matplotlib.colors import ListedColormap
from matplotlib_scalebar.scalebar import ScaleBar

import geopandas as gpd
import tifffile

mpl.rcParams['pdf.fonttype'] = 42


# --- Argparse arguments
parser = argparse.ArgumentParser(description="Plot TLS model overview for Figure 2I",
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

# Override sample_id
sample_id = "P12_LN_Vis"

print(f'Sample: {sample_id}')

sys.path.append(os.path.join(project_dir, 'scripts'))
from helper_functions import auto_crop_scanpy


# --- Load sample sheet

if args.get("sample_sheet"):
    sample_sheet_path = args["sample_sheet"]
else:
    sample_sheet_path = os.path.join(project_dir, 'config', 'sample_sheet_visium_lymph_nodes.csv')

sample_sheet = pd.read_csv(sample_sheet_path)
sample_lut = sample_sheet.set_index('section_id')


# --- Load adata object

adata = ad.read_h5ad(os.path.join(h5ad_dir, 'visium', f'{sample_id}.h5ad'))


# --- Load tissue annotations

annotation_dir = os.path.join(project_dir, 'data', 'pathology_annotations', 'visium_lymph_nodes')
tissue_annot_file = os.path.join(annotation_dir, f'{sample_id}_tissue_annotations.csv')

# Add the tissue annotations if available
tissue_annot_key = 'tissue_annot'

if os.path.exists(tissue_annot_file):
    print(f"Found tissue annotations for {sample_id}: {tissue_annot_file}")
    annotations = pd.read_csv(tissue_annot_file, index_col=0)
    # Map 'tissue' annotation to True/False
    adata.obs[tissue_annot_key] = annotations['annotation'] == 'tissue'
    # Subset the adata to just annotated spots
    adata = adata[adata.obs[tissue_annot_key] == True]
    print(f"Subset to {adata.n_obs} spots over tissue")
else:
    print(f"No tissue annotations found for {sample_id}")


# --- Calculate crop coords for optimal presentation

coords = auto_crop_scanpy(adata, border_size=0.2)
# Adjust the overview coords based on width_to_height_ratio
# Will be defined later after we load width_to_height_ratio
# This will be updated after width_to_height_ratio is set
# --- Set spot size for large capture area (11 mm)

spot_size = 1.8

# --- Load the TLS niche scores

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
        most_recent_dir_path = None
    return most_recent_dir_path

# Signature names mapping (variable name -> directory name)
signature_mapping = {
    'T_cell_zone': 'T_cell_zone_score',
    'B_cell_zone': 'B_cell_zone_score',
    'GC': 'GC_score'
}

# Note: For lymph nodes, there's a separate subdirectory
score_names = ['GC', 'B_cell_zone', 'T_cell_zone']
score_keys = {}

# Decide which number of top genes to use
top_genes = 20

# Load each signature score
for score_name in score_names:
    score_dir_name = signature_mapping[score_name]
    # Path to signature scores
    score_base_dir = os.path.join(res_dir, '03_TLS_niche_discovery', '04_map_niche_signatures', 
                                   '01_score_signatures', 'lymph_nodes', score_dir_name)
    
    # Find the most recent results
    most_recent_dir = find_most_recent(score_base_dir)
    
    if most_recent_dir is not None:
        csv_dir = os.path.join(most_recent_dir, 'csv_files')
        
        # Read in the results & add to anndata
        score_key = f'{score_name}_score_{top_genes}'
        score_file_name = f'{sample_id}_{score_dir_name}.csv'
        score_file_path = os.path.join(csv_dir, score_file_name)
        
        if os.path.exists(score_file_path):
            score_results = pd.read_csv(score_file_path, index_col=0)
            adata.obs[score_key] = score_results[score_key]
            score_keys[score_name] = score_key
            print(f"Loaded {score_name} score from {score_file_path}")
        else:
            print(f"Warning: Score file not found: {score_file_path}")
    else:
        print(f"Warning: No dated directory found for {score_name}")


# --- Load the GC annotations from GeoJSON

# Load the raw geojson annotations
anno_dir = os.path.join(project_dir, 'data', 'pathology_annotations', 'visium_lymph_nodes')

# Find geojson file for this sample
geojson_file = os.path.join(anno_dir, f'{sample_id}_GC_annotations.geojson')

if os.path.exists(geojson_file):
    print(f"Found geojson file: {geojson_file}")
    
    # Load the GeoJSON
    # This needs to be set so geopandas can read the file
    os.environ['PROJ_LIB'] = '/g/data/a56/conda/envs/usr/sl6147/share/proj'
    gdf = gpd.read_file(geojson_file)
    
    # Filter for GC annotations
    # If classification is a string, convert it to a dict and extract 'name'
    gdf = gdf[gdf['classification'].apply(lambda x: ast.literal_eval(x)['name'] == 'germinal_centre')]
    print(f"Found {len(gdf)} GC annotations")
else:
    print(f"Warning: GeoJSON file not found: {geojson_file}")
    gdf = None

# --- Prepare cmap per score

vmin = 0
vmax = 0.8

# Reverse order so GCs always end up on top
# Because areas with high GC score will also have high B cell zone score
score_names_r = score_names[::-1]

# Use function from cell2location to make colormap with increasing transparency gradient
# To be able to plot multiple scores in single plot
white_spacing = 20

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
}

color_dict_rgb = {}
for key, value in color_dict.items():
    rgb = tuple(int(c * 255) for c in to_rgb(value))
    color_dict_rgb[key] = rgb

cmap_dict = {}

for score_name in score_names:
    rgb = color_dict_rgb[score_name]
    cmap = create_colormap(rgb[0], rgb[1], rgb[2])
    cmap_dict[score_name] = cmap

# --- Load the full resolution image

# Find the original high-res image
image_dir = os.path.join(project_dir, 'data', 'images', 'he', 'aligned_images', 'visium')
image_file = f'{sample_id}.ome.tif'
image_file_path = os.path.join(image_dir, image_file)

# Load the original high-res scan using tifffile
img = tifffile.imread(image_file_path)

# Load the original high-res scan
img = tifffile.imread(image_file_path)

# Define key for adding to anndata
original_hires_key = 'original_hires'

# Calculate the scale factor
# Divide the height of the original hires image by the height of the downsampled hires image
size_fact = img.shape[0] / adata.uns['spatial'][sample_id]['images']['hires'].shape[0]

# Extract the downsampled hires scale factor
hires_fact = adata.uns['spatial'][sample_id]['scalefactors']['tissue_hires_scalef']

# Calculate the scale factor for the original hires image
fulres_fact = size_fact * hires_fact

# Add to the anndata object
adata.uns['spatial'][sample_id]['scalefactors'][f'tissue_{original_hires_key}_scalef'] = fulres_fact

# --- Calculate microns per pixel for scale bar

spot_diameter_fulres = adata.uns['spatial'][sample_id]['scalefactors']['spot_diameter_fullres']
spot_diameter = 55  # in microns
microns_per_pixel = spot_diameter / spot_diameter_fulres

# Extract the hires factor
hires_factor = adata.uns['spatial'][sample_id]['scalefactors']['tissue_hires_scalef']

# --- Define crop areas for zoomed views

# Coordinates from original script (in full resolution pixels)
crop_coords = {
    'Primary follicle': [30862, 30900],
    'Secondary follicle': [35077, 18986],
}

# Width to height ratio for crop areas (1 = square, 2 = width is 2x height, etc.)
width_to_height_ratio = 1.6

# Base dimension (used as height)
height = 3612
width = int(height * width_to_height_ratio)

# Calculate offset to center the crop
# If ratio != 1, we need to adjust x coordinate to keep crop centered
width_offset = int((width - height) / 2)

height_hires = height * hires_factor
width_hires = width * hires_factor

# Adjust the overview coords based on width_to_height_ratio
# coords format is (xmin, xmax, ymin, ymax)
if width_to_height_ratio != 1:
    # Calculate current dimensions
    coords_width = coords[1] - coords[0]
    coords_height = coords[3] - coords[2]
    # Calculate center point
    center_x = (coords[0] + coords[1]) / 2
    center_y = (coords[2] + coords[3]) / 2
    # Use the smaller dimension as height, then calculate new width
    new_height = min(coords_width, coords_height)
    new_width = new_height * width_to_height_ratio
    # Recalculate coords centered on the same point
    coords = (
        int(center_x - new_width / 2),
        int(center_x + new_width / 2),
        int(center_y - new_height / 2),
        int(center_y + new_height / 2)
    )

# --- Helper functions for plotting

def move_plots_x_y(axs, x_offset=0.02, y_offset=0.02, remove_title=True):
    """Move plots so they don't overlap, moving both in x & y direction"""
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

def crop_image(img, coords):
    """Crop the image with safe coordinates"""
    # Clip the coordinates to ensure they are within the image boundaries
    # Convert to int to ensure they can be used as slice indices
    x_min = int(max(coords[0], 0))
    x_max = int(min(coords[1], img.shape[1]))
    y_min = int(max(coords[2], 0))
    y_max = int(min(coords[3], img.shape[0]))
    # Crop the image using the clipped coordinates
    crop_img = img[y_min:y_max, x_min:x_max]
    return crop_img, (x_min, x_max, y_min, y_max)

def crop_and_subset_anndata(
    adata,
    adata_crop,
    img,
    sample_id,
    original_hires_key,
    zoom_coord_fulres,
    width_fulres,
    height_fulres,
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
        zoom_coord_fulres[1] + height_fulres
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
    # Extract spatial coordinates
    adata_crop.obs["pxl_col_in_fullres"] = adata_crop.obsm["spatial"][:, 0] * fulres_fact
    adata_crop.obs["pxl_row_in_fullres"] = adata_crop.obsm["spatial"][:, 1] * fulres_fact
    # Subset the AnnData object based on the crop region
    mask = (
        (adata_crop.obs["pxl_col_in_fullres"] >= xmin) & 
        (adata_crop.obs["pxl_col_in_fullres"] < xmax) &
        (adata_crop.obs["pxl_row_in_fullres"] >= ymin) & 
        (adata_crop.obs["pxl_row_in_fullres"] < ymax)
    )
    adata_sub = adata_crop[mask]
    return adata_sub, crop_img, (xmin, ymin, width_crop, height_crop)

# --- Create the figure

ncols = 3
nrows = 1 + len(crop_coords)

# Set the text size
sc.set_figure_params(fontsize=30)

# Set the scalebar size
scalebar_val = 1000  # in microns
scalebar_units = 'µm'

# Make figure
fig, axs = plt.subplots(nrows, ncols, figsize=(5, 5))

# --- Row 1: Overview plots

# Panel 1: H&E
ax = axs[0][0]
sc.pl.spatial(adata, color=None, size=spot_size, ax=ax, bw=False, alpha_img=1,
             crop_coord=coords, title='H&E', show=False)
ax.set_title(ax.get_title(), fontweight='bold')

scalebar = ScaleBar(microns_per_pixel/hires_factor, units=scalebar_units, location='lower right', 
                    box_alpha=0.5, fixed_value=scalebar_val, frameon=True, scale_loc="none", 
                    pad=0.3, height_fraction=0.02, font_properties={'size': 20})
scalebar.label = '1 mm'
ax.add_artist(scalebar)

# Add rectangles to mark the zoom areas
for area, crop_coord in crop_coords.items():
    # Adjust the coordinates to center the crop when width != height
    # Original crop_coord is top-left for a square, we adjust x to keep center the same
    adjusted_x = crop_coord[0] - width_offset
    crop_x = adjusted_x * hires_factor
    crop_y = crop_coord[1] * hires_factor
    rect = mpl.patches.Rectangle((crop_x, crop_y), width_hires, height_hires, 
                                  linewidth=2, edgecolor='black', facecolor='none')
    ax.add_patch(rect)

# Panel 2: GC annotations
ax = axs[0][1]
sc.pl.spatial(adata, color=None, size=spot_size, ax=ax, bw=False, alpha_img=1,
             crop_coord=coords, title='GC annotations', show=False)
ax.set_title(ax.get_title(), fontweight='bold')

scalebar = ScaleBar(microns_per_pixel/hires_factor, units=scalebar_units, location='lower right', 
                    box_alpha=0.5, fixed_value=scalebar_val, frameon=False, scale_loc="none", 
                    pad=0.3, height_fraction=0.02)
ax.add_artist(scalebar)

# Plot each polygon if GC annotations are available
if gdf is not None:
    for poly in gdf.geometry:
        if poly.geom_type == 'Polygon':
            x, y = poly.exterior.xy
            # Multiply by hires_factor to adjust for the image resolution
            x = [coord * hires_factor for coord in x]
            y = [coord * hires_factor for coord in y]
            # Plot the polygon
            coords_pol = list(zip(x, y))
            # Add filled polygon with edge fully opaque
            patch = mpl.patches.Polygon(
                coords_pol,
                closed=True,
                facecolor='cyan',
                edgecolor='cyan',
                linewidth=0.5,
            )
            facecolor = patch.get_facecolor()
            patch.set_facecolor((facecolor[0], facecolor[1], facecolor[2], 0.3))  # Semi-transparent fill
            ax.add_patch(patch)

# Add rectangles to mark the zoom areas
for area, crop_coord in crop_coords.items():
    adjusted_x = crop_coord[0] - width_offset
    crop_x = adjusted_x * hires_factor
    crop_y = crop_coord[1] * hires_factor
    rect = mpl.patches.Rectangle((crop_x, crop_y), width_hires, height_hires, 
                                  linewidth=2, edgecolor='black', facecolor='none')
    ax.add_patch(rect)

# Panel 3: Niche signatures (overlaid)
ax = axs[0][-1]

for n, score_name in enumerate(score_names_r):
    if score_name in score_keys:
        if n == 0:
            alpha_img = 0.7
        else:
            alpha_img = 0
        score_key = score_keys[score_name]
        cmap = cmap_dict[score_name]
        sc.pl.spatial(adata, color=score_key, size=spot_size, ax=ax, bw=False, alpha_img=alpha_img, cmap=cmap,
                     crop_coord=coords, title='Niche signatures', vmin=vmin, vmax=vmax, colorbar_loc=None)

ax.set_title(ax.get_title(), fontweight='bold')

# Add scalebar
scalebar = ScaleBar(microns_per_pixel/hires_factor, units=scalebar_units, location='lower right', 
                    box_alpha=0.5, fixed_value=scalebar_val, frameon=False, scale_loc="none", 
                    pad=0.3, height_fraction=0.02)
ax.add_artist(scalebar)

# Add rectangles to mark the zoom areas
for area, crop_coord in crop_coords.items():
    adjusted_x = crop_coord[0] - width_offset
    crop_x = adjusted_x * hires_factor
    crop_y = crop_coord[1] * hires_factor
    rect = mpl.patches.Rectangle((crop_x, crop_y), width_hires, height_hires, 
                                  linewidth=2, edgecolor='black', facecolor='none')
    ax.add_patch(rect)

# --- Rows 2+: Zoomed areas

for i, (area, crop_coord) in enumerate(crop_coords.items()):
    # Create fresh copy of adata for cropping
    adata_crop = adata.copy()
    
    # Panel 1: H&E zoomed
    ax = axs[i + 1][0]
    # Adjust crop coordinate to center when width != height
    adjusted_crop_coord = [crop_coord[0] - width_offset, crop_coord[1]]
    # Crop the H&E image
    adata_sub, crop_img, (xmin, ymin, width_crop, height_crop) = crop_and_subset_anndata(
        adata=adata,
        adata_crop=adata_crop,
        img=img,
        sample_id=sample_id,
        original_hires_key=original_hires_key,
        zoom_coord_fulres=adjusted_crop_coord,
        width_fulres=width,
        height_fulres=height,
        fulres_fact=fulres_fact
    )
    sc.pl.spatial(adata_sub, color=None, ax=ax, bw=False, alpha_img=1,
                 title='', img_key=original_hires_key, crop_coord=(0, width_crop, 0, height_crop))
    ax.set_title(ax.get_title(), fontweight='bold')
    scalebar = ScaleBar(microns_per_pixel, units=scalebar_units, location='lower right', box_alpha=0.5, 
                        fixed_value=100, frameon=True, scale_loc="top", pad=0.3, height_fraction=0.02, 
                        font_properties={'size': 20})
    ax.add_artist(scalebar)
    
    # Panel 2: GC annotations zoomed
    ax = axs[i + 1][1]
    sc.pl.spatial(adata_sub, color=None, ax=ax, bw=False, alpha_img=1,
                 title='', img_key=original_hires_key, crop_coord=(0, width_crop, 0, height_crop))
    ax.set_title(ax.get_title(), fontweight='bold')
    
    # Plot each polygon (if available)
    if gdf is not None:
        for poly in gdf.geometry:
            if poly.geom_type == 'Polygon':
                x, y = poly.exterior.xy
                # Transform coordinates to match the cropped image
                x = [coord - adjusted_crop_coord[0] for coord in x]
                y = [coord - adjusted_crop_coord[1] for coord in y]
                # Check whether any of the coordinates are outside the crop area
                if (max(x) < xmin or min(x) > xmin + width_crop or
                    max(y) < ymin or min(y) > ymin + height_crop):
                    continue
                # Plot the polygon
                coords_pol = list(zip(x, y))
                # Add filled polygon
                patch = mpl.patches.Polygon(
                    coords_pol,
                    closed=True,
                    facecolor='cyan',
                    edgecolor='cyan',
                    linewidth=3,
                )
                facecolor = patch.get_facecolor()
                patch.set_facecolor((facecolor[0], facecolor[1], facecolor[2], 0.2))  # Semi-transparent fill
                ax.add_patch(patch)
    
    scalebar = ScaleBar(microns_per_pixel, units=scalebar_units, location='lower right', box_alpha=0.5, 
                        fixed_value=100, frameon=False, scale_loc="none", pad=0.3, height_fraction=0.02)
    ax.add_artist(scalebar)
    
    # Panel 3: Niche signatures zoomed
    ax = axs[i + 1][-1]
    for n, score_name in enumerate(score_names_r):
        if score_name in score_keys:
            if n == 0:
                alpha_img = 0.7
            else:
                alpha_img = 0
            score_key = score_keys[score_name]
            cmap = cmap_dict[score_name]
            sc.pl.spatial(adata_sub, color=score_key, ax=ax, bw=False, alpha_img=alpha_img, cmap=cmap,
                         crop_coord=(0, width_crop, 0, height_crop), title='', vmin=vmin, vmax=vmax,
                         colorbar_loc=None, img_key=original_hires_key)
    
    scalebar = ScaleBar(microns_per_pixel, units=scalebar_units, location='lower right', box_alpha=0.5, 
                        fixed_value=100, frameon=False, scale_loc="none", pad=0.3, height_fraction=0.02)
    ax.add_artist(scalebar)
    ax.set_title(ax.get_title(), fontweight='bold')

# --- Format the figure

# Remove the x & y labels of all axs
for ax in axs.flat:
    ax.set_xlabel('')
    ax.set_ylabel('')

# Resize the plots
for ax in axs.flat:
    ax_pos = ax.get_position()
    # Set the width to be equal for all plots
    ax.set_position([ax_pos.x0, ax_pos.y0, 0.8, 0.8])

# Set area names as row labels
axs[0, 0].set_ylabel('Lymph node', fontweight='bold', labelpad=10)

for i, area in enumerate(crop_coords.keys()):
    axs[i+1, 0].set_ylabel(f'{area}', fontweight='bold', labelpad=10)

move_plots_x_y(axs, x_offset=0.05, y_offset=0.05, remove_title=True)

# --- Add colorbars

# Dict to update the score names for display
score_names_alt = {
    'T_cell_zone': 'T cell zone',
    'B_cell_zone': 'B cell zone',
    'GC': 'GC'
}

# Add colorbars for each score
target_ax = axs[-1, -1]  # Use the last subplot for colorbars

for i, (label, cmap) in enumerate(cmap_dict.items()):
    norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
    spacing = 0.1  # vertical space between colorbars
    start_y = -0.1  # starting y position (below axis)
    y_pos = start_y - i * spacing  # stack downward
    cbar_width = 0.4  # width of the colorbar
    cbar_height = 0.05
    x_pos = 0.5 - cbar_width / 2  # center the colorbar horizontally
    cbar_ax = target_ax.inset_axes([x_pos, y_pos, cbar_width, cbar_height])
    cbar = fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap),
                        cax=cbar_ax, orientation='horizontal')
    # Set the ylabel as title
    cbar_ax.set_ylabel(f'{score_names_alt[label]}', fontsize=23, labelpad=10, 
                       rotation=0, va='center', ha='right')
    cbar.ax.tick_params(labelsize=20)
    # Limit to 2 xticks
    cbar.ax.set_xticks([0, vmax])
    cbar.ax.set_xticklabels([0, vmax])
    # Remove the xticks if not the bottom colorbar
    if i < len(cmap_dict) - 1:
        cbar.ax.set_xticklabels([])

# --- Add title

# Calculate the suptitle position
title_x = (axs[0][0].get_position().x0 + axs[0][-1].get_position().x1) / 2
title_y = axs[0][0].get_position().y1 + 0.3

plt.suptitle(f'{sample_id} - GC annotations + TLS niche signatures',
             y=title_y, x=title_x, fontweight='bold')

# --- Save the plot

# Create output directory
out_dir = os.path.join(res_dir, 'figures')
os.makedirs(out_dir, exist_ok=True)

filename = 'Fig_2I.pdf'

plt.savefig(os.path.join(out_dir, filename), bbox_inches='tight')
plt.close('all')

print(f"Figure saved to {os.path.join(out_dir, filename)}")
print("Script successfully completed")
