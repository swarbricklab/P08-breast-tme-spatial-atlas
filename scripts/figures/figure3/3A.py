#!/usr/bin/env python3

# --- Aim
# Plot HEV signature at varying resolutions (Figure 3A)
# Shows HEV signature scores across different bin sizes with zoomed-in panels
# This script is set up for P29_HD

# --- Define environment
# env: vis-hd

# --- Load packages
import argparse
import os
import sys

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patheffects as path_effects
from matplotlib.colors import to_rgb, LinearSegmentedColormap

import numpy as np
import pandas as pd
import scanpy as sc
import tifffile
from matplotlib_scalebar.scalebar import ScaleBar

mpl.rcParams['pdf.fonttype'] = 42

sc.set_figure_params(facecolor="white", figsize=(8, 8))
sc.settings.verbosity = 3


# --- Argparse arguments
parser = argparse.ArgumentParser(description="Plot HEV signature at varying resolutions for Figure 3A",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument("--res_dir", help="Path to results directory")
parser.add_argument("--h5ad_dir", help="Path to directory for .h5ad output files")
parser.add_argument("--project_dir", help="Project directory (root)")
parser.add_argument("--sample_id", default="P29_HD", help="Sample shown in the figure")
parser.add_argument("--sample_sheet", help="Path to sample sheet CSV file", default=None)

args = vars(parser.parse_args())

# Arguments
res_dir = args["res_dir"]
h5ad_dir = args["h5ad_dir"]
sample_id = args["sample_id"]
project_dir = args["project_dir"]

print(f'Sample: {sample_id}')

sys.path.append(os.path.join(project_dir, 'scripts'))
from helper_functions import auto_crop_scanpy


# --- Extract some info from the sample sheet
# Use sample_sheet argument if provided, else default
if args.get("sample_sheet"):
    sample_sheet_path = args["sample_sheet"]
else:
    sample_sheet_path = os.path.join(project_dir, 'config', 'sample_sheet_visium_hd.csv')

sample_sheet = pd.read_csv(sample_sheet_path)
tissue_id = sample_sheet[sample_sheet['section_id'] == sample_id]['tissue_id'].values[0]


# --- Load the adata objects

# Run for all resolutions
resolutions = ['square_048um', 'square_016um', 'square_008um']

adata_dict = {}

for resolution in resolutions:
    adata_dict[resolution] = sc.read_h5ad(os.path.join(h5ad_dir, 'visium_hd', resolution, f'{sample_id}.h5ad'))


# --- Calculate crop coords for optimal presentation

coords = {}
for resolution in resolutions:
    coords[resolution] = auto_crop_scanpy(adata_dict[resolution], border_size=0.1)


# --- Load HEV scoring results

HEV_score_name = 'HEVs'
out_name = f'{HEV_score_name}_score'
top_genes = 20

for resolution in resolutions:
    HEV_res_dir = os.path.join(res_dir, '04_vascular_microenvironment', '01_identify_HEV_niches', '01_score_signature', out_name, resolution)
    most_recent_dir_path = HEV_res_dir
    
    # Read results and add to adata
    HEV_file_name = f'{sample_id}_{out_name}.csv'
    HEV_results = pd.read_csv(os.path.join(most_recent_dir_path, 'csv_files', HEV_file_name), index_col=0)
    HEV_key = f'{HEV_score_name}_score_{top_genes}'
    
    # Check if the key exists; otherwise, find the next highest available key
    available_keys = list(HEV_results.columns)
    if HEV_key not in available_keys:
        possible_keys = [k for k in available_keys if k.startswith(f'{HEV_score_name}_score_')]
        possible_numbers = sorted(
            [int(k.split('_')[-1]) for k in possible_keys if k.split('_')[-1].isdigit()],
            reverse=True
        )
        # Find the next highest available number
        top_genes = next((n for n in possible_numbers if n < top_genes), possible_numbers[-1] if possible_numbers else top_genes)
        HEV_key = f'{HEV_score_name}_score_{top_genes}'
    
    adata_dict[resolution].obs[HEV_key] = HEV_results[HEV_key]


# --- Import the original high-resolution scan

resolution = 'square_048um'
adata = adata_dict[resolution]

# Find the original high-res image
image_dir = os.path.join(project_dir, 'data', 'images', 'he', 'aligned_images', 'visium_hd')
image_file = f'{sample_id}.ome.tif'
image_file_path = os.path.join(image_dir, image_file)

# Load the original high-res scan using tifffile
img = tifffile.imread(image_file_path)

# Add the scan to the anndata object for easy plotting
original_hires_key = 'original_hires'

# Calculate the scale factor
# Divide the height of the original hires image by the height of the downsampled hires image
size_fact = img.shape[0] / adata.uns['spatial'][sample_id]['images']['hires'].shape[0]

# Extract the downsampled hires scale factor
hires_scalef = adata.uns['spatial'][sample_id]['scalefactors']['tissue_hires_scalef']

# Calculate the scale factor for the original hires image
fulres_fact = size_fact * hires_scalef

# Get the microns per pixel for the scalebar
microns_per_pixel = adata.uns['spatial'][sample_id]['scalefactors']['microns_per_pixel']


# --- Create colormap for HEV score

def create_colormap(R, G, B):
    N = 256
    vals = np.ones((N, 4))
    vals[:, 0] = R / 256
    vals[:, 1] = G / 256
    vals[:, 2] = B / 256
    # Create gradient with increasing transparency
    vals[:, 3] = np.linspace(0, 1, N)
    cmap = LinearSegmentedColormap.from_list("custom", vals)
    cmap.set_under(color=(1, 1, 1, 0))
    return cmap

HEV_color = mpl.colors.cnames["indigo"]
rgb = tuple(int(c * 255) for c in to_rgb(HEV_color))
cmap = create_colormap(rgb[0], rgb[1], rgb[2])


# --- Function to crop and subset anndata

def crop_and_subset_anndata(adata, adata_crop, img, sample_id, original_hires_key,
                            zoom_coord_fulres, width_fulres, fulres_fact):
    """Crop image and subset adata to region of interest"""
    # Calculate crop coordinates for the image (in pixel space)
    coords_group = [
        zoom_coord_fulres[0],
        zoom_coord_fulres[0] + width_fulres,
        zoom_coord_fulres[1],
        zoom_coord_fulres[1] + width_fulres
    ]
    
    # Clip coordinates to image boundaries
    x_min = max(coords_group[0], 0)
    x_max = min(coords_group[1], img.shape[1])
    y_min = max(coords_group[2], 0)
    y_max = min(coords_group[3], img.shape[0])
    
    # Crop the image
    crop_img = img[y_min:y_max, x_min:x_max]
    
    # Calculate the height and width of the cropped image
    height_crop = y_max - y_min
    width_crop = x_max - x_min
    
    # Add cropped image to anndata
    adata_crop.uns['spatial'][sample_id]['images'][original_hires_key] = crop_img
    adata_crop.uns['spatial'][sample_id]['scalefactors'][f'tissue_{original_hires_key}_scalef'] = fulres_fact
    
    # Shift the spot coordinates
    spot_coords = adata.obsm['spatial']
    shifted_coords = spot_coords - np.array([x_min, y_min])
    adata_crop.obsm['spatial'] = shifted_coords
    
    # Extract spatial coordinates in pixel space for filtering
    adata_crop.obs["pxl_col_in_fullres"] = adata_crop.obsm["spatial"][:, 0] * fulres_fact
    adata_crop.obs["pxl_row_in_fullres"] = adata_crop.obsm["spatial"][:, 1] * fulres_fact
    
    # Subset to spots within the crop (with 10% grace region)
    grace_region = 0.1
    xmin_grace = -int(width_crop * grace_region)
    xmax_grace = width_crop + int(width_crop * grace_region)
    ymin_grace = -int(height_crop * grace_region)
    ymax_grace = height_crop + int(height_crop * grace_region)
    
    mask = (
        (adata_crop.obs["pxl_col_in_fullres"] >= xmin_grace) & 
        (adata_crop.obs["pxl_col_in_fullres"] < xmax_grace) &
        (adata_crop.obs["pxl_row_in_fullres"] >= ymin_grace) & 
        (adata_crop.obs["pxl_row_in_fullres"] < ymax_grace)
    )
    
    adata_sub = adata_crop[mask].copy()
    
    return adata_sub, crop_img, (0, 0, width_crop, height_crop)


# --- Function to move plots

def move_plots_x_y(axs, x_offset=0.02, y_offset=0.02, remove_title=True):
    """Move plots to avoid overlap"""
    for i, axs_row in enumerate(axs):
        for j, ax in enumerate(axs_row):
            if i != 0:
                if remove_title:
                    ax.set_title('')
            # Adjust the positions of the subplots
            ax_pos = ax.get_position()
            # Shift the subplot to the desired location
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


# --- Define zoom areas

# Main zoom area coordinates (full resolution pixels)
zoom_coord_fulres = [37305, 54452]  # Top left corner
width_fulres = 8745  # Width in pixels

# Convert to hires coords for rectangle
zoom_coord_hires = [int(zoom_coord_fulres[0] * hires_scalef), 
                    int(zoom_coord_fulres[1] * hires_scalef)]
width_hires = int(width_fulres * hires_scalef)

# Sub-crop areas within the main zoom
subcrop_coords = {
    1: [37889, 56494],
    2: [42877, 58344],
}

width_fulres_sub = 2500

# Individual HEV coordinates
HEV_coords = {
    1: [38486, 58182],
    2: [44288, 58553],
}

width_HEV = 696

# Create crop dict
adata_crop_dict = {}
for resolution in resolutions:
    adata_crop_dict[resolution] = adata_dict[resolution].copy()


# --- Prepare the figure

# Set parameters
sc.set_figure_params(fontsize=20)
scalebar_val = 300  # microns
scalebar_units = 'µm'

vmin = 0
vmax = 0.8

ncols = 1 + len(resolutions)
nrows = 3

fig, axs = plt.subplots(nrows, ncols, figsize=(5, 5))

# --- First plot H&E zoomed in area

ax = axs[0, 0]

resolution = 'square_048um'
adata = adata_dict[resolution]
adata_crop = adata_crop_dict[resolution]

adata_sub, crop_img, (xmin, ymin, width_crop, height_crop) = crop_and_subset_anndata(
    adata, adata_crop, img, sample_id, original_hires_key,
    zoom_coord_fulres, width_fulres, fulres_fact
)

sc.pl.spatial(adata_sub, color=None, bw=False, alpha_img=1, 
              crop_coord=(0, width_crop, 0, height_crop), 
              title='H&E', na_in_legend=False, ax=ax, img_key=original_hires_key)
ax.set_title(ax.get_title(), fontweight='bold')

# Add rectangles to mark the sub-zoom areas
for area, coords in subcrop_coords.items():
    crop_x = coords[0] - zoom_coord_fulres[0]
    crop_y = coords[1] - zoom_coord_fulres[1]
    rect = mpl.patches.Rectangle((crop_x, crop_y), width_fulres_sub, width_fulres_sub, 
                                 linewidth=2, edgecolor='black', facecolor='none')
    ax.add_patch(rect)
    # Annotate the area
    ax.text(crop_x + width_fulres_sub / 2, crop_y - width_fulres_sub * 0.05, area, 
            color='black', fontsize=20, ha='center', va='bottom', fontweight='bold',
            path_effects=[path_effects.Stroke(linewidth=2, foreground='white'), 
                         path_effects.Normal()])

# Add scalebar
scalebar = ScaleBar(microns_per_pixel, units=scalebar_units, location='lower right', 
                   box_alpha=0.5, fixed_value=scalebar_val, frameon=True, 
                   scale_loc="top", pad=0.3, height_fraction=0.02)
ax.add_artist(scalebar)


# --- Plot HEV signature at varying resolutions

for n, resolution in enumerate(resolutions):
    print(f'Plotting resolution {resolution} ({n+1}/{len(resolutions)})')
    adata = adata_dict[resolution]
    adata_crop = adata_crop_dict[resolution]
    
    adata_sub, crop_img, (xmin, ymin, width_crop, height_crop) = crop_and_subset_anndata(
        adata, adata_crop, img, sample_id, original_hires_key,
        zoom_coord_fulres, width_fulres, fulres_fact
    )
    
    ax = axs[0, n + 1]
    
    sc.pl.spatial(adata_sub, color=f'{HEV_score_name}_score_{top_genes}', 
                 ax=ax, bw=False, alpha_img=0.7, cmap=cmap,
                 crop_coord=(0, width_crop, 0, height_crop), 
                 title='', vmin=vmin, vmax=vmax, colorbar_loc=None,
                 img_key=original_hires_key)
    
    # Set resolution label
    resolution_label = resolution.replace('square_', '').replace('um', 'μm').lstrip('0')
    resolution_label = f'{resolution_label} bins'
    ax.set_title(resolution_label, fontweight='bold')
    
    # Add scalebar
    scalebar = ScaleBar(microns_per_pixel, units=scalebar_units, 
                       location='lower right', box_alpha=0.5, 
                       fixed_value=scalebar_val, frameon=False, 
                       scale_loc='none', pad=0.3, height_fraction=0.02,
                       font_properties={'weight': 'bold'})
    ax.add_artist(scalebar)
    
    # Add rectangles to mark zoom areas
    for area, coords in subcrop_coords.items():
        crop_x = coords[0] - zoom_coord_fulres[0]
        crop_y = coords[1] - zoom_coord_fulres[1]
        rect = mpl.patches.Rectangle((crop_x, crop_y), width_fulres_sub, width_fulres_sub,
                                     linewidth=2, edgecolor='black', facecolor='none')
        ax.add_patch(rect)
        ax.text(crop_x + width_fulres_sub / 2, crop_y - width_fulres_sub * 0.05, area,
                color='black', fontsize=20, ha='center', va='bottom', fontweight='bold',
                path_effects=[path_effects.Stroke(linewidth=2, foreground='white'),
                             path_effects.Normal()])


# --- Plot HEV signature at zoomed-in areas

for i, (area, coords) in enumerate(subcrop_coords.items()):
    print(f'Plotting HEV signature for area {area}')
    
    # First plot H&E
    resolution = 'square_048um'
    adata = adata_dict[resolution]
    adata_crop = adata_crop_dict[resolution]
    
    adata_sub, crop_img, (xmin, ymin, width_crop, height_crop) = crop_and_subset_anndata(
        adata, adata_crop, img, sample_id, original_hires_key,
        coords, width_fulres_sub, fulres_fact
    )
    
    ax = axs[i + 1, 0]
    sc.pl.spatial(adata_sub, color=None, bw=False, alpha_img=1,
                 crop_coord=(0, width_crop, 0, height_crop),
                 title='', na_in_legend=False, ax=ax, img_key=original_hires_key)
    
    # Add scalebar
    scalebar = ScaleBar(microns_per_pixel, units=scalebar_units,
                       location='lower right', box_alpha=0.5, fixed_value=100,
                       frameon=True, scale_loc='top', pad=0.3, height_fraction=0.02)
    ax.add_artist(scalebar)
    
    # Then plot HEV signature at each resolution
    for n, resolution in enumerate(resolutions):
        adata = adata_dict[resolution]
        adata_crop = adata_crop_dict[resolution]
        
        adata_sub, crop_img, (xmin, ymin, width_crop, height_crop) = crop_and_subset_anndata(
            adata, adata_crop, img, sample_id, original_hires_key,
            coords, width_fulres_sub, fulres_fact
        )
        
        ax = axs[i + 1, n + 1]
        sc.pl.spatial(adata_sub, color=f'{HEV_score_name}_score_{top_genes}',
                     ax=ax, bw=False, alpha_img=0.7, cmap=cmap,
                     crop_coord=(0, width_crop, 0, height_crop),
                     title='', vmin=vmin, vmax=vmax, colorbar_loc=None,
                     img_key=original_hires_key)
        
        # Add scalebar
        scalebar = ScaleBar(microns_per_pixel, units=scalebar_units,
                           location='lower right', box_alpha=0.5, fixed_value=100,
                           frameon=False, scale_loc='none', pad=0.3, height_fraction=0.02)
        ax.add_artist(scalebar)


# --- Remove labels and resize plots

for ax in axs.flat:
    ax.set_xlabel('')
    ax.set_ylabel('')

for ax in axs.flat:
    ax_pos = ax.get_position()
    ax.set_position([ax_pos.x0, ax_pos.y0, 0.8, 0.8])

move_plots_x_y(axs, x_offset=0.05, y_offset=0.05, remove_title=True)


# --- Add colorbar

ax = axs[0, -1]
norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
cbar_ax = ax.inset_axes([1.1, 0.5, 0.3, 0.05])
cbar = fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap),
                   cax=cbar_ax, orientation='horizontal')
cbar.ax.set_ylabel('HEVs', rotation=0, ha='left', va='center', fontsize=18)
cbar.ax.yaxis.set_label_position("right")
cbar.ax.yaxis.set_label_coords(1.05, 1.05)
cbar.ax.tick_params(labelsize=15)


# --- Set ylabels for rows

for i in range(nrows):
    if i > 0:
        axs[i, 0].set_ylabel(f'{i}', fontweight='bold', rotation=0, labelpad=20)


# --- Add extra zoomed HEV panels

for i, (area, coords) in enumerate(HEV_coords.items()):
    print(f'Plotting HEV {area}')
    new_ax = fig.add_subplot(nrows, ncols, (nrows - 1) * ncols + i + 2)
    
    # Get target ax for positioning
    ax = axs[i + 1, -1]
    new_ax.set_position([ax.get_position().x1 + 0.05, ax.get_position().y0, 0.8, 0.8])
    
    # Add rectangle on previous ax
    target_crop_coords = subcrop_coords[area]
    HEV_x = coords[0] - target_crop_coords[0]
    HEV_y = coords[1] - target_crop_coords[1]
    rect = mpl.patches.Rectangle((HEV_x, HEV_y), width_HEV, width_HEV,
                                 linewidth=2, edgecolor='black', facecolor='none')
    ax.add_patch(rect)
    
    # Crop and plot HEV
    resolution = 'square_048um'
    adata = adata_dict[resolution]
    adata_crop = adata_crop_dict[resolution]
    
    adata_sub, crop_img, (xmin, ymin, width_crop, height_crop) = crop_and_subset_anndata(
        adata, adata_crop, img, sample_id, original_hires_key,
        coords, width_HEV, fulres_fact
    )
    
    sc.pl.spatial(adata_sub, color=None, bw=False, alpha_img=1,
                 crop_coord=(0, width_crop, 0, height_crop),
                 title='', na_in_legend=False, ax=new_ax, img_key=original_hires_key)
    
    if i == 0:
        new_ax.set_title('HEVs', fontweight='bold')
    
    new_ax.set_xlabel('')
    new_ax.set_ylabel('')
    
    # Add scalebar
    scalebar = ScaleBar(microns_per_pixel, units=scalebar_units,
                       location='lower right', box_alpha=0.5, fixed_value=20,
                       frameon=True, scale_loc='top', pad=0.3, height_fraction=0.02)
    new_ax.add_artist(scalebar)
    
    # Draw connector polygon
    trans_main = ax.transData.transform
    trans_zoom = new_ax.transAxes.transform
    inv = fig.transFigure.inverted().transform
    
    p1 = (HEV_x + width_HEV, HEV_y + width_HEV)  # top-right
    p2 = (HEV_x + width_HEV, HEV_y)  # bottom-right
    q1 = (0, 0)  # bottom-left of zoom
    q2 = (0, 1)  # top-left of zoom
    
    p1_fig = inv(trans_main(p1))
    p2_fig = inv(trans_main(p2))
    q1_fig = inv(trans_zoom(q1))
    q2_fig = inv(trans_zoom(q2))
    
    connector = mpl.patches.Polygon(
        [p1_fig, p2_fig, q2_fig, q1_fig],
        closed=True, edgecolor='black', facecolor=(1, 1, 1, 0.4),
        transform=fig.transFigure, linewidth=1
    )
    fig.patches.append(connector)


# --- Add header label

start_box = axs[0, 1].get_position()
end_box = axs[0, -1].get_position()

x0 = start_box.x0
x1 = end_box.x1
width = x1 - x0
y0 = start_box.y1 + 0.09
height = 0.015

highlight_rect = mpl.patches.Rectangle(
    (x0, y0), width, height,
    transform=fig.transFigure,
    facecolor='black', edgecolor='black',
    linewidth=0, alpha=1
)
fig.patches.append(highlight_rect)

x_mid = (start_box.x0 + end_box.x1) / 2
y_text = y0 + height + 0.01

fig.text(x_mid, y_text, "HEV signature", ha='center', va='bottom',
         fontsize=20, color='black', fontweight='bold')


# --- Add title

title_x = (axs[0][0].get_position().x0 + axs[0][-1].get_position().x1) / 2
title_y = axs[0][0].get_position().y1 + 0.3

plt.suptitle(f'{sample_id} - HEV signature at various resolutions (scalebar={scalebar_val}{scalebar_units}, top_genes={top_genes})',
             y=title_y, x=title_x, fontweight='bold')


# --- Save plot

out_dir = os.path.join(res_dir, 'figures')
os.makedirs(out_dir, exist_ok=True)

filename = 'Fig_3A.pdf'
plt.savefig(os.path.join(out_dir, filename), bbox_inches='tight')

print(f"Saved plot to: {os.path.join(out_dir, filename)}")
print("Script successfully completed")
