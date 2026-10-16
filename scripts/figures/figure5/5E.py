#!/usr/bin/env python3

# --- Aim
# Plot Fig. 5e: spatial transcript distribution of SPP1 in HEV regions

# --- Define environment
# env: spatdat

# --- Load packages

import argparse
import os
import shutil
import json

import anndata as ad
import scanpy as sc
import pandas as pd

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.lines import Line2D
from matplotlib_scalebar.scalebar import ScaleBar

import geopandas as gpd
import PIL
import spatialdata as sd
import tifffile

mpl.rcParams['pdf.fonttype'] = 42

sc.set_figure_params(figsize=(8, 8), fontsize=20)


# --- Argparse arguments

parser = argparse.ArgumentParser(description="Fig. 5e: SPP1 transcripts in HEV regions",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

# Required
parser.add_argument("--h5ad_dir", help="Path to directory for .h5ad output files")
parser.add_argument("--sample_id", default="P29_Xe5k", help="Sample shown in the figure")
parser.add_argument("--project_dir", help="Path to the root of the project directory")
parser.add_argument("--res_dir", help="Path to results directory")
parser.add_argument("--sample_sheet", help="Path to sample sheet CSV file", default=None)

args = vars(parser.parse_args())


# --- Arguments

h5ad_dir=args["h5ad_dir"]
sample_id=args["sample_id"]

res_dir = args['res_dir']
project_dir = args['project_dir']

# --- Extract some info from the sample sheet
# Use sample_sheet argument if provided, else default
if args.get("sample_sheet"):
    sample_sheet_path = args["sample_sheet"]
else:
    sample_sheet_path = os.path.join(project_dir, 'config', 'sample_sheet_xenium.csv')

sample_sheet = pd.read_csv(sample_sheet_path)

tissue_id = sample_sheet[sample_sheet['section_id'] == sample_id]['tissue_id'].values[0]

# Also load tissue metadata
metadata_dir = os.path.join(project_dir, 'resources',)
tissue_metadata_file = 'tissue_metadata.csv'
tissue_metadata = pd.read_csv(os.path.join(metadata_dir, tissue_metadata_file))

subtype = tissue_metadata[tissue_metadata['tissue_id'] == str(tissue_id)]['clinical_subtype'].values[0]


# --- Load data

adata = ad.read_h5ad(os.path.join(h5ad_dir, 'xenium', f'{sample_id}.h5ad'))


# --- Import the spatialdata object

spat_dat_dir = os.path.join(project_dir, 'data', 'spatialdata_objects', 'xenium')

filename = f'{sample_id}.zarr'

# Might need to unzip first
if not os.path.exists(os.path.join(spat_dat_dir, filename)):
    print(f"Unzipping {filename}")
    shutil.unpack_archive(os.path.join(spat_dat_dir, filename.replace('zarr', 'zip')), os.path.join(spat_dat_dir, filename))

sdata = sd.read_zarr(os.path.join(spat_dat_dir, filename))

# Filter transcripts to just transcripts of interest to save time & memory
transcripts_of_interest = ['SPP1']

# Assuming the gene column is 'feature_name'
filtered_transcripts = sdata.points['transcripts'][
    sdata.points['transcripts']['feature_name'].isin(transcripts_of_interest)
]

# Replace the transcripts in the sdata object with the filtered ones
sdata.points['transcripts'] = filtered_transcripts

# Transform the sdata to the image space
transformed_sdata = sd.transform(sdata, to_coordinate_system='image_space')


# --- Load the fulres H&E image

image_dir = os.path.join(project_dir, 'data', 'images', 'he', 'aligned_images', 'xenium')
image_file = f'{sample_id}.ome.tif'

# Need to set this or it won't let me import my image
PIL.Image.MAX_IMAGE_PIXELS = 17963679744
img = tifffile.imread(os.path.join(image_dir, image_file))

# Add the scan to the anndata object for easy plotting
original_hires_key = 'original_hires'
adata.uns['spatial'][sample_id]['images'][original_hires_key] = img

# -- Calculate the scale factor
# Divide the height of the original hires image by the height of the downsampled hires image
size_fact = img.shape[0]/adata.uns['spatial'][sample_id]['images']['hires'].shape[0]

# Extract the downsampled hires scale factor
hires_fact = adata.uns['spatial'][sample_id]['scalefactors']['tissue_hires_scalef']

# Calculate the scale factor for the original hires image and add to the anndata object
adata.uns['spatial'][sample_id]['scalefactors'][f'tissue_{original_hires_key}_scalef'] = size_fact * hires_fact


# --- Load the areas of interest

filename = f'{sample_id}_cropped_HE_TLS_SPP1.geojson'

resources_dir = os.path.join(project_dir, 'resources', 'areas_of_interest', 'xenium')
areas = gpd.read_file(os.path.join(resources_dir, filename))

# Extract a dict with x0, x1, y0, y1 for each area of interest
area_coords = {}

for idx, row in areas.iterrows():
    # Parse the 'classification' string into a dictionary
    classification = json.loads(row['classification'])
    area_name = classification['name']
    bounds = row['geometry'].bounds  # returns (minx, miny, maxx, maxy)
    area_coords[area_name] = (int(bounds[0]), int(bounds[2]), int(bounds[1]), int(bounds[3]))

# Make sure the overview area is first
area_coords = {k: area_coords[k] for k in sorted(area_coords.keys(), key=lambda x: (x != 'overview', x))}


# --- Load transcript colors

color_dir = os.path.join(project_dir, 'resources', 'plotting_colors')
color_file = 'transcript_colors.csv'

transcript_colors = pd.read_csv(os.path.join(color_dir, color_file), index_col=0)

# Make a dict
tx_color_dict = transcript_colors['color'].to_dict()


# --- Functions for plotting

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

# Make a pandas df from filtered_transcripts
filtered_transcripts = transformed_sdata.points['transcripts']
filtered_transcripts_df = filtered_transcripts.compute()


# --- Plot 1: Overview of SPP1 inside TLS

microns_per_pixel = adata.uns['spatial'][sample_id]['scalefactors']['microns_per_pixel']
scalebar_units = 'µm'

nrows = 2
ncols = len(area_coords)-1

fig, axs = plt.subplots(nrows, ncols, figsize=(5, 5))

scalebar_val = 20

# Subset to the non overview area
area_coords_subset = {k: v for i, (k, v) in enumerate(area_coords.items()) if i != 0}

# On the first row, plot crops of the H&E image
for n, (area, coords) in enumerate(area_coords_subset.items()):
    ax = axs[0, n]
    # Crop the image
    crop_img, adjusted_coords = crop_image(img, coords)
    ax.imshow(crop_img)
    # Add the scalebar
    scalebar = ScaleBar(microns_per_pixel, units=scalebar_units, location='lower right', box_alpha=0.5, fixed_value=scalebar_val,
                        frameon=False, scale_loc='none', pad=0.3, height_fraction=0.02, font_properties={"size": 18})
    ax.add_artist(scalebar)
    # Remove axis ticks
    ax.set_xticks([])
    ax.set_yticks([])

# On the second row, plot SPP1 transcripts
for n, (area, coords) in enumerate(area_coords_subset.items()):
    ax = axs[1, n]
    # Crop the H&E image
    crop_img, adjusted_coords = crop_image(img, coords)
    ax.imshow(crop_img, alpha=0.5)
    # Subset the transcripts based on the crop region
    mask = (
        (filtered_transcripts_df["x"] >= coords[0]) & 
        (filtered_transcripts_df["x"] < coords[1]) &
        (filtered_transcripts_df["y"] >= coords[2]) &
        (filtered_transcripts_df["y"] < coords[3]) &
        (filtered_transcripts_df["feature_name"] == 'SPP1')
    )
    transcripts_sub = filtered_transcripts_df[mask]
    # Adjust the coordinates for the cropped image
    transcripts_sub['x_adj'] = transcripts_sub['x'] - adjusted_coords[0]
    transcripts_sub['y_adj'] = transcripts_sub['y'] - adjusted_coords[2]
    # Plot the transcripts, coloring by gene
    if not transcripts_sub.empty:
        # Set s based on crop_width
        crop_width = coords[1] - coords[0]
        s = int(0.001 * (crop_width ** 2))
        # Create a darker version of the color for the edge
        color = tx_color_dict.get('SPP1', 'black')
        darker_color = tuple(max(c - 0.3, 0) for c in mcolors.to_rgb(color))
        ax.scatter(transcripts_sub['x_adj'], transcripts_sub['y_adj'], s=s, color=tx_color_dict['SPP1'], label='SPP1', alpha=1, linewidth=1, edgecolor=darker_color)
    # Reset the limits
    width_crop = coords[1] - coords[0]
    ax.set_xlim(0, width_crop)
    ax.set_ylim(width_crop, 0)  # Invert y-axis
    # Disable grid & remove axis ticks
    ax.grid(False)
    # Add the scalebar
    scalebar = ScaleBar(microns_per_pixel, units=scalebar_units, location='lower right', box_alpha=0.5, fixed_value=scalebar_val,
                        frameon=False, scale_loc='none', pad=0.3, height_fraction=0.02, font_properties={"size": 18})
    ax.add_artist(scalebar)

# Remove the x & y labels of all axs
for ax in axs.flat:
    ax.set_xlabel('')
    ax.set_ylabel('')
    ax.set_xticks([])
    ax.set_yticks([])

# Resize the plots
for ax in axs.flat:
    ax_pos = ax.get_position()
    # Set the width to be equal for all plots
    ax.set_position([ax_pos.x0, ax_pos.y0, 0.8, 0.8])

move_plots_x_y(axs, x_offset=0.05, y_offset=0.05, remove_title=True)

# Add numbers as titles
for n, ax in enumerate(axs[0, :]):
    ax.set_title(str(n+1), fontweight='bold', fontsize=25)

# Add SPP1 legend to the last plot
legend_items = ['SPP1']
legend_handles = []
for i, item in enumerate(legend_items):
    color = tx_color_dict.get(item, 'black')
    darker_color = tuple(max(c - 0.3, 0) for c in mcolors.to_rgb(color))
    legend_handles.append(Line2D([0], [0], marker='o', color='w', label=item, markerfacecolor=color, markersize=15, linewidth=2, markeredgecolor=darker_color))

axs[1, -1].legend(handles=legend_handles, loc='center left', bbox_to_anchor=(1, 0.5), fontsize=25, frameon=False, title_fontsize=18,
                handlelength=0.7)

# Position the main title
title_x = (axs[0,0].get_position().x0 + axs[0,-1].get_position().x1) / 2
title_y = axs[0,0].get_position().y1 + 0.4
plt.suptitle(f'{sample_id} ({subtype}) - SPP1 inside TLS (scalebar_val={scalebar_val}{scalebar_units})', x=title_x, y=title_y, fontsize=30, fontweight='bold')


# --- Save plot

out_dir = os.path.join(res_dir, 'figures')
os.makedirs(out_dir, exist_ok=True)

filename = 'Fig_5E.pdf'
plt.savefig(os.path.join(out_dir, filename),  bbox_inches='tight')


# --- Finish

print("Script successfully completed")