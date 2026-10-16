#!/usr/bin/env python3

# --- Aim
# Plot TLS model overview with zoomed areas (Figure 2A)
# Shows TLS model score on full tissue + 2 zoomed areas with H&E and pathology annotations
# This script is set up for sample P24_Vis


# --- Define environment
# env: spatial-pilot

# --- Load packages
import argparse
import os
import json
import sys

import anndata as ad
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib_scalebar.scalebar import ScaleBar
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
import scanpy as sc
import tifffile

from mycolorpy import colorlist as mcp

sc.set_figure_params(facecolor="white", figsize=(8, 8))
sc.settings.verbosity = 3

mpl.rcParams['pdf.fonttype'] = 42


# --- Argparse arguments
parser = argparse.ArgumentParser(description="Plot TLS model overview for Figure 2A",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument("--res_dir", help="Path to results directory")
parser.add_argument("--h5ad_dir", help="Path to directory for .h5ad output files")
parser.add_argument("--project_dir", help="Project directory (root)")
parser.add_argument("--sample_id", default="P24_Vis", help="Sample shown in the figure")
parser.add_argument("--sample_sheet", help="Path to sample sheet CSV file", default=None)

args = vars(parser.parse_args())

# Arguments
project_dir = args["project_dir"]
h5ad_dir = args["h5ad_dir"]
res_dir = args["res_dir"]
sample_id = args["sample_id"]

print(f'Sample: {sample_id}')

sys.path.append(os.path.join(project_dir, 'scripts'))
from helper_functions import auto_crop_scanpy, generate_color_variations


# --- Extract some info from the sample sheet
if args.get("sample_sheet"):
    sample_sheet_path = args["sample_sheet"]
else:
    sample_sheet_path = os.path.join(project_dir, 'config', 'sample_sheet_visium.csv')

sample_sheet = pd.read_csv(sample_sheet_path)
tissue_id = sample_sheet[sample_sheet['section_id'] == sample_id]['tissue_id'].values[0]


# --- Also get some tissue metadata

resource_dir = os.path.join(project_dir, 'resources')
metadata_file = 'tissue_metadata.csv'

tissue_metadata = pd.read_csv(os.path.join(resource_dir, metadata_file))

# Merge the two DataFrames on the common column 'tissue_id'
merged_df = pd.merge(sample_sheet, tissue_metadata, on='tissue_id')
subtype = merged_df[merged_df['section_id'] == sample_id]['clinical_subtype'].values[0]


# --- Load the adata object

adata = ad.read_h5ad(os.path.join(h5ad_dir, 'visium', f'{sample_id}.h5ad'))


# --- Load TLS model results

TLS_res_dir = os.path.join(res_dir, '03_TLS_niche_discovery', '01_identify_TLS', '01_run_TLS_model')

most_recent_dir_path = TLS_res_dir

# Read in the results & add to anndata
if most_recent_dir_path:
    TLS_file_name = f'{sample_id}_TLS_model.csv'
    TLS_results = pd.read_csv(os.path.join(most_recent_dir_path, TLS_file_name), index_col=0)
    adata.obs['TLS_model'] = TLS_results
    print(f"Loaded TLS model results from: {os.path.join(most_recent_dir_path, TLS_file_name)}")
else:
    print("Warning: TLS model results not found!")


# --- Calculate crop coords for optimal presentation

coords = auto_crop_scanpy(adata, border_size=0.1)
spot_size = 1.5


# --- Import the original high-resolution scan

# Define key for adding to anndata
original_hires_key = 'original_hires'

# Find the original high-res image
image_dir = os.path.join(project_dir, 'data', 'images', 'he', 'aligned_images', 'visium')
image_file = f'{sample_id}.ome.tif'
image_file_path = os.path.join(image_dir, image_file)

# Load the original high-res scan using tifffile
img = tifffile.imread(image_file_path)

# Add the scan to the anndata object for easy plotting
adata.uns['spatial'][sample_id]['images'][original_hires_key] = img

# Calculate the scale factor
# Divide the height of the original hires image by the height of the downsampled hires image
size_fact = img.shape[0] / adata.uns['spatial'][sample_id]['images']['hires'].shape[0]

# Extract the downsampled hires scale factor
hires_factor = adata.uns['spatial'][sample_id]['scalefactors']['tissue_hires_scalef']

# Calculate the scale factor for the original hires image and add to the anndata object
fulres_fact = size_fact * hires_factor
adata.uns['spatial'][sample_id]['scalefactors'][f'tissue_{original_hires_key}_scalef'] = fulres_fact


# --- Calculate microns per pixel for scale bar

spot_diameter_fulres = adata.uns['spatial'][sample_id]['scalefactors']['spot_diameter_fullres']
spot_diameter = 55  # in microns
microns_per_pixel = spot_diameter / spot_diameter_fulres


# --- Load the pathology annotations

annotation_dir = os.path.join(project_dir, 'data', 'pathology_annotations', 'visium')
filename = f'{sample_id}_pathology_annotations.csv'

annotations = pd.read_csv(os.path.join(annotation_dir, filename), index_col=0)
adata.obs['path_annotations'] = annotations[['single_layer']]


# --- Load the desired categories

annotation_dir_base = os.path.join(project_dir, 'resources', 'pathology_annotations')
filename = 'pathology_categories_per_layer.json'
path_to_file = os.path.join(annotation_dir_base, filename)

with open(path_to_file, 'r') as json_file:
    annotations_per_layer = json.load(json_file)


# --- Define a color_dict

def create_plotting_colors(category_dict, ordered_categories):
    # Get some pre-made colors
    colors = sc.plotting.palettes.default_20
    # Remove red to be more colorblind-friendly
    colors = colors[0:3] + colors[4:]
    # Get some color variations of the original colors (so variations of the same annotation are a similar color)
    color_variations = {}
    for n, (key, values) in enumerate(category_dict.items()):
        color_variations[key] = generate_color_variations(colors[n], len(values))
    # Create a new dictionary mapping values from dict1 to values from dict2
    color_dict = {val: color_variations[key][idx] for key, values in category_dict.items() for idx, val in enumerate(values)}
    return color_dict


# Extract a flat list of ordered categories
flat_ordered_categories = [item for sublist in annotations_per_layer.values() for item in sublist]

# Generate colors
color_dict = create_plotting_colors(annotations_per_layer, flat_ordered_categories)


# --- Functions for plotting

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
    adata_crop.obs["pxl_row_in_fullres"] = adata_crop.obsm["spatial"][:, 1] * fulres_fact  # Y-coordinates
    # Subset the AnnData object based on the crop region
    mask = (
        (adata_crop.obs["pxl_col_in_fullres"] >= xmin) & 
        (adata_crop.obs["pxl_col_in_fullres"] < xmax) &
        (adata_crop.obs["pxl_row_in_fullres"] >= ymin) & 
        (adata_crop.obs["pxl_row_in_fullres"] < ymax)
    )
    adata_sub = adata_crop[mask]
    return adata_sub, crop_img, (xmin, ymin, width_crop, height_crop)


# --- Define custom crop coordinates for zoomed in areas

zoom_coords = {
    "1": {
        "top_left": (9800, 4008),
        "width": 2707
    },
    "2": {
        "top_left": (5852, 14878),
        "width": 3394
    }
}


# --- Define custom colormap for plotting TLS model results

cmap = "YlOrRd"
color_list = mcp.gen_color(cmap=cmap, n=11)
color_list_extracted = ['lightgray'] + color_list[2:]

n_bins = 100
cmap_name = 'my_list'
cmap_adj = LinearSegmentedColormap.from_list(cmap_name, color_list_extracted, N=n_bins)

cmap_adj.set_under(color='lightgray', alpha=0.7)


# --- Set plotting parameters

sc.set_figure_params(fontsize=16)

# Set the scalebar size
scalebar_val = 200  # in microns
scalebar_units = 'µm'

# Set vmax & vmin
vmin = 0
vmax = 0.05

nrows = 1
ncols = 1

fig, ax = plt.subplots(nrows, ncols, figsize=(5, 5))


# --- Plot main TLS model

sc.pl.spatial(adata, color='TLS_model', size=spot_size, cmap=cmap_adj, bw=False, alpha_img=1, crop_coord=coords, title=f'', vmin=vmin, ax=ax,
                colorbar_loc=None, vmax=vmax)

# Get the axs
cur_axs = fig.get_axes()
main_ax = cur_axs[0]

# Resize the main ax
main_ax_pos = main_ax.get_position()
main_ax_width = 1
main_ax.set_position([main_ax_pos.x0, main_ax_pos.y0, main_ax_width, main_ax_width])
main_ax_pos = main_ax.get_position()
offset = 0.02

# Remove ax labels
main_ax.set_xlabel('')
main_ax.set_ylabel('')

# Add rectangles to mark the zoom areas
for area, crop_coord in zoom_coords.items():
    # Adjust the coordinates for the now cropped image
    crop_x = crop_coord['top_left'][0] * hires_factor
    crop_y = crop_coord['top_left'][1] * hires_factor
    rect = mpl.patches.Rectangle((crop_x, crop_y), crop_coord['width'] * hires_factor, crop_coord['width'] * hires_factor, linewidth=2, edgecolor='black', facecolor='none')
    main_ax.add_patch(rect)


# --- Add zoom axes

zoom_axs = []

zoom_ax = fig.add_axes([0.6, 0.6, 0.3, 0.3])
zoom_axs.append(zoom_ax)

# Position the zoom ax
mid_offset = offset * 2
ax_height = main_ax_pos.height/2 - mid_offset/2 + offset*2
start_offset = main_ax_pos.width*0.2
ax_pos_y = main_ax_pos.y1 - ax_height + offset*2
zoom_ax.set_position([main_ax_pos.x1 - start_offset, ax_pos_y, ax_height, ax_height])

# Crop the adata for first zoom
adata_crop = adata.copy()
crop_coord = zoom_coords['1']['top_left']
width = zoom_coords['1']['width']

adata_sub, crop_img, (xmin, ymin, width_crop, height_crop) = crop_and_subset_anndata(
    adata=adata,
    adata_crop=adata_crop,
    img=img,
    sample_id=sample_id,
    original_hires_key=original_hires_key,
    zoom_coord_fulres=crop_coord,
    width_fulres=width,
    fulres_fact=fulres_fact
)

ax = zoom_ax

sc.pl.spatial(adata_sub, color='TLS_model', ax=ax, bw=False, alpha_img=1, cmap=cmap_adj, vmin=vmin, vmax=vmax, colorbar_loc=None,
                title=f'', img_key=original_hires_key, crop_coord=(xmin, width_crop, ymin, height_crop))

ax.set_title('TLS model', fontweight='bold')

scalebar = ScaleBar(microns_per_pixel, units=scalebar_units, location='lower right', box_alpha=0.5, fixed_value=scalebar_val,
                    frameon=False, scale_loc="none", pad=0.3, height_fraction=0.02)
ax.add_artist(scalebar)

# Add second zoom ax
zoom_ax = fig.add_axes([0.6, 0.6, 0.3, 0.3])
zoom_axs.append(zoom_ax)

# Position the zoom ax
zoom_ax.set_position([main_ax_pos.x1 - start_offset, main_ax_pos.y0 - offset*2, ax_height, ax_height])

# Crop the adata for second zoom
adata_crop = adata.copy()
crop_coord = zoom_coords['2']['top_left']
width = zoom_coords['2']['width']

adata_sub, crop_img, (xmin, ymin, width_crop, height_crop) = crop_and_subset_anndata(
    adata=adata,
    adata_crop=adata_crop,
    img=img,
    sample_id=sample_id,
    original_hires_key=original_hires_key,
    zoom_coord_fulres=crop_coord,
    width_fulres=width,
    fulres_fact=fulres_fact
)

ax = zoom_ax

sc.pl.spatial(adata_sub, color='TLS_model', ax=ax, bw=False, alpha_img=1, cmap=cmap_adj, vmin=vmin, vmax=vmax, colorbar_loc=None,
                title=f'', img_key=original_hires_key, crop_coord=(xmin, width_crop, ymin, height_crop))

scalebar = ScaleBar(microns_per_pixel, units=scalebar_units, location='lower right', box_alpha=0.5, fixed_value=scalebar_val,
                    frameon=False, scale_loc="none", pad=0.3, height_fraction=0.02)
ax.add_artist(scalebar)

# Add & position the colorbar
ax_pos = ax.get_position()
norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
cbar_width = ax_pos.width * 0.4
cbar_height = ax_pos.height * 0.05
x0 = (ax_pos.x0 + ax_pos.width/2) - cbar_width/2
y0 = ax_pos.y0 - offset - cbar_height
cbar_ax = fig.add_axes([x0, y0, cbar_width, cbar_height])
cbar = fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap_adj), 
                    cax=cbar_ax, orientation='horizontal')

# Set xticks at vmin & vmax only
cbar.ax.set_xticks([vmin, vmax])

# Set xticklabel fontsize
cbar.ax.tick_params(labelsize=12)


# --- Add H&E and pathology annotation panels

# Keep track of annotations plotted
annotations_plotted = []

for n, zoom_ax in enumerate(zoom_axs):
    # Add H&E panel
    zoom_ax_pos = zoom_ax.get_position()
    new_zoom_ax = fig.add_axes([zoom_ax_pos.x1+offset, zoom_ax_pos.y0, zoom_ax_pos.width, zoom_ax_pos.height])
    
    # Crop the adata
    adata_crop = adata.copy()
    crop_coord = zoom_coords[str(n+1)]['top_left']
    width = zoom_coords[str(n+1)]['width']
    
    adata_sub, crop_img, (xmin, ymin, width_crop, height_crop) = crop_and_subset_anndata(
        adata=adata,
        adata_crop=adata_crop,
        img=img,
        sample_id=sample_id,
        original_hires_key=original_hires_key,
        zoom_coord_fulres=crop_coord,
        width_fulres=width,
        fulres_fact=fulres_fact
    )
    
    # Plot the H&E
    sc.pl.spatial(adata_sub, color=None, cmap=cmap_adj, bw=False, alpha_img=1, crop_coord=(xmin, width_crop, ymin, height_crop), 
                    title=f'', vmin=vmin, ax=new_zoom_ax, colorbar_loc=None, vmax=vmax, img_key=original_hires_key)
    if n == 0:
        new_zoom_ax.set_title(f'H&E', fontweight='bold')
    
    # Add scalebar
    scalebar = ScaleBar(microns_per_pixel, units=scalebar_units, location='lower right', box_alpha=0.5, fixed_value=scalebar_val,
                        frameon=False, scale_loc="none", pad=0.3, height_fraction=0.02)
    new_zoom_ax.add_artist(scalebar)
    
    # Add pathology annotations panel
    zoom_ax_pos = new_zoom_ax.get_position()
    new_zoom_ax = fig.add_axes([zoom_ax_pos.x1+offset, zoom_ax_pos.y0, zoom_ax_pos.width, zoom_ax_pos.height])
    sc.pl.spatial(adata_sub, color='path_annotations', ax=new_zoom_ax, bw=False, alpha_img=1, palette=color_dict, legend_loc=None,
                    title=f'', img_key=original_hires_key, crop_coord=(xmin, width_crop, ymin, height_crop))
    if n == 0:
        new_zoom_ax.set_title('Pathology\nannotations', fontweight='bold')
    
    # Add scalebar
    scalebar = ScaleBar(microns_per_pixel, units=scalebar_units, location='lower right', box_alpha=0.5, fixed_value=scalebar_val,
                        frameon=False, scale_loc="none", pad=0.3, height_fraction=0.02)
    new_zoom_ax.add_artist(scalebar)
    
    plotted_anno = adata_sub.obs['path_annotations'].dropna().unique()
    # Add the annotations to the list of plotted annotations
    annotations_plotted += list(plotted_anno)


# --- Remove x & y labels for all axes

axs = fig.get_axes()
for ax in axs:
    ax.set_xlabel('')
    ax.set_ylabel('')


# --- Add connecting lines between zoom areas and main plot

# Draw connecting lines from rect (right corners) to zoom plot (left corners)
# Use transforms to map data to figure coordinates
trans_main = main_ax.transData.transform

for n, zoom_ax in enumerate(zoom_axs):
    trans_zoom = zoom_ax.transAxes.transform
    inv = fig.transFigure.inverted().transform
    # Get the coordinates of the rectangle corners in data coords
    HEV_x = zoom_coords[str(n+1)]['top_left'][0] * hires_factor
    HEV_y = zoom_coords[str(n+1)]['top_left'][1] * hires_factor
    width_HEV = zoom_coords[str(n+1)]['width'] * hires_factor
    # Rect corners in data coords (main ax)
    p1 = (HEV_x + width_HEV, HEV_y + width_HEV)  # top-right
    p2 = (HEV_x + width_HEV, HEV_y)  # bottom-right
    # Zoom corners in axes coords (zoom ax)
    q1 = (0, 0)  # bottom-left of zoom ax
    q2 = (0, 1)  # top-left of zoom ax
    # Convert all to figure coords
    p1_fig = inv(trans_main(p1))
    p2_fig = inv(trans_main(p2))
    q1_fig = inv(trans_zoom(q1))
    q2_fig = inv(trans_zoom(q2))
    # Create a polygon from the four corners
    connector = mpl.patches.Polygon(
        [p1_fig, p2_fig, q2_fig, q1_fig],
        closed=True,
        edgecolor='black',
        facecolor=(1, 1, 1, 0.4),
        transform=fig.transFigure,
        linewidth=1
    )
    fig.patches.append(connector)


# --- Add annotations legend

annotations_plotted = list(set(annotations_plotted))  # Remove duplicates
# Sort according to the order in the annotations_per_layer dict
annotations_plotted = sorted(annotations_plotted, key=lambda x: flat_ordered_categories.index(x) if x in flat_ordered_categories else float('inf'))

# Create a legend for the pathology annotations
handles = []
for cat in annotations_plotted:
    if cat in color_dict:
        color = color_dict[cat]
        handle = Line2D([0], [0], marker='o', color='none', markerfacecolor=color, label=cat, markersize=8, alpha=1, 
                        markeredgewidth=0)
        handles.append(handle)

# Clean the labels for the legend
legend_labels = [label.capitalize().replace('_', ' ') for label in annotations_plotted]

# Have some labels in uppercase
uppercase = ['tls', 'lvi']

for i, label in enumerate(legend_labels):
    if label.lower() in uppercase:
        legend_labels[i] = label.upper()
    if label == 'Immune other':
        legend_labels[i] = 'Lymphocytes\npresent'

target_ax = axs[-1]  # Use the last ax as the target for the legend

target_ax.legend(handles=handles, labels=legend_labels, loc='center left', bbox_to_anchor=(1, 1.05), fontsize=12,
           handlelength=0.7, frameon=False, ncols=1, columnspacing=0.5)


# --- Position title

title_y = main_ax_pos.y1 + offset * 10
title_x = main_ax_pos.x0 + axs[-1].get_position().x1 / 2

plt.suptitle(f'{sample_id} ({subtype}) - TLS model score (scalebar={scalebar_val}{scalebar_units})', fontsize=16, fontweight='bold', y=title_y, x=title_x)


# --- Save the output

out_dir = os.path.join(res_dir, 'figures')
os.makedirs(out_dir, exist_ok=True)

filename = 'Fig_2A.pdf'
plt.savefig(os.path.join(out_dir, filename), bbox_inches='tight')

print(f"Saved plot to: {os.path.join(out_dir, filename)}")
print("Script successfully completed")
