#!/usr/bin/env python3

# --- Aim
# Plot selected TLSs with H&E, pathology annotations, TLS model, niche signatures, and deconvolution (Figure 2F)
# Horizontal orientation: 4 TLSs × 6 panels (H&E, pathology, TLS model, niche signatures, B-cells, T-cells)

# --- Define environment
# env: vis-hd

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
from matplotlib.colors import to_rgb
from matplotlib.colors import ListedColormap
import numpy as np
import pandas as pd
import scanpy as sc
import tifffile
from PIL import Image

from mycolorpy import colorlist as mcp

sc.set_figure_params(facecolor="white", figsize=(8, 8))
sc.settings.verbosity = 3

mpl.rcParams['pdf.fonttype'] = 42


# --- Argparse arguments
parser = argparse.ArgumentParser(description="Plot selected TLSs for Figure 2F",
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

print(f'Project: {project_dir}')

sys.path.append(os.path.join(project_dir, 'scripts'))
from helper_functions import generate_color_variations


# --- Load sample sheet
if args.get("sample_sheet"):
    sample_sheet_path = args["sample_sheet"]
else:
    sample_sheet_path = os.path.join(project_dir, 'config', 'sample_sheet_visium.csv')

sample_sheet = pd.read_csv(sample_sheet_path)
print(f'Sample sheet: {sample_sheet_path}')
print(f'Number of samples: {len(sample_sheet)}')


# --- Define TLS coordinates
# Each TLS has: sample_id (new_sample_id from sample sheet), top_left coordinates, and width
# Coordinates for the TLS subselection ['1', '4', '5', '6']

TLS_coords = {
    "1": {
        "new_sample_id": "P22_Vis",
        "top_left": (3865, 11311),
        "width": 2296
    },
    "2": {
        "new_sample_id": "P24_Vis",
        "top_left": (14981, 10509),
        "width": 1416
    },
    "3": {
        "new_sample_id": "P29_Vis_Rep1",
        "top_left": (23089, 8565),
        "width": 1883
    }
}

# Subselection for horizontal plot (3 TLSs)
subselection = ['1', '2', '3']
TLS_coords_sub = {k: TLS_coords[k] for k in subselection}


# --- Load pathology categories
annotation_dir_base = os.path.join(project_dir, 'resources', 'pathology_annotations')
filename = 'pathology_categories_per_layer.json'
path_to_file = os.path.join(annotation_dir_base, filename)

with open(path_to_file, 'r') as json_file:
    annotations_per_layer = json.load(json_file)

# Extract annotation layers
annotation_layers = list(annotations_per_layer.keys())

# Create ordered list of categories
ordered_labels_per_layer = annotations_per_layer
ordered_columns_major = []
for layer in annotation_layers:
    ordered_columns_major.extend(ordered_labels_per_layer[layer])


# --- Define color dictionary for pathology annotations
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


color_dict_path = create_plotting_colors(annotations_per_layer, ordered_columns_major)


# --- Define colormaps for TLS model and deconvolution
# TLS model colormap (YlOrRd with lightgray threshold)
cmap = "YlOrRd"
color_list = mcp.gen_color(cmap=cmap, n=11)
color_list_extracted = ['lightgray'] + color_list[2:]

n_bins = 100
cmap_name = 'my_list'
cmap_adj = LinearSegmentedColormap.from_list(cmap_name, color_list_extracted, N=n_bins)
cmap_adj.set_under(color='lightgray', alpha=0.7)

# Deconvolution colormap (BuPu, no lightgray)
cmap_decon = "BuPu"
color_list_decon = mcp.gen_color(cmap=cmap_decon, n=11)
color_list_decon = color_list_decon[1:]  # Remove first color
cmap_adj_decon = LinearSegmentedColormap.from_list('decon_cmap', color_list_decon, N=n_bins)
cmap_adj_decon.set_under(color=color_list_decon[0], alpha=0.4)


# --- Define colormaps for niche signatures
# Define score names
score_names_r = ['T_cell_zone', 'B_cell_zone', 'GC']
top_genes = 20  # Number of top genes used for scoring

# Create color dictionaries for each score
# Use function from cell2location
white_spacing=20

# Use function from cell2location to make colormap with increasing transparency gradient
# To be able to plot multiple scores in single plot
def create_colormap(R, G, B):
    spacing = int(white_spacing * 2.55)
    N = 255
    M = 3
    alphas = np.concatenate([[0] * spacing * M,
                                np.linspace(0, 1.0, (N - spacing) * M)])
    vals = np.ones((N * M, 4))
    #         vals[:, 0] = np.linspace(1, R / 255, N * M)
    #         vals[:, 1] = np.linspace(1, G / 255, N * M)
    #         vals[:, 2] = np.linspace(1, B / 255, N * M)
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

for score_name in score_names_r:
    rgb = color_dict_rgb[score_name]
    cmap = create_colormap(rgb[0], rgb[1], rgb[2])
    cmap_dict[score_name] = cmap



# --- Function to crop image
def crop_image(img, coords):
    # Clip the coordinates to ensure they are within the image boundaries
    x_min = max(coords[0], 0)
    x_max = min(coords[1], img.shape[1])
    y_min = max(coords[2], 0)
    y_max = min(coords[3], img.shape[0])
    # Crop the image using the clipped coordinates
    crop_img = img[y_min:y_max, x_min:x_max]
    return crop_img, (x_min, x_max, y_min, y_max)


# --- Function to crop & subset adata
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


# --- Function to move plots
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


# --- Adjust the crop coords
# So that the plots will be 2:1 aspect ratio
# This is done by adjusting the width of the crop to be 2x the height

aspect_ratio = 1.5  # Set aspect ratio width:height
height_adj = 1.3 # Adjust height to be a bit larger than originally defined

for TLS, coords in TLS_coords.items():
    sample = coords['new_sample_id']
    top_left = coords['top_left']
    width = coords['width']
    new_height = int(width * height_adj)
    # Adjust the top left coordinates to center the crop
    center_x = top_left[0] + width // 2
    center_y = top_left[1] + width // 2
    width_adjusted = int(new_height * aspect_ratio)
    top_left_adjusted_x = int(center_x - width_adjusted // 2)
    top_left_adjusted_y = int(center_y - new_height // 2)
    top_left_adjusted = (top_left_adjusted_x, top_left_adjusted_y)
    TLS_coords[TLS]['top_left'] = top_left_adjusted
    TLS_coords[TLS]['width'] = new_height


# --- Load fulres images & create cropped anndatas
cropped_adata_dict = {}

# Define key for adding to anndata
original_hires_key = 'original_hires'

# Load images
image_dir = os.path.join(project_dir, 'data', 'images', 'he', 'aligned_images', 'visium')

for n, TLS in enumerate(TLS_coords):
    coords = TLS_coords[TLS]
    new_sample_id = coords['new_sample_id']
    
    print(f'Processing TLS {TLS} - Sample: {new_sample_id}')
    
    # Load the h5ad file
    adata = ad.read_h5ad(os.path.join(h5ad_dir, 'visium', f'{new_sample_id}.h5ad'))
    
    # Load pathology annotations
    annotation_dir = os.path.join(project_dir, 'data', 'pathology_annotations', 'visium')
    filename = f'{new_sample_id}_pathology_annotations.csv'
    annotations = pd.read_csv(os.path.join(annotation_dir, filename), index_col=0)
    adata.obs['single_layer'] = annotations[['single_layer']]
    
    # Load TLS model results
    if 'TLS_model' not in adata.obs.columns:
        TLS_res_dir = os.path.join(res_dir, '03_TLS_niche_discovery', '01_identify_TLS', '01_run_TLS_model')
        if os.path.exists(TLS_res_dir):
            TLS_file = os.path.join(TLS_res_dir, f'{new_sample_id}_TLS_model.csv')
            if os.path.exists(TLS_file):
                adata.obs['TLS_model'] = pd.read_csv(TLS_file, index_col=0)
                print(f"  Loaded TLS model results from: {TLS_file}")
            else:
                print(f"  Warning: TLS model file not found: {TLS_file}")
        else:
            print(f"  Warning: TLS model results directory not found: {TLS_res_dir}")
    
    # Load signature scoring results
    # Map script variable names to file directory names
    signature_mapping = {
        'T_cell_zone': 'T_cell_zone_score',
        'B_cell_zone': 'B_cell_zone_score',
        'GC': 'GC_score'
    }
    top_genes = 20
    
    for script_name, file_name in signature_mapping.items():
        score_col = f'{script_name}_score_{top_genes}'
        if score_col not in adata.obs.columns:
            score_res_dir = os.path.join(res_dir, '03_TLS_niche_discovery', '04_map_niche_signatures',
                                        '01_score_signatures', file_name)
            if os.path.exists(score_res_dir):
                score_file = os.path.join(score_res_dir, 'csv_files',
                                        f'{new_sample_id}_{file_name}.csv')
                if os.path.exists(score_file):
                    score_results = pd.read_csv(score_file, index_col=0)
                    # Map file column name to script column name
                    file_col = f'{file_name}_{top_genes}'
                    if file_col in score_results.columns:
                        adata.obs[score_col] = score_results[file_col]
                        print(f"  Loaded {script_name} signature score from: {score_file}")
                    else:
                        print(f"  Warning: Column {file_col} not found in {score_file}")
                else:
                    print(f"  Warning: Signature score file not found: {score_file}")
            else:
                print(f"  Warning: Signature scoring directory not found: {score_res_dir}")
    
    # Load deconvolution results
    decon_cols_missing = ['B-cells', 'T-cells']
    decon_cols_missing = [col for col in decon_cols_missing if col not in adata.obs.columns]
    if decon_cols_missing:
        decon_file = os.path.join(res_dir, '02_cell_type_mapping', 'deconvolution', 'visium', 
                                 'major_level', f'{new_sample_id}_major_q05_cell_abundance_w_sf.csv')
        if os.path.exists(decon_file):
            decon_results = pd.read_csv(decon_file, index_col=0)
            for col in decon_results.columns:
                adata.obs[col] = decon_results[col]
            print(f"  Loaded deconvolution results from: {decon_file}")
        else:
            print(f"  Warning: Deconvolution file not found: {decon_file}")
    
    # Load the original high-res image using tifffile
    image_file = f'{new_sample_id}.ome.tif'
    image_file_path = os.path.join(image_dir, image_file)
    img = tifffile.imread(image_file_path)
    
    # Add the scan to the anndata object
    adata.uns['spatial'][new_sample_id]['images'][original_hires_key] = img
    
    # Calculate the scale factor
    size_fact = img.shape[0] / adata.uns['spatial'][new_sample_id]['images']['hires'].shape[0]
    hires_factor = adata.uns['spatial'][new_sample_id]['scalefactors']['tissue_hires_scalef']
    fulres_fact = size_fact * hires_factor
    adata.uns['spatial'][new_sample_id]['scalefactors'][f'tissue_{original_hires_key}_scalef'] = fulres_fact
    
    # Crop the adata
    top_left = coords['top_left']
    width = coords['width']
    
    adata_crop = adata.copy()
    adata_sub, crop_img, (xmin, ymin, width_crop, height_crop) = crop_and_subset_anndata(
        adata=adata,
        adata_crop=adata_crop,
        img=img,
        sample_id=new_sample_id,
        original_hires_key=original_hires_key,
        zoom_coord_fulres=top_left,
        width_fulres=width,
        fulres_fact=fulres_fact,
        aspect_ratio=aspect_ratio  # Adjusted to 2:1 aspect ratio
    )
    
    # Store the cropped adata
    cropped_adata_dict[TLS] = adata_sub.copy()


# --- Define plotting parameters
# Set scalebar
scalebar_val = 200  # in microns
scalebar_units = 'µm'

# Set fontsize
sc.set_figure_params(fontsize=20)

# Set microns per pixel for scale bar calculation
microns_per_pixel = 0.4426
# Visium spot diameter in microns (constant)
vis_spot_diameter = 55

# Deconvolution cell types to plot
decon_to_plot = ['B-cells', 'T-cells']

# Set vmin and vmax
vmin = 0
vmax = 0.8  # For niche signatures
vmax_tls = 0.05  # For TLS model
vmax_decon = 'p99'  # Use p99 for deconvolution results
vmin_decon = 0


# --- Plot: Horizontal orientation with selected TLSs

ncols = 4 + len(decon_to_plot)  # H&E, pathology, TLS model, niche signatures, + decon
nrows = len(TLS_coords_sub)  # Number of TLSs in the subselection

# Keep track of which annotations plotted for the legend
annotations_plotted = []

fig, axs = plt.subplots(nrows, ncols, figsize=(5, 5))

for n, TLS in enumerate(TLS_coords_sub):
    adata = cropped_adata_dict[TLS]
    # Get the crop coordinates
    coords = TLS_coords_sub[TLS]
    new_sample_id = coords['new_sample_id']
    top_left = coords['top_left']
    width = coords['width']
    
    # Plot the cropped H&E image
    ax = axs[n][0]
    sc.pl.spatial(adata, color=None, ax=ax, bw=False, alpha_img=1, title=f'',
                    crop_coord=(0, int(width*aspect_ratio), 0, width), img_key=original_hires_key)
    scalebar = ScaleBar(microns_per_pixel, units=scalebar_units, location='lower right', box_alpha=0.5, fixed_value=scalebar_val,
                    frameon=False, scale_loc='none', pad=0.3, height_fraction=0.04, font_properties={"size": 15})
    ax.add_artist(scalebar)
    
    # Plot the pathology annotations
    ax = axs[n][1]
    sc.pl.spatial(adata, color='single_layer', ax=ax, bw=False, alpha_img=0.7, palette=color_dict_path, 
                    spot_size=microns_per_pixel*vis_spot_diameter*6, 
                    crop_coord=(0, int(width*aspect_ratio), 0, width), title=f'', img_key=original_hires_key, legend_loc=None)
    annotations_plotted.extend(adata.obs['single_layer'].dropna().unique())
    
    # Plot the TLS model
    ax = axs[n][2]
    sc.pl.spatial(adata, color='TLS_model', ax=ax, bw=False, alpha_img=0.7, cmap=cmap_adj, 
                    spot_size=microns_per_pixel*vis_spot_diameter*6, 
                    crop_coord=(0, int(width*aspect_ratio), 0, width), title=f'', vmin=vmin, vmax=vmax_tls, 
                    colorbar_loc=None, img_key=original_hires_key)
    
    # Plot the niche signatures (overlay all three)
    ax = axs[n][3]
    for i, score_name in enumerate(score_names_r):
        if i == 0:
            alpha_img = 0.7
        else: 
            alpha_img = 0
        cmap = cmap_dict[score_name]
        sc.pl.spatial(adata, color=f'{score_name}_score_{top_genes}', ax=ax, bw=False, alpha_img=alpha_img, cmap=cmap, 
                      spot_size=microns_per_pixel*vis_spot_diameter*6,
                      crop_coord=(0, int(width*aspect_ratio), 0, width), title=f'', vmin=vmin, vmax=vmax, 
                      colorbar_loc=None, img_key=original_hires_key)
    
    # Plot the deconvolution results
    for i, decon_name in enumerate(decon_to_plot):
        ax = axs[n][4 + i]
        sc.pl.spatial(adata, color=decon_name, ax=ax, bw=False, alpha_img=0.7, cmap=cmap_adj_decon, 
                      spot_size=microns_per_pixel*vis_spot_diameter*6,
                      crop_coord=(0, int(width*aspect_ratio), 0, width), title=f'', vmin=vmin_decon, vmax=vmax_decon, 
                      colorbar_loc=None, img_key=original_hires_key)


# Remove the x & y labels of all axs
for ax in axs.flat:
    ax.set_xlabel('')
    ax.set_ylabel('')

# Resize the plots
for ax in axs.flat:
    ax_pos = ax.get_position()
    # Set the width to be equal for all plots
    ax.set_position([ax_pos.x0, ax_pos.y0, 0.4*aspect_ratio, 0.4])

move_plots_x_y(axs, x_offset=0.03, y_offset=0.03, remove_title=True)

# Set titles
titles = ['H&E', 'Pathology\nannotations', 'TLS model', 'Niche\nsignatures'] + decon_to_plot

for i, ax in enumerate(axs[0, :]):
    ax.set_title(titles[i], fontweight='bold')

# Add a pathology legend under the second column last plot
target_ax = axs[-1][1]
annotations_plotted = list(set(annotations_plotted))  # Remove duplicates

# Sort according to plot_order
annotations_plotted = [cat for cat in ordered_columns_major if cat in annotations_plotted]

# Clean the labels for the legend
annotations_plotted_adj = [label.capitalize().replace('_', ' ') for label in annotations_plotted]

# Have some labels in uppercase
uppercase = ['tls', 'lvi']

for i, label in enumerate(annotations_plotted_adj):
    if label.lower() in uppercase:
        annotations_plotted_adj[i] = label.upper()
    if label == 'Immune other':
        annotations_plotted_adj[i] = 'Lymphocytes\npresent'

# Create a legend for the pathology annotations
handles = []
for i, cat in enumerate(annotations_plotted):
    if cat in color_dict_path:
        color = color_dict_path[cat]
        handle = Line2D([0], [0], marker='o', color='none', markerfacecolor=color, label=annotations_plotted_adj[i], 
                        markersize=10, alpha=1, markeredgewidth=0)
        handles.append(handle)

# Set how many cols
n_cols = 2

# Position the legend
target_ax.legend(handles=handles, loc='upper center', bbox_to_anchor=(0.5, -0.03), fontsize=15,
           handlelength=0.7, frameon=False, ncol=n_cols)

# Add colorbars for each score
target_ax = axs[-1, 3]  # Use the last subplot for colorbars

# Dict to update the score names
score_names_alt = {
    'T_cell_zone': 'T cell zone',
    'B_cell_zone': 'B cell zone',
    'GC': 'GC'
}

for i, (label, cmap) in enumerate(cmap_dict.items()):
    norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
    spacing = 0.27  # vertical space between colorbars
    start_y = -0.3  # starting y position (below axis)
    y_pos = start_y - i * spacing  # stack downward
    cbar_width = 0.4  # width of the colorbar
    cbar_height = 0.1
    x_pos = 0.5 - cbar_width / 2  # center the colorbar horizontally
    cbar_ax = target_ax.inset_axes([x_pos, y_pos, cbar_width, cbar_height])
    cbar = fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap),
                        cax=cbar_ax, orientation='horizontal')
    cbar.ax.set_title(score_names_alt[label], fontsize=15)
    cbar.ax.tick_params(labelsize=15)
    # Limit to 2 xticks
    cbar.ax.set_xticks([0, vmax])
    # Remove the xticks if not the bottom colorbar
    if i < len(cmap_dict) - 1:
        cbar.ax.set_xticklabels([])

# Also add cbar for TLS model
target_ax = axs[-1, 2]  # Use the TLS model subplot for the colorbar
norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax_tls)
cbar_ax = target_ax.inset_axes([x_pos, start_y, cbar_width, cbar_height])
cbar = fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap_adj), 
                    cax=cbar_ax, orientation='horizontal')
cbar.ax.set_title('TLS model', fontsize=15)
cbar.ax.tick_params(labelsize=15)

# Limit to 2 xticks
cbar.ax.set_xticks([0, vmax_tls])

# And add a cbar for decon
target_ax = axs[-1, 4]  # Use the first deconvolution subplot for the colorbar
norm = mpl.colors.Normalize(vmin=0, vmax=1)
cbar_ax = target_ax.inset_axes([x_pos, start_y, cbar_width, cbar_height])
cbar = fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap_adj_decon), 
                    cax=cbar_ax, orientation='horizontal')
cbar.ax.set_title('Abundance', fontsize=15)
cbar.ax.tick_params(labelsize=15)

# Set only 2 labels: low and high
cbar.ax.set_xticks([0, 1])
cbar.ax.set_xticklabels(['Low', 'High'])

# Add the same cbar to the second ax
target_ax = axs[-1, 5]  # Use the second deconvolution subplot for the colorbar
cbar_ax = target_ax.inset_axes([x_pos, start_y, cbar_width, cbar_height])
cbar = fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap_adj_decon), 
                    cax=cbar_ax, orientation='horizontal')
cbar.ax.set_title('Abundance', fontsize=15)
cbar.ax.tick_params(labelsize=15)

# Set only 2 labels: low and high
cbar.ax.set_xticks([0, 1])
cbar.ax.set_xticklabels(['Low', 'High'])

# Add a rectangle marking deconvolution results
left_ax = axs[0][-2]
right_ax = axs[0][-1]

# Draw a line above the first to third axs in the first row
# Get subplot positions in figure coordinates
start_box = left_ax.get_position()
end_box = right_ax.get_position()

# Define rectangle dimensions
x0 = start_box.x0
x1 = end_box.x1
width_rect = abs(x1 - x0)

# Define y0 position
fig = plt.gcf()  # Get the current figure
fig.canvas.draw()
bbox = left_ax.get_tightbbox(fig.canvas.get_renderer())
inv_trans = fig.transFigure.inverted()
bbox_fig = inv_trans.transform(bbox)
top_edge = bbox_fig[1][1]  # Get the top edge of the first row

y0 = top_edge + 0.02  # bottom of the rectangle (just above the subplot)
height_rect = 0.01            # set height manually

# Create rectangle
highlight_rect = mpl.patches.Rectangle(
    (x0, y0), width_rect, height_rect,
    transform=fig.transFigure,
    facecolor='black',
    edgecolor='black',
    linewidth=0,
    alpha=1
)

# Add to figure
fig.patches.append(highlight_rect)

# Calculate midpoint y
y_text = y0 + height_rect + 0.01  # Slightly above the rectangle
x_text = x0 + width_rect / 2

# Add text label
fig.text(
    x_text, y_text,
    "Deconvolution",
    ha='center', va='bottom',
    color='black', fontweight='bold',
)

# Calculate position for title
title_x = axs[0, 0].get_position().x0 + (axs[0, -1].get_position().x1 - axs[0, 0].get_position().x0) / 2
title_y = axs[0, 0].get_position().y1 + 0.3  # Slightly above the first row

plt.suptitle(f'TLSs of interest: H&E, pathology annotations, TLS model & niche signatures (scalebar={scalebar_val} {scalebar_units})', 
             fontsize=24, fontweight='bold', y=title_y, x=title_x, ha='center')


# --- Save the output
out_dir = os.path.join(res_dir, 'figures')
os.makedirs(out_dir, exist_ok=True)

filename = 'Fig_2F.pdf'
plt.savefig(os.path.join(out_dir, filename), bbox_inches='tight')

print(f"Saved plot to: {os.path.join(out_dir, filename)}")
print("Script successfully completed")
