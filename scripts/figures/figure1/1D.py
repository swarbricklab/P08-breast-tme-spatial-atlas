#!/usr/bin/env python3

# --- Aim
# Make spatial overviews of 5 selected samples for Figure 1D
# Shows H&E images and pathology annotations for one sample per clinical subtype
# This creates Figure 1D

# --- Define environment
# env: spatial-pilot

# --- Load packages

import os
import sys
import json
import argparse

import scanpy as sc
import scanpy.plotting.palettes as sc_palettes
import anndata as ad
import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib_scalebar.scalebar import ScaleBar
from matplotlib.lines import Line2D

mpl.rcParams['pdf.fonttype'] = 42

# Custom functions
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))
from helper_functions import auto_crop_scanpy, generate_color_variations


# --- Arguments

# --- Argparse arguments
parser = argparse.ArgumentParser(description="Import Visium unfiltered data & plot",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument("-s", "--samples",
    nargs="*",  # 0 or more values expected => creates a list
    type=str,
    help="List of samples")
parser.add_argument("-r", "--res_dir",
                    help="Path to results directory (optional)")
parser.add_argument("--h5ad_dir", help="Path to directory for .h5ad output files") # Not necessary for this script
parser.add_argument("--project_dir", help="Project directory (root)")
parser.add_argument("--sample_sheet", help="Path to sample sheet CSV file")

args = vars(parser.parse_args())

# Arguments
h5ad_dir = args["h5ad_dir"]
res_dir = args["res_dir"]
sample_list = args["samples"]
project_dir = args["project_dir"]

print(f'Samples in python: {sample_list}')

sys.path.append(os.path.join(project_dir, 'scripts', 'figures'))
import dotplot_functions

# --- Load sample sheet

# Use sample_sheet argument if provided, else default
if args.get("sample_sheet"):
    sample_sheet_path = args["sample_sheet"]
else:
    sample_dir = os.path.join(project_dir, 'config')
    sample_sheet_file = 'sample_sheet_visium.csv'
    sample_sheet_path = os.path.join(sample_dir, sample_sheet_file)

sample_sheet = pd.read_csv(sample_sheet_path)
sample_lut = sample_sheet.set_index('section_id')


# --- Load tissue metadata

resource_dir = os.path.join(project_dir, 'resources')
metadata_file = 'tissue_metadata.csv'

tissue_metadata = pd.read_csv(os.path.join(resource_dir, metadata_file))

# Merge the two DataFrames on the common column 'tissue_id'
merged_df = pd.merge(sample_sheet, tissue_metadata, on='tissue_id')

# Create a sample - subtype dictionary
subtype_dict = dict(zip(merged_df['section_id'], merged_df['clinical_subtype']))

# --- Define samples of interest
# One sample per clinical subtype (all 6.5 mm capture area)

sample_list = ['P03_Vis', 'P15_Vis', 'P20_Vis', 'P23_Vis', 'P31_Vis']

# --- Load data

adata_dict = {}

for sample in sample_list:
    adata_dict[sample] = ad.read_h5ad(os.path.join(h5ad_dir, 'visium', f'{sample}.h5ad'))

# --- Calculate crop coords for optimal presentation

coords = {}

for sample in adata_dict.keys():
    coords[sample] = auto_crop_scanpy(adata_dict[sample], border_size=0.1)

# --- Set spot size for 6.5 mm capture area

spot_size = {}

for sample in adata_dict.keys():
    if sample_lut.loc[sample, 'slide_size'] == '11 mm':
        spot_size[sample] = 1.8
    else:
        spot_size[sample] = 1.6

# --- Load the pathology annotations

merged_column = 'single_layer'

pathology_dir = os.path.join(project_dir, 'data', 'pathology_annotations', 'visium')

for sample in sample_list:
    filename = f'{sample}_pathology_annotations.csv'
    path_to_file = os.path.join(pathology_dir, filename)
    if os.path.exists(path_to_file):   
        pathology_anno = pd.read_csv(path_to_file, index_col=0)
        adata_dict[sample].obs[merged_column] = pathology_anno[merged_column]
        print(f"Loaded pathology annotations for {sample}")
    else:
        # If we have no annotations: record that none of the spots have been annotated
        adata_dict[sample].obs[merged_column] = None
        print(f"Warning: No pathology annotations found for {sample}")

# --- Load the desired categories

annotation_layers = ['epithelial', 'immune', 'stroma', 'additional_features']

# Load the desired category order
annotation_dir_base = os.path.join(project_dir, 'resources', 'pathology_annotations')
filename = 'pathology_categories_per_layer.json'
path_to_file = os.path.join(annotation_dir_base, filename)

with open(path_to_file, 'r') as json_file:
    annotations_of_interest_dict = json.load(json_file)

# Collect all actual labels from the data (including those not in JSON)
all_data_labels = set()
for sample in sample_list:
    if sample in adata_dict:
        if merged_column in adata_dict[sample].obs.columns:
            all_data_labels.update(adata_dict[sample].obs[merged_column].dropna().unique())

# Get all unique annotation labels per layer that are present in the data
# Since single_layer is a merged column, we need to determine which layer each label belongs to
present_labels_per_layer = {}

for layer in annotation_layers:
    labels = set()
    # First, add labels from the JSON that are in the data
    desired_order = annotations_of_interest_dict[layer]
    labels.update([label for label in desired_order if label in all_data_labels])
    # Then add labels with the layer prefix (e.g., epithelial_mixed)
    labels.update([label for label in all_data_labels if label.startswith(f'{layer}_')])
    present_labels_per_layer[layer] = labels

# Now, for each layer, order the present labels according to annotations_of_interest_dict
# (same logic as 1D.py)
ordered_labels_per_layer = {}
for layer in annotation_layers:
    desired_order = annotations_of_interest_dict[layer]
    present_labels = present_labels_per_layer[layer]
    # Keep only those in present_labels, in the order of desired_order
    ordered_labels = [label for label in desired_order if label in present_labels]
    # Add any remaining labels that are in present_labels but not in desired_order
    remaining_labels = [label for label in present_labels if label not in ordered_labels]
    ordered_labels += remaining_labels
    ordered_labels_per_layer[layer] = ordered_labels

print("Ordered labels per layer:", ordered_labels_per_layer)

# --- Prepare for plotting

# Sort samples by subtype order
subtype_order = ['Luminal A', 'Luminal B', 'Luminal B (HER2+)', 'HER2+', 'TNBC']
sorted_samples = sorted(sample_list, key=lambda x: subtype_order.index(subtype_dict[x]) if subtype_dict[x] in subtype_order else float('inf'))

print(f"Sorted samples: {sorted_samples}")

# --- Define color generation function (same as 1D.py for consistency)

def create_plotting_colors(category_dict, ordered_categories, mixed=True, other=True, major_categories=False):
    # Get some colors
    if len(category_dict) < (21-2):
        colors = sc_palettes.default_20
        # Remove red (already using this for the mixed category)
        colors = colors[0:3] + colors[4:]
    else:
        colors = sc_palettes.godsnot_102[:len(category_dict.keys())]
    top_cat_colors = {}
    # If it's major categories: keep the colors consistent (epithelial blue, immune orange etc.)
    if major_categories:
        # Get the keys of the annotations of interest dict - this is the order to keep
        key_list = list(annotations_of_interest_dict.keys())
        for n, cat in enumerate(category_dict):
            i = key_list.index(cat)
            top_cat_colors[cat] = colors[i]
    else:
        for n, cat in enumerate(category_dict):
            top_cat_colors[cat] = colors[n]
    # Get some color variations of the original colors (so variations of the same annotation are a similar color)
    color_variations = {}
    for n, (key, values) in enumerate(category_dict.items()):
        if major_categories:
            i = key_list.index(key)
            color = colors[i]
        else:
            color = colors[n]
        num_cols = len(values)
        if num_cols > 4:
            step_size=0.1
        else:
            step_size=0.2
        # Import color variation function from dotplot_functions
        variations = generate_color_variations(color, num_cols, step_size=step_size)
        color_variations[key] = variations
    # Create a new dictionary mapping values from dict1 to values from dict2
    color_dict = {val: color_variations[key][idx] for key, values in category_dict.items() for idx, val in enumerate(values)}
    # Extract values from the dictionary based on the keys
    plotting_colors = [color_dict[key] for key in ordered_categories]
    return color_dict, plotting_colors, top_cat_colors

# --- Get colors per annotation

# Extract a flat list of ordered categories
ordered_columns_major = [item for sublist in ordered_labels_per_layer.values() for item in sublist]
print(f'Ordered columns for major categories: {ordered_columns_major}')

color_dict, plotting_colors, top_cat_colors = create_plotting_colors(
    ordered_labels_per_layer, ordered_columns_major, 
    mixed=False, major_categories=True
)

# --- Function for moving plots

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

# --- Plot: Spatial overview of pathology annotations

sc.set_figure_params(fontsize=25)

# Set Visium spot size
spot_diameter_micron = 55  # in microns

ncols = len(sorted_samples)
nrows = 2

# Set the value for the scalebar
scalebar_val = 1000  # in microns
scalebar_units = 'µm'

# Keep track of which annotations have been plotted
annotations_plotted = []

# Initialize the figure
fig, axs = plt.subplots(nrows, ncols, figsize=(5, 5))

for n, sample in enumerate(sorted_samples):
    # First plot the H&E image on the top row
    ax = axs[0][n]
    adata = adata_dict[sample]
    # Get the subtype for the sample
    subtype = subtype_dict[sample]
    # Get the microns per pixel for the scalebar
    spot_diameter_fulres = adata.uns['spatial'][sample]['scalefactors']['spot_diameter_fullres']
    microns_per_pixel = spot_diameter_micron / spot_diameter_fulres
    # Now divide by hires scalefactor
    scalefactor = adata.uns['spatial'][sample]['scalefactors']['tissue_hires_scalef']
    microns_per_pixel /= scalefactor
    # Plot the H&E image
    sc.pl.spatial(adata_dict[sample], color=None, size=spot_size[sample], ax=ax, alpha_img=1, 
                crop_coord=coords[sample], title=f'{sample}\n({subtype})', show=False)
    ax.set_title(ax.get_title(), fontweight='bold')
    # Add the scalebar
    # For the first plot, have a label
    if n == 0:
        scale_loc = 'none'
        frameon = True
    else:
        scale_loc = 'none'
        frameon = False
    scalebar = ScaleBar(microns_per_pixel, units=scalebar_units, location='lower right', box_alpha=0.5, 
                        fixed_value=scalebar_val, frameon=frameon, scale_loc=scale_loc, pad=0.3, 
                        height_fraction=0.02, font_properties={"size": 18})
    if n == 0:
        scalebar.label = '1 mm'
    ax.add_artist(scalebar)
    # Now plot the pathology annotations on the second row
    ax = axs[1][n]
    sc.pl.spatial(adata_dict[sample], color=merged_column, size=spot_size[sample], ax=ax,
                  crop_coord=coords[sample], title='', palette=color_dict, alpha_img=0.7,
                  legend_loc=None, show=False)
    # Add annotations plotted
    annotations_plotted += list(adata_dict[sample].obs[merged_column].dropna().unique())

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

# --- Add a pathology legend to the second row last column plot

target_ax = axs[1][-1]
annotations_plotted = list(set(annotations_plotted))  # Remove duplicates

# Sort according to plot_order
annotations_plotted = [cat for cat in ordered_columns_major if cat in annotations_plotted]

# Create a legend for the pathology annotations
handles = []
for cat in annotations_plotted:
    if cat in color_dict:
        color = color_dict[cat]
        handle = Line2D([0], [0], marker='o', color='none', markerfacecolor=color, label=cat, 
                        markersize=12, alpha=1, markeredgewidth=0)
        handles.append(handle)

# Clean the labels for the legend
legend_labels = [label.capitalize().replace('_', ' ') for label in annotations_plotted]

# Have some labels in uppercase
uppercase = ['tls', 'lvi']

for i, label in enumerate(legend_labels):
    if label.lower() in uppercase:
        legend_labels[i] = label.upper()

target_ax.legend(handles=handles, labels=legend_labels, loc='center left', bbox_to_anchor=(1.02, 1), 
                 fontsize=20, handlelength=0.7, frameon=False)

# --- Set ylabel for the first column

labels = ['H&E', 'Pathology\nannotations']

for i in range(nrows):
    axs[i][0].set_ylabel(labels[i], fontsize=25, fontweight='bold')

# --- Add title

# Calculate the suptitle position
title_x = (axs[0][0].get_position().x0 + axs[0][-1].get_position().x1) / 2
title_y = axs[0][0].get_position().y1 + 0.3  # Adjust y position slightly above the first row

plt.suptitle(f'Selected samples - Pathology annotations overview (scalebar={scalebar_val}{scalebar_units})',
             y=title_y, x=title_x, fontweight='bold')

# --- Save the plot

# Create output directory
out_dir = os.path.join(res_dir, 'figures')
os.makedirs(out_dir, exist_ok=True)

filename = "Fig_1D.pdf"

# Save with high resolution
plt.savefig(os.path.join(out_dir, filename), bbox_inches='tight')
plt.close('all')

print(f"Figure saved to {os.path.join(out_dir, filename)}")
print("Script successfully completed")
