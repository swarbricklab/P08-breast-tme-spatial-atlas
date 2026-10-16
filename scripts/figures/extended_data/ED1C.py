#!/usr/bin/env python3

# --- Aim
# Make a presence matrix showing which pathology annotation categories
# are present in which samples, split by annotation layer
# This creates Extended Data Fig. 1c

# --- Define environment
# env: spatial-pilot

# --- Load packages

import os
import sys
import json
import argparse
import re

import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.patches as mpatches
from matplotlib.colors import to_rgb, to_hex
import matplotlib.colors as mcolors

mpl.rcParams['pdf.fonttype'] = 42

# --- Argparse arguments

parser = argparse.ArgumentParser(description="Plot pathology annotation presence matrix",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument("-s", "--samples",
    nargs="*",  # 0 or more values expected => creates a list
    type=str,
    help="List of samples")
parser.add_argument("-r", "--res_dir",
                    help="Path to results directory (optional)")
parser.add_argument("--h5ad_dir", help="Path to directory for .h5ad output files")
parser.add_argument("--project_dir", help="Project directory (root)")
parser.add_argument("--sample_sheet", help="Path to sample sheet CSV file")

args = vars(parser.parse_args())

# Arguments
h5ad_dir = args["h5ad_dir"]
res_dir = args["res_dir"]
sample_list = args["samples"]
project_dir = args["project_dir"]

print(f'Samples in python: {sample_list}')


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


# --- Load pathology annotations

path_layers = ['epithelial', 'immune', 'stroma', 'additional_features']

pathology_dir = os.path.join(project_dir, 'data', 'pathology_annotations', 'visium')

annotation_dict = {}

for sample in sample_list:
    filename = f'{sample}_pathology_annotations.csv'
    path_to_file = os.path.join(pathology_dir, filename)
    if os.path.exists(path_to_file):
        pathology_anno = pd.read_csv(path_to_file, index_col=0)
        # Read annotations directly from the individual layer columns
        annotation_dict[sample] = pathology_anno[path_layers].copy()
        print(f"Loaded pathology annotations for {sample}")
    else:
        print(f"Warning: No pathology annotations found for {sample}")


# --- Extract a presence / absence matrix

def get_presence_matrix_and_value_lists(sample_dfs):
    # To collect unique values per column
    unique_values_per_column = {}
    # First pass: gather all unique values per column across all samples
    for df in sample_dfs.values():
        for col in df.columns:
            values = df[col].dropna().unique()
            unique_values_per_column.setdefault(col, set()).update(values)
    # Sort and freeze the unique values per column
    for col in unique_values_per_column:
        unique_values_per_column[col] = sorted(unique_values_per_column[col])
    # Prepare final columns
    final_columns = []
    for col, values in unique_values_per_column.items():
        for val in values:
            final_columns.append(f"{col}::{val}")
    # Build binary matrix per sample
    presence_data = {}
    for sample_id, df in sample_dfs.items():
        presence_row = dict.fromkeys(final_columns, 0)
        for col in df.columns:
            values = df[col].dropna().unique()
            for val in values:
                key = f"{col}::{val}"
                if key in presence_row:
                    presence_row[key] = 1
        presence_data[sample_id] = presence_row
    # Final DataFrame
    presence_df = pd.DataFrame.from_dict(presence_data, orient="index")
    presence_df.index.name = "sample_id"
    return presence_df, unique_values_per_column

presence_df, unique_vals = get_presence_matrix_and_value_lists(annotation_dict)


# --- Read in colors per subtype

color_dir = os.path.join(project_dir, 'resources', 'plotting_colors')
color_file = 'SCSubtype_colors.csv'

color_lut = pd.read_csv(os.path.join(color_dir, color_file))
color_lut = color_lut.set_index('subtype')
color_dict = color_lut['color_adj'].to_dict()

# Make a rowcolors df
row_df = pd.DataFrame.from_dict(subtype_dict, orient='index', columns=['subtype'])
row_df['subtype'] = row_df['subtype'].map(color_dict)

# Also add capture area size
row_df['capture_area'] = sample_lut['slide_size']
color_dict['6.5 mm'] = '#ebf4fb'
color_dict['11 mm'] = '#295687'

row_df['capture_area'] = row_df['capture_area'].map(color_dict)


# --- Sort the samples by subtype

# Define a plotting order for subtypes
subtype_order = ['Luminal A', 'Luminal B', 'Luminal B (HER2+)', 'HER2+', 'TNBC']
all_subtypes = list(set(subtype_dict[sample] for sample in sample_list if sample in subtype_dict))
sorted_subtypes = sorted(all_subtypes, key=lambda x: subtype_order.index(x) if x in subtype_order else float('inf'))

# Dictionary to store the lists of samples per type
samples_by_subtype = {}

# Group sample IDs by their type
for sample in sample_list:
    subtype = subtype_dict[sample]
    if subtype not in samples_by_subtype:
        samples_by_subtype[subtype] = []
    samples_by_subtype[subtype].append(sample)

# Sort the lists of sample IDs per type
for subtype in samples_by_subtype:
    samples_by_subtype[subtype].sort()

# Extract sample order by subtype
sorted_samples = []
for key in sorted_subtypes:
    sorted_samples += samples_by_subtype[key]

# Re-order the df
presence_df = presence_df.loc[sorted_samples]


# --- Sort the pathology annotations

annotation_dir_base = os.path.join(project_dir, 'resources', 'pathology_annotations')

# Load the desired categories per layer (to get the desired order)
filename = 'pathology_categories_per_layer.json'
path_to_file = os.path.join(annotation_dir_base, filename)

with open(path_to_file, 'r') as json_file:
    annotations_of_interest_dict = json.load(json_file)

ordered_major_cat_dict = {}

for path_layer in path_layers:
    sorted_keys = annotations_of_interest_dict[path_layer]
    # Sort values in unique_vals[path_layer] according to sorted_keys
    ordered_list = sorted(
        unique_vals[path_layer],
        key=lambda x: sorted_keys.index(x) if x in sorted_keys else float('inf')
    )
    # Check for missing values
    missing = [val for val in unique_vals[path_layer] if val not in sorted_keys]
    ordered_list += missing
    # Add the sorted list to the dictionary
    ordered_major_cat_dict[path_layer] = ordered_list

def sort_df_columns(df, ordered_major_cat_dict):
    def column_sort_key(col):
        if "::" not in col:
            return (float('inf'), col)  # Push unmatched columns to the end
        key, value = col.split("::", 1)
        if key in ordered_major_cat_dict:
            try:
                return (list(ordered_major_cat_dict.keys()).index(key), 
                        ordered_major_cat_dict[key].index(value))
            except ValueError:
                return (list(ordered_major_cat_dict.keys()).index(key), float('inf'))
        return (float('inf'), col)
    sorted_cols = sorted(df.columns, key=column_sort_key)
    return df[sorted_cols]

presence_df = sort_df_columns(presence_df, ordered_major_cat_dict)


# --- Load colors for pathology layers

color_file = 'path_colors.csv'
color_lut_path = pd.read_csv(os.path.join(color_dir, color_file), index_col=0)

for path_layer in path_layers:
    color_dict[path_layer] = color_lut_path.loc[path_layer, 'color_adj']


# --- Function to clean labels

def clean_label(label, remove_prefix=False):
    if "::" in label:
        prefix, suffix = label.split("::", 1)
        prefix = prefix.strip()
        suffix = suffix.strip()
        if suffix.startswith(prefix):
            suffix = suffix[len(prefix):].lstrip("_ ").strip()
        label = f"{prefix}: {suffix}"
        if remove_prefix:
            label = suffix.capitalize()
    label = label.replace("_", " ").strip()
    label = re.sub(r'\s+', ' ', label)  # collapse multiple spaces
    label = label.capitalize()
    if label.lower() in ['tls', 'lvi']:
        label = label.upper()
    return label


# --- Helper function to lighten colors

def lighten_color(color, amount=0.5):
    r, g, b = to_rgb(color)
    return to_hex([1 - (1 - c) * (1 - amount) for c in (r, g, b)])


# --- Prepare data for plotting

# Split the DataFrame into separate DataFrames for each layer
layer_dfs = {}
for layer in path_layers:
    layer_cols = [col for col in presence_df.columns if col.startswith(layer + "::")]
    layer_dfs[layer] = presence_df[layer_cols]

# Extract a colormatrix for each layer
layer_color_matrices = {}
for layer, df in layer_dfs.items():
    df = df.T  # Transpose for correct orientation
    # Convert hex color values to RGBA for the heatmap
    color_matrix = np.empty(df.shape + (4,), dtype=float)  # RGBA
    for row_idx, row in enumerate(df.columns):
        for col_idx, col in enumerate(df.index):
            val = df.loc[col, row]
            if val == 1:
                hex_color = color_dict.get(layer, "#cccccc")
            else:
                hex_color = lighten_color(color_dict.get(layer, "#cccccc"), amount=0.8)
            color_matrix[col_idx, row_idx] = mcolors.to_rgba(hex_color)
    layer_color_matrices[layer] = color_matrix


# --- Prepare the rowdf for plotting

row_df = row_df.T
row_df = row_df[presence_df.index]  # Reorder row_df to match presence_df index

color_matrix_list = [
    np.array([to_rgb(color) for color in row_df[column]])  # shape: (n_samples, 3)
    for column in row_df.columns
]

rgb_matrix = np.stack(color_matrix_list, axis=0)  # shape: (n_categories, n_samples, 3)
rgb_matrix = np.transpose(rgb_matrix, (1, 0, 2))  # shape: (n_samples, n_categories, 3)
rgb_matrix = np.flip(rgb_matrix, axis=0)  # flips the categories


# --- Plotting functions

def beautify_heatmap(ax, title=None, x_rot=90, spine_width=1.5, spine_color='black',
                     label_size=10, fontsize_title=13, space_per_square=0.07,
                     width_to_height_ratio=1, color_offset=0.01):
    # Remove axis labels
    ax.set_ylabel('')
    ax.set_xlabel('')
    # Rotate x tick labels
    ax.tick_params(axis='x', labelrotation=x_rot, labelsize=label_size)
    ax.tick_params(axis='y', labelsize=label_size)
    # Set spines
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(spine_width)
        spine.set_edgecolor(spine_color)
    # Adjust axis position based on number of ticks
    xlabels = ax.get_xticklabels()
    ylabels = ax.get_yticklabels()
    fig = ax.figure
    ax_pos = ax.get_position()
    new_width = space_per_square * width_to_height_ratio * len(xlabels)
    new_height = space_per_square * len(ylabels)
    ax.set_position((ax_pos.x0, ax_pos.y0, new_width, new_height))
    # Move ylabels to the right of the plot
    ax.yaxis.set_label_position("right")
    ax.yaxis.tick_right()
    # Make sure ylabels are horizontal
    ax.tick_params(axis='y', rotation=0)
    # Set title
    if title:
        ax.set_title(title, fontsize=fontsize_title, fontweight='bold')
    return ax

def position_axs_y(axs, offset=0.05):
    for n, ax in enumerate(axs):
        if n == 0:
            pass
        else:
            ax_pos = ax.get_position()
            ax_pos_prev = axs[n-1].get_position()
            ax.set_position((ax_pos.x0, ax_pos_prev.y0 - ax_pos.height - offset, ax_pos.width, ax_pos.height))

def position_legend(axs):
    """
    Given an array (or list) of axes, find the rightmost right_edge and the highest top_edge
    in figure coordinates, for robust legend placement.
    Returns (max_right_edge, max_top_edge)
    """
    # Accept a single axis as well
    if not isinstance(axs, (list, np.ndarray)):
        axs = [axs]
    fig = axs[0].figure
    fig.canvas.draw()
    right_edges = []
    top_edges = []
    for ax in axs:
        bbox = ax.get_tightbbox(fig.canvas.get_renderer())
        inv_trans = fig.transFigure.inverted()
        bbox_fig = inv_trans.transform(bbox)
        right_edges.append(bbox_fig[1][0])
        top_edges.append(bbox_fig[1][1])  # y1 is the top in figure coordinates
    max_right_edge = max(right_edges)
    max_top_edge = max(top_edges)
    return max_right_edge, max_top_edge

def add_color_legend(fig, labels, color_dict, title='', x=0, y=0, offset=0.01, spacing_per_label=0.05, loc='center left'):
    # Create legend handles (colored boxes with labels)
    handles = [
        mpatches.Patch(facecolor=color_dict[label], label=label, edgecolor='black', linewidth=0.8)
        for label in labels if label in color_dict
    ]
    # Estimate height needed
    est_height = spacing_per_label * len(handles)
    bbox_anchor = (x + offset, y)
    # Add the legend
    legend = fig.legend(
        handles=handles,
        title=title,
        loc=loc,
        bbox_to_anchor=bbox_anchor,
        borderaxespad=0,
        frameon=False,
        handlelength=0.8,
        bbox_transform=fig.transFigure,
    )
    legend.get_title().set_fontweight('bold')
    legend._legend_box.align = "left"
    # Return updated vertical position for next stacked legend
    return legend


# --- Create the plot

# Some parameters
spine_width = 1
offset = 0.01
space_per_square = 0.06
width_to_height_ratio = 0.8

# Plotting
fig, axs = plt.subplots(nrows=len(path_layers), ncols=1, figsize=(4, 4))

for i, layer in enumerate(path_layers):
    ax = axs[i]
    matrix = layer_color_matrices[layer]
    ax.imshow(matrix, aspect="auto",)
    # Draw grid lines
    data = layer_dfs[layer]
    ax.set_xticks(np.arange(-0.5, len(data.index), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(data.columns), 1), minor=True)
    ax.grid(which="minor", color="white", linestyle='-', linewidth=1)
    ax.tick_params(which="minor", bottom=False, left=False)
    # Set the yticklabels
    yticklabels = data.columns.tolist()
    cleaned_labels = [clean_label(label, remove_prefix=True) for label in yticklabels]
    # Set yticks accordingly
    ax.set_yticks(np.arange(len(cleaned_labels)))
    ax.set_yticklabels(cleaned_labels)
    # Set the xticklabels
    xticklabels = data.index.tolist()
    ax.set_xticks(np.arange(len(xticklabels)))
    ax.set_xticklabels(xticklabels, rotation=90)
    # Beautify the heatmap
    beautify_heatmap(ax, space_per_square=space_per_square, color_offset=0.02, spine_width=spine_width, width_to_height_ratio=width_to_height_ratio,)
    axs[i].set_ylabel('')
    axs[i].set_xlabel('')
    # Set a ylabel on the left of the plot
    ax.yaxis.set_label_position('left')
    ax.set_ylabel(layer.capitalize().replace('_', ' '), fontweight='bold')
    ax.tick_params(axis='y', which='minor', left=False, right=False)
    # Remove xlabels for all but the last plot
    if i < len(path_layers) - 1:
        axs[i].set_xticklabels([])
        # Also remove the xticks
        axs[i].tick_params(axis='x', which='both', bottom=False, top=False)

# Position the plots
position_axs_y(axs, offset=offset)

# Create an ax for the subtype color matrix
ax_pos = axs[0].get_position()

new_ax = plt.gcf().add_axes([ax_pos.x0, ax_pos.y1 + offset, ax_pos.width, ax_pos.height/len(axs[0].get_yticks())*len(row_df.index)],)
new_ax.imshow(rgb_matrix, aspect="auto",)
new_ax.set_xticks([])
new_ax.set_yticks([])

# Draw grid lines
new_ax.set_xticks(np.arange(-0.5, len(row_df.columns), 1), minor=True)
new_ax.set_yticks(np.arange(-0.5, len(row_df.index), 1), minor=True)
new_ax.grid(which="minor", color="white", linestyle='-', linewidth=1)
new_ax.tick_params(which="minor", bottom=False, left=False)

for spine in new_ax.spines.values():
    spine.set_visible(True)
    spine.set_linewidth(spine_width)

new_ax.set_ylim(-0.5, len(row_df.index)-0.5)  # Adjust y-limits to match the number of rows

# Add yticklabel on the right side
# Turn on the right y-axis ticks
new_ax.yaxis.tick_right()

new_ax.set_yticks(np.arange(len(row_df.index)))
new_ax.set_yticklabels(['Capture area', 'Clinical subtype'],)

# Remove all minor ticks
new_ax.tick_params(axis='y', which='minor', left=False, right=False)

# Add the new_ax to the axs list for legend positioning
fig = plt.gcf()
current_axs = fig.get_axes()

# Add the color legend
right_edge, top_edge = position_legend(current_axs)

right_edge += offset  # Add some space to the right edge for the legend

first_legend = add_color_legend(fig, labels=sorted_subtypes, color_dict=color_dict, title='Clinical subtype',
                                 x=right_edge, y=top_edge, offset=0, spacing_per_label=0.05, loc='upper left')

# Get the bounding box of the first legend in figure coordinates
fig.canvas.draw()  # Ensure the legend is drawn
legend_bbox = first_legend.get_window_extent(fig.canvas.get_renderer())
inv_trans = fig.transFigure.inverted()
legend_bbox_fig = inv_trans.transform(legend_bbox)
legend_bottom = legend_bbox_fig[0][1]  # y0 is the bottom in figure coordinates

# Add a small vertical offset if desired
vertical_spacing = 0
second_legend_y = legend_bottom - vertical_spacing

# Add the second legend just below the first
second_legend = add_color_legend(fig, x=right_edge, y=second_legend_y, labels=['6.5 mm', '11 mm'], color_dict=color_dict, title='Capture area',
                 offset=0, spacing_per_label=0.05, loc='upper left')

# Get the bounding box of the first legend in figure coordinates
fig.canvas.draw()  # Ensure the legend is drawn
legend_bbox = second_legend.get_window_extent(fig.canvas.get_renderer())
inv_trans = fig.transFigure.inverted()
legend_bbox_fig = inv_trans.transform(legend_bbox)
legend_bottom = legend_bbox_fig[0][1]  # y0 is the bottom in figure coordinates

# Add a small vertical offset if desired
vertical_spacing = 0
third_legend_y = legend_bottom - vertical_spacing

# Add the third legend just below the second
cleaned_path_layers = [clean_label(layer, remove_prefix=True) for layer in path_layers]

for label in path_layers:
    cleaned_label = clean_label(label, remove_prefix=True)
    color_dict[cleaned_label] = color_dict.get(label, "#cccccc")  # Ensure all layers have a color

third_legend = add_color_legend(fig, x=right_edge, y=third_legend_y, labels=cleaned_path_layers, color_dict=color_dict, title='Annotation layer',
                 offset=0, spacing_per_label=0.05, loc='upper left')

# Set a title for the top ax
new_ax.set_title('Pathology annotations', fontweight='bold')

# Calculate position for plot title
ax_pos = new_ax.get_position()
title_x = (ax_pos.x0 + ax_pos.x1) / 2
title_y = new_ax.get_position().y1 + 0.2

plt.suptitle('Overview of available pathology annotations per sample', y=title_y, x=title_x, fontsize=12, fontweight='bold')


# --- Save the plot

# Create output directory
out_dir = os.path.join(res_dir, 'figures', 'extended_data')
os.makedirs(out_dir, exist_ok=True)

filename = "Fig_ED1C.pdf"

# Save with high resolution
plt.savefig(os.path.join(out_dir, filename), bbox_inches='tight')
plt.close('all')

print(f"Figure saved to {os.path.join(out_dir, filename)}")
print("Script successfully completed")
