#!/usr/bin/env python3

# --- Aim
# Overview matrix of which data modalities were profiled for each patient,
# shaded by the number of tissue sections profiled per modality (0 / 1 / 2).
# Section counts are derived directly from the per-modality sample sheets
# (one row per profiled section), so the cohort composition is the single
# source of truth and no separate resource file is required.
# This creates Fig. 1c

# --- Define environment
# env: spatial-pilot

# --- Load packages

import argparse
import os

import matplotlib as mpl
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgb
import numpy as np

import pandas as pd
import seaborn as sns

mpl.rcParams['pdf.fonttype'] = 42


# --- Argparse arguments

parser = argparse.ArgumentParser(description="Plot the data-modalities-per-patient matrix (Fig. 1c)",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument("-s", "--samples",
    nargs="*",  # 0 or more values expected => creates a list
    type=str,
    help="List of samples (unused; the matrix always shows the full cohort)")
parser.add_argument("-r", "--res_dir", help="Path to results directory (optional)")
parser.add_argument("--h5ad_dir", help="Path to directory for .h5ad output files")  # Not necessary for this script
parser.add_argument("--project_dir", help="Project directory (root)")
parser.add_argument("--sample_sheet", help="Path to sample sheet CSV file (unused)")

args = vars(parser.parse_args())

# Arguments
res_dir = args["res_dir"]
project_dir = args["project_dir"]

config_dir = os.path.join(project_dir, 'config')
resource_dir = os.path.join(project_dir, 'resources')


# --- Clinical subtype per patient

tissue_metadata = pd.read_csv(os.path.join(resource_dir, 'tissue_metadata.csv'))

# One row per patient (primary tumors; lymph-node rows are a separate tissue)
primary = tissue_metadata[tissue_metadata['tissue_type'] == 'Primary tumor']
subtype_dict = dict(zip(primary['donor_id'], primary['clinical_subtype']))
sample_list = list(primary['donor_id'])


# --- Count sections profiled per patient per modality

# Each row of a per-modality sample sheet is one profiled section, so counting
# rows per patient gives the number of sections (0 / 1 / 2) shown in the matrix.
visium = pd.read_csv(os.path.join(config_dir, 'sample_sheet_visium.csv'))
visium_hd = pd.read_csv(os.path.join(config_dir, 'sample_sheet_visium_hd.csv'))
lymph_node = pd.read_csv(os.path.join(config_dir, 'sample_sheet_visium_lymph_nodes.csv'))
xenium = pd.read_csv(os.path.join(config_dir, 'sample_sheet_xenium.csv'))

presence_df = pd.DataFrame(index=sample_list)
presence_df['visium_6.5mm'] = visium[visium['slide_size'] == '6.5 mm'].groupby('donor_id').size()
presence_df['visium_11mm'] = visium[visium['slide_size'] == '11 mm'].groupby('donor_id').size()
presence_df['visium_hd'] = visium_hd.groupby('donor_id').size()
presence_df['xenium_prime'] = xenium.groupby('donor_id').size()
presence_df['lymph_node'] = lymph_node.groupby('donor_id').size()
presence_df = presence_df.fillna(0).astype(int)


# --- Read in colors per subtype

color_dir = os.path.join(project_dir, 'resources', 'plotting_colors')
color_file = 'SCSubtype_colors.csv'

color_lut = pd.read_csv(os.path.join(color_dir, color_file))
color_lut = color_lut.set_index('subtype')
color_dict = color_lut['color_adj'].to_dict()

# Make a rowcolors df (clinical subtype only)
row_df = pd.DataFrame.from_dict(subtype_dict, orient='index', columns=['subtype'])
row_df['subtype'] = row_df['subtype'].map(color_dict)


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


# --- Function for beautifying heatmap

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

# Position legend function
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

# Function to add legend
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
        bbox_transform=fig.transFigure,  # <-- This is the key fix!
    )
    legend.get_title().set_fontweight('bold')
    legend._legend_box.align = "left"
    # Return updated vertical position for next stacked legend
    return legend

# Function to position the axs in the y direction
def position_axs_y(axs, offset=0.05):
    for n, ax in enumerate(axs):
        if n == 0:
            ax_pos_prev = ax.get_position()
        else:
            ax_pos = ax.get_position()
            ax.set_position((ax_pos.x0, ax_pos_prev.y0 - ax_pos.height - offset, ax_pos.width, ax_pos.height))
            ax_pos_prev = ax.get_position()


# --- Prepare for saving

out_dir = os.path.join(res_dir, 'figures')
os.makedirs(out_dir, exist_ok=True)


# --- Plot: presence matrix recording the number of sections profiled per sample

# Make sure to drop any duplicate indexes
presence_df = presence_df.loc[~presence_df.index.duplicated(keep='first')]

# Prepare for splitting the modalities for the plot
# Define the modality split (row order as in the figure: lymph node last)
modality_split = {
    'Visium SD': ['visium_6.5mm', 'visium_11mm'],
    'Visium HD': ['visium_hd'],
    'Xenium prime': ['xenium_prime'],
    'Lymph nodes': ['lymph_node'],
}

# Also make a dict for the ylabels
modality_labels = {
    'visium_6.5mm': 'Visium (6.5 mm)',
    'visium_11mm': 'Visium (11 mm)',
    'visium_hd': 'Visium HD',
    'xenium_prime': 'Xenium Prime 5K',
    'lymph_node': 'Matching lymph node',
}

data = presence_df.T

# Define row splits (row indices you want in each subplot)
row_splits = []

for modality, modalities in modality_split.items():
    # Keep the modalities present in the data, in the order defined above
    indices = [m for m in modalities if m in data.index]
    if indices:
        row_splits.append(indices)

# Prepare a colormatrix from the row_df
row_df = row_df.loc[data.columns]

# Create a color matrix for the clinical subtype strip
color_matrix = np.array([to_rgb(color) for color in row_df['subtype']]).reshape(1, -1, 3)
rgb_matrix = color_matrix

# Some parameters
spine_width = 1
offset = 0.01

# Plotting
fig, axs = plt.subplots(nrows=len(row_splits), ncols=1, figsize=(4, 4))

for i, rows in enumerate(row_splits):
    ax = axs[i]
    sns.heatmap(data.loc[rows], ax=ax, cmap='Blues', vmin=-0.1, vmax=1.1,
                cbar=False, xticklabels=True, yticklabels=True, lw=1, linecolor='white',)
    beautify_heatmap(ax, space_per_square=0.06, color_offset=0.02, spine_width=spine_width, width_to_height_ratio=0.6,)
    axs[i].set_ylabel('')
    axs[i].set_xlabel('')
    # Remove xlabels for all but the last plot
    if i < len(row_splits) - 1:
        axs[i].set_xticklabels([])
        # Also remove the xticks
        axs[i].tick_params(axis='x', which='both', bottom=False, top=False)

# Position the plots
position_axs_y(axs, offset=offset)

# Update the ylabels
for ax in axs:
    yticklabels = [tick.get_text() for tick in ax.get_yticklabels()]
    cleaned_labels = [modality_labels.get(label, label) for label in yticklabels]
    ax.set_yticklabels(cleaned_labels)

# Add the subtype color matrix on top
# Get the position of the first heatmap axis
ax_pos = axs[0].get_position()  # Get position of the first heatmap axis

# Create an ax for the subtype color matrix
new_ax = plt.gcf().add_axes([ax_pos.x0, ax_pos.y1 + offset, ax_pos.width, ax_pos.height/len(axs[0].get_yticks())*len(row_df.columns)],)
new_ax.imshow(rgb_matrix, aspect="auto",)
new_ax.set_xticks([])
new_ax.set_yticks([])

# Draw grid lines
new_ax.set_xticks(np.arange(-0.5, len(data.columns), 1), minor=True)
new_ax.set_yticks(np.arange(-0.5, len(data.index), 1), minor=True)
new_ax.grid(which="minor", color="white", linestyle='-', linewidth=1)
new_ax.tick_params(which="minor", bottom=False, left=False)

for spine in new_ax.spines.values():
    spine.set_visible(True)
    spine.set_linewidth(spine_width)

new_ax.set_ylim(-0.5, len(row_df.columns)-0.5)  # Adjust y-limits to match the number of rows

# Add yticklabel on the right side
# Turn on the right y-axis ticks
new_ax.yaxis.tick_right()

new_ax.set_yticks(np.arange(len(row_df.columns)))
new_ax.set_yticklabels(['Clinical subtype'], fontsize=10)

# Remove all minor ticks
new_ax.tick_params(axis='y', which='minor', left=False, right=False)

# Add the new_ax to the axs list for legend positioning
fig = plt.gcf()
current_axs = fig.get_axes()

# Add the color legend
right_edge, top_edge = position_legend(current_axs)

right_edge += offset  # Add some space to the right edge for the legend

add_color_legend(fig, labels=sorted_subtypes, color_dict=color_dict, title='Clinical subtype',
                 x=right_edge, y=top_edge, offset=0, spacing_per_label=0.05, loc='upper left')

# Add title
new_ax.set_title('Data modalities', fontsize=10, fontweight='bold')

# Calculate position for plot title
title_x = (ax_pos.x0 + ax_pos.x1) / 2
title_y = new_ax.get_position().y1 + 0.1

plt.suptitle('Overview of data modalities per tissue', y=title_y, x=title_x, fontsize=12, fontweight='bold')


# --- Save the output

filename = 'Fig_1C.pdf'
plt.savefig(os.path.join(out_dir, filename), bbox_inches='tight')
plt.close('all')


# --- Finish

print(f"Figure saved to {os.path.join(out_dir, filename)}")
print("Script successfully completed")
