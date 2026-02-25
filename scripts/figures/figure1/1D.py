#!/usr/bin/env python3

# --- Aim
# Plot Figure 1D: barplot of sample composition by pathology annotations & deconvolution results
# Shows stacked bars for each sample with proportions of pathology categories and cell types
# This creates Figure 1D

# --- Define environment
# env: spatial-pilot

# --- Load packages

import json
import os
from collections import defaultdict
import argparse

import pandas as pd
import numpy as np

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.lines as mlines
import scanpy.plotting.palettes as sc_palettes

# Custom functions
from sophiesfunctions.save_multi_image import save_multi_image
from sophiesfunctions.color_functions import generate_color_variations

mpl.rcParams['pdf.fonttype'] = 42


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
parser.add_argument("--path_combo_of_interest", help="Pathology layer combinations to be plotted (e.g. epithelial_immune)")
parser.add_argument("--sample_sheet", help="Path to sample sheet CSV file")

args = vars(parser.parse_args())

# Arguments
h5ad_dir = args["h5ad_dir"]
res_dir = args["res_dir"]
sample_list = args["samples"]
project_dir = args["project_dir"]

print(f'Samples in python: {sample_list}')


# --- Extract some info from the sample sheet

# Use sample_sheet argument if provided, else default
if args.get("sample_sheet"):
    sample_sheet_path = args["sample_sheet"]
else:
    sample_dir = os.path.join(project_dir, 'config')
    sample_sheet_file = 'sample_sheet_visium.csv'
    sample_sheet_path = os.path.join(sample_dir, sample_sheet_file)

sample_sheet = pd.read_csv(sample_sheet_path)


# --- Also get some tissue metadata

resource_dir = os.path.join(project_dir, 'resources')
metadata_file = 'tissue_metadata.csv'

tissue_metadata = pd.read_csv(os.path.join(resource_dir, metadata_file))

# Merge the two DataFrames on the common column 'tissue_id'
merged_df = pd.merge(sample_sheet, tissue_metadata, on='tissue_id')

# Create a sample - subtype dictionary
subtype_dict = dict(zip(merged_df['section_id'], merged_df['clinical_subtype']))


# --- Load the annotations

# Define annotation path
reworked_anno_dir = os.path.join(project_dir, 'data', 'pathology_annotations', 'visium')

anno_dict = {}

for sample in sample_list:
    annotation_file = f'{sample}_pathology_annotations.csv'
    annotation_path = os.path.join(reworked_anno_dir, annotation_file)
    if os.path.exists(annotation_path):
        single_layer = pd.read_csv(annotation_path, index_col=0)
        anno_dict[sample] = single_layer
    else:
        print(f'No annotations found for {sample} - skipping')


# --- Load the desired categories
# To get the desired order
annotation_dir_base = os.path.join(project_dir, 'resources', 'pathology_annotations')
filename = 'pathology_categories_per_layer.json'
path_to_file = os.path.join(annotation_dir_base, filename)

with open(path_to_file, 'r') as json_file:
                annotations_of_interest_dict = json.load(json_file)

# Get all unique annotation labels per layer that are present in the data
present_labels_per_layer = {}
all_layers = ['epithelial', 'immune', 'stroma', 'additional_features']

for layer in all_layers:
    labels = set()
    for sample in sample_list:
        labels.update(anno_dict[sample][layer].dropna().unique())
    present_labels_per_layer[layer] = labels

# Now, for each layer, order the present labels according to annotations_of_interest_dict
ordered_labels_per_layer = {}
for layer in all_layers:
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

# Split the adatas by subtype for plotting
# Take care to only include the adatas that have annotations for all layers
subtype_dict = {key: value for key, value in subtype_dict.items() if key in sample_list}

all_subtypes = set(subtype_dict.values())

# Define a desired order for plotting the subtypes
subtype_order = ['Luminal A', 'Luminal B (ER+)', 'Luminal B (HER2+)', 'HER2+', 'TNBC']
sorted_subtypes = sorted(all_subtypes, key=lambda x: subtype_order.index(x) if x in subtype_order else float('inf'))

ncols = len(sorted_subtypes)


# --- Read in colors per subtype

color_dir = os.path.join(project_dir, 'resources', 'plotting_colors')
color_file = 'SCSubtype_colors.csv'

color_lut = pd.read_csv(os.path.join(color_dir, color_file))
color_lut = color_lut.set_index('subtype')


# --- Prepare for plotting the collapsed annotations

# Collect counts for each sample
overwrite_counts = []
sample_names = []

for sample in sample_list:
    counts = anno_dict[sample]['single_layer'].dropna().value_counts()
    overwrite_counts.append(counts)
    sample_names.append(sample)

# Combine into a single DataFrame
overwrite_counts_df = pd.DataFrame(overwrite_counts, index=sample_names).fillna(0).astype(int)

# If you want to order columns by all possible labels (e.g. from ordered_labels_per_layer or a flat list)
all_labels = sorted(set().union(*[df['single_layer'].dropna().unique() for df in anno_dict.values()]))
overwrite_counts_df = overwrite_counts_df.reindex(columns=all_labels, fill_value=0)

print(overwrite_counts_df)

# Convert to proportions
overwrite_proportions_df = overwrite_counts_df.div(overwrite_counts_df.sum(axis=1), axis=0)


# --- Functions

# -- Step 1: prepare for plotting: gridspec and prepare the figure
def prepare_barplot(barplot_dict_annotation_major, ncols=ncols, groupby='sample_id', all_subtypes=sorted_subtypes, height=5, width=10, wspace=0.1):
    # Calculate the width of each subtype subplot based on how many samples per subtype
    grid_spec = []
    for subtype in all_subtypes:
        grid_spec.append(len(barplot_dict_annotation_major[subtype].index))
    fig, axs = plt.subplots(ncols=ncols, figsize=(width, height), sharey=False, gridspec_kw={'width_ratios': grid_spec, 'wspace': wspace})
    # If ncols = 1: need to reshape to make axs iterable
    if ncols == 1:
        axs = np.reshape(axs, 1)
    return axs, grid_spec

# -- Step 2: actual barplot

# Define a function for creating barplots
def plot_stacked_barplot(barplot_dict, axs, grid_spec, title=None, all_subtypes=sorted_subtypes, offset=0.02, ncols_legend=1, ylabel='', stacked=True,
                        plot_subtype_cols=True, col_height=0.03, color_dict=None, col_width=0.03):
    for n, (subtype, ax) in enumerate(zip(all_subtypes, axs)):
        barplot_dict[subtype].sort_index().plot(kind='bar', stacked=stacked, legend=False, ax=ax, width=0.9, grid=False,
                                                        color=color_dict, linewidth=0.5)
        for i, bar in enumerate(ax.patches):  # Adjust the colors of the bar
            # Get the RGBA values of the fill color
            fill_color = bar.get_facecolor()
            fill_color = list(fill_color)
            darkened_color = fill_color.copy()
            for j in range(len(darkened_color)-1):
                darkened_color[j] *= 0.8
            fill_color[-1] = 0.8 # make the fill color a bit transparent
            bar.set_color(tuple(fill_color))
            bar.set_edgecolor(tuple(darkened_color))
        if title is None:
            ax.set_title(subtype, fontweight='bold')
            # If subtype is Luminal B with (xxx), put an enter before the brackets
            if subtype.startswith('Luminal B ('):
                ax.set_title(subtype.replace(' (', '\n('), fontweight='bold')
        else:
            ax.set_title(title, fontweight='bold')
        # move x-axis tick labels to vertical
        xlabels = ax.get_xticklabels()
        ax.set_xticklabels(xlabels, rotation=90, ha='center')
        # remove the x-axis label
        ax.set_xlabel('')
        # Set the ylabel
        ax.set_ylabel(ylabel)
        if n > 0:
            ax.set_ylabel('')
                # Adjust the positions of the subplots
        ax_pos = axs[n].get_position()
        # Shift the subplot to the desired location + set the width so it's equal for all plots
        if n == 0:
            axs[n].set_position([ax_pos.x0, ax_pos.y0, col_width*grid_spec[n], ax_pos.height])
        else:
            new_pos = ax_pos_prev.x0 + ax_pos_prev.width + offset
            axs[n].set_position([new_pos, ax_pos.y0, col_width*grid_spec[n], ax_pos.height])
        ax_pos_prev = axs[n].get_position()
    # Add a legend
    plt.legend(title_fontproperties={'weight':'bold'}, bbox_to_anchor=(1, 0.5), loc='center left', ncol=ncols_legend,
                            frameon=False, alignment='left', handlelength=0.7)
    # Add subtype colors at the top
    if plot_subtype_cols:
        for n, subtype in enumerate(all_subtypes):
            # Add subtype color at top & bottom
            # Add colors for subtypes
            color = color_lut.loc[subtype, 'color']
            width = grid_spec[n] + 1
            xlim = axs[n].get_xlim()
            # Increase ylim a bit
            axs[n].set_ylim(0, 1.02 * axs[n].get_ylim()[1])
            ylim = axs[n].get_ylim()
            range_y = ylim[1] - ylim[0]
            height = col_height * range_y
            # Markings at top
            rect = patches.Rectangle((xlim[0], ylim[1]-height), width=width, height=height, facecolor=color, alpha=0.7)
            axs[n].add_patch(rect)

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
        variations = generate_color_variations(color, num_cols, step_size=step_size)
        color_variations[key] = variations
    # Create a new dictionary mapping values from dict1 to values from dict2
    color_dict = {val: color_variations[key][idx] for key, values in category_dict.items() for idx, val in enumerate(values)}
    # Extract values from the dictionary based on the keys
    plotting_colors = [color_dict[key] for key in ordered_categories]
    return color_dict, plotting_colors, top_cat_colors

# Extract a flat list of ordered categories
flat_ordered_categories = [item for sublist in ordered_labels_per_layer.values() for item in sublist]

# Re-extract colors
color_dict, plotting_colors, top_cat_colors = create_plotting_colors(ordered_labels_per_layer, flat_ordered_categories, mixed=False, major_categories=True)

# Reverse the order for plotting
flat_ordered_categories = flat_ordered_categories[::-1]

# Create a list of handles and labels for the legend
handles = []
labels = []
for cat, color in color_dict.items():
    # Make sure not to include non-epithelial
    if cat != 'other':
        handles.append(mlines.Line2D([], [], color=color, marker='o', alpha=1, lw=0, mew=0, markersize=10))
        labels.append(cat)

# Make a barplotdict per subtype from the overwrite_proportions_df
barplot_dict_annotation_major = {}

for subtype in sorted_subtypes:
    # Get a list of samples belonging to this subtype
    samples = [sample for sample, sub in subtype_dict.items() if sub == subtype]
    # Order the list
    samples = sorted(samples)
    # Extract a df
    df = overwrite_proportions_df.loc[samples]
    # Order the columns
    df = df[flat_ordered_categories]
    barplot_dict_annotation_major[subtype] = df

# Set fontsize
plt.rcParams.update({'font.size': 12})


# --- Also make the same plots for decon results

level = 'major'
decon_dir = os.path.join(project_dir, 'results', '02_cell_type_mapping', 'deconvolution', 'visium', f'{level}_level')


# --- Read in colors for the barplot

# Also extract a row_color_df
# Load miniatlas classifications
classifications_dir = os.path.join(project_dir, 'resources', 'cell_types', 'miniatlas')
classifications_file = 'miniatlas_subset_classifications.csv'

classifications = pd.read_csv(os.path.join(classifications_dir, classifications_file))

# Read in the results & add to anndata

decon_dict = {}

for sample in sample_list:
    decon_filename = f'{sample}_{level}_q05_cell_abundance_w_sf.csv'
    # Load decon results
    decon_results = pd.read_csv(os.path.join(decon_dir, decon_filename), index_col=0)
    decon_dict[sample] = decon_results

# Calculate per sample proportions of cells & extract single df for barplot
decon_proportions_dict = {}
for sample in sample_list:
    decon_df = decon_dict[sample]
    # Calculate total counts per cell type
    total_counts = decon_df.sum(axis=0)
    # Calculate proportions from total counts
    total_cells = total_counts.sum()
    proportions = total_counts / total_cells
    # Store proportions in a dictionary
    decon_proportions_dict[sample] = proportions
    # Convert the dictionary to a DataFrame
    decon_proportions_df = pd.DataFrame(decon_proportions_dict).T

classifications.index = classifications[f'miniatlas_{level}']

# Filter for cell types that are in the decon results
classifications = classifications[classifications[f'miniatlas_{level}'].isin(decon_proportions_df.columns)]

cell_types = decon_proportions_df.columns.tolist()
lineage_colors_file = 'miniatlas_lineage_colors.csv'
miniatlas_major_colors_file = f'miniatlas_{level}_colors.csv'

lineage_colors = pd.read_csv(os.path.join(color_dir, lineage_colors_file))
miniatlas_major_colors = pd.read_csv(os.path.join(color_dir, miniatlas_major_colors_file))

# Make some lookup tables
color_dict_decon = dict(zip(miniatlas_major_colors.iloc[:,0], miniatlas_major_colors['color']))


# --- Define a plotting order by lineage

# Desired lineage order
lineage_order = ['Epithelial', 'Immune', 'Stromal']

# Step 1: Convert lineage column to categorical with the predefined order
classifications['lineage'] = pd.Categorical(classifications['lineage'], categories=lineage_order, ordered=True)

# Step 2: Sort by lineage and category
df_sorted = classifications.reset_index(drop=True).sort_values(by=['lineage', f'miniatlas_{level}'])

# Drop any NaN values
df_sorted = df_sorted.dropna(subset=[f'miniatlas_{level}'])

# Step 3: Group by lineage and get sorted unique categories within each lineage
grouped = df_sorted.groupby('lineage')[f'miniatlas_{level}'].unique()

# Step 4: Sort categories alphabetically within each lineage
sorted_categories = {lineage: sorted(categories) for lineage, categories in grouped.items()}

# Flatten the sorted categories into a single unnested list
unnested_categories = [category for categories in sorted_categories.values() for category in categories]

# But remove any items from unnested_categories that are not in the df
unnested_categories_filtered = [cat for cat in unnested_categories if cat in decon_proportions_df.columns]

# Add any missing cell types at the end
for cell_type in decon_proportions_df.columns:
    if cell_type not in unnested_categories_filtered:
        unnested_categories_filtered.append(cell_type)

# If the cell type is not in the lookup table, use a default color
# Figure out how many cell types are not in the color_dict_decon
not_in_dict = [cell_type for cell_type in unnested_categories_filtered if cell_type not in color_dict_decon]

for n, cell_type in enumerate(not_in_dict):
    if len(not_in_dict) < 21:
        # Use a default color from the sc_palettes
        color_dict_decon[cell_type] = sc_palettes.default_20[n]
    else:
        # Use larger palette
        color_dict_decon[cell_type] = sc_palettes.godsnot_102[n]


# --- Prepare for barplot decon results

# Make a barplotdict per subtype from the overwrite_proportions_df
barplot_dict_decon = {}

for subtype in sorted_subtypes:
    # Get a list of samples belonging to this subtype
    samples = [sample for sample, sub in subtype_dict.items() if sub == subtype]
    # Only keep samples present in decon_proportions_df
    available_samples = [s for s in samples if s in decon_proportions_df.index]
    if not available_samples:
        print(f"Warning: No decon results for subtype {subtype}")
        continue
    # Order the list
    available_samples = sorted(available_samples)
    # Extract a df
    df = decon_proportions_df.loc[available_samples]
    # Order the columns according to unnested_categories
    df = df[unnested_categories_filtered[::-1]]
    barplot_dict_decon[subtype] = df

# Figure out which cells are epithelial
epithelial = classifications[classifications['lineage'] == 'Epithelial'][f'miniatlas_{level}'].unique().tolist()

# Figure out the index of epithelial in the unnested categories
epithelial_index = [i for i, cat in enumerate(unnested_categories_filtered) if cat in epithelial]

# Figure out which cells are non-epithelial
non_epithelial = classifications[classifications['lineage'] != 'Epithelial'][f'miniatlas_{level}'].unique().tolist()

# Filter the present cell types
non_epithelial = [cat for cat in unnested_categories_filtered if cat in non_epithelial]

# Make a new barplot dict for non-epithelial cells
barplot_dict_decon_non_epithelial = {}

for subtype in sorted_subtypes:
    # Get a list of samples belonging to this subtype
    samples = [sample for sample, sub in subtype_dict.items() if sub == subtype]
    # Only keep samples present in decon_proportions_df
    available_samples = [s for s in samples if s in decon_proportions_df.index]
    if not available_samples:
        print(f"Warning: No decon results for subtype {subtype}")
        continue
    # Order the list
    available_samples = sorted(available_samples)
    # Extract a df
    df = decon_proportions_df.loc[available_samples]
    # Order the columns according to non_epithelial
    df = df[non_epithelial[::-1]]
    # Recalculate proportions
    df = df.div(df.sum(axis=1), axis=0)
    # Add to the barplot dict
    barplot_dict_decon_non_epithelial[subtype] = df


# --- Plot 1: Stacked barplots showing sample composition by pathology annotations & decon

# -- First the pathology annotations

col_width = 0.025

# Plotting
axs, grid_spec = prepare_barplot(barplot_dict_annotation_major, height=2.8,)

plot_stacked_barplot(axs=axs, barplot_dict=barplot_dict_annotation_major, grid_spec=grid_spec, ylabel='Proportion of spots', 
                            ncols_legend=2, offset=0.005, col_height=0.05, color_dict=color_dict, col_width=col_width)

# Reverse the legend
legend = axs[-1].get_legend()
handles = legend.legendHandles[::-1]
labels = flat_ordered_categories[::-1]
# Clean the labels for the legend
labels = [label.capitalize().replace('_', ' ') for label in labels]

# Have some labels in uppercase
uppercase = ['tls', 'lvi']

for i, label in enumerate(labels):
    if label.lower() in uppercase:
        labels[i] = label.upper()

# Automatically calculate ncols for legend
# Want 8 per column
ncols_legend = max(1, len(labels) // 8)

axs[-1].legend(handles, labels, title_fontproperties={'weight':'bold'}, bbox_to_anchor=(1, 0.5), loc='center left', ncol=ncols_legend,
                            frameon=False, alignment='left', handlelength=0.7)

# Center the title
title_x = (axs[0].get_position().x0 + axs[-1].get_position().x1)/2

# Adjust the margins of the plos
for ax in axs:
    xlim = ax.get_xlim()
    ax.set_xlim(xlim[0] + 0.1, xlim[1] - 0.1)

# Remove yticks after first
for ax in axs:
    if ax != axs[0]:
        ax.set_yticks([])
    # Also remove the xticklabels
    ax.set_xticklabels([])

# --- Then the decon results

offset = 0.05
height_fraction = 0.95

fig = plt.gcf()  # Get the current figure

# Create new axs under the previous ones
new_axs = []

for n, ax in enumerate(axs):
    # Get the position of the current ax
    pos = ax.get_position()
    # Create a new ax below the current one
    new_ax = fig.add_axes([pos.x0, pos.y0 - pos.height*height_fraction - offset, pos.width, pos.height*height_fraction])
    new_axs.append(new_ax)

plot_stacked_barplot(axs=new_axs, barplot_dict=barplot_dict_decon, grid_spec=grid_spec, ylabel='Estimated frequencies', 
                            ncols_legend=1, offset=0.005, col_height=0.05, color_dict=color_dict_decon, col_width=col_width,
                            plot_subtype_cols=False)

# Hash the epithelial cell types
for ax in new_axs:
    # Group patches by x-position (i.e., each stacked bar group)
    bars_by_x = defaultdict(list)
    for patch in ax.patches:
        x = patch.get_x()
        bars_by_x[x].append(patch)
    # For each stack (x-position), style the top 2 bars
    for x_val, bars in bars_by_x.items():
        # Sort bars from top to bottom by y (top = higher y)
        bars_sorted = sorted(bars, key=lambda p: p.get_y() + p.get_height(), reverse=True)
        for n, patch in enumerate(bars_sorted):  # top 2 bars
            if n in epithelial_index:
                # patch.set_linewidth(1.5)
                patch.set_hatch('//')

# Adjust the margins of the plos
for ax in new_axs:
    xlim = ax.get_xlim()
    ax.set_xlim(xlim[0] + 0.1, xlim[1] - 0.1)
    ax.set_ylim(0, 1.02)  # Increase ylim a bit
    # And remove the title
    ax.set_title('')
    # And remove the xticklabels
    ax.set_xticklabels([])

# Remove yticks after first
for ax in new_axs:
    if ax != new_axs[0]:
        ax.set_yticks([])

# Reverse the legend
legend = new_axs[-1].get_legend()
handles = legend.legendHandles[::-1]
labels = unnested_categories_filtered

# Automatically calculate ncols for legend
# Want 8 per column
ncols_legend = max(1, len(labels) // 8)

new_axs[-1].legend(handles, labels, title_fontproperties={'weight':'bold'}, bbox_to_anchor=(1, 0.5), loc='center left', ncol=ncols_legend,
                            frameon=False, alignment='left', handlelength=0.7,)

# --- Then the decon results for just the TME

# Create new axs under the previous ones
new_axs_tme = []

for n, ax in enumerate(new_axs):
    # Get the position of the current ax
    pos = ax.get_position()
    # Create a new ax below the current one
    new_ax = fig.add_axes([pos.x0, pos.y0 - pos.height - offset, pos.width, pos.height])
    new_axs_tme.append(new_ax)

plot_stacked_barplot(axs=new_axs_tme, barplot_dict=barplot_dict_decon_non_epithelial, grid_spec=grid_spec, ylabel='Estimated frequencies', 
                            ncols_legend=1, offset=0.005, col_height=0.05, color_dict=color_dict_decon, col_width=col_width,
                            plot_subtype_cols=False)

# Adjust the margins of the plos
for ax in new_axs_tme:
    xlim = ax.get_xlim()
    ax.set_xlim(xlim[0] + 0.1, xlim[1] - 0.1)
    ax.set_ylim(0, 1.02)  # Increase ylim a bit
    # And remove the title
    ax.set_title('')

# Remove yticks after first
for ax in new_axs_tme:
    if ax != new_axs_tme[0]:
        ax.set_yticks([])

# Reverse the legend
legend = new_axs_tme[-1].get_legend()
handles = legend.legendHandles[::-1]
labels = non_epithelial

# Automatically calculate ncols for legend
# Want 8 per column
ncols_legend = max(1, len(labels) // 8)

new_axs_tme[-1].legend(handles, labels, title_fontproperties={'weight':'bold'}, bbox_to_anchor=(1, 0.5), loc='center left', ncol=ncols_legend,
                            frameon=False, alignment='left', handlelength=0.7)

plt.suptitle(f'All samples - {all_layers} annotations major categories collapsed to single layer', y=1.2, x=title_x, fontweight='bold')

# Draw a line from the bottom of the last row to the top of the second row
# Just to the left of the first column of plots
top_ax = new_axs[0]
bottom_ax = new_axs_tme[0]

# Draw a line above the first to third axs in the first row
# Get subplot positions in figure coordinates
start_box = top_ax.get_position()
end_box = bottom_ax.get_position()

# Define rectangle dimensions
y1 = start_box.y1
y0 = end_box.y0
height = abs(y1 - y0)

# Define x0 position (left side of the rectangle)
fig = plt.gcf()  # Get the current figure
fig.canvas.draw()
bbox = top_ax.get_tightbbox(fig.canvas.get_renderer())
inv_trans = fig.transFigure.inverted()
bbox_fig = inv_trans.transform(bbox)
left_edge = bbox_fig[0][0]  # x0 is the left in figure coordinates

x0 = left_edge - 0.02  # bottom of the rectangle (just above the subplot)
width = 0.005            # set width manually

# Create rectangle
highlight_rect = mpl.patches.Rectangle(
    (x0, y0), width, height,
    transform=fig.transFigure,
    facecolor='black',
    edgecolor='black',
    linewidth=0,
    alpha=1
)

# Add to figure
fig.patches.append(highlight_rect)

# Calculate midpoint y
y_mid = y0 + height / 2
x_text = x0 - 0.005 # slightly above the line

# Add vertical text
fig.text(
    x_text, y_mid,
    "Deconvolution",
    ha='right', va='center',
    color='black', fontweight='bold',
    rotation=90,
)

# Same for the top row
ax = axs[0]

# Define rectangle dimensions
y1 = ax.get_position().y1
y0 = ax.get_position().y0
height = abs(y1 - y0)

# Create rectangle
highlight_rect = mpl.patches.Rectangle(
    (x0, y0), width, height,
    transform=fig.transFigure,
    facecolor='black',
    edgecolor='black',
    linewidth=0,
    alpha=1
)

# Add to figure
fig.patches.append(highlight_rect)

# Calculate midpoint y
y_mid = y0 + height / 2
x_text = x0 - 0.005 # slightly above the line

# Add vertical text
fig.text(
    x_text, y_mid,
    "Pathology annotations",
    ha='right', va='center',
    color='black', fontweight='bold',
    rotation=90,
)


# --- Save the plots

out_dir = os.path.join(res_dir, 'figures')
os.makedirs(out_dir, exist_ok=True)

filename = f'Fig_1D.pdf'

save_multi_image(os.path.join(out_dir, filename))


# --- Finish

print("Script successfully completed")
