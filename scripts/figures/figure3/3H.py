#!/usr/bin/env python3

# --- Aim
# After calculating composition vs distance per sample - plot the average composition across samples


# --- Define environment
# env: vis-hd

# --- Load packages
import argparse
import math
import os
from datetime import datetime

import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import scanpy.plotting.palettes as sc_palettes

# Import custom functions
from sophiesfunctions.color_functions import generate_color_variations

mpl.rcParams['pdf.fonttype'] = 42


# --- Argparse arguments
parser = argparse.ArgumentParser(description="Import Visium unfiltered data & plot",
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
sample_list = args["samples"]

print(f'Samples in python: {sample_list}')


# --- Extract some info from the sample sheet

# Use sample_sheet argument if provided, else default
if args.get("sample_sheet"):
    sample_sheet_path = args["sample_sheet"]
else:
    sample_dir = os.path.join(project_dir, 'config')
    sample_sheet_file = 'sample_sheet_visium_hd.csv'
    sample_sheet_path = os.path.join(sample_dir, sample_sheet_file)

sample_sheet = pd.read_csv(sample_sheet_path)
sample_lut = sample_sheet.set_index('section_id')


# --- Also get some tissue metadata

resource_dir = os.path.join(project_dir, 'resources')
metadata_file = 'tissue_metadata.csv'

tissue_metadata = pd.read_csv(os.path.join(resource_dir, metadata_file))

# Merge the two DataFrames on the common column 'tissue_id'
merged_df = pd.merge(sample_sheet, tissue_metadata, on='tissue_id')
# Create a sample - subtype dictionary
subtype_dict = dict(zip(merged_df['section_id'], merged_df['clinical_subtype']))


# --- Load miniatlas classifications

decon_level = 'minor'

# To make it easier to visualize trends by lineage
class_dir = os.path.join(project_dir, 'resources', 'cell_types', 'miniatlas')
class_file = 'miniatlas_subset_classifications.csv'

class_df = pd.read_csv(os.path.join(class_dir, class_file))

# Get endothelial classes
endothelial = list(class_df[class_df['miniatlas_major'] == 'Endothelial'][f'miniatlas_{decon_level}'].unique())
endothelial = sorted(endothelial)

# And a dict for the endothelial names
endo_name_dict = {
    'CXCL12': 'Arterial',
    'RGS5': 'Capillary',
    'ACKR1': 'Venous',
    'LYVE1': 'Lymphatic',
}


# --- Load composition results

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
    return most_recent_dir_path

# Define which endothelial type want to look at
endo_of_interest = 'Endothelial ACKR1'
out_name = endo_of_interest.split(' ')[-1]

# Specify path
TLS_group_dir = os.path.join(project_dir, 'results', '04_vascular_microenvironment', '03_compare_HEV_to_other_venous', '01_compare_microenvironment', 'csv_files')
most_recent_dir_path = find_most_recent(TLS_group_dir)

# Grouped by whether overlap with HEVs
group_types = ['no_overlap', 'overlap']

# Load per sample
avg_composition = {}

for cell_type in group_types:
    res_dict = {}
    for sample in sample_list:
        filename = f'{sample}_{out_name}_{cell_type}_composition_vs_distance.csv'
        res_dict[sample] = pd.read_csv(os.path.join(most_recent_dir_path, filename), index_col=0)
        # Also grab the index order to make sure it stays the same
        index_order = res_dict[sample].index
    # Compute the mean DataFrame
    average_composition_df = pd.concat(res_dict.values(), axis=0).groupby(level=0).mean()
    cell_types = list(average_composition_df.columns)
    # Reorder distance bins
    average_composition_df = average_composition_df.loc[index_order, cell_types]
    # After averaging:
    average_composition_df = average_composition_df.div(average_composition_df.sum(axis=1), axis=0)
    avg_composition[cell_type] = average_composition_df

average_composition_df.sum(axis=1)


# --- Prepare for plotting

plt.rcParams.update({'font.size': 8})


# --- Get some colors for the pathology annotations

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
all_cats = list(average_composition_df.columns)

colors = sc_palettes.default_20[:len(annotation_labels)]

top_cat_colors = {}
for n, cat in enumerate(annotation_labels):
    top_cat_colors[cat] = colors[n]

# Get some color variations of the original colors (so variations of the same annotation are a similar color)
color_dict = {}

# And also define a plot order
plot_order = []

for top_cat, labels in annotation_labels.items():
    base_color = top_cat_colors[top_cat]
    # Check how many of the labels are found in the object
    num_cols = sum(item in set(all_cats) for item in labels)
    # Get the matching items
    matching_items = [item for item in labels if item in set(all_cats)]
    plot_order += matching_items
    if num_cols > 4:
        step_size=0.1
    else:
        step_size=0.2
    variations = generate_color_variations(base_color, num_cols, step_size=step_size)
    for n, item in enumerate(matching_items):
        color_dict[item] = variations[n]


# --- Plot 1: Cell type composition across distance

# Disable the grid
plt.rcParams['axes.grid'] = False

nrows = 1
ncols = len(group_types)

fig, axs = plt.subplots(ncols=ncols, nrows=nrows, figsize=(5, 5), sharey=True)

for n, cell_type in enumerate(group_types):
    ax = axs[n]
    composition = avg_composition[cell_type][plot_order[::-1]]
    composition.plot.area(ax=ax, color=color_dict, alpha=0.8, linewidth=0)
    if n == (ncols-1):
        # Move legend
        legend = ax.legend()
        legend.set_title('')
        # Automatically calculate how many columns
        ncol_legend = math.ceil(len(plot_order)/10)
        ax.legend(handlelength=0.8, loc='center left', bbox_to_anchor=(1, 0.5), alignment='left', frameon=False, markerscale=2, reverse=True, ncols=ncol_legend)
    else:
        # Remove legend otherwise
        ax.get_legend().remove()
    # Make sure all xticklabels are displayed
    ax.set_xticklabels(ax.get_xticklabels(), rotation=90)  # Rotate x-tick labels if needed
    ax.set_xticks(range(len(composition)))  # Ensure all x-ticks are placed
    ax.set_xticklabels(composition.index)  # Set the x-tick labels to match the DataFrame index
    # Adjust lims
    ax.set_ylim(0, 1)
    ax.set_xlim(0, ax.get_xlim()[1]-0.75)
    # Set xlabel
    ax.set_xlabel('Distance (µm)')
    ax.set_ylabel('Proportion')
    # Set title
    title = f'{endo_name_dict[out_name]}'
    if cell_type == 'overlap':
        title = f'HEV'
    ax.set_title(title, fontweight='bold')

# Resize the plots
plot_width = 0.25

for ax in axs:
    ax_pos = ax.get_position()
    ax.set_position([ax_pos.x0, ax_pos.y0, plot_width, plot_width])
    # And annotate only every second xtick
    for n, label in enumerate(ax.get_xticklabels()):
        if n % 2 != 0:
            label.set_visible(False)

# Adjust positioning
for n, ax in enumerate(axs):
    ax_pos = ax.get_position()
    # Shift the subplot to the desired location + set the width so it's equal for all plots
    if n == 0:
        ax.set_position([ax_pos.x0, ax_pos.y0, ax_pos.width, ax_pos.height])
    else:
        new_pos = ax_pos_prev.x1 + 0.02
        ax.set_position([new_pos, ax_pos.y0, ax_pos.width, ax_pos.height])
    ax_pos_prev = ax.get_position()

# Adjust the legend texts
legend = ax.get_legend()
for text in legend.get_texts():
    # Capitalize the first word 
    text.set_text(text.get_text().replace('_', ' ').capitalize())
    # Leave LVI capitalized
    if text.get_text() == 'Lvi':
        text.set_text('LVI')

# Position the title
title_x = (axs[0].get_position().x0 + axs[-1].get_position().x1)/2
title_y = axs[0].get_position().y1 + 0.1

plt.suptitle(f'All samples - Average pathology annotations across {endo_of_interest} vessel distance bins (overlapping with HEV vs not overlapping)', y=title_y, fontweight='bold', x=title_x)


# --- Save output

out_dir = os.path.join(res_dir, 'figures')
os.makedirs(out_dir, exist_ok=True)

filename = 'Fig_3H.pdf'
plt.savefig(os.path.join(out_dir, filename), bbox_inches='tight')


# --- Finish

print("Script successfully completed")
