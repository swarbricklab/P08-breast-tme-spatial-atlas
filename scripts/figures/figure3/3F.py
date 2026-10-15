#!/usr/bin/env python3

# --- Aim
# Plot Figure 3F: PVL (perivascular-like) subtype composition across distance from vessels
# Shows how different PVL subtypes associate with different vessel types based on proximity

# --- Define environment
# env: vis-hd

# --- Load packages
import argparse
import os

import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import sys
import numpy as np

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

# Import function script as a module
sys.path.append(os.path.join(project_dir, 'scripts', 'figures',))
import dotplot_functions


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

# Sort according to desired order
plot_order = ['CXCL12', 'RGS5', 'ACKR1', 'LYVE1']

# plot_order contains substrings, endothelial contains full names
sorted_endothelial = []
for key in plot_order:
    matches = [e for e in endothelial if key in e]
    sorted_endothelial.extend(matches)

# Optionally, check for any endothelial types not matched (to append at the end)
remaining = [e for e in endothelial if e not in sorted_endothelial]
sorted_endothelial.extend(remaining)

endothelial = sorted_endothelial

# Define decon level to use
decon_level = 'subset'

# Also load the colors for the endothelial types
colors_file = os.path.join(project_dir, 'resources', 'plotting_colors', 'endothelial_colors.csv')
colors_df = pd.read_csv(colors_file)

# Make a color_dict
color_dict = colors_df.set_index('endothelial_type')['color'].to_dict()
    
# Get PVL classes for this decon level
PVL_classes = list(class_df[class_df['miniatlas_major'] == 'PVL'][f'miniatlas_{decon_level}'].unique())
PVL_classes = sorted(PVL_classes)


# --- Load composition results

# Specify path
TLS_group_dir = os.path.join(res_dir, '04_vascular_microenvironment', '02_identify_all_vasculature', '04_assess_PVL_association', 'csv_files')
most_recent_dir_path = TLS_group_dir

# Load per sample
res_dict = {}

for sample in sample_list:
    filename = f'{sample}_PVL_composition_vs_distance.csv.gz'
    res_dict[sample] = pd.read_csv(os.path.join(most_recent_dir_path, filename), index_col=0)

# Extract a single df
res_df = pd.concat(
    [df.assign(sample_id=key) for key, df in res_dict.items()],
    ignore_index=True
)


# --- Prepare for plotting

plt.rcParams.update({'font.size': 10})

# Prepend dist_bin_ to endothelial keys
endothelial_keys = [f'dist_bin_{cell_type}' for cell_type in endothelial]

# Enable gridlines behind the plot for all plots
sns.set_style("ticks",{'axes.grid' : True})
plt.rcParams['axes.axisbelow'] = True


# --- Define a dict to update the PVL names to something more interpretable

# And edit the names on top
name_dict = {
    'PVL Differentiated s3': 'PVL s3\n(Contractile)',
    'PVL Immature s1': 'PVL s1\n(Inflammatory)',
    'PVL_Immature s2': 'PVL s2\n(Hypoxia response)',
}

# And a dict for the endothelial names
endo_name_dict = {
    'Endothelial CXCL12': 'Arterial (CXCL12+)',
    'Endothelial RGS5': 'Capillary (RGS5+)',
    'Endothelial ACKR1': 'Venous (ACKR1+)',
    'Endothelial Lymphatic LYVE1': 'Lymphatic (LYVE1+)',
}

# --- Plot 5: Similar plots as above put plotting each endothelial subtype on its own axis

# Also limit the bins to make the plots smaller
bin_limit = 200

nrows = len(endothelial)
ncols = len(PVL_classes)

fig, axs = plt.subplots(ncols=ncols, nrows=nrows, figsize=(2, 2.2), squeeze=False)

for n, PVL_class in enumerate(PVL_classes):
    # Extract data to plot
    cell_df = res_df[endothelial_keys + [PVL_class]]
    # Melt to long format: one row per (bin, subtype, PVL value)
    df_long = pd.DataFrame()
    for endo_type in endothelial_keys:
        temp = cell_df[[endo_type, PVL_class]].copy()
        temp = temp.rename(columns={endo_type: 'distance_bin', PVL_class: 'PVL_value'})
        temp['endothelial_subtype'] = endo_type.replace('dist_bin_', '')
        df_long = pd.concat([df_long, temp], axis=0)
    # Drop rows with missing bins or PVL values
    df_long = df_long.dropna(subset=['distance_bin', 'PVL_value'])
    # Order bins numerically
    bin_order = sorted(df_long['distance_bin'].unique(), key=lambda x: int(x.split('-')[0]))
    df_long['distance_bin'] = pd.Categorical(df_long['distance_bin'], categories=bin_order, ordered=True)
    # Limit the bins
    bin_order = [b for b in bin_order if int(b.split('-')[1]) <= bin_limit]
    df_long = df_long[df_long['distance_bin'].isin(bin_order)]
    # Plot each subtype separately
    ymax = df_long.groupby(['distance_bin', 'endothelial_subtype'])['PVL_value'].mean().max() * 1.1  # Initialize ymax for the current PVL_class
    ymin = df_long.groupby(['distance_bin', 'endothelial_subtype'])['PVL_value'].mean().min() * 0.9  # Initialize ymin for the current PVL_class
    for i, endo_subtype in enumerate(endothelial):
        ax = axs[i, n]  # Access the correct subplot
        sub_df = df_long[df_long['endothelial_subtype'] == endo_subtype]
        # Calculate mean per bin
        sub_df = sub_df.groupby('distance_bin')['PVL_value'].mean().reset_index()
        sub_df.set_index('distance_bin').plot.area(ax=ax, color=color_dict[endo_subtype], legend=False, alpha=0.5)
        # Make sure all xticklabels are displayed
        ax.set_xticklabels(ax.get_xticklabels(), rotation=90)  # Rotate x-tick labels if needed
        ax.set_xticks(range(len(bin_order)))  # Ensure all x-ticks are placed
        ax.set_xticklabels(bin_order)  # Set the x-tick labels to match the DataFrame index
        if PVL_class in name_dict:
            PVL_class_name = name_dict[PVL_class]
        else:
            PVL_class_name = PVL_class
        ax.set_title(f'{PVL_class_name}', fontweight='bold', fontsize=10)
        # # Rotate xtick labels 90
        # ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
        # Adjust labels
        ax.set_xlabel('Distance (µm)')
        ax.set_ylabel('')
        # Adjust x-limits to fit the bins
        ax.set_xlim(0, len(bin_order) - 1)
        # Adjust lims
        ax.set_ylim(ymin, ymax)

# Adjust positioning & size
for n, ax in enumerate(axs.ravel()):
    ax_pos = ax.get_position()
    ax.set_position([ax_pos.x0, ax_pos.y0, len(bin_order)*0.08, ax_pos.height])

dotplot_functions.move_plots_x_y(axs, x_offset=0.25, y_offset=0.05, remove_ylabels=False)

# Move legend
ax = axs[0,-1]

# Adjust the endothelial names according to the dict
dotplot_functions.add_legend(ax=ax, ordered_list=endothelial, color_lut=color_dict, title='Vessel type', 
                            loc='upper left', y_offset=0, x_offset=1)

# Make the facecolor of the legend handles more transparent
for handle in ax.get_legend().legend_handles:
    facecolor = handle.get_facecolor()
    handle.set_facecolor(handle.get_facecolor()[:3] + (0.5,))  # Add transparency to the facecolor
    # Set the edgecolor to be the same as the facecolor
    handle.set_edgecolor(facecolor)
    # Set the linewidth to 2
    handle.set_linewidth(1.25)

# Also adjust the legend texts according to the dict
legend = ax.get_legend()
for text in legend.get_texts():
    if text.get_text() in endo_name_dict:
        text.set_text(endo_name_dict[text.get_text()])

# Get the left-most x-position among all subplots
lefts = [ax.get_position().x0 for ax in axs.flat]
left_min = min(lefts)

# Add shared y-label just to the left of all subplots
fig = plt.gcf()
ax = axs[0, 0]  # First subplot

# Draw the figure to ensure ticklabels are rendered
fig.canvas.draw()

# Get the renderer
renderer = fig.canvas.get_renderer()

# Get the bounding boxes of the yticklabels in figure coordinates
yticklabels = ax.get_yticklabels()
bboxes = [label.get_window_extent(renderer=renderer).transformed(fig.transFigure.inverted()) for label in yticklabels if label.get_text()]

# Find the minimum x0 (leftmost edge) among all yticklabels
if bboxes:
    min_x = min(bbox.x0 for bbox in bboxes)
    # Add a small offset to the left
    label_x = min_x - 0.02
else:
    # Fallback if no labels
    label_x = ax.get_position().x0 - 0.03

# Add the shared y-label at the calculated position
fig.text(label_x, 0.5, 'Average proportion', va='center', ha='right', rotation='vertical')

# Calculate x to center
title_x = (axs[0][0].get_position().x0 + axs[0][-1].get_position().x1)/2

plt.suptitle(f'All samples - PVL composition across distance bins (decon_level={decon_level})',
            y=1.2, fontweight='bold', x=title_x)


# --- Save the output

out_dir = os.path.join(res_dir, 'figures')
os.makedirs(out_dir, exist_ok=True)

filename = 'Fig_3F.pdf'
plt.savefig(os.path.join(out_dir, filename), bbox_inches='tight')

print("Script successfully completed")
