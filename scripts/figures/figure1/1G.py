#!/usr/bin/env python3

# --- Aim
# Plot Figure 1G: B-cell to T-cell ratios per clinical subtype
# Shows ratios across TLS/aggregate annotated spots
# This creates Figure 1G


# --- Define environment
# env: spatial-pilot

# --- Load packages

import argparse
import os
import pandas as pd
import numpy as np

import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches


# --- Argparse arguments

parser = argparse.ArgumentParser(description="Import Visium data & do basic QC",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument("--data_dir", help="Path to raw data")
parser.add_argument("-s", "--samples",
    nargs="*",  # 0 or more values expected => creates a list
    type=str,
    help="List of samples")
parser.add_argument("-r", "--res_dir",
                    help="Path to results directory (optional)")
parser.add_argument("--h5ad_dir", help="Path to directory for .h5ad output files") # Not necessary for this script
parser.add_argument("--gene_file", help="Path to a csv file with genes, genes stored in first column") # Not necessary for this script
parser.add_argument("--use_data", help="Which data to use (SCT or log_norm)") # Not necessary for this script
parser.add_argument("--project_dir", help="Project directory (root)")
parser.add_argument("--sample_sheet", help="Path to sample sheet CSV file")

args = vars(parser.parse_args())

# Arguments
data_dir = args["data_dir"]
h5ad_dir = args["h5ad_dir"]
sample_list = args["samples"]
project_dir = args['project_dir']


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


# --- Load the decon results

level = 'major'
decon_dir = os.path.join(project_dir, 'results', '02_cell_type_mapping', 'deconvolution', 'visium', f'{level}_level')

decon_dict = {}

for sample in sample_list:
    decon_filename = f'{sample}_{level}_q05_cell_abundance_w_sf.csv'
    # Load decon results
    decon_results = pd.read_csv(os.path.join(decon_dir, decon_filename), index_col=0)
    decon_dict[sample] = decon_results

all_cell_types = list(decon_dict[sample].columns)


# --- Also load pathology annotations

# Define annotation path
reworked_anno_dir = os.path.join(project_dir, 'data', 'pathology_annotations', 'visium')

anno_df_dict = {}

for sample in sample_list:
    annotation_file = f'{sample}_pathology_annotations.csv'
    annotation_path = os.path.join(reworked_anno_dir, annotation_file)
    if os.path.exists(annotation_path):
        path_anno = pd.read_csv(annotation_path, index_col=0)
        anno_df_dict[sample] = path_anno
        decon_dict[sample][path_anno.columns] = path_anno
    else:
        print(f'No annotations found for {sample} - skipping')


# --- Sort the samples per subtype

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

print(samples_by_subtype)

# Define a plotting order for subtypes
subtype_order = ['Luminal A', 'Luminal B (ER+)', 'Luminal B (HER2+)', 'HER2+', 'TNBC']
sorted_subtypes = sorted(samples_by_subtype.keys(), key=lambda x: subtype_order.index(x) if x in subtype_order else float('inf'))


# --- Calculate B/T ratios in TLS and aggregate annotated spots

annotations_of_interest = ['tls', 'aggregate']
annotations_of_interest = ['tls', 'aggregate']

b_t_ratio_annotated_dict = {}
n_annotated_spots_dict = {}

for sample, df in decon_dict.items():
    # Filter the DataFrame for specific annotations
    annotated_df = df[df['immune'].isin(annotations_of_interest)]
    n_annotated_spots = annotated_df.shape[0]
    if not annotated_df.empty:
        # Calculate the ratio of B cells to T cells
        total_b = annotated_df[['B-cells']].sum(axis=0).values[0]
        total_t = annotated_df[['T-cells']].sum(axis=0).values[0]
        b_t_ratio = total_b / total_t
        # Log2 transform the ratio
        if total_b > 0:
            b_t_ratio = np.log2(b_t_ratio)
        else:
            b_t_ratio = float('inf')
    else:
        b_t_ratio = float('nan')  # No annotated spots available
    # Store the results in the dictionaries
    b_t_ratio_annotated_dict[sample] = b_t_ratio
    n_annotated_spots_dict[sample] = n_annotated_spots/df.shape[0]

# Convert the dictionaries into DataFrames
b_t_ratio_annotated_df = pd.DataFrame.from_dict(b_t_ratio_annotated_dict, orient='index', columns=['B_T_ratio_annotated'])
n_annotated_spots_df = pd.DataFrame.from_dict(n_annotated_spots_dict, orient='index', columns=['n_annotated_spots'])

# Add subtype info
b_t_ratio_annotated_df['subtype'] = b_t_ratio_annotated_df.index.map(subtype_dict)

# Merge the two DataFrames
b_t_ratio_annotated_df = b_t_ratio_annotated_df.merge(n_annotated_spots_df, left_index=True, right_index=True)


# --- Prepare for plotting

sns.set_style("ticks")


# --- Get colors per subtype

color_dir = os.path.join(project_dir, 'resources', 'plotting_colors')
color_file = 'SCSubtype_colors.csv'

color_lut = pd.read_csv(os.path.join(color_dir, color_file))
color_lut = color_lut.set_index('subtype')

# Extract a color dict
color_dict = color_lut['color'].to_dict()


# --- Prepare for making some boxplots

# Properties for boxplot
whisker_props = {'color': 'black'}
median_props = {'color': 'black'}

alpha = 0.7

# Enable gridlines behind the plot for all plots
sns.set_style("ticks",{'axes.grid' : True})
plt.rcParams['axes.axisbelow'] = True

# Function to avoid repetition
def plot_seaborn_boxplot(df, ax, x, y, hue, plot_order, stripplot=True, title=None, gap=0, legend=False,
                                annotation_size=10, spot_size=1, box_width=0.2, palette=None, width=0.8,
                                ylabel=None, height=0.4):
    if stripplot:
        showfliers = False
    else:
        showfliers = True
    sns.boxplot(df, x=x, y=y, showfliers=showfliers, linewidth=1, whiskerprops=whisker_props, 
                medianprops=median_props, width=width, ax=ax, palette=palette,
                saturation=1, legend=legend, hue=hue, fliersize=1, order=plot_order, gap=gap)
    # Plot the individual dots on the boxplot
        # Plot dots if you want
    if stripplot:
        sns.stripplot(data=df, x=x, y=y, color="black", order=plot_order, jitter=True, size=spot_size, ax=ax)
    for i, bar in enumerate(ax.patches):  # Adjust the colors of the bar
        # Get the RGBA values of the fill color
        fill_color = bar.get_facecolor()
        fill_color = list(fill_color)
        darkened_color = fill_color.copy()
        for j in range(len(darkened_color)-1):
            darkened_color[j] *= 0.8
        fill_color[-1] = alpha # make the fill color a bit transparent
        bar.set_color(tuple(fill_color))
        bar.set_edgecolor(tuple(darkened_color))
    # Add title
    ax.set_title(title, fontweight='bold')
    # move x-axis tick labels to the left
    xlabels = ax.get_xticklabels()
    ax.set_xticklabels(xlabels, rotation=45, ha='right')
    # remove the x-axis label
    ax.set_xlabel('')
    # set ylabel
    if ylabel is None:
        ax.set_xlabel(y)
    else:
        ax.set_ylabel(ylabel)
    # Adjust the positions of the subplots
    ax_pos = ax.get_position()
    # Set the width so it's equal for all plots
    ax.set_position([ax_pos.x0, ax_pos.y0, box_width*len(plot_order), height])

# Function to move axs if overlap
def move_axs(axs, offset=0.01):
    for n, ax in enumerate(axs):
        ax_pos = ax.get_position()
        # Shift the subplot to the desired location
        if n == 0:
            ax.set_position([ax_pos.x0, ax_pos.y0, ax_pos.width, ax_pos.height])
        else:
            new_pos = ax_pos_prev.x0 + ax_pos_prev.width + offset
            ax.set_position([new_pos, ax_pos.y0, ax_pos.width, ax_pos.height])
        ax_pos_prev = ax.get_position()

# Function to increase the distance between rows of subplots
def increase_row_spacing(axs_by_row_dict, offset=0.05):
    """
    Adjust the vertical position of rows of subplots by adding an offset between rows.
    
    Parameters:
    - axs_by_row_dict: A dictionary where each key is a row number and each value is a list of axes in that row.
    - offset: The amount of vertical space to add between rows (default is 0.05).
    """
    for row_num, axs in axs_by_row_dict.items():
        for ax in axs:
            ax_pos = ax.get_position()
            if row_num == 0:
                # Keep the first row at the original position
                ax.set_position([ax_pos.x0, ax_pos.y0, ax_pos.width, ax_pos.height])
            else:
                # Shift the row downward by the cumulative offset
                new_y0 = ax_pos_prev.y0 - ax_pos.height - offset
                ax.set_position([ax_pos.x0, new_y0, ax_pos.width, ax_pos.height])
        ax_pos_prev = ax.get_position()

# Function to add count bars at top of plot
def add_count_bars(df, x, y, target_ax, alpha=0.7):
    fig = plt.gcf()
    # Position to target_ax
    target_ax_pos = target_ax.get_position()
    # Create a new axis above the target axis
    count_ax_height = target_ax_pos.height * 0.2
    count_ax = fig.add_axes([target_ax_pos.x0, target_ax_pos.y1 + 0.01, target_ax_pos.width, count_ax_height])
    # Plot the count bars
    sns.barplot(x=df[x], y=df[y], hue=df[x], ax=count_ax, color='lightgray', alpha=0.7, width=0.8,
                edgecolor='black', linewidth=1, order=sorted_subtypes, err_kws={'color': 'black', 'linewidth': 1}, capsize=0.3,
                palette=color_dict, saturation=1)
    for i, bar in enumerate(count_ax.patches):  # Adjust the colors of the bar
        # Get the RGBA values of the fill color
        fill_color = bar.get_facecolor()
        fill_color = list(fill_color)
        darkened_color = fill_color.copy()
        for j in range(len(darkened_color)-1):
            darkened_color[j] *= 0.8
        fill_color[-1] = alpha # make the fill color a bit transparent
        bar.set_color(tuple(fill_color))
        bar.set_edgecolor(tuple(darkened_color))
    # Get the title
    title = target_ax.get_title()
    count_ax.set_title(title, fontweight='bold')
    target_ax.set_title('')  # Remove the title from the target axis
    count_ax.set_ylabel('Prop.\nspots', fontsize=10)
    # Remove the xticks
    count_ax.set_xticks([])
    count_ax.set_xlabel('')
    # Increase the ylim a bit
    ylim = count_ax.get_ylim()
    count_ax.set_ylim(0, ylim[1] * 1.1)


# --- Plot 1: Horizontal boxplot of B/T cell ratio in annotated spots

# Columns to plot
title_dict = {
    'B_T_ratio_annotated': 'B/T cell ratio\n(TLS/aggregates)'
}

fig = plt.figure(figsize=(2, 2))
ax = plt.subplot(1, 1, 1)

# Create horizontal boxplot
sns.boxplot(b_t_ratio_annotated_df, y='subtype', x='B_T_ratio_annotated', showfliers=False, linewidth=1, 
            whiskerprops=whisker_props, medianprops=median_props, width=0.8, ax=ax, 
            palette=color_dict, saturation=1, legend=False, hue='subtype', order=sorted_subtypes, 
            orient='h', gap=0)

# Add stripplot
sns.stripplot(data=b_t_ratio_annotated_df, y='subtype', x='B_T_ratio_annotated', color="black", 
              order=sorted_subtypes, jitter=True, size=3, ax=ax)

# Adjust colors
for i, bar in enumerate(ax.patches):
    fill_color = bar.get_facecolor()
    fill_color = list(fill_color)
    darkened_color = fill_color.copy()
    for j in range(len(darkened_color)-1):
        darkened_color[j] *= 0.8
    fill_color[-1] = alpha
    bar.set_color(tuple(fill_color))
    bar.set_edgecolor(tuple(darkened_color))

# Set labels
ax.set_ylabel('')
ax.set_xlabel('log2(B/T cell ratio)')
ax.set_title(f'B/T cell ratio', fontweight='bold', fontsize=10)

# Draw vertical line at x=0
ax.axvline(0, color='lightgray', linestyle='-', linewidth=4, zorder=0)

# Adjust axis position
ax_pos = ax.get_position()
ax.set_position([ax_pos.x0, ax_pos.y0, 0.7, 0.6])

# Add horizontal count bars to the right
ax_pos = ax.get_position()
count_ax_width = ax_pos.width * 0.4
count_ax_h = fig.add_axes([ax_pos.x1 + 0.05, ax_pos.y0, count_ax_width, ax_pos.height])
count_ax_h.set_title(f'Prevalence', fontweight='bold', fontsize=10)

# Plot the count bars horizontally
sns.barplot(y=b_t_ratio_annotated_df['subtype'], x=b_t_ratio_annotated_df['n_annotated_spots'], hue=b_t_ratio_annotated_df['subtype'], 
            ax=count_ax_h, width=0.8, edgecolor='black', linewidth=1, order=sorted_subtypes, 
            palette=color_dict, saturation=1, orient='h', err_kws={'color': 'black', 'linewidth': 1}, capsize=0.3)

# Adjust colors
for i, bar in enumerate(count_ax_h.patches):
    fill_color = bar.get_facecolor()
    fill_color = list(fill_color)
    darkened_color = fill_color.copy()
    for j in range(len(darkened_color)-1):
        darkened_color[j] *= 0.8
    fill_color[-1] = alpha
    bar.set_color(tuple(fill_color))
    bar.set_edgecolor(tuple(darkened_color))

count_ax_h.set_xlabel('Prop. spots', fontsize=10)
count_ax_h.set_ylabel('')
# Remove yticklabels
count_ax_h.set_yticklabels([])

# Increase the xlim a bit
xlim = count_ax_h.get_xlim()
count_ax_h.set_xlim(0, xlim[1] * 1.1)

# Set a shared title centered above the 2 axs
title_x = (ax.get_position().x0 + count_ax_h.get_position().x1) / 2
title_y = ax.get_position().y1 + 0.15
fig.text(title_x, title_y, 'TLS/aggregates', 
         ha='center', fontsize=12, fontweight='bold')

# Add overall title
title_x = (ax.get_position().x0 + count_ax_h.get_position().x1) / 2
title_y = ax.get_position().y1 + 0.4
plt.suptitle(f'All samples - B/T cell proportion in {annotations_of_interest}-annotated spots (annotation_level={level})', 
            y=title_y, x=title_x, fontweight='bold', fontsize=12)


# --- Save the output

out_dir = os.path.join(project_dir, 'results', 'figures',)
os.makedirs(out_dir, exist_ok=True)

filename = 'Fig_1G.pdf'
plt.savefig(os.path.join(out_dir, filename), bbox_inches='tight')


# --- Finish

print("Script successfully completed")