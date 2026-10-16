#!/usr/bin/env python3

# --- Aim
# Plot Figure 2B: ridgeplot showing TLS model scores per pathology annotation
# Shows distribution of TLS scores across different tissue compartments (expect highest in TLS/aggregate regions)

# --- Define environment
# env: spatial-pilot

# --- Load packages
import argparse
import os

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
from mycolorpy import colorlist as mcp

import seaborn as sns
from matplotlib import patheffects

mpl.rcParams['pdf.fonttype'] = 42  # Use Type 42 (TrueType) fonts in PDF files


# --- Argparse arguments
parser = argparse.ArgumentParser(description="Fig. 2b: TLS model scores per pathology annotation",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument("-s", "--samples",
    nargs="*",  # 0 or more values expected => creates a list
    type=str,
    help="List of samples")
parser.add_argument("-r", "--res_dir", help="Path to results directory")
parser.add_argument("--h5ad_dir", help="Path to directory for .h5ad output files")
parser.add_argument("--project_dir", help="Project directory (root)")
parser.add_argument("--TLS_threshold", help="Threshold for selecting TLS-high spots (TLS model)", type=float)
parser.add_argument("--GC_threshold", help="Threshold for selecting GC-high spots (GC score)", type=float)
parser.add_argument("--sample_sheet", help="Path to sample sheet CSV file", default=None)

args = vars(parser.parse_args())

# Arguments
project_dir = args["project_dir"]
h5ad_dir = args["h5ad_dir"]
res_dir = args["res_dir"]
sample_list = args["samples"]
TLS_threshold=args["TLS_threshold"]
GC_threshold=args["GC_threshold"]

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
sample_lut = sample_sheet.set_index('section_id')


# --- Also get some tissue metadata

resource_dir = os.path.join(project_dir, 'resources')
metadata_file = 'tissue_metadata.csv'

tissue_metadata = pd.read_csv(os.path.join(resource_dir, metadata_file))

# Merge the two DataFrames on the common column 'tissue_id'
merged_df = pd.merge(sample_sheet, tissue_metadata, on='tissue_id')

# Create a sample - subtype dictionary
subtype_dict = dict(zip(merged_df['section_id'], merged_df['clinical_subtype']))


# --- Read in the immune path annotations

path_layer = 'immune'

# Define annotation path
reworked_anno_dir = os.path.join(project_dir, 'data', 'pathology_annotations', 'visium')

immune_annotations = {}

for sample in sample_list:
    annotation_file = f'{sample}_pathology_annotations.csv'
    annotation_path = os.path.join(reworked_anno_dir, annotation_file)
    if os.path.exists(annotation_path):
        path_anno = pd.read_csv(annotation_path, index_col=0)
        immune_annotations[sample] = path_anno[[path_layer]]
    else:
        print(f'No annotations found for {sample} - skipping')


# --- Load the TLS model results for each sample

TLS_res_dir = os.path.join(res_dir, '03_TLS_niche_discovery', '01_identify_TLS', '01_run_TLS_model')

most_recent_dir_path = TLS_res_dir

# Read in the results & add to anndata
for sample in sample_list:
    TLS_file_name = f'{sample}_TLS_model.csv'
    TLS_results = pd.read_csv(os.path.join(most_recent_dir_path, TLS_file_name), index_col=0)
    immune_annotations[sample]['TLS_model'] = TLS_results


# --- Prepare for plotting the result per subtype: define a plotting order based on subtype

# Extract all unique subtypes from the metadata
all_subtypes = set(subtype_dict[sample] for sample in sample_list if sample in subtype_dict)

# Sort according to preferred plotting order
subtype_order = ['Luminal A', 'Luminal B', 'Luminal B (HER2+)', 'HER2+', 'TNBC']
sorted_subtypes = sorted(all_subtypes, key=lambda x: subtype_order.index(x) if x in subtype_order else float('inf'))


# --- Extract a single df for all samples

for sample in sample_list:
    if sample in immune_annotations:
        # Add the sample_id to the annotations
        immune_annotations[sample]['sample_id'] = sample
    else:
        print(f'No immune annotations found for {sample}')

# Concatenate all the annotations into a single DataFrame
TLS_path_anno = pd.concat(immune_annotations.values(), ignore_index=True)


# --- Check that immune categories match the annotations of interest

# Define annotations of interest
annotations_of_interest = ['tls', 'aggregate', 'lymphocytes_present']
all_immune_cat = list(TLS_path_anno['immune'].dropna().unique())

# Check all items are present in annotations_of_interest
all_items_in_list = all(item in annotations_of_interest for item in all_immune_cat)
items_not_in_list2 = [item for item in all_immune_cat if item not in annotations_of_interest]

if all_items_in_list:
    print(f'All good: annotations of interest ({annotations_of_interest}) match immune categories ({all_immune_cat})')
else:
    print(f'WARNING: Not all categories are in annotations of interest ({annotations_of_interest})')
    print(f'Missing categories ({items_not_in_list2}) will be replaced with immune_other')
    # Replace anything that isn't in annotations of interest with immune_other
    for sample in sample_list:
        immune_anno = immune_annotations[sample]
        if 'immune' in immune_anno.columns:
            # Replace items not in annotations_of_interest with 'immune_other'
            immune_anno['immune'] = np.where(immune_anno['immune'].isin(items_not_in_list2), 'immune_other', immune_anno['immune'])
            # Also update the all_adata object
            if sample in TLS_path_anno.index:
                TLS_path_anno.loc[TLS_path_anno['sample_id'] == sample, 'immune'] = np.where(TLS_path_anno['immune'].isin(items_not_in_list2), 'immune_other', TLS_path_anno['immune'])

# Now sort the immune categories that are present
all_immune_cat = list(TLS_path_anno['immune'].dropna().unique())
plotting_order = sorted(all_immune_cat, key=lambda x: annotations_of_interest.index(x) if x in annotations_of_interest else float('inf'))


# --- Extract dataframe

# Mark all the NaNs as non_immune
TLS_path_anno['immune'] = TLS_path_anno['immune'].fillna('non_immune')

plotting_order.append('non_immune')


# --- Get some colors for plotting the annotations

# Prepare color dictionary for plotting
cmap = "YlOrRd"
anno_color_list = mcp.gen_color(cmap=cmap, n=(len(annotations_of_interest) + 3))[2:-1] # Don't want to get the lightest color (hard to see)

# Reverse the list
anno_color_list = anno_color_list[::-1]

# Make a new color dict
color_lut_anno = {}
for n, annotation in enumerate(annotations_of_interest):
    color_lut_anno[annotation] = anno_color_list[n]

# Add other
color_lut_anno['non_immune'] = 'lightgray'


# --- Some plotting functions to avoid repetition

# Add labels
def label(label, ax, fontsize=8, y=0.08):
    effects = [patheffects.Stroke(linewidth=2, foreground='white'), patheffects.Normal()]
    ax.text(0, y, label, color='black', fontsize=fontsize,
            ha="left", va="center", transform=ax.transAxes, path_effects=effects)

# Function to add some lines
def add_lines(ax, indicate_threshold=True, add_text=True, ymax=0.5, fontsize=8, text_height=0.4):
    # Make some vertical lines
    # Indicate 0
    ax.axvline(0, ymin=0, ymax=ymax, linestyle='--', linewidth=0.8, c='gray', zorder=-1)
    if indicate_threshold:
        # Indicate the TLS threshold
        ax.axvline(TLS_threshold, ymin=0, ymax=ymax, linestyle='-', linewidth=0.8, c='red', zorder=-1, alpha=0.6)
        if add_text:
            # Create white outline effect
            effects = [patheffects.Stroke(linewidth=2, foreground='white'), patheffects.Normal()]
            ax.text(x=TLS_threshold+0.02*ax.get_xlim()[1], y=text_height*ax.get_ylim()[1], s=f'{TLS_threshold}', fontsize=fontsize,
                    path_effects=effects)

def make_ridgeplot(df, x, plot_axs, group_key, group_order, xlabel=None, title='GC score per annotation', y_label=0.05, adjust=0.5, cut=3,
                   add_numbers=False, fontsize=8, title_y=1):
    for n, group in enumerate(group_order):
        # Subset df
        df_sub = df[df[group_key] == group]
        # Get ax
        ax = plot_axs[n]
        # Make distribution plot for that subtype
        color = color_lut_anno[group]
        sns.kdeplot(df_sub, x=x, hue=group_key, fill=True, ax=ax, palette=[color], legend=False,
                    alpha=0.6, cut=cut)
        # Modify ax
        # Hide the y-axis and its labels
        ax.yaxis.set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        # Get the current x-tick labels
        ax.tick_params(axis='x', labelsize=8)
        if xlabel is None:
            xlabel = x
        ax.set_xlabel(xlabel, fontsize=9)
        if add_numbers:
            total = len(df_sub[x].dropna())
            reformat = reformat_number(total)
            text = f'n = {reformat}'
            effects = [patheffects.Stroke(linewidth=2, foreground='white'), patheffects.Normal()]
            ax.text(1, y_label, text, color='black', fontsize=fontsize,
                    ha="right", va="center", transform=ax.transAxes, path_effects=effects)
    # Make them overlap a bit
    plt.subplots_adjust(hspace=-adjust)
    # Plot the subtype label
    for n, subtype in enumerate(group_order):
        label(subtype, plot_axs[n], y=y_label, fontsize=fontsize)
    # Add title
    plot_axs[0].set_title(title, fontweight='bold', fontsize=9, y=title_y)


# Function to reformat numbers
def reformat_number(number):
    if number >= 1000000:
        # Divide the number by 1,000,000 and format it to one decimal place
        formatted_number = "{:.1f}m".format(number / 1000000)
    elif number >= 1000:
        # Divide the number by 1000 and format it to one decimal place
        formatted_number = "{:.1f}k".format(number / 1000)
    else:
        formatted_number = str(round(number, 1))
    return formatted_number

# Calculate percentage above threshold
# Group by 'type' and calculate the percentage of values above the threshold
def calculate_percentage_above_threshold(group, threshold, key='TLS_model'):
    count_above_threshold = (group[key] > threshold).sum()
    total_count = group[key].count()
    percentage = (count_above_threshold / total_count) * 100
    return percentage

def add_shared_ylabel(ylabel, axs, fontsize=9):
    # Ensure the figure is drawn so positions are accurate
    plt.gcf().canvas.draw()
    # Get positions in figure coordinates
    top_pos = axs[0].get_position()
    bottom_pos = axs[-1].get_position()
    # Calculate the vertical midpoint
    y_mid = (top_pos.y1 + bottom_pos.y0) / 2
    # Choose an x-position just left of the axes (e.g., 0.02)
    x_label = top_pos.x0 - 0.02
    # Add the shared ylabel
    plt.gcf().text(x_label, y_mid, ylabel, va='center', ha='right', rotation='vertical',
                   fontsize=fontsize)


# --- Define a dict to update the labels

clean_labels = {
    'aggregate': 'Aggregate',
    'tls': 'TLS',
    'lymphocytes_present': 'Lymphocytes present',
    'non_immune': 'Non-immune',
}


# --- Plot 1a: Ridgeplot split by immune annotation type

# Plotting
sns.set_theme(style="white", rc={"axes.facecolor": (0, 0, 0, 0), 'axes.linewidth':1})

# Make axs
nrows = len(plotting_order)
ncols = 1

fig, axs = plt.subplots(ncols=ncols, nrows=nrows, figsize=(3, 0.8*nrows), sharex='all', sharey='all')

make_ridgeplot(TLS_path_anno, x='TLS_model', plot_axs=axs, group_key='immune', group_order=plotting_order,
            y_label=0.1, title=f"TLS score per annotation", add_numbers=True,
            title_y=1, adjust=0.2)

# Add some lines
for n, ax in enumerate(axs):
    add_lines(ax, indicate_threshold=False, ymax=0.8, text_height=0.7)

# Update the labels
for ax in axs:
    # Find the text objects in the axes
    for text in ax.texts:
        # Check if the text is in the clean_labels dictionary
        if text.get_text() in clean_labels:
            # Update the text to the cleaned label
            text.set_text(clean_labels[text.get_text()])

axs[-1].set_xlabel('TLS model score', fontsize=9)

# Add shared y-label
add_shared_ylabel('Pathology annotation', axs)

plt.suptitle(f'All samples - TLS score distribution per immune annotatation (pathologist annotation) - n = number of spots', y=1.1, fontweight='bold', fontsize=10)


# --- Save the output

output_dir = os.path.join(res_dir, 'figures')
os.makedirs(output_dir, exist_ok=True)

filename = 'Fig_2B.pdf'
filepath = os.path.join(output_dir, filename)
plt.savefig(filepath, bbox_inches='tight')

print("Script succesfully completed")
