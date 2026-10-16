#!/usr/bin/env python3

# --- Aim
# Plot GC score distribution per annotation as a ridgeplot for lymph node samples
# This creates Fig. 2j

# --- Define environment
# env: spatial-pilot

# --- Load packages
import argparse
import os

import matplotlib.pyplot as plt
from matplotlib import patheffects
import numpy as np
import pandas as pd
import seaborn as sns
from mycolorpy import colorlist as mcp


# --- Argparse arguments
parser = argparse.ArgumentParser(description="Plot GC score distribution per annotation for Fig. 2j",
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
res_dir = args["res_dir"]
sample_list = args["samples"]

print(f'Samples in python: {sample_list}')


# --- Load sample sheet

if args.get("sample_sheet"):
    sample_sheet_path = args["sample_sheet"]
else:
    sample_sheet_path = os.path.join(project_dir, 'config', 'sample_sheet_visium_lymph_nodes.csv')

sample_sheet = pd.read_csv(sample_sheet_path)
sample_lut = sample_sheet.set_index('section_id')


# --- Load GC annotations and scores, create plotting dataframe

# Only plot annotated GCs
annotations_of_interest = ['germinal_centre']
annotation_key = 'annotation'
GC_key = 'GC_score_20'

# Get paths
GC_annotations_dir = os.path.join(project_dir, 'data', 'pathology_annotations', 'visium_lymph_nodes')
GC_res_dir = os.path.join(res_dir, '03_TLS_niche_discovery', '04_map_niche_signatures',
                          '01_score_signatures', 'lymph_nodes', 'GC_score')
most_recent_scores_dir = GC_res_dir

# Load data per sample and merge
df_list = []

for sample in sample_list:
    # Load GC scores (this has all spots)
    scores_file = os.path.join(most_recent_scores_dir, 'csv_files', f'{sample}_GC_score.csv')
    if os.path.exists(scores_file):
        scores_df = pd.read_csv(scores_file, index_col=0)
    else:
        print(f'Warning: {scores_file} not found')
        continue
    
    # Load GC annotations (this only has annotated spots)
    annotations_file = os.path.join(GC_annotations_dir, f'{sample}_GC_annotations.csv')
    if os.path.exists(annotations_file):
        annotations_df = pd.read_csv(annotations_file, index_col=0)
        # Left merge to keep all spots from scores_df
        merged_df = pd.merge(scores_df[[GC_key]], annotations_df[[annotation_key]], 
                             left_index=True, right_index=True, how='left')
    else:
        print(f'Warning: {annotations_file} not found, all spots will be marked as no_annotation')
        merged_df = scores_df[[GC_key]].copy()
        merged_df[annotation_key] = None
    
    merged_df['sample_id'] = sample
    df_list.append(merged_df)

# Concatenate all samples
GC_plotting_df = pd.concat(df_list, ignore_index=True)

# Record spots with no annotation
na_fill = 'no_annotation'
GC_plotting_df[annotation_key] = GC_plotting_df[annotation_key].fillna(na_fill)

# Also overwrite annotations not in the annotations of interest
GC_plotting_df[annotation_key] = np.where(GC_plotting_df[annotation_key].isin(annotations_of_interest), 
                                          GC_plotting_df[annotation_key], na_fill)


# --- Prepare for plotting

plotting_order = annotations_of_interest + [na_fill]

# Define some colors
cmap = "BuPu"
anno_color_list = mcp.gen_color(cmap=cmap, n=(len(annotations_of_interest) + 3))[2:-1]

# Reverse the list
anno_color_list = anno_color_list[::-1]

# Make a new color dict
color_lut_anno = {}
for n, annotation in enumerate(annotations_of_interest):
    color_lut_anno[annotation] = anno_color_list[n]

# Add other
color_lut_anno[na_fill] = 'lightgray'


# --- Plotting functions

# Add labels
def label(label, ax, fontsize=8, y=0.08):
    effects = [patheffects.Stroke(linewidth=2, foreground='white'), patheffects.Normal()]
    ax.text(0, y, label, color='black', fontsize=fontsize,
            ha="left", va="center", transform=ax.transAxes, path_effects=effects)


# Function to add some lines
def add_lines(ax, indicate_threshold=True, add_text=True, ymax=0.5, fontsize=8, 
              text_height=0.4, threshold=0):
    # Make some vertical lines
    # Indicate 0
    ax.axvline(0, ymin=0, ymax=ymax, linestyle='--', linewidth=0.8, c='gray', zorder=-1)
    if indicate_threshold:
        # Indicate the GC threshold
        ax.axvline(threshold, ymin=0, ymax=ymax, linestyle='-', linewidth=0.8, c='red', 
                   zorder=-1, alpha=0.6)
        if add_text:
            # Create white outline effect
            effects = [patheffects.Stroke(linewidth=2, foreground='white'), patheffects.Normal()]
            ax.text(x=threshold+0.02*ax.get_xlim()[1], y=text_height*ax.get_ylim()[1], 
                    s=f'{threshold}', fontsize=fontsize, path_effects=effects)


def make_ridgeplot(df, x, plot_axs, group_key, group_order, xlabel=None, 
                   title='GC score per annotation', y_label=0.05, adjust=0.5, cut=3,
                   add_numbers=False, fontsize=8, title_y=1):
    for n, group in enumerate(group_order):
        # Subset df
        df_sub = df[df[group_key] == group]
        # Get ax
        ax = plot_axs[n]
        # Make distribution plot for that subtype
        color = color_lut_anno[group]
        sns.kdeplot(df_sub, x=x, hue=group_key, fill=True, ax=ax, palette=[color], 
                    legend=False, alpha=0.6, cut=cut)
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


# --- Plot: Ridgeplot split by annotation type

# Plotting
sns.set_theme(style="white", rc={"axes.facecolor": (0, 0, 0, 0), 'axes.linewidth':1})

# Make axs
nrows = len(plotting_order)
ncols = 1

fig, axs = plt.subplots(ncols=ncols, nrows=nrows, figsize=(3, 0.8*nrows), 
                        sharex='all', sharey='all')

make_ridgeplot(GC_plotting_df, x=GC_key, plot_axs=axs, group_key=annotation_key, 
               group_order=plotting_order, y_label=0.1, 
               title="GC score per annotation", add_numbers=True,
               title_y=1, adjust=0.2)

# Add some lines
for n, ax in enumerate(axs):
    add_lines(ax, indicate_threshold=False, ymax=0.8, text_height=0.7)

plt.suptitle('Lymph node samples - GC score distribution per annotation', 
             y=1.15, fontweight='bold', fontsize=10)


# --- Save the output

out_dir = os.path.join(res_dir, 'figures')
os.makedirs(out_dir, exist_ok=True)

filename = 'Fig_2J.pdf'
plt.savefig(os.path.join(out_dir, filename), bbox_inches='tight')

print(f'Figure saved to {os.path.join(out_dir, filename)}')


# --- Finish

print("Script successfully completed")
