#!/usr/bin/env python3

# --- Aim
# Analyze PVL composition across distance bins from vessels.
# Assess tendencies for PVL association with different vessel types.

# --- Define environment
# env: vis-hd

# --- Load packages
import argparse
import os
from datetime import date, datetime

import anndata as ad
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import squidpy as sq
import seaborn as sns
from matplotlib.colors import ListedColormap
import sys
import re

sc.set_figure_params(facecolor="white", figsize=(8, 8))
sc.settings.verbosity = 3

mpl.rcParams['pdf.fonttype'] = 42


# --- Argparse arguments

parser = argparse.ArgumentParser(description="Import Visium data & do basic QC",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)


# Required
parser.add_argument("--res_dir", help="Path to results directory")
parser.add_argument("--h5ad_dir", help="Path to directory for .h5ad output files")
parser.add_argument("--project_dir", help="Project directory (root)")
parser.add_argument("--sample_id", help="Sample ID")
parser.add_argument("--sample_sheet", help="Path to sample sheet CSV", default=None)


args = vars(parser.parse_args())


# --- Arguments

res_dir = args["res_dir"]
h5ad_dir=args["h5ad_dir"]
sample_id=args["sample_id"]
project_dir=args["project_dir"]

sys.path.append(os.path.join(project_dir, 'scripts'))
from helper_functions import auto_crop_scanpy, auto_crop, save_multi_image

# Import function script as a module
sys.path.append(os.path.join(project_dir, 'scripts', '04_vascular_microenvironment'))
import plotting_functions


# --- Get some info from the sample sheet
# Use --sample_sheet if provided, otherwise default to config/sample_sheet_visium_hd.csv
if args["sample_sheet"] is not None:
    sample_sheet_path = args["sample_sheet"]
else:
    sample_sheet_path = os.path.join(project_dir, 'config', 'sample_sheet_visium_hd.csv')

sample_sheet = pd.read_csv(sample_sheet_path)
tissue_id = sample_sheet[sample_sheet['section_id'] == sample_id]['tissue_id'].values[0]


# --- Also get some tissue metadata

resource_dir = os.path.join(project_dir, 'resources')
metadata_file = 'tissue_metadata.csv'

tissue_metadata = pd.read_csv(os.path.join(resource_dir, metadata_file))

# Create a sample - subtype dictionary
subtype = tissue_metadata[tissue_metadata['tissue_id'] == str(tissue_id)]['clinical_subtype'].values[0]


# --- Load the adata object

use_data = 'log_norm'
# Need both the 48 & 16 res
base_res = 'square_016um'
decon_resolution = 'square_016um'
resolutions = [base_res]

if decon_resolution not in resolutions:
    resolutions.append(decon_resolution)

adata_dict = {}

for resolution in resolutions:
    adata_dict[resolution] = ad.read_h5ad(os.path.join(h5ad_dir, 'visium_hd', resolution, f'{sample_id}.h5ad'))


# --- Calculate crop coords for optimal presentation

coords_dict  = {}

for resolution in resolutions:
    coords = auto_crop_scanpy(adata_dict[resolution], border_size=0.1)
    coords_dict[resolution] = coords

# Also calculate some coords for squidpy
coords_dict_sq  = {}

for resolution in resolutions:
    coords = auto_crop(adata_dict[resolution], border_size=0.1)
    coords_dict_sq[resolution] = coords


# --- Load miniatlas classifications

decon_level = 'minor'

# To make it easier to visualize trends by lineage
class_dir = os.path.join(project_dir, 'resources', 'cell_types', 'miniatlas')
class_file = 'miniatlas_subset_classifications.csv'

class_df = pd.read_csv(os.path.join(class_dir, class_file))

# Get endothelial classes
endothelial = list(class_df[class_df['miniatlas_major'] == 'Endothelial'][f'miniatlas_{decon_level}'].unique())
endothelial = sorted(endothelial)

# Decide which decon level to use
decon_level = 'subset'

# Also get PVL classes
PVL_classes = list(class_df[class_df['miniatlas_major'] == 'PVL'][f'miniatlas_{decon_level}'].unique())
PVL_classes = sorted(PVL_classes)

# Also load the colors for the endothelial types
colors_file = os.path.join(project_dir, 'resources', 'plotting_colors', 'endothelial_colors.csv')
colors_df = pd.read_csv(colors_file)

# Make a color_dict
color_dict = colors_df.set_index('endothelial_type')['color'].to_dict()


# --- Load vessel groups

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

# Find most recent results
threshold = 0.5
threshold_adj = str(threshold).split('.')[-1]

# Specify path
TLS_group_dir = os.path.join(project_dir, 'results', '04_vascular_microenvironment', '02_identify_all_vasculature', '01_classify_vessels', 'csv_files')

most_recent_dir_path = find_most_recent(TLS_group_dir)

filename = f'{sample_id}_vessel_classification.csv'
vessel_classes = pd.read_csv(os.path.join(most_recent_dir_path, filename), index_col=0)
adata_dict[base_res].obs[vessel_classes.columns] = vessel_classes

# Collapse to single vessel class per spot
collapsed_col = vessel_classes[endothelial].apply(lambda row: row.first_valid_index() if row.first_valid_index() is not None else None, axis=1)
adata_dict[base_res].obs['vessel_class'] = collapsed_col


# --- Load the deconvolution results

resolution_numbers = re.findall(r'\d+', decon_resolution)[0]
resolution_numbers = int(resolution_numbers)

decon_dir = os.path.join(project_dir, 'results', '02_cell_type_mapping', 'deconvolution', 'visium_hd', f'{decon_level}_level', decon_resolution,)

decon_file = f'{sample_id}_{decon_level}_{resolution}_q05_cell_abundance_w_sf.csv'

decon_res = pd.read_csv(os.path.join(decon_dir, decon_file), index_col=0)

# Only load the PVL classes
adata_dict[decon_resolution].obs[PVL_classes] = decon_res[PVL_classes]


# --- And load the distance calculations

dist_dir = os.path.join(project_dir, 'results', '04_vascular_microenvironment', '02_identify_all_vasculature', '03_calculate_distances', decon_resolution, 'csv_files')
most_recent = find_most_recent(dist_dir)

distance_dict = {}

for cell_type in endothelial:
    dist_file = f'{sample_id}_distance_closest_{cell_type.replace(" ", "_")}_per_spot.csv'
    distance_df = pd.read_csv(os.path.join(most_recent, 'min_dist', dist_file), index_col=0)
    distance_dict[cell_type] = distance_df
    # Also add to adata
    adata_dict[base_res].obs[f'dist_{cell_type}'] = distance_df['nearest_016um_dist_um']


# --- Bin by distance

def bin_distances(df, column_name, bin_size, max_distance):
    # Define the bin edges
    bins = np.arange(0, max_distance + bin_size, bin_size)
    # Use pd.cut to assign each value in the column to a bin
    df['binned'] = pd.cut(df[column_name], bins=bins, right=False, labels=[f'{bins[i]}-{bins[i+1]}' for i in range(len(bins)-1)], include_lowest=True)
    return df['binned']

# Bin distances
bin_size = 25
column_name = 'nearest_016um_dist_um'
max_distance=400

for cell_type in endothelial:
    df_min_dist = distance_dict[cell_type]
    df_min_dist['binned'] = bin_distances(df_min_dist, column_name=column_name, bin_size=bin_size, max_distance=max_distance)
    # Add to adata
    adata_dict[decon_resolution].obs[f'dist_bin_{cell_type}'] = df_min_dist['binned']


# --- Plot 1: Plot distances after binning

sc.set_figure_params(fontsize=15)

ncols = len(endothelial)
nrows = 1

plt.figure(figsize=(8*ncols, 6*nrows))

for n, cell_type in enumerate(endothelial):
    # add a new subplot iteratively
    ax = plt.subplot(nrows, ncols, n + 1)
    # Add legend only to the last one
    if n == (ncols-1):
        legend_loc = 'right margin'
    else:
        legend_loc = None
    adata_sub = adata_dict[decon_resolution][adata_dict[decon_resolution].obs[f'dist_bin_{cell_type}'].notna()]
    sq.pl.spatial_scatter(adata_sub, color=f'dist_bin_{cell_type}', img_alpha=0.7,
                crop_coord=coords_dict_sq[decon_resolution], title=f'{sample_id} ({subtype})\nHEV classification', ax=ax,
                shape='square', vmax=400, legend_loc=legend_loc, palette='viridis')
    adata_sub = adata_dict[base_res][adata_dict[base_res].obs['vessel_class'] == cell_type]
    sq.pl.spatial_scatter(adata_sub, color='vessel_class', img_alpha=0, crop_coord=coords_dict_sq[base_res], 
                title=f'{sample_id} ({subtype})\n{cell_type}',
                palette=ListedColormap(['white']), shape='square', ax=ax, legend_loc=None)
    ax.set_title(ax.get_title(), fontweight='bold')

plt.suptitle(f'{sample_id} ({subtype}) - Distances after binning (resolution={decon_resolution})', y=1.1, fontweight='bold', x=0.5)


# --- Plot 2: Plot PVL composition across distance bins

sc.set_figure_params(fontsize=10)

# Prepend dist_bin_ to endothelial keys
endothelial_keys = [f'dist_bin_{cell_type}' for cell_type in endothelial]

# Enable gridlines behind the plot for all plots
plt.rcParams['axes.axisbelow'] = True

ncols = len(PVL_classes)

plt.figure(figsize=(2, 2.2),)

for n, PVL_class in enumerate(PVL_classes):
    # add a new subplot iteratively
    ax = plt.subplot(nrows, ncols, n + 1)
    # Extract data to plot
    cell_df = adata_dict[decon_resolution].obs[endothelial_keys + [PVL_class]]
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
    sns.lineplot(
        data=df_long,
        x='distance_bin',
        y='PVL_value',
        hue='endothelial_subtype',
        marker='o',
        errorbar='ci',
        ax=ax,
        palette=color_dict,
        legend=False,
        markeredgecolor=None,
        markersize=4,
    )
    ax.set_title(f'{PVL_class}', fontweight='bold')
    # Rotate xtick labels 90
    ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
    # Adjust labels
    ax.set_xlabel('Distance (µm)')
    ax.set_ylabel('Proportion')
    # Adjust x-limits to fit the bins
    ax.set_xlim(0, ax.get_xlim()[1]-0.75)

axs = plt.gcf().get_axes()

# Adjust positioning & size
for n, ax in enumerate(axs):
    ax_pos = ax.get_position()
    # Shift the subplot to the desired location + set the width so it's equal for all plots
    if n == 0:
        ax.set_position([ax_pos.x0, ax_pos.y0, len(bin_order)*0.08, ax_pos.height])
    else:
        new_pos = ax_pos_prev.x1 + 0.4
        ax.set_position([new_pos, ax_pos.y0, len(bin_order)*0.08, ax_pos.height])
    ax_pos_prev = ax.get_position()

# Move legend
plotting_functions.add_legend(ax=ax, ordered_list=endothelial, color_lut=color_dict, title='Vessel type', 
                               loc='center left', y_offset=0.5, x_offset=1)

# Calculate x to center
title_x = (axs[0].get_position().x0 + axs[-1].get_position().x1)/2

plt.suptitle(f'{sample_id} ({subtype}) - PVL composition across distance bins (threshold={threshold}, decon_res={decon_resolution}, decon_level={decon_level})',
             y=1.2, fontweight='bold', x=title_x)


# --- Plot 3: Same as above but shared ymax

plt.figure(figsize=(2, 2.2),)

ymax = 0

for n, PVL_class in enumerate(PVL_classes):
    # add a new subplot iteratively
    ax = plt.subplot(nrows, ncols, n + 1)
    # Extract data to plot
    cell_df = adata_dict[decon_resolution].obs[endothelial_keys + [PVL_class]]
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
    sns.lineplot(
        data=df_long,
        x='distance_bin',
        y='PVL_value',
        hue='endothelial_subtype',
        marker='o',
        errorbar='ci',
        ax=ax,
        palette=color_dict,
        legend=False,
        markeredgecolor=None,
        markersize=4,
    )
    ax.set_title(f'{PVL_class}', fontweight='bold')
    # Rotate xtick labels 90
    ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
    # Adjust labels
    ax.set_xlabel('Distance (µm)')
    ax.set_ylabel('Proportion')
    if ax.get_ylim()[1] > ymax:
        ymax = ax.get_ylim()[1]
    # Adjust x-limits to fit the bins
    ax.set_xlim(0, ax.get_xlim()[1]-0.75)

axs = plt.gcf().get_axes()

for ax in axs:
    # Set the same y-limits for all subplots
    ax.set_ylim(0, ymax)

# Adjust positioning & size
for n, ax in enumerate(axs):
    ax_pos = ax.get_position()
    # Shift the subplot to the desired location + set the width so it's equal for all plots
    if n == 0:
        ax.set_position([ax_pos.x0, ax_pos.y0, len(bin_order)*0.08, ax_pos.height])
    else:
        new_pos = ax_pos_prev.x1 + 0.4
        ax.set_position([new_pos, ax_pos.y0, len(bin_order)*0.08, ax_pos.height])
    ax_pos_prev = ax.get_position()

# Move legend
plotting_functions.add_legend(ax=ax, ordered_list=endothelial, color_lut=color_dict, title='Vessel type', 
                               loc='center left', y_offset=0.5, x_offset=1)

# Calculate x to center
title_x = (axs[0].get_position().x0 + axs[-1].get_position().x1)/2

plt.suptitle(f'{sample_id} ({subtype}) - PVL composition across distance bins - same ymax (threshold={threshold}, decon_res={decon_resolution}, decon_level={decon_level})',
             y=1.2, fontweight='bold', x=title_x)


# --- Plot 4: Same as above but coloring area under plot

plt.figure(figsize=(2, 2.2),)

for n, PVL_class in enumerate(PVL_classes):
    # add a new subplot iteratively
    ax = plt.subplot(nrows, ncols, n + 1)
    # Extract data to plot
    cell_df = adata_dict[decon_resolution].obs[endothelial_keys + [PVL_class]]
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
    # Plot each subtype separately
    ymax = df_long.groupby(['distance_bin', 'endothelial_subtype'])['PVL_value'].mean().max() * 1.1  # Initialize ymax for the current PVL_class
    for endo_subtype in endothelial:
        sub_df = df_long[df_long['endothelial_subtype'] == endo_subtype]
        # Calculate mean per bin
        sub_df = sub_df.groupby('distance_bin')['PVL_value'].mean().reset_index()
        sub_df.set_index('distance_bin').plot.area(ax=ax, color=color_dict[endo_subtype], legend=False, alpha=0.2)
    sns.lineplot(
        data=df_long,
        x='distance_bin',
        y='PVL_value',
        hue='endothelial_subtype',
        marker='o',
        errorbar=None,
        ax=ax,
        palette=color_dict,
        legend=False,
        markeredgecolor=None,
        markersize=0,
    )
    # Make sure all xticklabels are displayed
    ax.set_xticklabels(ax.get_xticklabels(), rotation=90)  # Rotate x-tick labels if needed
    ax.set_xticks(range(len(bin_order)))  # Ensure all x-ticks are placed
    ax.set_xticklabels(bin_order)  # Set the x-tick labels to match the DataFrame index
    # Adjust lims
    # ymax = max(composition_collapsed[cell_type]) * 1.1
    ax.set_ylim(0, ymax)
    ax.set_xlim(0, ax.get_xlim()[1]-0.75)
    ax.set_title(f'{PVL_class}', fontweight='bold')
    # Rotate xtick labels 90
    ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
    # Adjust labels
    ax.set_xlabel('Distance (µm)')
    ax.set_ylabel('Proportion')
    # Adjust x-limits to fit the bins
    ax.set_xlim(0, ax.get_xlim()[1]-0.75)

axs = plt.gcf().get_axes()

# Adjust positioning & size
for n, ax in enumerate(axs):
    ax_pos = ax.get_position()
    # Shift the subplot to the desired location + set the width so it's equal for all plots
    if n == 0:
        ax.set_position([ax_pos.x0, ax_pos.y0, len(bin_order)*0.08, ax_pos.height])
    else:
        new_pos = ax_pos_prev.x1 + 0.4
        ax.set_position([new_pos, ax_pos.y0, len(bin_order)*0.08, ax_pos.height])
    ax_pos_prev = ax.get_position()

# Move legend
plotting_functions.add_legend(ax=ax, ordered_list=endothelial, color_lut=color_dict, title='Vessel type', 
                               loc='center left', y_offset=0.5, x_offset=1)

# Calculate x to center
title_x = (axs[0].get_position().x0 + axs[-1].get_position().x1)/2

plt.suptitle(f'{sample_id} ({subtype}) - PVL composition across distance bins - same ymax (threshold={threshold}, decon_res={decon_resolution}, decon_level={decon_level})',
             y=1.2, fontweight='bold', x=title_x)


# --- Plot 5: Same as above but shared ymax

plt.figure(figsize=(2, 2.2),)

for n, PVL_class in enumerate(PVL_classes):
    # add a new subplot iteratively
    ax = plt.subplot(nrows, ncols, n + 1)
    # Extract data to plot
    cell_df = adata_dict[decon_resolution].obs[endothelial_keys + [PVL_class]]
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
    # Plot each subtype separately
    ymax = df_long.groupby(['distance_bin', 'endothelial_subtype'])['PVL_value'].mean().max() * 1.1  # Initialize ymax for the current PVL_class
    for endo_subtype in endothelial:
        sub_df = df_long[df_long['endothelial_subtype'] == endo_subtype]
        # Calculate mean per bin
        sub_df = sub_df.groupby('distance_bin')['PVL_value'].mean().reset_index()
        sub_df.set_index('distance_bin').plot.area(ax=ax, color=color_dict[endo_subtype], legend=False, alpha=0.2)
    sns.lineplot(
        data=df_long,
        x='distance_bin',
        y='PVL_value',
        hue='endothelial_subtype',
        marker='o',
        errorbar=None,
        ax=ax,
        palette=color_dict,
        legend=False,
        markeredgecolor=None,
        markersize=0,
    )
    # Make sure all xticklabels are displayed
    ax.set_xticklabels(ax.get_xticklabels(), rotation=90)  # Rotate x-tick labels if needed
    ax.set_xticks(range(len(bin_order)))  # Ensure all x-ticks are placed
    ax.set_xticklabels(bin_order)  # Set the x-tick labels to match the DataFrame index
    # Adjust lims
    # ymax = max(composition_collapsed[cell_type]) * 1.1
    ax.set_ylim(0, ymax)
    ax.set_xlim(0, ax.get_xlim()[1]-0.75)
    ax.set_title(f'{PVL_class}', fontweight='bold')
    # Rotate xtick labels 90
    ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
    # Adjust labels
    ax.set_xlabel('Distance (µm)')
    ax.set_ylabel('Proportion')
    # Adjust x-limits to fit the bins
    ax.set_xlim(0, ax.get_xlim()[1]-0.75)

axs = plt.gcf().get_axes()

ymax = 0

for ax in axs:
    if ax.get_ylim()[1] > ymax:
        ymax = ax.get_ylim()[1]

# Adjust positioning & size
for n, ax in enumerate(axs):
    ax_pos = ax.get_position()
    # Shift the subplot to the desired location + set the width so it's equal for all plots
    if n == 0:
        ax.set_position([ax_pos.x0, ax_pos.y0, len(bin_order)*0.08, ax_pos.height])
    else:
        new_pos = ax_pos_prev.x1 + 0.4
        ax.set_position([new_pos, ax_pos.y0, len(bin_order)*0.08, ax_pos.height])
    ax_pos_prev = ax.get_position()
    ax.set_ylim(0, ymax)  # Set the same y-limits for all subplots

# Move legend
plotting_functions.add_legend(ax=ax, ordered_list=endothelial, color_lut=color_dict, title='Vessel type', 
                               loc='center left', y_offset=0.5, x_offset=1)

# Calculate x to center
title_x = (axs[0].get_position().x0 + axs[-1].get_position().x1)/2

plt.suptitle(f'{sample_id} ({subtype}) - PVL composition across distance bins - same ymax (threshold={threshold}, decon_res={decon_resolution}, decon_level={decon_level})',
             y=1.2, fontweight='bold', x=title_x)


# --- Save output

qc_dir = os.path.join(project_dir, 'results', '04_vascular_microenvironment', '02_identify_all_vasculature', '04_assess_PVL_association')

today = date.today()
today = today.strftime("%Y%m%d")

out_dir = os.path.join(qc_dir)

os.makedirs(out_dir, exist_ok=True)

# Save plots
out_dir_plots = os.path.join(out_dir, 'plots', today)
os.makedirs(out_dir_plots, exist_ok=True)

# Save plots
filename = f'{sample_id}_all_vessels_PVL_composition_vs_distance.pdf'

save_multi_image(os.path.join(out_dir_plots, filename))
plt.close('all')

# Also save csv file for convenience
out_dir_csv = os.path.join(out_dir, 'csv_files', today)
os.makedirs(out_dir_csv, exist_ok=True)

cell_df = adata_dict[decon_resolution].obs[endothelial_keys + PVL_classes]
# Save as gzipped csv
cell_df.to_csv(os.path.join(out_dir_csv, f'{sample_id}_PVL_composition_vs_distance.csv.gz'), index=True, compression='gzip')


# --- Finish

print("Script successfully completed")
