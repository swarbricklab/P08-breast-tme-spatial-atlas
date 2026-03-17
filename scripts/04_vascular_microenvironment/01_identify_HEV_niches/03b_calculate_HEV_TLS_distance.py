#!/usr/bin/env python3

# --- Aim
# Calculate the distance from each HEV group to the nearest TLS area.
# Classify HEVs as intra-TLS, peri-TLS, or non-TLS based on spatial proximity.

# --- Define environment
# env: vis-hd

# --- Load packages
import argparse
import os
import sys
from datetime import date, datetime

import anndata as ad
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from matplotlib.colors import ListedColormap, to_rgb
from scipy.spatial.distance import cdist

sc.set_figure_params(facecolor="white", figsize=(8, 8))
sc.settings.verbosity = 3


# --- Argparse arguments

parser = argparse.ArgumentParser(description="Import Visium data & do basic QC",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

# Required

parser.add_argument("--res_dir", help="Path to results directory")
parser.add_argument("--h5ad_dir", help="Path to directory for .h5ad output files")
parser.add_argument("--project_dir", help="Project directory (root)")
parser.add_argument("--sample_id", help="Sample ID")
parser.add_argument("--sample_sheet", help="Path to sample sheet CSV file", default=None)

args = vars(parser.parse_args())


# --- Arguments

res_dir = args["res_dir"]
h5ad_dir=args["h5ad_dir"]
sample_id=args["sample_id"]
project_dir=args["project_dir"]
sys.path.append(os.path.join(project_dir, 'scripts'))
from helper_functions import auto_crop_scanpy, save_multi_image


# --- Get some info from the sample sheet
# Use sample_sheet argument if provided, else default
if args.get("sample_sheet"):
    sample_sheet_path = args["sample_sheet"]
else:
    sample_dir = os.path.join(project_dir, 'config')
    sample_sheet_file = 'sample_sheet_visium_hd.csv'
    sample_sheet_path = os.path.join(sample_dir, sample_sheet_file)

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
resolutions = ['square_048um', 'square_016um']

adata_dict = {}

for resolution in resolutions:
    adata_dict[resolution] = ad.read_h5ad(os.path.join(h5ad_dir, 'visium_hd', resolution, f'{sample_id}.h5ad'))


# --- Calculate crop coords for optimal presentation
coords_dict  = {}

for resolution in resolutions:
    coords = auto_crop_scanpy(adata_dict[resolution], border_size=0.1)
    coords_dict[resolution] = coords


# --- Load TLS zone groups
# Function to find the most recent results
def find_most_recent(dir):
    subdirs = [d for d in os.listdir(dir) if os.path.isdir(os.path.join(dir, d))]
    date_format = '%Y%m%d'
    valid_subdirs = []
    for subdir in subdirs:
        try:
            datetime.strptime(subdir, date_format)
            valid_subdirs.append(subdir)
        except ValueError:
            pass
    if valid_subdirs:
        most_recent_dir = max(valid_subdirs, key=lambda x: datetime.strptime(x, date_format))
        most_recent_dir_path = os.path.join(dir, most_recent_dir)
        print("Most recent results:", most_recent_dir_path)
    else:
        print("No valid directories found.")
        most_recent_dir_path = dir
    return most_recent_dir_path

# Find most recent results
threshold = 0.5
threshold_adj = str(threshold).split('.')[-1]

TLS_group_dir = os.path.join(res_dir, '04_vascular_microenvironment', '01_identify_HEV_niches', '03a_identify_TLS_areas', 'square_048um', 'csv_files')

most_recent_dir_path = find_most_recent(TLS_group_dir)

# Read in the results & add to anndata
TLS_file_name = f'{sample_id}_TLS_zones_groups.csv'
TLS_results = pd.read_csv(os.path.join(most_recent_dir_path, TLS_file_name), index_col=0)

# Decide minimum number of spots wish to consider
min_spots = 2
group_key = 'group'

# Find groups that are too small
too_small = TLS_results[group_key].value_counts()[TLS_results[group_key].value_counts() < min_spots].index

# Make a temporary df
temp_df = TLS_results[[group_key]]

# Set small groups as too_small
too_small_value = 'too_small'
temp_df[group_key] = np.where(temp_df[group_key].isin(too_small), too_small_value, temp_df[group_key])

# Save singletons as separat df
too_small = temp_df[temp_df[group_key] == too_small_value]

# Remove the singletons for now
temp_df = temp_df[temp_df[group_key] != too_small_value]

# Update group names based on size
group_counts = temp_df[group_key].value_counts()
sorted_groups = group_counts.sort_values(ascending=False).index
group_map = {group: f"group{i+1}" for i, group in enumerate(sorted_groups)}

temp_df[group_key] = temp_df[group_key].map(group_map)

# Add results to a new column
adata_dict['square_048um'].obs['group_too_small_removed'] = temp_df[group_key]

group_key = 'group_too_small_removed'


# --- Also load the HEV groups 

threshold_hev = 0.5
threshold_hev_adj = str(threshold_hev).split('.')[-1]

# Use the same directory structure as in 04_vascular_microenvironment/01_identify_HEV_niches/03a_identify_TLS_areas.py
TLS_group_dir = os.path.join(res_dir, '04_vascular_microenvironment', '01_identify_HEV_niches', '02b_expand_HEV_groups', 'square_016um', f'threshold_{threshold_hev_adj}')
TLS_group_dir = os.path.join(TLS_group_dir, 'csv_files')

most_recent_dir_path = find_most_recent(TLS_group_dir)

# Read in the results & add to anndata
score_of_i = 'HEVs'
TLS_file_name = f'{sample_id}_{score_of_i}_groups.csv'
TLS_results = pd.read_csv(os.path.join(most_recent_dir_path, TLS_file_name), index_col=0)

group_key = 'group_too_small_removed'
group_key_source = 'expanded_groups'

adata_dict['square_016um'].obs[group_key] = TLS_results[group_key_source]


# --- And load the scoring results

clusters_of_i = {
    'B_cell_zone': 'b_cell_zone_top_genes_filtered.csv',
    'T_cell_zone': 't_cell_zone_top_genes_filtered.csv',
    'HEVs': 'hevs_filtered_endothelial_genes.csv',
}   

# Define the TLS zone names
TLS_zone_names = ['T_cell_zone', 'B_cell_zone']

score_keys_dict = {}

for score_name, filename in clusters_of_i.items():
    out_name = f'{score_name}_score'
    if score_name in TLS_zone_names:
        resolution = 'square_048um'
    else:
        resolution = 'square_016um'
    GC_res_dir = os.path.join(res_dir, '04_vascular_microenvironment', '01_identify_HEV_niches', '01_score_signature', f'{score_name}_score', resolution)
    most_recent_dir_path = find_most_recent(GC_res_dir)
    # Read results and add to adata
    GC_file_name = f'{sample_id}_{out_name}.csv'
    GC_results = pd.read_csv(os.path.join(most_recent_dir_path, 'csv_files', GC_file_name), index_col=0)
    # Choose number of genes used for scoring
    top_gene = 20
    GC_key = f'{score_name}_score_{top_gene}'
    score_keys_dict[score_name] = GC_key
    adata_dict[resolution].obs[GC_key] = GC_results[GC_key]


# Also make colormaps per score
# Use function from cell
white_spacing=20

# Use function from cell2location to make colormap with increasing transparency gradient
# To be able to plot multiple scores in single plot
def create_colormap(R, G, B):
    spacing = int(white_spacing * 2.55)
    N = 255
    M = 3
    alphas = np.concatenate([[0] * spacing * M,
                                np.linspace(0, 1.0, (N - spacing) * M)])
    vals = np.ones((N * M, 4))
    for i, color in enumerate([R, G, B]):
        vals[:, i] = color / 255
    vals[:, 3] = alphas
    return ListedColormap(vals)

# Define some colors per score
color_dict = {
    'GC': mpl.colors.cnames["cyan"],
    'B_cell_zone': mpl.colors.cnames["cornflowerblue"],
    'T_cell_zone': mpl.colors.cnames["crimson"],
    'HEVs': mpl.colors.cnames["indigo"],    
}

color_dict_rgb = {}
for key, value in color_dict.items():
    rgb = tuple(int(c * 255) for c in to_rgb(value))
    color_dict_rgb[key] = rgb

cmap_dict = {}

for score_name in clusters_of_i:
    rgb = color_dict_rgb[score_name]
    cmap = create_colormap(rgb[0], rgb[1], rgb[2])
    cmap_dict[score_name] = cmap


# --- Extract 48 µm groups (TLS zones)

adata = adata_dict['square_048um']

coords = adata.obsm['spatial']  # Extract the spot coords
groups = adata.obs[group_key]  # Group labels (48 µm)

# Convert to DataFrame for easier manipulation
df = pd.DataFrame(coords, columns=["x", "y"])
df["group"] = groups.values
df = df.dropna(subset='group')

# Extract coordinates for each group
grouped_coords_48 = {
    group: df[df["group"] == group][["x", "y"]].to_numpy()
    for group in df["group"].unique()
}


# --- Extract 16 µm groups (HEVs)

adata = adata_dict['square_016um']

coords = adata.obsm['spatial']  # Extract the spot coords
groups = adata.obs[group_key]  # Group labels (48 µm)

# Convert to DataFrame for easier manipulation
df = pd.DataFrame(coords, columns=["x", "y"])
df["group"] = groups.values
df = df.dropna(subset='group')

# Extract coordinates for each group
grouped_coords_16 = {
    group: df[df["group"] == group][["x", "y"]].to_numpy()
    for group in df["group"].unique()
}


# --- Calculate distance for each HEV to closest TLS group

# Store the results (16um group, 48um group, min distance)
results = []

# Calculate pairwise distances and find the minimum distance for each 16um group
for group_16, coords_16 in grouped_coords_16.items():
    min_distance = np.inf  # Initialize with a very large value
    closest_group_48 = None    
    # Compare each 16um group to each 48um group
    for group_48, coords_48 in grouped_coords_48.items():
        # Calculate pairwise distances between each 16um point and 48um points
        distances = cdist(coords_16, coords_48)
        # Find the smallest distance for this pair of groups
        min_dist_to_group_48 = np.min(distances)
        # If this is the smallest distance found, update the results
        if min_dist_to_group_48 < min_distance:
            min_distance = min_dist_to_group_48
            closest_group_48 = group_48
    # Append the result for this 16um group
    results.append({
        'HEV_group_16um': group_16,
        'TLS_zones_group_48um': closest_group_48,
        'min_distance': min_distance
    })

# Convert results to a DataFrame for easier analysis
results_df = pd.DataFrame(results)


# --- Divide HEVs into intra-TLS, peri-TLS, non-TLS
# Based on distance to closest TLS

# First get conversion values pixel to um
adata = adata_dict['square_048um']

# Get a scalefactor for converting pixel coords to um
spot_diameter_pixel = adata.uns['spatial'][sample_id]['scalefactors']['spot_diameter_fullres']

spot_diameter_um = 48 # Square size is whatever the resolution is 
scalefactor_um = spot_diameter_um/spot_diameter_pixel

# Add um distances
results_df['min_dist_um'] = results_df['min_distance'] * scalefactor_um

# Now classify HEVs based on how many spots from closest TLS
peri_TLS_bound = 4
intra_TLS_bound = 1

results_df['HEV_classification'] = 'non_TLS'

results_df['HEV_classification'] = np.where(results_df['min_distance'] < spot_diameter_pixel*peri_TLS_bound,
         'peri_TLS',
         results_df['HEV_classification'])

results_df['HEV_classification'] = np.where(results_df['min_distance'] < spot_diameter_pixel*intra_TLS_bound,
         'intra_TLS',
         results_df['HEV_classification'])

TLS_HEV_cats = ['intra_TLS', 'peri_TLS', 'non_TLS']
results_df['HEV_classification'] = pd.Categorical(results_df['HEV_classification'], ordered=True, categories=TLS_HEV_cats)

# --- Add results to adata

adata = adata_dict['square_016um']
obs_df = adata.obs

merge_df = pd.merge(obs_df, results_df, how='left', left_on=group_key, right_on='HEV_group_16um')

# Re-add obsnames
merge_df.index = obs_df.index

adata.obs = merge_df


# --- Plot 1: Plot HEV group classifications on adata

# Define some colors
color_dict_HEV = {
    'non_TLS': 'royalblue',
    'peri_TLS': 'darkorange',
    'intra_TLS': 'indigo'
}

# Plotting
ncols = 3
nrows = 1

sc.set_figure_params(fontsize=15)

fig, axs = plt.subplots(nrows, ncols, figsize=(8*ncols, 6.5*nrows))

# -- First plot TLS groups & HEV classifications
n = 0
ax = axs[n]

sc.pl.spatial(adata_dict['square_048um'], color=group_key, bw=False, alpha_img=0, crop_coord=coords_dict['square_048um'], title=f'{sample_id} ({subtype})\n > {threshold}',
                palette=['lightblue'], ax=ax, legend_loc=None)

sc.pl.spatial(adata, color='HEV_classification', bw=False, alpha_img=0, crop_coord=coords_dict['square_016um'], title=f'{sample_id} ({subtype})\nHEV classification',
                palette=color_dict_HEV, ax=ax)

ax.set_title(ax.get_title(), fontweight='bold')


# -- Same but with the H&E image in the background
n += 1
ax = axs[n]

sc.pl.spatial(adata, color='HEV_classification', bw=False, alpha_img=0.7, crop_coord=coords_dict['square_016um'], title=f'{sample_id} ({subtype})\nHEV classification',
                palette=color_dict_HEV, ax=ax)

ax.set_title(ax.get_title(), fontweight='bold')

# -- Then the score results
n += 1
ax = axs[n]

vmax = 0.8

# Then all scores in single plot
# Use scanpy because squidpy ignores info abt alpha in cmap
score_names = list(clusters_of_i.keys())
score_names_r = score_names[::-1]

for i, score_name in enumerate(score_names_r):
    if i == 0:
        alpha_img = 0.7
    else: 
        alpha_img = 0
    if score_name in TLS_zone_names:
        resolution = 'square_048um'
    else:
        resolution = 'square_016um'
    cmap = cmap_dict[score_name]
    sc.pl.spatial(adata_dict[resolution], color=score_keys_dict[score_name], alpha_img=alpha_img, cmap=cmap,
        crop_coord=coords_dict[resolution], title=f'{sample_id} ({subtype})\nTLS zones', ax=ax, vmin=0, vmax=vmax, colorbar_loc=None,)

ax.set_title(ax.get_title(), fontweight='bold')

# Add colorbars per score
for i, (label, cmap) in enumerate(cmap_dict.items()):
    norm = mpl.colors.Normalize(vmin=0, vmax=vmax)
    spacing = 0.13
    mid_y = 0.5
    num_cbars = len(cmap_dict)
    y_pos = mid_y + (i - (num_cbars - 1) / 2) * spacing  # Distribute around the midpoint
    cbar_ax = ax.inset_axes([1.1, y_pos, 0.3, 0.05])
    cbar = fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), 
                        cax=cbar_ax, orientation='horizontal')
    cbar.ax.set_ylabel(f'{score_keys_dict[label]}', rotation=0, ha='left', va='center', fontsize=18)
    cbar.ax.yaxis.set_label_position("right")
    cbar.ax.yaxis.set_label_coords(1.05, 1.05)
    cbar.ax.tick_params(labelsize=15)

plt.suptitle(f'{sample_id} ({subtype}) - TLS score + TLS_zone scores (top_gene={top_gene}, threshold={threshold}, resolution={resolution}, min_spots_TLS={min_spots})', y=1.1, fontweight='bold', x=0.5)


# --- Some functions for plotting

# Properties for boxplot
whisker_props = {'color': 'black'}
median_props = {'color': 'black'}

alpha = 0.7

# Function to avoid repetition
def plot_seaborn_boxplot(df, ax, x, y, hue, plot_order, stripplot=True, title=None, gap=0, legend=False,
                                annotation_size=10, spot_size=1, box_width=0.2, palette=None, width=0.8, title_size=10,
                                hue_order=None, dodge=True, highlight=None, rot_xtick=True, jitter=True):
    if stripplot:
        showfliers = False
    else:
        showfliers = True
    sns.boxplot(df, x=x, y=y, showfliers=showfliers, linewidth=1, whiskerprops=whisker_props, 
                medianprops=median_props, width=width, ax=ax, palette=palette,
                saturation=1, legend=legend, hue=hue, fliersize=1, order=plot_order, gap=gap,
                hue_order=hue_order)
    # Plot the individual dots on the boxplot
    # Plot dots if you want
    if stripplot:
        sns.stripplot(data=df, x=x, y=y, color="black", order=plot_order, size=spot_size, ax=ax, dodge=dodge, hue=hue, legend=False,
                    hue_order=hue_order, jitter=jitter)
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
    ax.set_title(title, fontweight='bold', fontsize=title_size)
    # move x-axis tick labels to the left
    xlabels = ax.get_xticklabels()
    if rot_xtick:
        ax.set_xticklabels(xlabels, rotation=45, ha='right')
    else:
        ax.set_xticklabels(xlabels, rotation=90)
    # remove the x-axis label
    ax.set_xlabel('')
    # Adjust the positions of the subplots
    ax_pos = ax.get_position()
    # Set the width so it's equal for all plots
    ax.set_position([ax_pos.x0, ax_pos.y0, box_width*len(plot_order), ax_pos.height])

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

# Define a function for annotation the bars
def annotate_bars(ax, df, annotation_size=8, margin=1.05):
    for index, value in enumerate(df):
        text = reformat_number(int(value))
        ax.annotate(text, xy=(index, value), ha='center', va='bottom', fontsize=annotation_size)
    # Adjust ymax to accommodate highest value
    ymax = ax.get_ylim()[1]
    ax.set_ylim(0, ymax * margin)

# Define a function for adding count bars at the top
def add_count_bars(df, target_ax, groupby, color_dict, plot_order, n, margin=1.05, width=0.9, alpha=0.7, height_div=4):
    # Make a new ax
    # add a new subplot iteratively
    ax = plt.subplot(nrows, ncols, n + 1, zorder=-1)
    value_counts = df[groupby].value_counts()
    # Reindex the value_counts series
    sorted_value_counts = value_counts.reindex(plot_order, fill_value=0)
    # Extract palette
    colors = [color_dict[k] for k in plot_order if k in color_dict]
    # Extract colors using list comprehension, using .get to handle missing keys
    sorted_value_counts.plot(kind='bar', legend=False, ax=ax, width=width, grid=False, color=colors, linewidth=1)
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
    # Remove the title from the target ax
    title = target_ax.get_title()
    target_ax.set_title('')
    ax.set_title(title, fontweight='bold')
    # Set the spines to be invisible
    for spine in ax.spines.values():
        spine.set_visible(False)
    # Remove xticks, yticks, and labels
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    # Optional: Remove axis labels if they were set
    ax.set_xlabel('')
    ax.set_ylabel('')
    # Move plot to correct position
    target_pos = target_ax.get_position()
    ax.set_position([target_pos.x0, target_pos.y1, target_pos.width, target_pos.height/height_div])
    # Set xlim to be same
    target_xlim = target_ax.get_xlim()
    ax.set_xlim(target_xlim)
    # Annotate the bars
    annotate_bars(ax, sorted_value_counts, margin=margin)


# --- Plot 2: HEV distance to TLS distribution

sc.set_figure_params(fontsize=10)

ncols = 3
nrows = 1

xval = 'min_dist_um'

# Plotting
fig, axs = plt.subplots(nrows, ncols, figsize=(4*ncols, 2.5*nrows))

# -- Histogram showing distance distribution
n=0
ax = axs[n]

sns.histplot(results_df, x=xval, bins=20, kde=True, ax=ax, alpha=0.8, linewidth=0)

# Change xlab
ax.set_xlabel('Distance (µm)')

# Set grid behind plot
ax.set_axisbelow(True)

# Add title
ax.set_title('HEV distance to TLS\nDistribution', fontweight='bold')

# -- Same but with lines indicating thresholds
# Calculate in um
peri_TLS = peri_TLS_bound * spot_diameter_um
intra_TLS = intra_TLS_bound * spot_diameter_um

positions = [intra_TLS, peri_TLS]

n+=1
ax = axs[n]

sns.histplot(results_df, x=xval, bins=20, kde=True, ax=ax, alpha=0.8, linewidth=0)

# Change xlab
ax.set_xlabel('Distance (µm)')

# Set grid behind plot
ax.set_axisbelow(True)

# Add title
ax.set_title('HEV distance to TLS\nDistribution', fontweight='bold')

# Add some lines indicating where the boundaries are
for pos in positions:
    ax.axvline(x=pos, color='red', alpha=0.5,)

# -- Then a boxplot showing how the split is

n+=1
ax = axs[n]

plot_seaborn_boxplot(results_df, ax=ax, x='HEV_classification', y='min_dist_um', hue='HEV_classification',
                     plot_order=TLS_HEV_cats, palette=color_dict_HEV, dodge=False, hue_order=TLS_HEV_cats, title='HEV distance to TLS\nClassification',
                     rot_xtick=True, spot_size=1.5, jitter=0.15)

ax.set_ylabel('Distance (µm)')

# Add some lines indicating where the boundaries are
for pos in positions:
    ax.axhline(y=pos, color='red', alpha=0.5, zorder=-1)

# Annotate each hline with its y-value
xmax = ax.get_xlim()[1]

for y in positions:
    ax.annotate(f"{y:.1f}", 
                xy=(xmax, y),  # Pointing to (max x, y)
                xytext=(xmax*1.05, y),  # Offset the text slightly
                arrowprops=dict(arrowstyle="-", color="black", shrinkA=0, shrinkB=0, linewidth=0.8),
                va='center', ha='left', fontsize=9)


plt.tight_layout()

add_count_bars(results_df, target_ax=ax, groupby='HEV_classification', color_dict=color_dict_HEV, plot_order=TLS_HEV_cats, n=n, width=0.8,
               height_div=6, margin=1.2)

title_x = (axs[-1].get_position().x1 - axs[0].get_position().x0)/2

plt.suptitle(f'{sample_id} ({subtype}) - Relationship #HEVs vs TLS size (threshold={threshold}, min_spots_TLS={min_spots})', y=1.2, fontweight='bold', x=title_x)


# --- Plot number of HEVs per TLS group

# Filter out what classified as non-TLS
results_df_sub = results_df[results_df['HEV_classification'] != 'non_TLS']

TLS_zone_df = pd.DataFrame(results_df_sub['TLS_zones_group_48um'].value_counts().reset_index())
TLS_zone_df.columns = ['TLS_zones_group_48um', 'HEV_count']

# Calculate TLS zone size as well
TLS_zone_size_df = pd.DataFrame(adata_dict['square_048um'].obs[group_key].dropna().value_counts().reset_index())
TLS_zone_size_df.columns = ['TLS_zones_group_48um', 'TLS_size']

TLS_zone_df = TLS_zone_df.merge(TLS_zone_size_df, on='TLS_zones_group_48um')


# --- Plot 2: Scatterplot TLS zone size vs # HEVs

ncols = 2
nrows = 1

xval = 'HEV_count'
yval = 'TLS_size'

fig, axs = plt.subplots(nrows, ncols, figsize=(3*ncols, 3*nrows))

# -- Plot Relationship TLS size vs number of HEVs

n = 0
ax = axs[n]

sns.scatterplot(TLS_zone_df, x=xval, y=yval, ax=ax, linewidth=0, alpha=0.8,)

# Calculate coordinates for the text (pearson correlation value)
x_coordinate = ax.get_xlim()[0] + 0.05 * (ax.get_xlim()[1] - ax.get_xlim()[0])
y_coordinate = ax.get_ylim()[1] - 0.05 * (ax.get_ylim()[1] - ax.get_ylim()[0])
pearson_corr_val = TLS_zone_df[[xval, yval]].corr(method='pearson')
pearson_corr_val = pearson_corr_val.loc[xval, yval].round(2)
text = ax.text(x_coordinate, y_coordinate, f'r = {pearson_corr_val}', verticalalignment='top', horizontalalignment='left', fontsize=10)

# Set grid behind plot
ax.set_axisbelow(True)

# Add title
ax.set_title('# HEVs vs TLS size', fontweight='bold')

# -- Histogram of number of HEVs per TLS

n+=1
ax = axs[n]

sns.histplot(TLS_zone_df, x=xval, bins=20, kde=True, ax=ax, alpha=0.8, linewidth=0)

# Set grid behind plot
ax.set_axisbelow(True)

# Add title
ax.set_title('# HEVs per TLS', fontweight='bold')

plt.suptitle(f'{sample_id} ({subtype}) - Relationship #HEVs vs TLS size (threshold={threshold}, min_spots_TLS={min_spots})', y=1.05, fontweight='bold', x=0.5)

plt.tight_layout()


# --- Plot 3: HEV size vs TLS classification

# Get the HEV size and add to the df
HEV_size = adata_dict['square_016um'].obs[group_key].value_counts().to_dict()
results_df['HEV_size'] = results_df['HEV_group_16um'].map(HEV_size)

# Plotting
ncols = 1
nrows = 1

fig, ax = plt.subplots(nrows, ncols, figsize=(2, 2))

plot_seaborn_boxplot(results_df, ax=ax, x='HEV_classification', y='HEV_size', palette=color_dict_HEV, plot_order=TLS_HEV_cats, hue='HEV_classification',
                     dodge=False, spot_size=2, stripplot=False, title='HEV size per classification',)

ax.set_ylabel('HEV size\n(# spots 16 um)')

add_count_bars(results_df, target_ax=ax, groupby='HEV_classification', plot_order=TLS_HEV_cats, color_dict=color_dict_HEV, n=0, width=0.7,
               height_div=6, margin=1.4)

plt.suptitle(f'{sample_id} ({subtype}) - Relationship HEV size vs TLS classification (threshold={threshold})', y=1.4, fontweight='bold', x=0.5)


# --- Save output

qc_dir = os.path.join(project_dir, 'results', '04_vascular_microenvironment', '01_identify_HEV_niches', '03b_calculate_HEV_TLS_distance')

today = date.today()
today = today.strftime("%Y%m%d")

threshold_adj = str(threshold).split('.')[-1]
out_dir = os.path.join(qc_dir, 'square_016um', f'threshold_{threshold_adj}')

os.makedirs(out_dir, exist_ok=True)

# Save plots
out_dir_plots = os.path.join(out_dir, 'plots', today)
os.makedirs(out_dir_plots, exist_ok=True)

# Also save the groupings as a .csv
out_dir_csv = os.path.join(out_dir, 'csv_files', today)
os.makedirs(out_dir_csv, exist_ok=True)

# Save plots
filename = f'{sample_id}_{score_of_i}_TLS_classification.pdf'

save_multi_image(os.path.join(out_dir_plots, filename))
plt.close('all')

# And save the distance calculations + classification
filename = f'{sample_id}_HEV_closest_TLS_distance_classifications.csv'
results_df.to_csv(os.path.join(out_dir_csv, filename), index=False)

# Also save the number of HEVs per TLS
filename = f'{sample_id}_HEV_per_TLS_group.csv'
TLS_zone_df.to_csv(os.path.join(out_dir_csv, filename), index=False)


# --- Finish

print("Script successfully completed")
