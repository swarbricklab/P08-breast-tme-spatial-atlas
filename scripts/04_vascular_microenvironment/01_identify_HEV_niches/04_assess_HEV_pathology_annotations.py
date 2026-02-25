#!/usr/bin/env python3

# --- Aim
# Compare HEV group locations to pathology annotations.
# For each sample, determine the dominant pathology annotation associated with each HEV group.
# (Note: pathology annotations were done on 8um resolution, HEVs on 16um, so give each HEV group the dominant annotation of the 8um spots that fall within its area)

# --- Define environment
# env: vis-hd

# --- Load packages
import argparse
import math
import os
import sys
from datetime import date, datetime

import anndata as ad
import matplotlib as mpl
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
import scanpy.plotting.palettes as sc_palettes

import seaborn as sns

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
parser.add_argument("--sample_sheet", help="Path to sample sheet CSV file", default=None)

parser.add_argument("--min_count", help="Min count for QC", type=int)
parser.add_argument("--min_gene", help="Min gene for QC", type=int)
parser.add_argument("--min_spots", help="Minimal spots per gene", type=int)

args = vars(parser.parse_args())


# --- Arguments

res_dir = args["res_dir"]
h5ad_dir=args["h5ad_dir"]
sample_id=args["sample_id"]
project_dir=args["project_dir"]
sys.path.append(os.path.join(project_dir, 'scripts'))
from helper_functions import auto_crop_scanpy, save_multi_image, generate_color_variations



# --- Get some info from the sample sheet
# Use sample_sheet argument if provided, else default
if args.get("sample_sheet"):
    sample_sheet_path = args["sample_sheet"]
else:
    sample_dir = os.path.join(project_dir, 'config')
    sample_sheet_file = 'sample_sheet.csv'
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
# Path annotations were made on 8µm resolution
resolutions = ['square_016um', 'square_008um']

adata_dict = {}

for resolution in resolutions:
    adata_dict[resolution] = ad.read_h5ad(os.path.join(h5ad_dir, 'visium_hd', resolution, f'{sample_id}.h5ad'))


# --- Calculate crop coords for optimal presentation
coords_dict  = {}

for resolution in resolutions:
    coords = auto_crop_scanpy(adata_dict[resolution], border_size=0.1)
    coords_dict[resolution] = coords


# --- Load the HEV groups 

resolution = 'square_016um'

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

# Decide whether to use the expanded groups
TLS_group_dir = os.path.join(res_dir, '04_vascular_microenvironment', '01_identify_HEV_niches', '02b_expand_HEV_groups', resolution, f'threshold_{threshold_adj}',)
TLS_group_dir = os.path.join(TLS_group_dir, 'csv_files')

most_recent_dir_path = find_most_recent(TLS_group_dir)

# Read in the results & add to anndata
score_of_i = 'HEVs'
TLS_file_name = f'{sample_id}_{score_of_i}_groups.csv'
TLS_results = pd.read_csv(os.path.join(most_recent_dir_path, TLS_file_name), index_col=0)

group_key = 'group_too_small_removed'
group_key_source = 'expanded_groups'

adata_dict['square_016um'].obs[group_key] = TLS_results[group_key_source]

# --- Also load the HEV TLS classifications

TLS_group_dir = os.path.join(res_dir, '04_vascular_microenvironment', '01_identify_HEV_niches', '03b_calculate_HEV_TLS_distance', 'square_016um', f'threshold_{threshold_adj}', 'csv_files')
most_recent_dir_path = find_most_recent(TLS_group_dir)

# Read in the results & add to anndata
TLS_file_name = f'{sample_id}_HEV_closest_TLS_distance_classifications.csv'
HEV_TLS_classifications = pd.read_csv(os.path.join(most_recent_dir_path, TLS_file_name))

# Add the HEV TLS classification to adata
class_dict = HEV_TLS_classifications.set_index('HEV_group_16um')['HEV_classification'].to_dict()

adata_dict['square_016um'].obs['HEV_TLS_classification'] = adata_dict['square_016um'].obs[group_key].map(class_dict)

# Define some colors
color_dict_HEV = {
    'non_TLS': 'royalblue',
    'peri_TLS': 'darkorange',
    'intra_TLS': 'indigo'
}

# --- Load the reworked annotations

annotation_dir = os.path.join(project_dir, 'data', 'pathology_annotations', 'visium_hd')
filename = f'{sample_id}_pathology_annotations.csv'

annotation_key = 'path_anno'

annotations = pd.read_csv(os.path.join(annotation_dir, filename), index_col=0)
adata_dict['square_008um'].obs[annotation_key] = annotations


# --- Get some colors for the pathology annotations

adata = adata_dict['square_008um']

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
all_cats = list(adata.obs[annotation_key].dropna().unique())

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

# Set the plot_order
adata.obs[annotation_key] = pd.Categorical(adata.obs[annotation_key], categories=plot_order, ordered=True)


# --- Plot 1: Spatial plots pathology annotations & HEV groups

# Plotting
sc.set_figure_params(fontsize=20)

ncols = 4
nrows = 1

# Make the figure
plt.figure(figsize=(8*ncols, 7*nrows))
n=0
ax = plt.subplot(nrows, ncols, n + 1)

# First plot the H&E
sc.pl.spatial(adata_dict['square_016um'], color=None, bw=False, alpha_img=1, crop_coord=coords_dict['square_016um'], title=f'{sample_id} - H&E', ax=ax)
ax.set_title(ax.get_title(), fontweight='bold')

# Then plot the reworked annotations
n+= 1
ax = plt.subplot(nrows, ncols, n + 1)

sc.pl.spatial(adata, color=annotation_key, title=f'Pathology annotations', crop_coord=coords_dict['square_008um'], ax=ax, palette=color_dict,
              na_in_legend=False, alpha_img=0.7)
ax.set_title(ax.get_title(), fontweight='bold')

# Then plot pathology annotations + HEV groups on top
n+= 1
ax = plt.subplot(nrows, ncols, n + 1)

sc.pl.spatial(adata, color=annotation_key, title=f'Pathology annotations', crop_coord=coords_dict['square_008um'], ax=ax, palette=color_dict,
              na_in_legend=False, alpha_img=0.7)

sc.pl.spatial(adata_dict['square_016um'], color=group_key, title=f'Pathology annotations + HEVs', crop_coord=coords_dict['square_016um'], ax=ax, palette=['indigo'],
              na_in_legend=False, alpha_img=0, legend_loc=None)
ax.set_title(ax.get_title(), fontweight='bold')


# Then plot pathology annotations + TLS classification

n+= 1
ax = plt.subplot(nrows, ncols, n + 1)

sc.pl.spatial(adata, color=annotation_key, title=f'Pathology annotations', crop_coord=coords_dict['square_008um'], ax=ax, palette=color_dict,
              na_in_legend=False, alpha_img=0.7, legend_loc=None)

sc.pl.spatial(adata_dict['square_016um'], color='HEV_TLS_classification', title=f'HEV TLS classification', crop_coord=coords_dict['square_016um'], ax=ax, palette=color_dict_HEV,
              na_in_legend=False, alpha_img=0,)
ax.set_title(ax.get_title(), fontweight='bold')


plt.subplots_adjust(wspace=0.7)

plt.suptitle(f'{sample_id} ({subtype}) - Pathology annotations + HEV groups (resolution={resolution})', y=1.1, fontweight='bold')


# --- Find the dominant path annotations per HEV group


def get_pathology_counts_per_group(adata_008um, adata_016um, group_key="groups", pathology_key="pathology", sample_id=None):
    """
    Computes the count of each pathology annotation for each 16µm group using corresponding 8µm spots.

    Parameters:
    - adata_008um: AnnData object at 8µm resolution (contains pathology annotations)
    - adata_016um: AnnData object at 16µm resolution (contains spot groupings)
    - group_key: Column in `adata_016um.obs` defining spot groups (e.g., "group1", "group2")
    - pathology_key: Column in `adata_008um.obs` defining pathology annotations
    - sample_id: Sample ID to extract `microns_per_pixel` from `adata_016um.uns['spatial']`

    Returns:
    - DataFrame where each row is a 16µm group and columns are counts of pathology annotations.
    """
    # Ensure spatial coordinates exist
    if "spatial" not in adata_008um.obsm or "spatial" not in adata_016um.obsm:
        raise ValueError("Spatial coordinates not found in AnnData objects.")
    # Convert spatial coordinates to DataFrame
    df_008 = pd.DataFrame(adata_008um.obsm["spatial"], index=adata_008um.obs_names, columns=["x", "y"])
    df_016 = pd.DataFrame(adata_016um.obsm["spatial"], index=adata_016um.obs_names, columns=["x", "y"])
    # Add pathology annotations to df_008
    df_008["pathology"] = adata_008um.obs[pathology_key]
    # Add group annotations to df_016
    df_016["group"] = adata_016um.obs[group_key]
    # Drop NaN groups (unlabeled spots)
    df_016 = df_016.dropna(subset=["group"])
    # Get the microns per pixel
    microns_per_pixel = adata_016um.uns['spatial'][sample_id]['scalefactors']['microns_per_pixel']
    # Dictionary to store dominant pathology per group
    group_to_pathology = {}
    # Calculate the diagonal using pythagoras
    spot_dist_pixel=16
    c = pow(spot_dist_pixel/2, 2) + pow(spot_dist_pixel/2, 2)
    c = math.sqrt(c)
    # Convert to pixels
    diagonal = c/microns_per_pixel
    # Dictionary to store pathology counts per group
    group_pathology_counts = {}
    # Iterate over unique groups
    for group, group_spots in df_016.groupby("group"):
        # Get spatial coordinates of all spots in this group
        group_coords = group_spots[["x", "y"]].values
        # Initiate empty df
        corresponding_8um_spots = pd.DataFrame()
        # Find 8µm spots that fall within each 16µm spot's bounding box
        for x_16, y_16 in group_coords:
            # Find all 8um spots that are inside this 16um spot (2x2 grid)
            mask = (
                (df_008["x"] >= x_16 - diagonal) & (df_008["x"] < x_16 + diagonal) &
                (df_008["y"] >= y_16 - diagonal) & (df_008["y"] < y_16 + diagonal)
            )
            # Extract corresponding 8um spots
            corresponding_8um_spots = pd.concat([corresponding_8um_spots, df_008[mask]])
        # Count pathology annotations
        pathology_counts = corresponding_8um_spots["pathology"].value_counts().to_dict()
        # Store counts
        group_pathology_counts[group] = pathology_counts
    # Convert to DataFrame
    df_group_pathology = pd.DataFrame.from_dict(group_pathology_counts, orient="index").fillna(0).astype(int)
    # Reset index for a clean output
    df_group_pathology.index.name = "group"
    df_group_pathology.reset_index(inplace=True)
    return df_group_pathology


pathology_counts = get_pathology_counts_per_group(adata_008um=adata_dict['square_008um'], adata_016um=adata_dict['square_016um'],
                          group_key=group_key, pathology_key=annotation_key, sample_id=sample_id)

# Now summarize to just the most dominant category per group
dominant_pathology = pathology_counts.set_index("group").idxmax(axis=1)
dominant_count = pathology_counts.set_index("group").max(axis=1)

# Create a summary DataFrame
df_summary = pd.DataFrame({
    "group": pathology_counts["group"],
    "dominant_pathology": dominant_pathology.values,
    "count": dominant_count.values
})


# --- Plot 2: Plot the HEV dominant annotations spatially to check it makes sense

# Make the colors a bit lighter to see the HEVs easier
def lighten_color(color, factor=0.5):
    rgb = mcolors.to_rgb(color)  # Convert to RGB tuple
    light_rgb = [(1 - factor) * c + factor for c in rgb]  # Blend with white
    return light_rgb

color_dict_lighter = {}

for score_name, color in color_dict.items():
    color_dict_lighter[score_name] = lighten_color(color)

# Add the dominant pathology annotation back to adata
path_dict = df_summary.set_index('group')['dominant_pathology'].to_dict()
adata_dict['square_016um'].obs['dominant_path'] = adata_dict['square_016um'].obs[group_key].map(path_dict)

ncols = 2
nrows = 1

# Make the figure
plt.figure(figsize=(8*ncols, 7*nrows))
n=0
ax = plt.subplot(nrows, ncols, n + 1)

# First plot the H&E
sc.pl.spatial(adata_dict['square_016um'], color=None, bw=False, alpha_img=1, crop_coord=coords_dict['square_016um'], title=f'{sample_id} - H&E', ax=ax)
ax.set_title(ax.get_title(), fontweight='bold')

# Then plot pathology annotations + HEV groups on top
n+= 1
ax = plt.subplot(nrows, ncols, n + 1)

sc.pl.spatial(adata_dict['square_008um'], color=annotation_key, title=f'Pathology annotations', crop_coord=coords_dict['square_008um'], ax=ax, palette=color_dict_lighter,
              na_in_legend=False, alpha_img=0.7, legend_loc=None)

sc.pl.spatial(adata_dict['square_016um'], color='dominant_path', title=f'HEVs\Pathology annotation', crop_coord=coords_dict['square_016um'], ax=ax, palette=color_dict,
              na_in_legend=False, alpha_img=0,)

ax.set_title(ax.get_title(), fontweight='bold')


plt.suptitle(f'{sample_id} ({subtype}) - Dominant pathology annotation per HEV group', y=1.2, fontweight='bold')


# --- Prepare for plotting: define some functions

# Prepare for plotting
sc.set_figure_params(fontsize=10)
sns.set_style("ticks")

# Function for beautifying heatmap
def beautify_heatmap(g, col_colors=True, row_colors=True, cbar_title='Pearson correlation',
                     fontsize_title=13, label_size=12, cbar_title_pad=-60, cbar_pad=6, title=None, x_rot=90,
                     cbar_height_square=5, space_per_square=0.07, title_ax_n=None, cbar_width_square=0.1):
    # Get the current figure
    fig = plt.gcf()
    # Extract the heatmap ax
    heatmap_ax = g.ax_heatmap
    # -- Alter the main heatmap axis
    heatmap_ax.set_ylabel('')
    heatmap_ax.set_xlabel('')
    heatmap_ax.tick_params(axis='x', labelsize=label_size, rotation=x_rot)  # Set font size for xticklabels
    heatmap_ax.tick_params(axis='y', labelsize=label_size)  # Set font size for yticklabels
    heatmap_ax.set_facecolor('darkgray') # Set facecolor to gray so masked values appear in gray
    for spine in heatmap_ax.spines.values():
        spine.set_visible(True)
    heatmap_ax_pos = heatmap_ax.get_position()
    # --- Alter the top dendrogram (put the title above of the dendrogram)
    top_dendrogram_ax = g.ax_col_dendrogram
    if title_ax_n is None:
        heatmap_ax.set_title(f'{title}', fontsize=fontsize_title, fontweight='bold')
    else:
        fig.axes[title_ax_n].set_title(f'{title}', fontsize=fontsize_title, fontweight='bold')
    left_dendrogram_ax = g.ax_row_dendrogram
    # -- Adjust colorbar columns
    if col_colors:
        col_colors_ax = g.ax_col_colors
        for spine in col_colors_ax.spines.values():
            spine.set_visible(True)
        # Move column colors & dendrogram up
        top_dendrogram_ax_pos = top_dendrogram_ax.get_position()
        col_colors_ax_pos = col_colors_ax.get_position()
        col_colors_ax.set_position((heatmap_ax_pos.x0, heatmap_ax_pos.y1 + 0.01, col_colors_ax_pos.width, col_colors_ax_pos.height))
        top_dendrogram_ax.set_position((heatmap_ax_pos.x0, heatmap_ax_pos.y1 + col_colors_ax_pos.height + 0.01, top_dendrogram_ax_pos.width, top_dendrogram_ax_pos.height))
     # -- Adjust row color
    if row_colors:
        row_colors_ax = g.ax_row_colors
        for spine in row_colors_ax.spines.values():
            spine.set_visible(True)
    # Adjust heatmap ax based on number of samples
    xlabels = heatmap_ax.get_xticklabels()
    ylabels = heatmap_ax.get_yticklabels()
    new_width = space_per_square * len(xlabels)
    new_height = space_per_square * len(ylabels)
    heatmap_ax.set_position((heatmap_ax_pos.x0, heatmap_ax_pos.y0, new_width, new_height))
    # Adjust position other axes
    new_heatmap_pos = heatmap_ax.get_position()
    if col_colors:
        ylabels_col = col_colors_ax.get_yticklabels()
        new_height = space_per_square * len(ylabels_col) / 2
        col_colors_ax.set_position((new_heatmap_pos.x0, new_heatmap_pos.y1 + 0.01, new_heatmap_pos.width, new_height))
        new_col_pos = col_colors_ax.get_position()
        top_dendrogram_ax.set_position((new_col_pos.x0, new_col_pos.y1, new_col_pos.width, top_dendrogram_ax_pos.height))
    if row_colors:
        xlabels_row = row_colors_ax.get_xticklabels()
        new_width = space_per_square * len(xlabels_row) / 2
        row_colors_ax_pos = row_colors_ax.get_position()
        row_colors_ax.set_position((new_heatmap_pos.x0 - row_colors_ax_pos.width - 0.05, new_heatmap_pos.y0, new_width,
                                    new_heatmap_pos.height))
    left_dendrogram_ax_pos = left_dendrogram_ax.get_position()
    left_dendrogram_ax.set_position((new_heatmap_pos.x0 - left_dendrogram_ax_pos.width, new_heatmap_pos.y0, left_dendrogram_ax_pos.width,
                                     new_heatmap_pos.height))
    # -- Adjust colorbar axis
    cbar_ax = g.ax_cbar
    cbar_ax.set_ylabel(cbar_title, labelpad=10, fontsize=label_size)
    cbar_ax.yaxis.set_label_position("left")
    cbar_ax.tick_params(axis='y', labelsize=label_size)  # Set font size for yticklabels
    for spine in cbar_ax.spines.values():
        spine.set_visible(True)
    # Adjust position
    cbar_ax_pos = cbar_ax.get_position()
    heatmap_ax_pos = heatmap_ax.get_position()
    cbar_ax.set_position((heatmap_ax_pos.x1 + space_per_square*cbar_pad, heatmap_ax_pos.y0, space_per_square*cbar_width_square, space_per_square*cbar_height_square))
    return fig.axes


# --- Plot 3: Barplots showing classification

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

# Define a function for creating barplots
def plot_seaborn_barplot(df, plot_key, ax, groupby='sample', title=None, offset=0.02, annotate=False, annotation_size=8, plot_order=None,
                         color_dict=None, legend=True, bar_width=1, margin=1.03):
    sns.barplot(df, hue=groupby, legend=legend, x=groupby, y=plot_key, ax=ax, linewidth=1, palette=color_dict, 
                order=plot_order, saturation=1)
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
    # Set title
    ax.set_title(title, fontweight='bold')
    # move x-axis tick labels to the left
    xlabels = ax.get_xticklabels()
    ax.set_xticklabels(xlabels, rotation=45, ha='right')
    # remove the x-axis label
    ax.set_xlabel('')
    if n > 0:
        ax.set_ylabel('')
    # Annotate the bars
    xlabels = [label.get_text() for label in ax.get_xticklabels()]
    if annotate:
        for index in xlabels:
            value = df.set_index(groupby).loc[index, plot_key]
            reformat_value = reformat_number(value)
            ax.annotate(str(reformat_value), xy=(index, value), ha='center', va='bottom', fontsize=annotation_size)
    # Adjust the positions of the subplots
    ax_pos = axs[n].get_position()
    # Shift the subplot to the desired location + set the width so it's equal for all plots
    axs[n].set_position([ax_pos.x0, ax_pos.y0, bar_width*len(xlabels), ax_pos.height])
    # Adjust y-axis limits to accommodate the highest value
    if annotate:
        ymax = ax.get_ylim()[1]
        ax.set_ylim(0, ymax * margin)

# Properties for boxplot
whisker_props = {'color': 'black'}
median_props = {'color': 'black'}

alpha = 0.7

# Function to avoid repetition
def plot_seaborn_boxplot(df, ax, x, y, hue, plot_order=None, stripplot=True, title=None, gap=0, legend=False,
                                annotation_size=10, spot_size=1, box_width=0.2, palette=None, width=0.8, title_size=10,
                                hue_order=None, dodge=True, highlight=None):
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
                    hue_order=hue_order)
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
    ax.set_xticklabels(xlabels, rotation=45, ha='right')
    # remove the x-axis label
    ax.set_xlabel('')
    # Adjust the positions of the subplots
    ax_pos = ax.get_position()
    # Set the width so it's equal for all plots
    xlabels = [label.get_text() for label in ax.get_xticklabels()]
    ax.set_position([ax_pos.x0, ax_pos.y0, box_width*len(xlabels), ax_pos.height])

# Function to move plots so they don't overlap (only move in x direction)
def move_plots(axs, offset=0.02):
    for n, ax in enumerate(axs):
        # Adjust the positions of the subplots
        ax_pos = axs[n].get_position()
        # Shift the subplot to the desired location + set the width so it's equal for all plots
        if n == 0:
            axs[n].set_position([ax_pos.x0, ax_pos.y0, ax_pos.width, ax_pos.height])
        else:
            new_pos = ax_pos_prev.x0 + ax_pos_prev.width + offset
            axs[n].set_position([new_pos, ax_pos.y0, ax_pos.width, ax_pos.height])
        ax_pos_prev = axs[n].get_position()


# Add the TLS classifications
df_summary['HEV_TLS_classification'] = df_summary['group'].map(class_dict)

# Also add HEV size
size_dict = HEV_TLS_classifications.set_index('HEV_group_16um')['HEV_size'].to_dict()
df_summary['HEV_size'] = df_summary['group'].map(size_dict)

# Count number of occurences
plot_df = df_summary['dominant_pathology'].value_counts().to_frame()

# Ensure all columns exist
plot_df = plot_df.reindex(index=plot_order, fill_value=0)
plot_df = plot_df.reset_index()

# Plotting
nrows = 1
ncols = 2

fig, axs = plt.subplots(nrows, ncols, figsize=(3, 2))

n=0
ax = axs[n]

plot_seaborn_barplot(plot_df, plot_key='count', groupby='dominant_pathology', ax=ax, bar_width=0.1, annotate=True,
                     legend=False, title='HEV location', margin=1.1, plot_order=plot_order, color_dict=color_dict)

# Boxplot HEV size
n+=1
ax = axs[n]

plot_seaborn_boxplot(df_summary, ax=ax, x='dominant_pathology', y='HEV_size', hue='dominant_pathology',
                     palette=color_dict, box_width=0.1, stripplot=False, title='HEV size', plot_order=plot_order)

# Move plots so they don't overlap
move_plots(axs, offset=0.2)

# Set title
title_x = (axs[-1].get_position().x0 + axs[0].get_position().x1)/2

plt.suptitle(f'{sample_id} ({subtype}) - HEV location & size pathology annotations', y=1.2, x=title_x, fontweight='bold')


# --- Plot 4: Heatmap TLS classification + dominant pathology annotation

# Count occurrences of each pair
pair_counts = df_summary.groupby(['dominant_pathology', 'HEV_TLS_classification']).size().unstack(fill_value=0)

# Ensure all columns exist
pair_counts = pair_counts.reindex(index=plot_order, columns=list(color_dict_HEV.keys())[::-1], fill_value=0)

# Make heatmap
g = sns.clustermap(pair_counts, row_cluster=False, col_cluster=False, cmap='Blues', linewidths=0.5, figsize=(4, 4), annot=True,
               colors_ratio=(0.01, 0.1), linecolor='black', xticklabels=True, yticklabels=True,)

axs = beautify_heatmap(g, row_colors=False, col_colors=False, cbar_title=None, label_size=10, fontsize_title=11, title=f'{sample_id}',
                    title_ax_n=2, cbar_pad=0.1, cbar_height_square=len(pair_counts.index), cbar_width_square=0.3, space_per_square=0.1)

heatmap_ax = g.ax_heatmap

# Move y-axis labels to the left
heatmap_ax.yaxis.set_label_position("left")  # Move axis label
heatmap_ax.yaxis.tick_left()  # Move tick labels

for spine in heatmap_ax.spines.values():
    spine.set_linewidth(1)  # Change to your preferred thickness
    spine.set_visible(True)  # Change to your preferred thickness

heatmap_ax.set_title(f'{sample_id}\nHEV location', fontweight='bold')

# Add labels
heatmap_ax.set_xlabel('TLS classification')
heatmap_ax.set_ylabel('Pathology annotation')

# Calculate the midpoint to position the title
title_x = (heatmap_ax.get_position().x0 + heatmap_ax.get_position().x1)/2
title_y = heatmap_ax.get_position().y1 + 0.2

plt.suptitle(f'{sample_id} ({subtype}) - HEV location pathology annotations vs TLS distance classification', y=title_y, x=title_x, fontweight='bold')


# --- Plot 5: Stacked barplot pathology annotation per HEV TLS class

# Prepare the df for the plot
df_long = pair_counts.reset_index().melt(id_vars="dominant_pathology", var_name="HEV_TLS_classification", value_name="count")

# Calculate rowsums
rowsums = pd.DataFrame(pair_counts.T.sum(axis=1))

# Also make a proportion df
proportion_df = pair_counts.T.div(pair_counts.T.sum(axis=1), axis=0)

nrows = 1
ncols = 2

# Make a function for ease
def plot_composition(pivot_df, ax, ylabel=None, title=None, annotate=True, legend=True, margin=1.03, palette=None):
    pivot_df.plot(kind='bar', stacked=True, ax=ax, grid=False, width=0.8, linewidth=0.2, edgecolor='black', legend=legend,
                  color=palette)
    # Annotate the bars
    if annotate:
        for index, value in enumerate(rowsums[0]):
            ax.annotate(str(int(value)), xy=(index, value), ha='center', va='bottom', fontsize=8)
    # Rotate xlabels
    xlabels = ax.get_xticklabels()
    ax.set_xticklabels(xlabels, rotation=45, ha='right')
    # Set axis labels
    ax.set_xlabel('')
    ax.set_ylabel(ylabel)
    # Set title
    ax.set_title(title, fontweight='bold')
    if legend:
        # Move legend
        ax.legend(loc='upper left', bbox_to_anchor=(1.0, 1.0), fancybox=False, frameon=False, handlelength=1.5)
        # Adjust y-axis limits to accommodate the highest value
    if annotate:
        ymax = ax.get_ylim()[1]
        ax.set_ylim(0, ymax * margin)

sc.set_figure_params(fontsize=10)

# Get some colors
# Map colors to index order
colors = [color_dict[idx] for idx in pair_counts.index]


fig, axs = plt.subplots(nrows, ncols, figsize=(2*ncols, 2*nrows), gridspec_kw={'hspace':0.5})

plot_composition(pair_counts.T, ax=axs[0], ylabel='Count', title='HEVs', legend=False, margin=1.05, palette=colors)
plot_composition(proportion_df, ax=axs[1], ylabel='Proportion', title='HEVs', annotate=False, palette=colors)

plt.subplots_adjust(wspace=1)

plt.suptitle(f'{sample_id} ({subtype}) - HEV location - pathology annotations vs TLS distance classification', y=1.2, fontweight='bold')


# --- Save the output

# Put in a timestamped folder to avoid overwriting older plots
qc_dir = os.path.join(res_dir, '04_vascular_microenvironment', '01_identify_HEV_niches', '04_assess_HEV_pathology_annotations')
os.makedirs(qc_dir, exist_ok=True)

today = date.today()
today = today.strftime("%Y%m%d")

out_dir = os.path.join(qc_dir)
os.makedirs(out_dir, exist_ok=True)

# Save plots
out_dir_plots = os.path.join(out_dir, 'plots', today)
os.makedirs(out_dir_plots, exist_ok=True)

filename = f'{sample_id}_HEVs_per_pathology_annotation.pdf'

save_multi_image(os.path.join(out_dir_plots, filename))

# And save the csv files
out_dir_csv = os.path.join(out_dir, 'csv_files', today)
os.makedirs(out_dir_csv, exist_ok=True)

filename = f'{sample_id}_dominant_pathology_per_HEV_group.csv'
df_summary.to_csv(os.path.join(out_dir_csv, filename), index=None)

# Also save the number of path labels per HEV group
filename = f'{sample_id}_pathology_counts_per_HEV_group.csv'
pathology_counts.to_csv(os.path.join(out_dir_csv, filename), index=None)


# --- Finish

print("Script successfully completed")
