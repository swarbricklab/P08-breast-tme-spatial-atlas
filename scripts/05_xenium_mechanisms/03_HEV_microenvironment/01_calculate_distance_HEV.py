#!/usr/bin/env python3

# --- Aim
# Calculate the distance from each cell to the nearest HEV group in Xenium data.
# Enable spatial profiling of the HEV microenvironment and proximity-based enrichment analyses.
# Also enable classification by distance to TLS.

# --- Define environment
# env: xen-hev

# --- Load packages
import argparse
import os
import sys
from datetime import datetime

import anndata as ad
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

import pandas as pd
import numpy as np
import scanpy as sc

import seaborn as sns
from scipy.spatial import cKDTree

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

args = vars(parser.parse_args())


# --- Arguments

res_dir = args["res_dir"]
h5ad_dir=args["h5ad_dir"]
sample_id=args["sample_id"]
project_dir=args["project_dir"]

# Import custom functions
sys.path.append(os.path.join(project_dir, 'scripts'))
from helper_functions import auto_crop_xenium, save_multi_image


# --- Get some info from the sample sheet
# Use sample_sheet argument if provided, else default
if args.get("sample_sheet"):
    sample_sheet_path = args["sample_sheet"]
else:
    sample_dir = os.path.join(project_dir, 'config')
    sample_sheet_file = 'sample_sheet_xenium.csv'
    sample_sheet_path = os.path.join(sample_dir, sample_sheet_file)

sample_sheet = pd.read_csv(sample_sheet_path)
tissue_id = sample_sheet[sample_sheet['section_id'] == sample_id]['tissue_id'].values[0]


# --- Also get some tissue metadata

resource_dir = os.path.join(project_dir, 'resources')
metadata_file = 'tissue_metadata.csv'

tissue_metadata = pd.read_csv(os.path.join(resource_dir, metadata_file))

# Create a sample - subtype dictionary
subtype = tissue_metadata[tissue_metadata['tissue_id'] == str(tissue_id)]['clinical_subtype'].values[0]


# --- Load data

adata = ad.read_h5ad(os.path.join(h5ad_dir, 'xenium', f'{sample_id}.h5ad'))


# --- Calculate crop coordinates for optimal presentation

crop_coord = auto_crop_xenium(adata)


# --- Load the HEV niches

# Find most recent results
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
        return most_recent_dir_path
    else:
        print("No valid directories found.")

hev_dir = os.path.join(res_dir, '05_xenium_mechanisms', '02_identify_HEV_regions')
most_recent = find_most_recent(hev_dir)

hev_file = f'{sample_id}_HEV_niches_cell_to_niche.csv'
tmp_df = pd.read_csv(os.path.join(most_recent, 'csv_files', hev_file), index_col=0)
adata.obs['niche_id'] = tmp_df['niche_id']
adata.obs['HEV_score'] = tmp_df['HEV_score']

# Also load info whether HEV is near TLS or not
filename = f'{sample_id}_HEV_niches_info.csv'
tmp_df = pd.read_csv(os.path.join(most_recent, 'csv_files', filename),)

# Make a dict close to TLS vs not
near_tls_dict = dict(zip(tmp_df['niche_id'], tmp_df['close_to_TLS']))

# Map to adata
adata.obs['HEV_in_TLS'] = adata.obs['niche_id'].map(near_tls_dict)


# --- Also load the TLS annotations

tls_dir = os.path.join(res_dir, '05_xenium_mechanisms', '01_identify_TLS_regions', 'csv_files')

tls_file = f'{sample_id}_cell_TLS_annotations.csv'

tmp_df = pd.read_csv(os.path.join(tls_dir, tls_file), index_col=0)
adata.obs['in_TLS'] = tmp_df['in_TLS']


# --- Prepare for saving

out_dir = os.path.join(res_dir, '05_xenium_mechanisms', '03_HEV_microenvironment', '01_calculate_distance')

today = datetime.now().strftime('%Y%m%d')
out_dir_plots = os.path.join(out_dir, today, 'plots')
os.makedirs(out_dir_plots, exist_ok=True)

out_dir_csv = os.path.join(out_dir, today, 'csv_files')
os.makedirs(out_dir_csv, exist_ok=True)


# --- Plot 1: Plot the identified HEV niches

# Set fontsize
sc.set_figure_params(fontsize=20)

ncols = 3
fig, axs = plt.subplots(nrows=1, ncols=ncols, figsize=(9*ncols, 7))

# -- First plot just the H&E
n = 0
ax = axs[n]

# Create an emtpy adata to speed things up
adata_sub = adata[[]].copy()
sc.pl.spatial(adata_sub, color=None, title=f'{sample_id} ({subtype})',
            crop_coord=crop_coord, ax=ax, na_in_legend=False, show=False)

ax.set_title(ax.get_title(), fontweight='bold')

# -- Then plot the HEV niches on top

n = 1
ax = axs[n]

# Subset adata to only spots in niches
adata_sub = adata[adata.obs['niche_id'].notna()]
sc.pl.spatial(adata_sub, color='niche_id', title=f'HEV niches', palette=['blue'],
            crop_coord=crop_coord, ax=ax, na_in_legend=False, legend_loc=None, show=False,
            alpha_img=0.5, size=0.5)

ax.set_title(ax.get_title(), fontweight='bold')

# -- Then plot HEVs + TLS niches
n = 2
ax = axs[n]

# Subset adata to only TLS spots
adata_sub = adata[adata.obs['in_TLS'] == True]
adata_sub.obs['in_TLS'] = adata_sub.obs['in_TLS'].astype('category')
sc.pl.spatial(adata_sub, color='in_TLS', title=f'HEV niches in TLS', palette=['cornflowerblue'],
            crop_coord=crop_coord, ax=ax, na_in_legend=False, legend_loc=None, show=False,
            alpha_img=0.7, size=0.5)

adata_sub = adata[adata.obs['niche_id'].notna()]
sc.pl.spatial(adata_sub, color='niche_id', title=f'HEV + TLS niches', palette=['blue'],
            crop_coord=crop_coord, ax=ax, na_in_legend=False, legend_loc=None, show=False,
            alpha_img=0, size=0.5)       

ax.set_title(ax.get_title(), fontweight='bold')

# Add legend with TLS & HEV
# Make the markers circles
color_dict = {
    'TLS': 'cornflowerblue',
    'HEV niche': 'blue',
}

# Create custom legend handles with circles
handles = [
    Line2D([0], [0], marker='o', color='w', markerfacecolor=color, markersize=10, label=label)
    for label, color in color_dict.items()
]

# Remove all x & ylabels
for ax in axs:
    ax.set_xlabel('')
    ax.set_ylabel('')

# Adjust the positioning of the axes
for n, ax in enumerate(axs):
    ax_pos = ax.get_position()
    ax_width = ax_pos.width
    if n == 0:
        ax_pos_prev = ax_pos
    else:
        new_x = ax_pos_prev.x1 + ax_width*0.1
        ax.set_position([new_x, ax_pos.y0, ax_pos.width, ax_pos.height])
        ax_pos_prev = ax.get_position()

plt.legend(handles=handles, loc='center left', bbox_to_anchor=(1, 0.5), fontsize=16, frameon=False)

# Add a title
plt.suptitle(f'{sample_id} ({subtype}) - HEV niches (n={adata_sub.obs["niche_id"].nunique()}) identified based on HEV scoring', y=1.05, fontweight='bold')


# --- Extract the coordinates in um

coords = adata.obs[['x_centroid', 'y_centroid']].copy()
# Rename columns
coords.columns = ['x', 'y']
groups = adata.obs['niche_id']

# Convert to DataFrame for easier manipulation
df = pd.DataFrame(coords, columns=["x", "y"])
df["HEV_group"] = groups.values


# --- For each cell - calculate the distance to the closest HEV group

def calculate_distance(df, group_key="HEV_group", x_col="x", y_col="y"):
    """
    Calculate the closest distance for each row in the DataFrame to any group recorded in the `group_key` column.

    Parameters:
        df (pd.DataFrame): A DataFrame with coordinates (`x`, `y`) and group labels (`group_key`).
        group_key (str): The column name for group labels.
        x_col (str): The column name for x-coordinates.
        y_col (str): The column name for y-coordinates.

    Returns:
        pd.DataFrame: A DataFrame with the original data and additional columns for the nearest group and distance.
    """
    # Ensure the input DataFrame has the required columns
    if not all(col in df.columns for col in [x_col, y_col, group_key]):
        raise ValueError(f"The DataFrame must contain the columns: {x_col}, {y_col}, and {group_key}")
    # Separate rows with valid group labels
    group_df = df.dropna(subset=[group_key])
    coords_group = group_df[[x_col, y_col]].to_numpy()
    group_labels = group_df[group_key].to_numpy()
    # Build a KDTree for the group coordinates
    tree = cKDTree(coords_group)
    # Query the tree for each row in the DataFrame
    coords_query = df[[x_col, y_col]].to_numpy()
    distances, nearest_idx = tree.query(coords_query)
    # Get the nearest group for each row
    nearest_groups = group_labels[nearest_idx]
    # Add results to the DataFrame
    df["nearest_group"] = nearest_groups
    df["nearest_distance"] = distances
    return df

df_results = calculate_distance(df, group_key="HEV_group", x_col="x", y_col="y")

# Add distances to adata
adata.obs['dist_to_HEV_um'] = df_results['nearest_distance'].values


# --- Plot 2: Plot the distances

ncols = 3
fig, axs = plt.subplots(nrows=1, ncols=ncols, figsize=(9*ncols, 7))

# -- First plot just the H&E
n = 0
ax = axs[n]

# Create an emtpy adata to speed things up
adata_sub = adata[[]].copy()
sc.pl.spatial(adata_sub, color=None, size=0.5, title=f'{sample_id} ({subtype})',
            crop_coord=crop_coord, ax=ax, na_in_legend=False, show=False)

ax.set_title(ax.get_title(), fontweight='bold')

# -- Then plot distance to closest HEV

n = 1
ax = axs[n]

sc.pl.spatial(adata, color='dist_to_HEV_um', size=0.5, title=f'Distance to HEV', cmap='viridis',
            crop_coord=crop_coord, ax=ax,
            alpha_img=0.5, vmax=400)

ax.set_title(ax.get_title(), fontweight='bold')

# -- Then plot only distances for cells inside TLS

n = 2
ax = axs[n]

# Subset adata to only TLS spots
adata_sub = adata[adata.obs['in_TLS'] == True]
sc.pl.spatial(adata_sub, color='dist_to_HEV_um', size=0.5, title=f'Distance to HEV (in TLS)', cmap='viridis',
            crop_coord=crop_coord, ax=ax,
            alpha_img=0.5, vmax=400)

# Plot the HEVs on top
adata_sub = adata[adata.obs['niche_id'].notna()]
sc.pl.spatial(adata_sub, color='niche_id', size=0.5, title=f'Distance to HEV (in TLS)', palette=['cornflowerblue'],
            crop_coord=crop_coord, ax=ax, na_in_legend=False, legend_loc=None, show=False,
            alpha_img=0)

ax.set_title(ax.get_title(), fontweight='bold')

# Add a title
plt.suptitle(f'{sample_id} ({subtype}) - Distance to closest HEV niche (µm)', y=1.05, fontweight='bold')


# --- Now also calculate distances, splitting HEVs close to TLS vs not

# Add columns recording whether HEV is close to TLS or not
df['HEV_in_TLS'] = df['HEV_group'].map(near_tls_dict)

# Transfer the niche id if HEV is close to TLS
df['HEV_group_TLS'] = df.apply(lambda row: row['HEV_group'] if row['HEV_in_TLS'] else np.nan, axis=1)

# Add column for HEVs not close to TLS
df['HEV_group_not_TLS'] = df.apply(lambda row: row['HEV_group'] if not row['HEV_in_TLS'] else np.nan, axis=1)

dist_close_df = calculate_distance(df[['x', 'y', 'HEV_group_TLS']].copy(), group_key='HEV_group_TLS', x_col='x', y_col='y')
dist_far_df = calculate_distance(df[['x', 'y', 'HEV_group_not_TLS']].copy(), group_key='HEV_group_not_TLS', x_col='x', y_col='y')

# Add distances to adata
adata.obs['dist_to_HEV_close_to_TLS_um'] = dist_close_df['nearest_distance'].values
adata.obs['dist_to_HEV_not_close_to_TLS_um'] = dist_far_df['nearest_distance'].values


# --- Plot 3: Plot to check

ncols = 3
fig, axs = plt.subplots(nrows=1, ncols=ncols, figsize=(9*ncols, 7))

# -- First plot H&E + TLS + HEVs
n = 0
ax = axs[n]

# Subset adata to only TLS spots
adata_sub = adata[adata.obs['in_TLS'] == True]
adata_sub.obs['in_TLS'] = adata_sub.obs['in_TLS'].astype('category')
sc.pl.spatial(adata_sub, color='in_TLS', size=0.5, title=f'HEV niches in TLS', palette=['cornflowerblue'],
            crop_coord=crop_coord, ax=ax, na_in_legend=False, legend_loc=None, show=False,
            alpha_img=0.7)

adata_sub = adata[adata.obs['niche_id'].notna()]
sc.pl.spatial(adata_sub, color='niche_id', size=0.5, title=f'HEV + TLS niches', palette=['blue'],
            crop_coord=crop_coord, ax=ax, na_in_legend=False, legend_loc=None, show=False,
            alpha_img=0)       

ax.set_title(ax.get_title(), fontweight='bold')


# -- Then distances HEV close TLS

n = 1
ax = axs[n]

sc.pl.spatial(adata, color='dist_to_HEV_close_to_TLS_um', size=0.5, title=f'Distance to HEV\n(close to TLS)', cmap='viridis',
            crop_coord=crop_coord, ax=ax,
            alpha_img=0.5, vmax=400)

ax.set_title(ax.get_title(), fontweight='bold')


# -- Then distances HEV not close TLS

n = 2
ax = axs[n]

sc.pl.spatial(adata, color='dist_to_HEV_not_close_to_TLS_um', size=0.5, title=f'Distance to HEV\n(far from TLS)', cmap='viridis',
            crop_coord=crop_coord, ax=ax,
            alpha_img=0.5, vmax=400)

ax.set_title(ax.get_title(), fontweight='bold')

# Add a title
plt.suptitle(f'{sample_id} ({subtype}) - Distance to closest HEV niche (µm) - Separated by HEVs close to TLS vs not', y=1.05, fontweight='bold')


# --- Plot 4: Plot distribution of distances

sc.set_figure_params(fontsize=10)

# Enable gridlines behind the plot for all plots
plt.rcParams['axes.axisbelow'] = True

ncols = 3
nrows = 2

fig, axs = plt.subplots(nrows, ncols, figsize=(3*ncols, 3*nrows))

n=0
ax = axs[0][n]

sns.histplot(df_results, x=f'nearest_distance', bins=40, kde=False, ax=ax, alpha=0.8, linewidth=0.5, edgecolor='black')
ax.set_title('HEV distance', fontweight='bold')

# Now first 400 um
ax = axs[1][n]

sns.histplot(df_results[df_results['nearest_distance'] <= 400], x=f'nearest_distance', bins=40, kde=False, ax=ax, alpha=0.8, linewidth=0.5, edgecolor='black')
ax.set_title('HEV distance\n(first 400 um)', fontweight='bold')

# Now to HEVs close to TLS
n=1
ax = axs[0][n]

sns.histplot(dist_close_df, x=f'nearest_distance', bins=40, kde=False, ax=ax, alpha=0.8, linewidth=0.5, edgecolor='black')
ax.set_title('HEV distance\n(close to TLS)', fontweight='bold')

ax = axs[1][n]

sns.histplot(dist_close_df[dist_close_df['nearest_distance'] <= 400], x=f'nearest_distance', bins=40, kde=False, ax=ax, alpha=0.8, linewidth=0.5, edgecolor='black')
ax.set_title('HEV distance\n(close to TLS, first 400 um)', fontweight='bold')

# Now to HEVs not close to TLS
n=2
ax = axs[0][n]

sns.histplot(dist_far_df, x=f'nearest_distance', bins=40, kde=False, ax=ax, alpha=0.8, linewidth=0.5, edgecolor='black')
ax.set_title('HEV distance\n(far from TLS)', fontweight='bold')

ax = axs[1][n]
sns.histplot(dist_far_df[dist_far_df['nearest_distance'] <= 400], x=f'nearest_distance', bins=40, kde=False, ax=ax, alpha=0.8, linewidth=0.5, edgecolor='black')
ax.set_title('HEV distance\n(far from TLS, first 400 um)', fontweight='bold')

# Adjust yticklabels to show numbers in thousands
for ax in axs.flatten():
    yticks = ax.get_yticks()
    ax.set_yticklabels([f'{int(y/1000)}k' if y >= 1000 else f'{int(y)}' for y in yticks])

# Increase hspace
plt.subplots_adjust(hspace=0.5, wspace=0.4)

# Add a title
plt.suptitle(f'{sample_id} ({subtype}) - Distribution of distances to closest HEV niche (µm)', y=1.02, fontweight='bold')


# --- Save output

# Save plots
filename = f'{sample_id}_HEV_distances.pdf'

save_multi_image(os.path.join(out_dir_plots, filename))
plt.close('all')

# And save the distances
filename = f'{sample_id}_closest_HEV_distances.csv'
df_results[['x', 'y', 'HEV_group', 'nearest_group', 'nearest_distance']].to_csv(os.path.join(out_dir_csv, filename), index=True)

# Also save close vs far from TLS
filename = f'{sample_id}_HEV_distances_close_TLS.csv'
dist_close_df.to_csv(os.path.join(out_dir_csv, filename), index=True)

filename = f'{sample_id}_HEV_distances_far_from_TLS.csv'
dist_far_df.to_csv(os.path.join(out_dir_csv, filename), index=True)


# --- Finish

print("Script successfully completed")