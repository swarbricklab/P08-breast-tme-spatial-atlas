#!/usr/bin/env python3

# --- Aim
# Calculate distances from each vessel to surrounding spots.
# Enable analysis of microenvironment composition as a function of spatial proximity to vessels.

# --- Define environment
# env: vis-hd

# --- Load packages
import argparse
import os
from datetime import date, datetime

import anndata as ad
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
import squidpy as sq
import seaborn as sns
from matplotlib.colors import ListedColormap
from scipy.spatial import cKDTree
import re


# Import custom functions
from sophiesfunctions.auto_crop import auto_crop
from sophiesfunctions.save_multi_image import save_multi_image

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

parser.add_argument("--min_count", help="Min count for QC", type=int)
parser.add_argument("--min_gene", help="Min gene for QC", type=int)
parser.add_argument("--min_spots", help="Minimal spots per gene", type=int)

args = vars(parser.parse_args())


# --- Arguments

res_dir = args["res_dir"]
h5ad_dir=args["h5ad_dir"]
sample_id=args["sample_id"]
project_dir=args["project_dir"]



# --- Get some info from the sample sheet
# Use --sample_sheet if provided, otherwise default to config/sample_sheet.csv
if args["sample_sheet"] is not None:
    sample_sheet_path = args["sample_sheet"]
else:
    sample_sheet_path = os.path.join(project_dir, 'config', 'sample_sheet.csv')

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
target_res = 'square_008um'

resolutions = [base_res]
if target_res not in resolutions:
    resolutions.append(target_res)

adata_dict = {}

for resolution in resolutions:
    adata_dict[resolution] = ad.read_h5ad(os.path.join(h5ad_dir, 'visium_hd', resolution, f'{sample_id}.h5ad'))


# --- Load miniatlas classifications

decon_level = 'minor'

# To make it easier to visualize trends by lineage
class_dir = os.path.join(project_dir, 'resources', 'cell_types', 'miniatlas')
class_file = 'miniatlas_subset_classifications.csv'

class_df = pd.read_csv(os.path.join(class_dir, class_file))

# Get endothelial classes
endothelial = list(class_df[class_df['miniatlas_major'] == 'Endothelial'][f'miniatlas_{decon_level}'].unique())
endothelial = sorted(endothelial)


# --- Calculate crop coords for optimal presentation

# Also calculate some coords for squidpy
coords_dict_sq  = {}

for resolution in resolutions:
    coords = auto_crop(adata_dict[resolution], border_size=0.1)
    coords_dict_sq[resolution] = coords

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


# --- Prepare for saving

qc_dir = os.path.join(project_dir, 'results', '04_vascular_microenvironment', '02_identify_all_vasculature', '03_calculate_distances')

today = date.today()
today = today.strftime("%Y%m%d")

threshold_adj = str(threshold).split('.')[-1]
out_dir = os.path.join(qc_dir)

os.makedirs(out_dir, exist_ok=True)

# Save plots
out_dir_plots = os.path.join(out_dir, 'plots', today)
os.makedirs(out_dir_plots, exist_ok=True)

# Also save the groupings as a .csv
out_dir_csv = os.path.join(out_dir, 'csv_files', today)
os.makedirs(out_dir_csv, exist_ok=True)


# --- Figure out pixel to um scalefactor

resolution_numbers = re.findall(r'\d+', target_res)[0]
resolution_numbers = int(resolution_numbers)

spot_diameter_um = resolution_numbers # Square size is whatever the resolution is 

# Get a scalefactor for converting pixel coords to um
spot_diameter_pixel = adata_dict[target_res].uns['spatial'][sample_id]['scalefactors']['spot_diameter_fullres']

scalefactor_um = spot_diameter_um/spot_diameter_pixel


# --- Calculate distances per vessel group

# Function to calculate distances
def calculate_distance(base_res, target_res, group_key):
    # Extract shortened versions of resolution
    base_res_short = base_res.split('_')[1]
    target_res_short = target_res.split('_')[1]
    adata = adata_dict[base_res]
    coords = adata.obsm['spatial']  # Extract the spot coords
    groups = adata.obs[group_key]  # Group labels
    # Convert to DataFrame for easier manipulation
    df = pd.DataFrame(coords, columns=["x", "y"])
    df["vessel_group"] = groups.values
    # Add barcodes to the 16µm DataFrame
    df[f"barcode_{base_res_short}"] = adata.obs.index
    df = df.dropna(subset='vessel_group')
    # Extract coordinates and barcodes for each group
    grouped_data_16 = {
        group: df[df["vessel_group"] == group][["x", "y", f"barcode_{base_res_short}"]].reset_index(drop=True)
        for group in df["vessel_group"].unique()
    }
    # --- Extract decon_resolution spot coordinates
    adata_48um = adata_dict[target_res]
    df_48um = pd.DataFrame(adata_48um.obsm['spatial'], columns=["x", "y"])
    df_48um["barcode"] = adata_48um.obs.index
    # --- Compute the minimum distance to any 16 µm spot in each group
    results = []
    for group, group_df_16 in grouped_data_16.items():
        if group_df_16.empty:
            continue  # Skip empty groups
        coords_16um = group_df_16[["x", "y"]].to_numpy()
        barcodes_16um = group_df_16[f"barcode_{base_res_short}"].to_numpy()
        # Build KDTree for the current group
        tree = cKDTree(coords_16um)
        # Query the tree to find the closest 16 µm spot for each 48 µm spot
        min_distances, nearest_idx = tree.query(df_48um[['x', 'y']].values)
        # Get coordinates and barcodes of the closest 16 µm spot
        nearest_coords = coords_16um[nearest_idx]
        nearest_barcodes = barcodes_16um[nearest_idx]
        # Store results in a DataFrame, including barcodes
        df_temp = pd.DataFrame({
            f"barcode_{target_res_short}": df_48um["barcode"],  # Save barcode
            "vessel_group": group,
            f"x_{target_res_short}": df_48um["x"],
            f"y_{target_res_short}": df_48um["y"],
            f"nearest_{base_res_short}_x": nearest_coords[:, 0],
            f"nearest_{base_res_short}_y": nearest_coords[:, 1],
            f"nearest_{base_res_short}_barcode": nearest_barcodes,
            f"nearest_{base_res_short}_dist": min_distances
        })
        results.append(df_temp)
    # Combine results for all groups
    df_results = pd.concat(results, ignore_index=True)
    # --- Convert pixel distances to um
    df_results[f'nearest_{base_res_short}_dist_um'] = df_results[f'nearest_{base_res_short}_dist'] * scalefactor_um
    return df_results

# Save results in a dict
dist_dict = {}

for group_key in endothelial:
    print(f'Calculating distances for {group_key}')
    dist_dict[group_key] = calculate_distance(base_res, target_res, group_key)


# Also get a df with smallest distance to any vessel within a vessel type
dist_dict_min_dist = {}

base_res_short = base_res.split('_')[1]
target_res_short = target_res.split('_')[1]

for group_key in endothelial:
    print(f'Calculating minimal distances for {group_key}')
    df_results = dist_dict[group_key]
    df_min_dist = df_results.sort_values(f"nearest_{base_res_short}_dist_um").drop_duplicates(f"barcode_{target_res_short}", keep="first")
    dist_dict_min_dist[group_key] = df_min_dist


# Make some plots to check it worked
# --- Plot 1: All vessel groups together

ncols = len(endothelial)
nrows = 1

fig, axs = plt.subplots(ncols=ncols, nrows=nrows, figsize=(8*ncols, 6*nrows),)

for n, group_key in enumerate(endothelial):
    # Set ax
    ax = axs[n]
    # Add to adata for plotting
    df_min_dist = dist_dict_min_dist[group_key]
    adata_dict[target_res].obs['dist_to_HEV_um'] = df_min_dist.set_index(f'barcode_{target_res_short}')[f'nearest_{base_res_short}_dist_um']
    # Plot all HEV groups
    adata_sub = adata_dict[base_res][adata_dict[base_res].obs[group_key].notna()]
    sq.pl.spatial_scatter(adata_dict[target_res], color='dist_to_HEV_um', img_alpha=0.7,
                    crop_coord=coords_dict_sq[target_res], title=f'{sample_id} ({subtype})\nHEV classification', ax=ax,
                    shape='square', vmax=400)
    sq.pl.spatial_scatter(adata_sub, color=group_key, img_alpha=0,
                    crop_coord=coords_dict_sq[target_res], title=f'{sample_id} ({subtype})\n{group_key}\nDistance (µm)', ax=ax,
                    shape='square', palette=ListedColormap(['white']), legend_loc=None)
    ax.set_title(ax.get_title(), fontweight='bold')


plt.suptitle(f'{sample_id} ({subtype}) - vessel groups - distance (um) to each {target_res} spot (threshold={threshold})', y=1.15, fontweight='bold', x=0.5)


# --- Plot 2: Distribution of distance values

# Enable gridlines behind the plot for all plots
plt.rcParams['axes.axisbelow'] = True

nrows = 1
ncols = 4

for group_key in endothelial:

    df_min_dist = dist_dict_min_dist[group_key]

    fig, axs = plt.subplots(nrows, ncols, figsize=(4*ncols, 4*nrows))

    n=0
    ax = axs[n]

    sns.histplot(df_min_dist, x=f'nearest_{base_res_short}_dist_um', bins=40, kde=False, ax=ax, alpha=0.8, linewidth=0.5, edgecolor='black')
    ax.set_title(f'{group_key} distance', fontweight='bold')

    # Now first 400 um only
    n+=1
    ax = axs[n]

    df_min_dist_sub = df_min_dist[df_min_dist[f'nearest_{base_res_short}_dist_um'] < 400]
    sns.histplot(df_min_dist_sub, x=f'nearest_{base_res_short}_dist_um', bins=40, kde=False, ax=ax, alpha=0.8, linewidth=0.5, edgecolor='black')
    ax.set_title('<400 µm', fontweight='bold')

    # Now first 400 um only
    n+=1
    ax = axs[n]

    sns.histplot(df_min_dist_sub, x=f'nearest_{base_res_short}_dist_um', bins=16, kde=False, ax=ax, alpha=0.8, linewidth=0.5, edgecolor='black')
    ax.set_title('<400µm - 25 µm bins', fontweight='bold')

    plt.suptitle(f'{sample_id} ({subtype}) - {group_key} groups - distribution of distance measurements - {target_res} spot (threshold={threshold})', y=1.1, fontweight='bold', x=0.5)

    plt.tight_layout()

    # Now first 400 um only
    n+=1
    ax = axs[n]

    sns.histplot(df_min_dist_sub, x=f'nearest_{base_res_short}_dist_um', bins=25, kde=False, ax=ax, alpha=0.8, linewidth=0.5, edgecolor='black')
    ax.set_title('<400µm - 16 µm bins', fontweight='bold')

    plt.suptitle(f'{sample_id} ({subtype}) - {group_key} groups - distribution of distance measurements - {target_res} spot (threshold={threshold})', y=1.1, fontweight='bold', x=0.5)

    plt.tight_layout()


# --- Save output

# Save plots
filename = f'{sample_id}_vessel_distances.pdf'

save_multi_image(os.path.join(out_dir_plots, filename))
plt.close('all')

# And save the distances
# Actually don't save for the 008um resolution because the files are big and I doubt I will use them (don't have decon results)
if target_res != 'square_008um':
    os.makedirs(os.path.join(out_dir_csv, 'all_dist'), exist_ok=True)
    for group_key in endothelial:
        filename = f'{sample_id}_{group_key.replace(" ", "_")}_distances_all_spots.csv'
        df_results = dist_dict[group_key]
        # The files are fairly large so I'll subset these dfs to be up to 400 µm bc it's unlikely I'll use the higher distances
        if target_res != 'square_048um':
            df_results = df_results[df_results[f'nearest_{base_res_short}_dist_um'] < 400]
        df_results.to_csv(os.path.join(out_dir_csv, 'all_dist', filename), index=False)

# Also save the distance for each spot to closest HEV
os.makedirs(os.path.join(out_dir_csv, 'min_dist'), exist_ok=True)

for group_key in endothelial:
    filename = f'{sample_id}_distance_closest_{group_key.replace(" ", "_")}_per_spot.csv'
    df_min_dist = dist_dict_min_dist[group_key]
    df_min_dist.to_csv(os.path.join(out_dir_csv, 'min_dist', filename), index=False)


# --- Finish

print("Script successfully completed")
