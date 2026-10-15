#!/usr/bin/env python3

# --- Aim
# Identify differences in ligand-receptor interactions (LRI) among HEV groups using stlearn.
# Compare LRI patterns for HEVs near TLS versus those distant from TLS.

# --- Define environment
# env: vis-hd

# --- Load packages
import argparse
import os
import sys

import anndata as ad
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import stlearn as st
import re
import scanpy.plotting.palettes as sc_palettes
from collections.abc import Iterable
import math

sc.set_figure_params(facecolor="white", figsize=(8, 8))
sc.settings.verbosity = 3

mpl.rcParams['pdf.fonttype'] = 42


parser = argparse.ArgumentParser(description="Import Visium data & do basic QC",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

# Required
parser.add_argument("--res_dir", help="Path to results directory")
parser.add_argument("--h5ad_dir", help="Path to directory for .h5ad output files")
parser.add_argument("--project_dir", help="Project directory (root)")
parser.add_argument("--sample_id", help="Sample ID")
parser.add_argument("--sample_sheet", help="Path to sample sheet CSV file", default=None)
parser.add_argument("--bin_size", type=int, default=0,
                    help="Include bins within this distance (µm) of an HEV (0 = HEV bins only, 25 = HEV + perivascular bins)")

args = vars(parser.parse_args())


# --- Arguments

res_dir = args["res_dir"]
h5ad_dir=args["h5ad_dir"]
sample_id=args["sample_id"]
project_dir=args["project_dir"]
sys.path.append(os.path.join(project_dir, 'scripts'))
from helper_functions import auto_crop_scanpy, save_multi_image


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
resolution = 'square_016um'

adata = ad.read_h5ad(os.path.join(h5ad_dir, 'visium_hd', resolution, f'{sample_id}.h5ad'))


# --- Calculate crop coords for optimal presentation

coords = auto_crop_scanpy(adata, border_size=0.1)


# --- Load miniatlas classifications

# To make it easier to visualize trends by lineage
class_dir = os.path.join(project_dir, 'resources', 'cell_types', 'miniatlas')
class_file = 'miniatlas_subset_classifications.csv'

class_df = pd.read_csv(os.path.join(class_dir, class_file))


# --- Load HEV groups

threshold = 0.5
threshold_adj = str(threshold).split('.')[-1]

# Use the expanded groups
TLS_group_dir = os.path.join(res_dir, '04_vascular_microenvironment', '01_identify_HEV_niches', '02b_expand_HEV_groups', resolution, f'threshold_{threshold_adj}',)
TLS_group_dir = os.path.join(TLS_group_dir, 'csv_files')

most_recent_dir_path = TLS_group_dir

# Read in the results & add to anndata
score_of_i = 'HEVs'
TLS_file_name = f'{sample_id}_{score_of_i}_groups.csv'
TLS_results = pd.read_csv(os.path.join(most_recent_dir_path, TLS_file_name), index_col=0)

group_key = 'group_too_small_removed'
group_key_source = 'expanded_groups'

adata.obs[group_key] = TLS_results[group_key_source]


# --- Also load the HEV distance calculations

decon_resolution = resolution
dist_dir = os.path.join(res_dir, '04_vascular_microenvironment', '01_identify_HEV_niches', '05_calculate_all_HEV_distances', decon_resolution, 'csv_files')

most_recent = dist_dir
dist_file = f'{sample_id}_HEV_distances_all_spots.csv'

distance_df = pd.read_csv(os.path.join(most_recent, dist_file))

# For each HEV group, find all spots within bin_size µm (--bin_size; 0 = HEV bins only)
bin_size = args["bin_size"]
distance_col = 'nearest_016um_dist_um'
# Filter the distance DataFrame to only include spots within the specified bin size
distance_df = distance_df[distance_df[distance_col] <= bin_size]

# Make a dict for each barcode with all barcodes within the bin size
def create_spot_dict(df, group_key, distance_col):
    """Create a dictionary mapping each barcode to a list of barcodes within the specified distance."""
    spot_dict = {}
    for _, row in df.iterrows():
        barcode = row['barcode_016um']
        group = row[group_key]
        distance = row[distance_col]
        if group not in spot_dict:
            spot_dict[group] = []
        # Add the barcode to the list for this group
        spot_dict[group].append(barcode)
    return spot_dict

spot_dict = create_spot_dict(distance_df, 'HEV_group', distance_col)


# --- Load the HEV classifications

class_key = 'HEV_classification'

base_dir = os.path.join(res_dir, '04_vascular_microenvironment', '01_identify_HEV_niches', '03b_calculate_HEV_TLS_distance', resolution, f'threshold_{threshold_adj}',)
TLS_group_dir = os.path.join(base_dir, 'csv_files')

most_recent_dir_path = TLS_group_dir
filename = f'{sample_id}_HEV_closest_TLS_distance_classifications.csv'

HEV_TLS_class_df = pd.read_csv(os.path.join(most_recent_dir_path, filename), index_col=0)

# Add to adata
mapping = HEV_TLS_class_df[class_key].to_dict()
adata.obs[class_key] = adata.obs[group_key].map(mapping)


# --- Define some colordicts

# Define some colors
color_dict_HEV = {
    'non_TLS': 'cornflowerblue',
    'peri_TLS': 'darkorange',
    'intra_TLS': 'indigo',
}


# --- Plot 1: Spatial plots of HEVs & HEV classifications for sanity

# Plotting
sc.set_figure_params(fontsize=15)

ncols = 2
nrows = 1

fig, axs = plt.subplots(ncols=ncols, nrows=nrows, figsize=(7*ncols, 5*nrows))

# -- First the HEVs
n=0
ax = axs[n]

adata_sub = adata[adata.obs[group_key].notna()]

sc.pl.spatial(adata_sub, color=group_key, bw=False, alpha_img=0.5, crop_coord=coords, 
            title=f'{sample_id} ({subtype})\nHEVs',
            palette=['indigo'], ax=ax, legend_loc=None)

ax.set_title(ax.get_title(), fontweight='bold')

# -- Then plot TLS distance classification
n+=1
ax = axs[n]

adata_sub = adata[adata.obs[class_key].notna()]

sc.pl.spatial(adata_sub, color=class_key, bw=False, alpha_img=0.5, crop_coord=coords, 
            title=f'{sample_id} ({subtype})\nTLS distance classification',
            palette=color_dict_HEV, ax=ax,)

ax.set_title(ax.get_title(), fontweight='bold')

plt.suptitle(f'{sample_id} ({subtype}) - HEV classifications', y=1.2, x=0.5, fontweight='bold')


# --- For simplicity, adjust the grouping a bit

group_list = ['near_TLS_HEV', 'non_TLS_HEV']

class_of_i = 'grouping_adj'
adata.obs[class_of_i] = adata.obs[class_key]
adata.obs[class_of_i] = adata.obs[class_of_i].replace('intra_TLS', 'near_TLS_HEV')
adata.obs[class_of_i] = adata.obs[class_of_i].replace('peri_TLS', 'near_TLS_HEV')
adata.obs[class_of_i] = adata.obs[class_of_i].replace('non_TLS', 'non_TLS_HEV')

print(f'All classes: {adata.obs[class_of_i].dropna().unique()}')


# --- Plot 2: Check spots in bins around HEVs per group

def flatten_unique_values(d, keys):
    result = set()
    for key in keys:
        val = d.get(key)
        if isinstance(val, Iterable) and not isinstance(val, (str, bytes)):
            result.update(val)
        else:
            result.add(val)
    return list(result)

# Quick function to avoid repetition
def make_plot(ax, group_key, class_key, group_name, vmax=None):
    # Get the group_name barcodes
    all_groups = list(adata[adata.obs[class_key] == group_name].obs[group_key].dropna().unique())
    # spots in group
    all_barcodes = flatten_unique_values(spot_dict, all_groups)
    # Note which spots are in the group
    adata.obs['near_group'] = np.where(adata.obs_names.isin(all_barcodes), 'nearby_spots', None)
    # filter df and plot ticker on the new subplot axis
    adata_sub = adata[adata.obs['near_group'].notna()]
    sc.pl.spatial(adata_sub, color='near_group', alpha_img=0.7,
              crop_coord=coords, title=f'{sample_id} ({subtype})\nSpots {group_name}', ax=ax,
              palette=['lightskyblue'], vmax=vmax)
    # Subset adata to speed things up
    adata_sub = adata[adata.obs[group_key].isin(all_groups)]
    sc.pl.spatial(adata_sub, color=group_key, alpha_img=0,
            crop_coord=coords, title=f'{sample_id} ({subtype})\nSpots {group_name} (<{bin_size}um)', ax=ax,
            palette=['indigo'], legend_loc=None)
    ax.set_title(ax.get_title(), fontweight='bold')

# Groups to plot
ncols = len(group_list)
nrows = 1

fig, axs = plt.subplots(ncols=ncols, nrows=nrows, figsize=(8*ncols, 6*nrows),)

for n, group in enumerate(group_list):
    ax = axs[n]
    make_plot(ax, group_key=group_key, group_name=group, class_key=class_of_i, vmax=None)

plt.suptitle(f'{sample_id} ({subtype}) - Spots (resolution{decon_resolution}) within {bin_size} µm of HEVs (threshold={threshold})', y=1.1, fontweight='bold', x=0.5)


# --- Prep anndata for LR interactions analysis

categories = group_list + ['nearby_spots']

# Need to make a column with discrete labels
# Prep work
all_groups = list(adata.obs[group_key].dropna().unique())
# spots in group
all_barcodes = flatten_unique_values(spot_dict, all_groups)
# Note which spots are in the group
adata.obs['near_group'] = np.where(adata.obs_names.isin(all_barcodes), 'nearby_spots', None)

adata.obs['cluster_key'] = adata.obs['near_group']
adata.obs['cluster_key'] = np.where(adata.obs[class_of_i].notna(),
                                    adata.obs[class_of_i],
                                    adata.obs['cluster_key'])

adata.obs['cluster_key'] = adata.obs['cluster_key'].astype('category')
adata_sub = adata[adata.obs['cluster_key'].notna()]


# --- Calculate radius for stlearn to identify neighboring spots

# Get a scalefactor for converting pixel coords to um
spot_diameter_pixel = adata.uns['spatial'][sample_id]['scalefactors']['spot_diameter_fullres']

resolution_numbers = re.findall(r'\d+', resolution)[0]
resolution_numbers = int(resolution_numbers)

spot_diameter_um = resolution_numbers # Square size is whatever the resolution is 
scalefactor_um = spot_diameter_pixel/spot_diameter_um

# Vis hd spots have no distance between them, so distance from center to center is the square size
spot_dist_um = resolution_numbers
spot_dist_pixel = spot_dist_um * scalefactor_um

# But Vis hd spots are squares, I also want to catch the diagonal neighbors
# So need to adjust the radius for this (will not catch these neighbors if I use the length of the square's edge)

# Lucky I remember how to do pythagoras
c = pow(spot_dist_pixel/2, 2) + pow(spot_dist_pixel/2, 2)
c = math.sqrt(c)

radius_max = c*2

# Round up to whole pixel
radius_max = math.ceil(radius_max)

# Enlarge the neighbour radius slightly so all spots are captured; the exact
# radius occasionally misses spots near region boundaries.
radius_max = radius_max * 1.1

# stlearn expects the radius in hires-image coordinates, so rescale it
scalef = adata.uns['spatial'][sample_id]['scalefactors']['tissue_hires_scalef']
radius_max = radius_max*scalef


# --- Prepare an object for stlearn

# Make an object with just the normalized data (not log1p)
lr_obj = ad.AnnData(adata_sub.layers['norm'])
lr_obj.var = adata_sub.var
lr_obj.obs = adata_sub.obs
lr_obj.obsm = adata_sub.obsm
lr_obj.uns = adata_sub.uns

# Create a new object for stlearn so it doesn't get confused by the visium hd format
# Or maybe it's confused by scanpy but in any case it's confused
count_matrix = lr_obj.to_df()
spatial = lr_obj.obsm['spatial']
spatial = pd.DataFrame(spatial, columns=['imagecol', 'imagerow'])

data = st.create_stlearn(count=count_matrix, spatial=spatial, library_id=sample_id, scale=scalef,)

# Copy over uns
data.uns = adata_sub.uns


# --- Run LR analysis

# Loading the LR databases available within stlearn (from NATMI)
lrs = st.tl.cci.load_lrs(['connectomeDB2020_lit'], species='human')
print(len(lrs))

# Running the analysis #
st.tl.cci.run(data, lrs,
                  min_spots = 20, #Filter out any LR pairs with no scores for less than min_spots
                  distance=radius_max, # None defaults to spot+immediate neighbours; distance=0 for within-spot mode
                  n_pairs=10000, # Number of random pairs to generate; low as example, recommend ~10,000
                  n_cpus=None, # Number of CPUs for parallel. If None, detects & use all available.
                  )

lr_info = data.uns['lr_summary'] # A dataframe detailing the LR pairs ranked by number of significant spots.
print('\n', lr_info)

# Adjust the pvalues
st.tl.cci.adj_pvals(data, correct_axis='spot',
                   pval_adj_cutoff=0.05, adj_method='fdr_bh')

# Showing the rankings of the LR from a global and local perspective.
# Ranking based on number of significant hotspots.
# Make font smaller
sc.set_figure_params(fontsize=10)

st.pl.lr_summary(data, n_top=500)

plt.suptitle(f'{sample_id} ({subtype}) - LR HEVs - spots <{bin_size}µm HEV - top 500 LR pairs based on number of significant spots', y=1, fontweight='bold')


st.pl.lr_summary(data, n_top=50, figsize=(10,3))

plt.suptitle(f'{sample_id} ({subtype}) - LR HEVs - spots <{bin_size}µm HEV - top 50 LR pairs based on number of significant spots', y=1.3, fontweight='bold')


# --- Plot 3: Check stlearn's spatial neighbors calculations

# Disable the grid
plt.rcParams['axes.grid'] = False

# Step 1: Get spot coordinates
coords = data.obsm['spatial']

# Step 2: Convert neighbors to proper format
# Assuming it's a DataFrame:
neighbors_df = data.obsm['spot_neighbours']
neighbor_strs = neighbors_df['neighbour_indices']

# Convert string -> list of ints
neighbors = [list(map(int, s.split(','))) for s in neighbor_strs]

# Step 3: Build mapping from index to coordinate
# AnnData index must match neighbors_df index
spot_index_map = {spot: i for i, spot in enumerate(data.obs_names)}

# Step 4: Plot
fig, ax = plt.subplots(figsize=(4, 4))
ax.scatter(coords[:, 0], coords[:, 1], s=10, c='grey', alpha=0.5)

for spot_id, nbrs in zip(neighbors_df.index, neighbors):
    i = spot_index_map[spot_id]  # central spot index
    for j in nbrs:
        ax.plot([coords[i, 0], coords[j, 0]],
                [coords[i, 1], coords[j, 1]],
                color='blue', alpha=0.3, linewidth=0.5)

ax.set_aspect('equal')
ax.invert_yaxis()
ax.set_title("Spatial Neighbors (stlearn)")

plt.suptitle(f'{sample_id} ({subtype}) - LR HEVs - spots <{bin_size}µm HEV - stlearn spatial neighbors', y=1.1, fontweight='bold')


# --- Diagnostic plots

# Should show next to no correlation between the hotspots of the LR pairs and the expression level and frequency of expression
# If not: need larger number of permutations
st.pl.lr_diagnostics(data, figsize=(10,2.5))
plt.suptitle(f'{sample_id} ({subtype}) - LR HEVs - spots <{bin_size}µm HEV - Correlation LR pairs & expression level / frequency', y=1.1, fontweight='bold')

st.pl.lr_n_spots(data, n_top=50, figsize=(11, 3),
                    max_text=100)
plt.suptitle(f'{sample_id} ({subtype}) - LR HEVs - spots <{bin_size}µm HEV - top 50 LRs: number of spots significant vs non-significant', y=1.1, fontweight='bold')

st.pl.lr_n_spots(data, n_top=500, figsize=(11, 3),
                    max_text=100)
plt.suptitle(f'{sample_id} ({subtype}) - LR HEVs - spots <{bin_size}µm HEV - top 500 LRs: number of spots significant vs non-significant', y=1.1, fontweight='bold')


# --- LR statistics visualizations

best_lr = data.uns['lr_summary'].index.values[0] # Just choosing one of the top from lr_summary

# Need smaller dots for HD
spot_size = 0.5


# -- Plot: LR statistics

stats = ['lr_scores', 'p_vals', 'p_adjs', '-log10(p_adjs)']
fig, axes = plt.subplots(ncols=len(stats), figsize=(16,6))
for i, stat in enumerate(stats):
    st.pl.lr_result_plot(data, use_result=stat, use_lr=best_lr, show_color_bar=True, ax=axes[i],
                         size=spot_size)
    axes[i].set_title(f'{best_lr} {stat}')

plt.suptitle(f'{sample_id} ({subtype}) - LR HEVs - spots <{bin_size}µm HEV - Top LR ({best_lr}) statistics', y=1, fontweight='bold')


# -- Plot: -log10(padj)

fig, axes = plt.subplots(ncols=2, figsize=(8,6))
st.pl.lr_result_plot(data, use_result='-log10(p_adjs)', use_lr=best_lr, show_color_bar=True, ax=axes[0], size=spot_size)
st.pl.lr_result_plot(data, use_result='lr_sig_scores', use_lr=best_lr, show_color_bar=True, ax=axes[1], size=spot_size)
axes[0].set_title(f'{best_lr} -log10(p_adjs)')
axes[1].set_title(f'{best_lr} lr_sig_scores')

plt.suptitle(f'{sample_id} ({subtype}) - LR HEVs - spots <{bin_size}µm HEV - Top LR ({best_lr}) -log10(p_adj) + significant scores', y=1, fontweight='bold')


# --- These plots are more informative for LRs of interest
# # --- LR interpretation visualizations

# # All spots
# st.pl.lr_plot(data, best_lr, inner_size_prop=0.1, outer_mode='binary', pt_scale=spot_size,
#               use_label=None, show_image=True,
#               sig_spots=False,)

# # Significant spots only
# st.pl.lr_plot(data, best_lr, outer_size_prop=1, outer_mode='binary', pt_scale=spot_size,
#               use_label=None, show_image=True,
#               sig_spots=True)


# # --- Continuous LR expression for significant spots

# # All spots #
# st.pl.lr_plot(data, best_lr,
#               inner_size_prop=0.04, middle_size_prop=.07, outer_size_prop=.4,
#               outer_mode='continuous', pt_scale=60,
#               use_label=None, show_image=True,
#               sig_spots=False)

# # Only significant spots #
# st.pl.lr_plot(data, best_lr,
#               inner_size_prop=0.04, middle_size_prop=.07, outer_size_prop=.4,
#               outer_mode='continuous', pt_scale=60,
#               use_label=None, show_image=True,
#               sig_spots=True)


# # --- Binary LR coexpression with spot annotations

# data.obs['cell_type'] = adata.obs['cluster_key']

# st.pl.lr_plot(data, best_lr,
#               inner_size_prop=0.08, middle_size_prop=.3, outer_size_prop=.5,
#               outer_mode='binary', pt_scale=spot_size*2,
#               use_label='cell_type', show_image=True,
#               sig_spots=True)


# --- Predicting significant CCIs

data.obs['cell_type'] = adata_sub.obs['cluster_key']

# Running the counting of co-occurence of cell types and LR expression hotspots #
st.tl.cci.run_cci(data, 'cell_type', # Spot cell information either in data.obs or data.uns
                  min_spots=3, # Minimum number of spots for LR to be tested.
                  spot_mixtures=False, # If True will use the label transfer scores,
                                      # so spots can have multiple cell types if score>cell_prop_cutoff
                  cell_prop_cutoff=0.2, # Spot considered to have cell type if score>0.2
                  sig_spots=True, # Only consider neighbourhoods of spots which had significant LR scores.
                  n_perms=1000 # Permutations of cell information to get background, recommend ~1000
                 )


# --- Plot: Check interaction and cell type frequency correlation
# should show little to no correlation if the number of permutations is adequate; otherwise recommend increasing n_perms above

st.pl.cci_check(data, 'cell_type', figsize=(5,3), axis_text_size=10, cell_label_size=10, tick_size=10)

plt.suptitle(f'{sample_id} ({subtype}) - LR HEVs - spots <{bin_size}µm HEV - interaction & frequency correlation', y=1.1, fontweight='bold')


# --- Plot: CCI network

# Visualising the no. of interactions between cell types across all LR pairs #
pos_1 = st.pl.ccinet_plot(data, 'cell_type', return_pos=True)
plt.suptitle(f'{sample_id} ({subtype}) - LR HEVs - spots <{bin_size}µm HEV - CCI network', y=1.1, fontweight='bold')

# This one is also better if you have a pair of interest
# # Just examining the cell type interactions between selected pairs #
# lrs = data.uns['lr_summary'].index.values[0:3]
# for n, best_lr in enumerate(lrs[0:3]):
#     st.pl.ccinet_plot(data, 'cell_type', best_lr, min_counts=2,
#                          figsize=(10,7.5), pos=pos_1,
#                       )


# --- Plot: CCI chord-plot

# # Need to downgrade to matplotlib 3.7 for this to work
# st.pl.lr_chord_plot(data, 'cell_type')

# for lr in lrs:
#     st.pl.lr_chord_plot(data, 'cell_type', lr)


# --- Heatmap visualizations

# -- Plot 1: Top interactions in single heatmap
n_int = 30
n_cat = pow(len(categories), 2)

# This will automatically select the top interacting CCIs and their respective LRs #
st.pl.lr_cci_map(data, 'cell_type', lrs=None, min_total=100, figsize=(0.5*n_cat,0.25*n_int),
                 n_top_lrs=n_int)

plt.suptitle(f'{sample_id} ({subtype}) - LR HEVs - spots <{bin_size}µm HEV - top {n_int} interactions heatmap', y=1, fontweight='bold')


# -- Plot 2: Top interactions separate heatmaps

nrows = 4
ncols = 4

fig, axs = plt.subplots(nrows, ncols, figsize=(4*ncols, 3*nrows))
axs = axs.ravel()

lrs = data.uns['lr_summary'].index.values[0:nrows*ncols]
for n, lr in enumerate(lrs):
    ax = axs[n]
    st.pl.cci_map(data, 'cell_type', lr, ax=ax)
    ax.set_title(lr, fontweight='bold')
    # Adjust plot lims
    ymax = ax.get_ylim()[1]
    ax.set_ylim(ymax-len(categories)+0.5, ymax+0.5)
    xmax = ax.get_xlim()[1]
    ax.set_xlim(0, xmax+2)

plt.suptitle(f'{sample_id} ({subtype}) - LR HEVs - spots <{bin_size}µm HEV - top {nrows*ncols} interactions', y=1, fontweight='bold')

plt.tight_layout()


# --- Save the output

qc_dir = os.path.join(res_dir, '04_vascular_microenvironment', '04_ligand_receptor_analysis', '01_run_LR_analysis_stlearn',)

out_dir = os.path.join(qc_dir, f'{bin_size}um_bin')
os.makedirs(out_dir, exist_ok=True)

# Save plots
out_dir_plots = os.path.join(out_dir, 'plots')
os.makedirs(out_dir_plots, exist_ok=True)

# Also save the groupings as a .csv
out_dir_csv = os.path.join(out_dir, 'csv_files')
os.makedirs(out_dir_csv, exist_ok=True)

# Save plots
filename = f'{sample_id}_HEV_LR_analysis.pdf'
save_multi_image(os.path.join(out_dir_plots, filename))

# And save the dfs
filename = f'{sample_id}_LR_summary.csv'
lr_info = data.uns['lr_summary']
lr_info.to_csv(os.path.join(out_dir_csv, filename), index=True)

# -- And save the anndata object
# Save lrfeatures separately as a csv because it throws an error when saving as anndata
to_save = data.uns.get("lrfeatures", None)
filename = f'{sample_id}_HEV_LR_analysis_uns_lrfeatures.csv'
to_save.to_csv(os.path.join(out_dir_csv, filename))

# Remove the lrfeatures from the uns
data.uns.pop("lrfeatures", None)

filename = f'{sample_id}_HEV_LR_analysis_stlearn.h5ad'
ad.AnnData.write(data, filename=os.path.join(out_dir_csv, filename), compression='gzip')


# --- Finish

print("Script successfully completed")
