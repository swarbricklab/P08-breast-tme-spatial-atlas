#!/usr/bin/env python3

# --- Aim
# Automatically identify high endothelial venule (HEV) regions by scoring for HEV-enriched endothelial transcripts in Xenium data.
# Segment and annotate HEV-rich areas for spatial and functional analysis.

# --- Define environment
# env: xen-hev

# --- Load packages

import argparse
import os
import sys
import json

import anndata as ad
import scanpy as sc
import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.patches as mpatches

from scipy.spatial import cKDTree
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components
from scipy.spatial import distance_matrix

sc.set_figure_params(figsize=(8, 8), fontsize=20)


# --- Argparse arguments

parser = argparse.ArgumentParser(description="Import Visium data & save as .h5ad file",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

# Required

parser.add_argument("--raw_data_dir", help="Path to raw data")
parser.add_argument("--h5ad_dir", help="Path to directory for .h5ad output files")
parser.add_argument("--sample_id", help="Sample ID")
parser.add_argument("--project_dir", help="Path to the root of the project directory")
parser.add_argument("--outs_path", help="Path to spaceranger 'outs'")
parser.add_argument("--use_data", help="Which data to use (log_norm or SCT)")
parser.add_argument("--res_dir", help="Path to results directory")
parser.add_argument("--seed", help="Seed for reproducibility", type=int)
parser.add_argument("--sample_sheet", help="Path to sample sheet CSV file", default=None)


args = vars(parser.parse_args())


# --- Arguments

h5ad_dir=args["h5ad_dir"]
sample_id=args["sample_id"]

res_dir = args['res_dir']
project_dir = args['project_dir']

# Import custom functions
sys.path.append(os.path.join(project_dir, 'scripts'))
from helper_functions import auto_crop_xenium


# --- Extract some info from the sample sheet
# Use sample_sheet argument if provided, else default
if args.get("sample_sheet"):
    sample_sheet_path = args["sample_sheet"]
else:
    sample_dir = os.path.join(project_dir, 'config')
    sample_sheet_file = 'sample_sheet_xenium.csv'
    sample_sheet_path = os.path.join(sample_dir, sample_sheet_file)

sample_sheet = pd.read_csv(sample_sheet_path)
tissue_id = sample_sheet[sample_sheet['section_id'] == sample_id]['tissue_id'].values[0]

# Also load tissue metadata
metadata_dir = os.path.join(project_dir, 'resources',)
tissue_metadata_file = 'tissue_metadata.csv'
tissue_metadata = pd.read_csv(os.path.join(metadata_dir, tissue_metadata_file))

subtype = tissue_metadata[tissue_metadata['tissue_id'] == str(tissue_id)]['clinical_subtype'].values[0]


# --- Load data

adata = ad.read_h5ad(os.path.join(h5ad_dir, 'xenium', f'{sample_id}.h5ad'))


# --- Calculate crop coordinates for optimal presentation

crop_coord = auto_crop_xenium(adata)


# --- Load the tangram annotations

# Decide which annotation level to use
level = 'major'

annotation_dir = os.path.join(project_dir, 'results', '02_cell_type_mapping', 'label_transfer', 'xenium', f'{level}_level')
annotation_file = f'{sample_id}_{level}_pred_cell_types.csv'
annotations = pd.read_csv(os.path.join(annotation_dir, annotation_file), index_col=0)

# Add to anndata object
annotation_key = 'tangram_annotation'
probability_key = 'tangram_probability'

adata.obs[annotation_key] = annotations['pred_cell_type']
adata.obs[probability_key] = annotations['pred_cell_type_prob']

# Also add the subset annotations
level = 'subset'

annotation_dir = os.path.join(project_dir, 'results', '02_cell_type_mapping', 'label_transfer', 'xenium', f'{level}_level')
annotation_file = f'{sample_id}_{level}_pred_cell_types.csv'
annotations = pd.read_csv(os.path.join(annotation_dir, annotation_file), index_col=0)

adata.obs[level] = annotations['pred_cell_type']


# --- Load the TLS annotations

tls_dir = os.path.join(res_dir, '05_xenium_mechanisms', '01_identify_TLS_regions', 'csv_files')

tls_file = f'{sample_id}_cell_TLS_annotations.csv'
tmp_df = pd.read_csv(os.path.join(tls_dir, tls_file), index_col=0)

# Add to adata
adata.obs['in_TLS'] = tmp_df['in_TLS']

# Make categorical
adata.obs['in_TLS'] = adata.obs['in_TLS'].astype('category')

color_dict_TLS = {
    True: 'cornflowerblue',
    False: 'lightgrey'
}


# --- Load plotting colors

level = 'major'
color_dir = os.path.join(project_dir, 'resources', 'plotting_colors')
color_file = f'miniatlas_{level}_colors.csv'

color_dict = pd.read_csv(os.path.join(color_dir, color_file), index_col=0).to_dict()['color']

all_cell_types = adata.obs[annotation_key].astype('category').cat.categories

# Make sure all cell types have a color
for cell_type in all_cell_types:
    if cell_type not in color_dict.keys():
        print(f'Adding color for {cell_type}')
        color_dict[cell_type] = 'lightgrey'

# Also add subset colors
miniatlas_dir = os.path.join(project_dir, 'resources', 'cell_types', 'miniatlas')
miniatlas_file = 'miniatlas_subset_classifications.csv'

miniatlas_df = pd.read_csv(os.path.join(miniatlas_dir, miniatlas_file),)

# Add the other levels to the color_dict
level = 'subset'

for cell_type in miniatlas_df[f'miniatlas_{level}'].unique():
    if cell_type not in color_dict.keys():
        # Check that the value is not nan
        if pd.isna(cell_type):
            continue
        print(f'Adding color for {cell_type}')
    # Get the corresponding major cell type
    major_type = miniatlas_df[miniatlas_df[f'miniatlas_{level}'] == cell_type][f'miniatlas_major'].values
    # Find which value is in the color_dict
    major_type_in_dict = [m for m in major_type if m in color_dict.keys()]
    if len(major_type_in_dict) > 0:
        color_dict[cell_type] = color_dict[major_type_in_dict[0]]
    else:
        color_dict[cell_type] = 'lightgray'


# --- Plot 1: Spatial overview of B+T cells & TLS regions

# Make sure it's plotting w/ png backend
sc.set_figure_params(vector_friendly=True)

cell_types_of_interest = ['B-cells', 'T-cells']

ncols = 2
nrows = 1

plt.figure(figsize=(9*ncols, 7*nrows))

# add a new subplot iteratively
n=0
ax = plt.subplot(nrows, ncols, n + 1)

sc.pl.spatial(adata, color=[annotation_key], groups=cell_types_of_interest, size=0.5, title=f'{sample_id} ({subtype})',
            crop_coord=crop_coord, ax=ax, palette=color_dict, na_in_legend=False, alpha_img=0.3)

ax.set_title(f'{sample_id} ({subtype})\nT & B cells', fontweight='bold')

# Remove the labels
ax.set_xlabel('')
ax.set_ylabel('')

# Add another subplot
ax = plt.subplot(nrows, ncols, n + 2)

sc.pl.spatial(adata, color='in_TLS', size=0.5, title=f'Cells inside TLS',
            crop_coord=crop_coord, ax=ax, na_in_legend=False, palette=color_dict_TLS, alpha_img=0.3)

ax.set_title(f'Cells inside TLS', fontweight='bold')

# Remove the labels
ax.set_xlabel('')
ax.set_ylabel('')

plt.suptitle(f'{sample_id} ({subtype}) - T & B cells based on Tangram label transfer + TLS regions (annotation_level={level})', y=1.1, fontweight='bold')


# --- Identify HEVs by scoring

# Define the HEV transcripts
HEV_transcripts = ['TSPAN7', 'MEOX2', 'RAMP3', 'SPNS2', 'IL33']

sc.tl.score_genes(adata, gene_list=HEV_transcripts, score_name='HEV_score')

# Define a threshold for HEV score
threshold = 0.6
adata.obs['HEV_positive'] = adata.obs['HEV_score'] > threshold


# --- Plot 2: Plot HEV score + HEV high

plt.figure(figsize=(9*ncols, 7*nrows))

# add a new subplot iteratively
n=0
ax = plt.subplot(nrows, ncols, n + 1)

sc.pl.spatial(adata, color='HEV_score', size=0.5, title=f'{sample_id} ({subtype})',
            crop_coord=crop_coord, ax=ax, na_in_legend=False, vmin=0, vmax=0.8, cmap='viridis', alpha_img=0.3)

ax.set_title(f'HEV score', fontweight='bold')

# Remove the labels
ax.set_xlabel('')
ax.set_ylabel('')

# Add another subplot
ax = plt.subplot(nrows, ncols, n + 2)

adata.obs['HEV_positive'] = adata.obs['HEV_positive'].astype('category')
sc.pl.spatial(adata, color='HEV_positive', size=0.5, title=f'{sample_id} ({subtype})',
            crop_coord=crop_coord, ax=ax, na_in_legend=False, palette=['lightgray'], legend_loc=None, alpha_img=0.3)

adata_sub = adata[adata.obs['HEV_positive'] == True]
sc.pl.spatial(adata_sub, color='HEV_positive', spot_size=25, title=f'{sample_id} ({subtype})',
            crop_coord=crop_coord, ax=ax, na_in_legend=False, palette=['blue'], alpha_img=0.3)

ax.set_title(f'HEV score >{threshold}', fontweight='bold')

# Remove the labels
ax.set_xlabel('')
ax.set_ylabel('')

plt.suptitle(f'{sample_id} ({subtype}) - HEV scoring based on HEV transcripts {HEV_transcripts}', y=1.1, fontweight='bold')


# --- Create spatial groups of HEVs

def find_spatial_niches(
    adata,
    score_key: str,
    xy_keys=("x", "y"),
    r_um: float = 30.0,
    high_thr: float = 0.5,
    low_thr: float | None = None,
    min_neighbors: int = 1,
    min_core_size: int = 3,
    min_expanded_size: int = 5, 
    do_expand: bool = True,
):
    # --- Extract coordinates
    if isinstance(xy_keys, str) and xy_keys in adata.obsm:
        xy = adata.obsm[xy_keys]
        if xy.shape[1] < 2:
            raise ValueError(f"obsm['{xy_keys}'] must have at least 2 columns (x,y).")
        coords = np.asarray(xy[:, :2], dtype=float)
    else:
        xkey, ykey = xy_keys
        if xkey not in adata.obs or ykey not in adata.obs:
            raise KeyError(f"Expected '{xkey}' and '{ykey}' in adata.obs (or pass an obsm key).")
        coords = adata.obs[[xkey, ykey]].to_numpy(dtype=float)
    n = coords.shape[0]
    scores = np.asarray(adata.obs[score_key].values, dtype=float)
    # --- 1) High-scoring mask
    mask_high = scores >= high_thr
    adata.obs["sig_high"] = mask_high
    adata.obs["niche_core_id"] = 0
    adata.obs["niche_expanded_id"] = 0
    if not mask_high.any():
        return {
            "core_stats": pd.DataFrame(columns=["niche_core_id", "size", "x_centroid", "y_centroid"]),
            "core_membership": pd.DataFrame(columns=["cell_id", "niche_id"]),
            "expanded_membership": pd.DataFrame(columns=["cell_id", "niche_id"]),
        }
    # --- KD-tree on ALL cells for fast neighbor queries
    tree_all = cKDTree(coords)
    # --- 2) Lone-cell removal among high cells
    idx_high = np.where(mask_high)[0]
    coords_high = coords[idx_high]
    tree_high = cKDTree(coords_high)
    neigh_lists = tree_high.query_ball_tree(tree_high, r_um)
    has_neighbors = np.fromiter(
        (max(0, len(lst) - 1) >= min_neighbors for lst in neigh_lists),
        count=len(neigh_lists), dtype=bool
    )
    idx_core = idx_high[has_neighbors]
    if idx_core.size == 0:
        return {
            "core_stats": pd.DataFrame(columns=["niche_core_id", "size", "x_centroid", "y_centroid"]),
            "core_membership": pd.DataFrame(columns=["cell_id", "niche_id"]),
            "expanded_membership": pd.DataFrame(columns=["cell_id", "niche_id"]),
        }
    # --- 3) Connected components of core cells
    coords_core = coords[idx_core]
    tree_core = cKDTree(coords_core)
    core_lists = tree_core.query_ball_tree(tree_core, r_um)
    rows, cols = [], []
    for i, lst in enumerate(core_lists):
        for j in lst:
            if i != j:
                rows.append(i)
                cols.append(j)
    m = len(idx_core)
    adj = csr_matrix((np.ones(len(rows), dtype=np.uint8), (rows, cols)), shape=(m, m))
    adj = adj.maximum(adj.T)
    n_comp, labels = connected_components(csgraph=adj, directed=False, return_labels=True)
    comp_sizes = np.bincount(labels, minlength=n_comp)
    keep_comp = comp_sizes >= min_core_size
    comp_map = {c: (k + 1) for k, c in enumerate(np.where(keep_comp)[0])}
    core_ids = np.array([comp_map.get(c, 0) for c in labels], dtype=int)
    niche_core_id = np.zeros(n, dtype=int)
    niche_core_id[idx_core] = core_ids
    adata.obs["niche_core_id"] = niche_core_id
    core_stats = []
    for cid in sorted(set(niche_core_id) - {0}):
        gidx = np.where(niche_core_id == cid)[0]
        cx, cy = coords[gidx].mean(axis=0)
        core_stats.append({"niche_core_id": cid, "size": gidx.size, "x_centroid": cx, "y_centroid": cy})
    core_stats = pd.DataFrame(core_stats).sort_values("niche_core_id")
    # --- Iterative Expansion
    if do_expand and niche_core_id.max() > 0:
        niche_expanded_id = niche_core_id.copy()  # Initialize niche_expanded_id
        expanded = True
        while expanded:
            expanded = False
            if low_thr is None:
                mask_cand = niche_expanded_id == 0
            else:
                mask_cand = (scores >= low_thr) & (niche_expanded_id == 0)
            idx_cand = np.where(mask_cand)[0]
            if idx_cand.size > 0:
                core_local_to_label = niche_expanded_id[idx_core]
                cand_lists = tree_core.query_ball_point(coords[idx_cand], r_um)
                d_nn, i_nn = tree_core.query(coords[idx_cand], k=1)
                i_nn = np.atleast_1d(i_nn)
                chosen = np.zeros(idx_cand.size, dtype=int)
                for k, lst in enumerate(cand_lists):
                    if len(lst) == 0:
                        chosen[k] = 0
                    else:
                        labels_here = core_local_to_label[np.asarray(lst)]
                        labels_here = labels_here[labels_here > 0]
                        if labels_here.size == 0:
                            chosen[k] = 0
                        else:
                            vals, counts = np.unique(labels_here, return_counts=True)
                            if (counts == counts.max()).sum() > 1:
                                chosen[k] = core_local_to_label[i_nn[k]]
                            else:
                                chosen[k] = vals[np.argmax(counts)]
                # Check if any new cells were added
                if np.any(chosen > 0):
                    expanded = True
                    niche_expanded_id[idx_cand] = chosen
    # --- Merge expanded groups if any members are within r_um of each other (recursive)
    idx_exp = np.where(niche_expanded_id > 0)[0]
    coords_exp = coords[idx_exp]
    group_ids = niche_expanded_id[idx_exp]
    tree_exp = cKDTree(coords_exp)
    m = len(idx_exp)
    merged = True
    while merged and m > 0:
        rows, cols = [], []
        for i in range(m):
            neighbors = tree_exp.query_ball_point(coords_exp[i], r_um)
            for j in neighbors:
                if i != j and group_ids[i] != group_ids[j]:
                    rows.append(i)
                    cols.append(j)
        if not rows:
            merged = False
            break
        adj = csr_matrix((np.ones(len(rows), dtype=np.uint8), (rows, cols)), shape=(m, m))
        adj = adj.maximum(adj.T)
        n_comp, labels = connected_components(csgraph=adj, directed=False, return_labels=True)
        # Assign new merged group ids
        new_group_ids = np.zeros_like(group_ids)
        for comp in np.unique(labels):
            members = np.where(labels == comp)[0]
            merged_group_id = group_ids[members].min()
            new_group_ids[members] = merged_group_id
        # If no change, stop
        if np.array_equal(group_ids, new_group_ids):
            merged = False
        else:
            group_ids = new_group_ids
            niche_expanded_id[idx_exp] = group_ids
            adata.obs["niche_expanded_id"] = niche_expanded_id
            # Rebuild tree for next iteration
            tree_exp = cKDTree(coords_exp)
    # --- Filter expanded groups by minimum size
    expanded_sizes = pd.Series(niche_expanded_id).value_counts()
    valid_groups = expanded_sizes[expanded_sizes >= min_expanded_size].index
    niche_expanded_id = np.where(np.isin(niche_expanded_id, valid_groups), niche_expanded_id, 0)
    adata.obs["niche_expanded_id"] = niche_expanded_id
    # --- Build membership DataFrames
    core_membership = pd.DataFrame({
        "cell_id": adata.obs_names,
        "niche_id": niche_core_id
    }).query("niche_id > 0")
    expanded_membership = pd.DataFrame({
        "cell_id": adata.obs_names,
        "niche_id": niche_expanded_id
    }).query("niche_id > 0")
    return {
        "core_stats": core_stats,
        "core_membership": core_membership,
        "expanded_membership": expanded_membership
    }

# Define parameters
r_um = 30.0         # Radius (microns) for neighbor search and region expansion
high_thr = threshold      # Threshold for high HEV score (core cell inclusion)
low_thr = 0.2       # Threshold for low HEV score (expansion candidates)
min_neighbors = 0   # Minimum number of neighbors for a cell to be considered part of a core
min_core_size = 1   # Minimum number of core cells to form a spatial niche
min_expanded_size = 2  # Minimum number of cells in expanded niche to be kept
do_expand = True    # Whether to expand regions to include lower-scoring neighbors

# Run the function
niche_dict = find_spatial_niches(
    adata,
    score_key="HEV_score",
    xy_keys=("x_centroid", "y_centroid"),
    r_um=r_um,
    high_thr=high_thr,
    low_thr=low_thr,
    min_neighbors=min_neighbors,
    min_core_size=min_core_size,
    do_expand=do_expand,
)

expanded_membership = niche_dict['expanded_membership']

# Add niche ids to adata
adata.obs['niche_id'] = expanded_membership.set_index('cell_id')

# Prepend HEV if the niche id is not NaN
adata.obs.loc[~adata.obs['niche_id'].isna(), 'niche_id'] = 'HEV_' + adata.obs.loc[~adata.obs['niche_id'].isna(), 'niche_id'].astype(int).astype(str)

# Make categorical
adata.obs['niche_id'] = adata.obs['niche_id'].astype('category')


# --- Plot 3: Plot HEV spatial groups

ncols=2

# Get scanpy's palette
palette = sc.pl.palettes.godsnot_102

plt.figure(figsize=(9*ncols, 7*nrows))

# add a new subplot iteratively
n=0
ax = plt.subplot(nrows, ncols, n + 1)

sc.pl.spatial(adata, color='niche_id', size=0.5, title=f'HEV spatial groups',
            crop_coord=crop_coord, ax=ax, palette=palette, na_in_legend=False, legend_loc=None, alpha_img=0.3)

ax.set_title(ax.get_title(), fontweight='bold')

# Add another subplot
ax = plt.subplot(nrows, ncols, n + 2)

sc.pl.spatial(adata, color='in_TLS', size=0.5, title=f'HEV spatial groups',
            crop_coord=crop_coord, ax=ax, palette=color_dict_TLS, na_in_legend=False, legend_loc=None, alpha_img=0.3)

adata_sub = adata[adata.obs['niche_id'].notna()]
sc.pl.spatial(adata_sub, color='niche_id', size=0.5, title=f'HEVs + TLS regions',
            crop_coord=crop_coord, ax=ax, na_in_legend=False, palette=['blue'], legend_loc=None, alpha_img=0.3)

ax.set_title(ax.get_title(), fontweight='bold')

# Remove the labels
ax.set_xlabel('')
ax.set_ylabel('')

plt.suptitle(f'{sample_id} ({subtype}) - HEV identification based on scoring (high_thresh={high_thr}, low_thresh={low_thr}, dist={r_um}, min_neighbors={min_neighbors}, min_core_size={min_core_size}, min_expanded_size={min_expanded_size})', 
                y=1.1, fontweight='bold')


# --- Plot 4: Plot the cell type annotations in the HEV niches

lineage_order = ['Epithelial', 'Immune', 'Stromal']

# Sort based on lineage first (recorded in the miniatlas df)
def get_lineage(cell_type, level='major'):
    if cell_type in miniatlas_df[f'miniatlas_{level}'].values:
        lineage = miniatlas_df[miniatlas_df[f'miniatlas_{level}'] == cell_type]['lineage'].values[0]
    else:
        lineage = 'other'
    return lineage

major_cell_types_sorted = sorted(all_cell_types, key=lambda x: (lineage_order.index(get_lineage(x)) if get_lineage(x) in lineage_order else len(lineage_order), x))

# Enable grid behind plots
plt.rcParams['axes.axisbelow'] = True

# Set fontsize
sc.set_figure_params(fontsize=10)

# Count the cell annotations
tmp_df = adata.obs[['niche_id', 'subset']]

# Drop any NaN values in niche id
tmp_df = tmp_df[~tmp_df['niche_id'].isna()]

# Count cell types
count_df = tmp_df['subset'].value_counts()

# Make a bar plot
fig, ax = plt.subplots(figsize=(3, 2))

# Extract colors for the barplot
bar_colors = [color_dict.get(ct, 'lightgrey') for ct in count_df.index]
count_df.plot(kind='bar', color=bar_colors, ax=ax, width=0.9)

# Automatically adjust the ax size
ax_pos = ax.get_position()
n_cells = len(count_df)
ax.set_position([ax_pos.x0, ax_pos.y0, n_cells * 0.06, ax_pos.height])

# Annotate each bar
for i, v in enumerate(count_df):
    # If the number is more than a 1000, use a k
    if v >= 1000:
        text = f'{v/1000:.1f}k'
    else:
        text = str(v)
    ax.text(i, v + 0.03*ax.get_ylim()[1], text, color='black', ha='center', fontsize=8)

# Adjust the ylim a bit
ax.set_ylim(0, ax.get_ylim()[1]*1.1)

# Set a title
ax.set_title(f'{sample_id} ({subtype})\nCell labels in HEVs', fontweight='bold')

# Add a legend based on major_cell_types_sorted
handles = []
labels = []
for cell_type in major_cell_types_sorted:
    if cell_type in color_dict:
        handles.append(mpatches.Patch(facecolor=color_dict[cell_type], label=cell_type, edgecolor='black', linewidth=0.8))
        labels.append(cell_type)

ax.legend(handles=handles, labels=labels, bbox_to_anchor=(1, 0.5), loc='center left', handlelength=0.8, frameon=False)

# Position the suptitle
title_x = ax.get_position().x0 + ax.get_position().width / 2

plt.suptitle(f'{sample_id} ({subtype}) - Cell type annotations in HEV niches (annotation_level={level})', y=1.2, x=title_x, fontweight='bold')


# --- Calculate the distance to nearest TLS region

# For each HEV niche, calculate the distance to the nearest TLS region

# Get the coordinates of the TLS regions
tls_coords = adata[adata.obs['in_TLS'] == True].obs[['x_centroid', 'y_centroid']].to_numpy()

# Get the coordinates of the HEV niches
hev_coords = adata[adata.obs['niche_id'].notna()].obs[['x_centroid', 'y_centroid']].to_numpy()

# Calculate the distance matrix
dist_matrix = distance_matrix(hev_coords, tls_coords)

# Find the minimum distance for each HEV niche
min_distances = dist_matrix.min(axis=1)

# Make a df with niche id and distance
hev_niches = adata[adata.obs['niche_id'].notna()].obs[['niche_id']].copy()
hev_niches['min_distance_to_TLS'] = min_distances

# Add the distance to the adata
dist_dict = hev_niches.set_index('niche_id')['min_distance_to_TLS'].to_dict()
adata.obs['min_distance_to_TLS'] = adata.obs['niche_id'].map(dist_dict)

# Also classify HEV niches based on distance to TLS
distance_threshold = 100  # microns
hev_niches['close_to_TLS'] = hev_niches['min_distance_to_TLS'] <= distance_threshold


# --- Plot 5: Plot the distance to nearest TLS region

sc.set_figure_params(fontsize=20)

ncols=1

plt.figure(figsize=(9*ncols, 7*nrows))

# add a new subplot iteratively
n=0
ax = plt.subplot(nrows, ncols, n + 1)

sc.pl.spatial(adata, color='in_TLS', size=0.5, title=f'HEV spatial groups',
            crop_coord=crop_coord, ax=ax, palette=color_dict_TLS, na_in_legend=False, legend_loc=None, alpha_img=0.3)


adata_sub = adata[adata.obs['niche_id'].notna()]
vmax=200
sc.pl.spatial(adata_sub, color='min_distance_to_TLS', size=0.5, title=f'HEV distance to nearest TLS region',
            crop_coord=crop_coord, ax=ax, cmap='viridis', na_in_legend=False, legend_loc=None, vmax=vmax, alpha_img=0.3)

ax.set_title(ax.get_title(), fontweight='bold')

plt.suptitle(f'{sample_id} ({subtype}) - Distance of HEV niches to nearest TLS region (vmax={vmax}µm)', y=1.05, fontweight='bold')


# --- Wrangle the dfs a bit

# Add niche size as a column
niche_sizes = hev_niches['niche_id'].value_counts().to_dict()
hev_niches['niche_size'] = hev_niches['niche_id'].map(niche_sizes)

# Rename the niches from largest to smallest
niche_order = hev_niches['niche_id'].value_counts().index.tolist()
niche_rename_dict = {old_name: f'HEV_{i+1}' for i, old_name in enumerate(niche_order)}
hev_niches['niche_id'] = hev_niches['niche_id'].map(niche_rename_dict)

# Extract a df with cell to niche mapping
cell_to_niche = hev_niches[['niche_id']].copy()

# Also add the HEV score of each cell
# Even if not in the cell_to_niche df
cell_to_niche = pd.merge(cell_to_niche, adata.obs[['HEV_score']], left_index=True, right_index=True, how='outer')

# Extract a df with unique niches only
hev_niches = hev_niches.drop_duplicates(subset=['niche_id']).reset_index(drop=True)

# Sort by numerical part of the niche id
hev_niches['niche_id'] = hev_niches['niche_id'].str.extract('(\d+)').astype(int)
hev_niches = hev_niches.sort_values(by='niche_id')

hev_niches['niche_id'] = 'HEV_' + hev_niches['niche_id'].astype(str)


# --- Plot 6: Histogram of distances

ncols = 2
sc.set_figure_params(fontsize=10)

plt.figure(figsize=(2, 2))

# add a new subplot iteratively
n=0
ax = plt.subplot(nrows, ncols, n + 1)

# Plot histogram of distances
hev_niches['min_distance_to_TLS'].plot(kind='hist', bins=30, color='cornflowerblue', edgecolor='black', ax=ax)

ax.set_xlabel('Distance (µm)')
ax.set_ylabel('# HEV niches')
ax.set_title(f'Distance to TLS', fontweight='bold')

# Add a new subplot
ax = plt.subplot(nrows, ncols, n + 2)

# Barplot of close vs far
distance_counts = hev_niches['close_to_TLS'].value_counts().sort_index()
distance_colors = ['lightgrey', 'cornflowerblue']
distance_counts.index = ['Far (>100µm)', 'Close (≤100µm)']
distance_counts.plot(kind='bar', color=distance_colors, edgecolor='black', ax=ax)

ax.set_ylabel('# HEV niches')
ax.set_title(f'HEV niches close vs far from TLS', fontweight='bold')

# Annotate each bar
for i, v in enumerate(distance_counts):
    ax.text(i, v + 0.03*ax.get_ylim()[1], str(v), color='black', ha='center', fontsize=8)

# Adjust the ylim a bit
ax.set_ylim(0, ax.get_ylim()[1]*1.1)

# Resize the axs
axs = plt.gcf().get_axes()
for ax in axs:
    ax_pos = ax.get_position()
    ax.set_position([ax_pos.x0, ax_pos.y0, 1, ax_pos.height])

# Reposition the axs
for n, ax in enumerate(axs):
    ax_pos = ax.get_position()
    if n == 0:
        ax.set_position([ax_pos.x0, ax_pos.y0, ax_pos.width, ax_pos.height])
    else:
        new_x = ax_pos_prev.x1 + 0.4
        ax.set_position([new_x, ax_pos.y0, ax_pos.width, ax_pos.height])
    ax_pos_prev = ax.get_position()

# Add a suptitle
title_x = (axs[0].get_position().x0 + axs[-1].get_position().x1) / 2
title_y = max(ax.get_position().y1 for ax in axs) + 0.3

plt.suptitle(f'{sample_id} ({subtype}) - Distance of HEV niches to nearest TLS region', y=title_y, x=title_x, fontweight='bold')


# --- Save the output

today = pd.Timestamp.now().strftime('%Y%m%d')
out_dir = os.path.join(res_dir, '05_xenium_mechanisms', '02_identify_HEV_regions')

out_dir = os.path.join(out_dir, today)
os.makedirs(out_dir, exist_ok=True)

out_dir_plots = os.path.join(out_dir, 'plots')
os.makedirs(out_dir_plots, exist_ok=True)

out_dir_csv = os.path.join(out_dir, 'csv_files')
os.makedirs(out_dir_csv, exist_ok=True)

# Save the plots
filename = f'{sample_id}_identify_HEV_regions.pdf'

# Save with high resolution
filename_dir=os.path.join(out_dir_plots, filename)

pp = PdfPages(filename_dir)
fig_nums = plt.get_fignums()

figs = [plt.figure(n) for n in fig_nums]

for fig in figs:
    fig.savefig(pp, format='pdf', bbox_inches='tight', dpi=300)

pp.close()


# -- Save the HEV niches

# Save both dfs
filename = f'{sample_id}_HEV_niches_cell_to_niche.csv'
cell_to_niche.to_csv(os.path.join(out_dir_csv, filename))

filename = f'{sample_id}_HEV_niches_info.csv'
hev_niches.to_csv(os.path.join(out_dir_csv, filename), index=False)

# Also save the parameters as a json
params = {
    "r_um": r_um,
    "high_thr": high_thr,
    "low_thr": low_thr,
    "min_neighbors": min_neighbors,
    "min_core_size": min_core_size,
    "min_expanded_size": min_expanded_size,
    "distance_threshold": distance_threshold
}

filename = f'{sample_id}_HEV_niches_parameters.json'
with open(os.path.join(out_dir_csv, filename), 'w') as f:
    json.dump(params, f, indent=4)


# --- Finish

print("Script successfully completed")