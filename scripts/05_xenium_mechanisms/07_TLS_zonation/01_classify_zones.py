#!/usr/bin/env python3

# --- Aim
# Classify 50 um bins inside previously-identified TLS regions into B-cell / T-cell
# zones AND germinal-centre (GC) zones in a single pass (Figures S9C/S9D and S10G/S10H
# support). Approach (identical per-bin / pseudobulk maths to the source scripts):
#   1. Bin cells into a 50 um square grid (micron coordinates)
#   2. Pseudo-bulk raw counts per bin, then normalize_total + log1p
#   3. Restrict to bins that intersect the TLS regions from 05/01
#   4. Score each bin for the B cell zone, T cell zone and Germinal center niche
#      signatures with sc.tl.score_genes (top 20 genes per niche, by rank)
#   5. Classify bins:
#        B/T threshold 0.3  (B>0.3 & T<0.3 -> B zone; T>0.3 & B<0.3 -> T zone;
#                            both>0.3 -> mixed; else unclassified)
#        GC threshold 0.8   (score>0.8 -> GC)
#   6. Group GC bins into discrete GCs (8-connectivity, merging across a single
#      empty bin)
#   7. Propagate the bin labels to the cells contained in each bin


# --- Define environment
# env: xen-hev

# --- Load packages

import argparse
import os
import sys

import anndata as ad
import scanpy as sc
import pandas as pd
import numpy as np
import scipy.sparse as sp

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle
from matplotlib.cm import ScalarMappable

import geopandas as gpd

from scipy.ndimage import label as cc_label
from scipy.spatial import cKDTree

mpl.rcParams['pdf.fonttype'] = 42
sc.set_figure_params(facecolor="white", figsize=(8, 8), fontsize=20)
# Dense spatial layers (50 um bin squares, per-cell scatter) are drawn with
# rasterized=True; this DPI keeps them crisp while the PDF stays vector for text/axes.
mpl.rcParams['savefig.dpi'] = 200


# --- Argparse arguments

parser = argparse.ArgumentParser(description="Classify TLS B/T cell zones and GC zones per 50um bin",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

# Runner contract (submit_array_jobs_nci.sh -> run_single_sample_job.sh)
parser.add_argument("--h5ad_dir", help="Path to directory for .h5ad input files")
parser.add_argument("--sample_id", help="Sample (section_id) to process")
parser.add_argument("--project_dir", help="Path to the root of the project directory")
parser.add_argument("--res_dir", help="Path to results directory")
parser.add_argument("--sample_sheet", help="Path to sample sheet CSV file", default=None)

# Analysis parameters (Methods-fixed defaults; exposed but never required)
parser.add_argument("--zone_set", choices=['all', 'bt', 'gc'], default='all',
                    help="Which zones to classify (default 'all' does both B/T and GC)")
parser.add_argument("--grid_um", type=float, default=50.0, help="Bin size in microns")
parser.add_argument("--top_genes", type=int, default=20, help="Top N signature genes (by rank) per niche")
parser.add_argument("--bt_threshold", type=float, default=0.3, help="Score threshold for the B/T zone call")
parser.add_argument("--gc_threshold", type=float, default=0.8, help="Score threshold for the GC call")
parser.add_argument("--cc_connectivity", type=int, default=8, choices=[4, 8],
                    help="4- or 8-connectivity for discrete GC connected components")
parser.add_argument("--merge_gap_bins", type=int, default=1,
                    help="Merge GC components separated by <= this many empty bins (0 disables)")
parser.add_argument("--min_gc_bins", type=int, default=1,
                    help="Minimum bins per connected component to keep as a discrete GC")
parser.add_argument("--seed", type=int, default=42, help="Random seed for sc.tl.score_genes")

args = vars(parser.parse_args())

h5ad_dir        = args["h5ad_dir"]
sample_id       = args["sample_id"]
project_dir     = args["project_dir"]
res_dir         = args["res_dir"]
zone_set        = args["zone_set"]
grid_um         = args["grid_um"]
top_genes       = args["top_genes"]
bt_threshold    = args["bt_threshold"]
gc_threshold    = args["gc_threshold"]
cc_connectivity = args["cc_connectivity"]
merge_gap_bins  = args["merge_gap_bins"]
min_gc_bins     = args["min_gc_bins"]
seed            = args["seed"]

do_bt = zone_set in ('all', 'bt')
do_gc = zone_set in ('all', 'gc')

# Local helper functions
sys.path.append(os.path.join(project_dir, 'scripts'))
from helper_functions import save_multi_image


# --- Sample metadata

if args.get("sample_sheet"):
    sample_sheet_path = args["sample_sheet"]
else:
    sample_sheet_path = os.path.join(project_dir, 'config', 'sample_sheet_xenium.csv')

sample_sheet    = pd.read_csv(sample_sheet_path)
tissue_metadata = pd.read_csv(os.path.join(project_dir, 'resources', 'tissue_metadata.csv'))

tissue_id = sample_sheet.loc[sample_sheet['section_id'] == sample_id, 'tissue_id'].values[0]
subtype   = tissue_metadata.loc[tissue_metadata['tissue_id'].astype(str) == str(tissue_id),
                                'clinical_subtype'].values[0]


# --- Load data

adata = ad.read_h5ad(os.path.join(h5ad_dir, 'xenium', f'{sample_id}.h5ad'))

# Raw counts live in layers['filtered_counts'] for Xenium Prime 5K (.X is log-norm)
raw_counts = adata.layers['filtered_counts']
if not sp.issparse(raw_counts):
    raw_counts = sp.csr_matrix(raw_counts)
else:
    raw_counts = raw_counts.tocsr()

# Micron coordinates (x_centroid/y_centroid); obsm['spatial'] is image/pixel space in P08
xy = adata.obs[['x_centroid', 'y_centroid']].values.astype(float)


def square_crop(x, y, border=0.1):
    """Square crop (x0, x1, y0, y1) in micron space around all cells."""
    x0, x1 = float(x.min()), float(x.max())
    y0, y1 = float(y.min()), float(y.max())
    xr, yr = x1 - x0, y1 - y0
    cx, cy = x0 + xr / 2, y0 + yr / 2
    half = max(xr, yr) / 2
    b = (border / 2) * max(xr, yr)
    return (cx - half - b, cx + half + b, cy - half - b, cy + half + b)


crop_coord = square_crop(xy[:, 0], xy[:, 1])


# --- Load TLS regions (from existing 05/01_identify_TLS_regions)

tls_regions_dir = os.path.join(res_dir, '05_xenium_mechanisms', '01_identify_TLS_regions', 'csv_files')

tls_geojson = os.path.join(tls_regions_dir, f'{sample_id}_TLS_regions.geojson')
tls_gdf = gpd.read_file(tls_geojson)
print(f'Loaded {len(tls_gdf)} TLS regions for {sample_id}')

cell_tls_csv = os.path.join(tls_regions_dir, f'{sample_id}_cell_TLS_annotations.csv')
cell_tls = pd.read_csv(cell_tls_csv, index_col=0)
adata.obs['in_TLS'] = cell_tls['in_TLS'].reindex(adata.obs_names).astype(str)
adata.obs['TLS_id'] = cell_tls['TLS_id'].reindex(adata.obs_names).fillna(0).astype(int)


# --- Output paths

out_base  = os.path.join(res_dir, '05_xenium_mechanisms', '07_TLS_zonation', '01_classify_zones')
plots_dir = os.path.join(out_base, 'plots')
csv_dir   = os.path.join(out_base, 'csv_files')
os.makedirs(plots_dir, exist_ok=True)
os.makedirs(csv_dir, exist_ok=True)


# --- Bin cells into a 50 um grid

xmin, ymin = xy.min(axis=0)
xmax, ymax = xy.max(axis=0)
nx = int(np.ceil((xmax - xmin) / grid_um))
ny = int(np.ceil((ymax - ymin) / grid_um))

ix = np.minimum(((xy[:, 0] - xmin) / grid_um).astype(int), nx - 1)
iy = np.minimum(((xy[:, 1] - ymin) / grid_um).astype(int), ny - 1)

cell_bin_flat = iy * nx + ix
n_bins_full = nx * ny
n_cells = adata.n_obs


# --- Pseudo-bulk: sum raw counts per bin

agg = sp.csr_matrix(
    (np.ones(n_cells, dtype=np.float32),
     (cell_bin_flat, np.arange(n_cells))),
    shape=(n_bins_full, n_cells),
)

bin_counts = agg @ raw_counts.astype(np.float32)
bin_cell_count = np.asarray(agg.sum(axis=1)).ravel().astype(int)

nonempty = np.where(bin_cell_count > 0)[0]
bin_counts = bin_counts[nonempty]
bin_cell_count = bin_cell_count[nonempty]

bin_iy = nonempty // nx
bin_ix = nonempty %  nx
bin_x = xmin + (bin_ix + 0.5) * grid_um
bin_y = ymin + (bin_iy + 0.5) * grid_um

print(f'Pseudo-bulk: {len(nonempty)} non-empty bins of {n_bins_full} total ({n_cells} cells)')


# --- Restrict to bins inside TLS regions (using bin centroid)

bin_points = gpd.GeoDataFrame(
    {'bin_idx': np.arange(len(nonempty))},
    geometry=gpd.points_from_xy(bin_x, bin_y),
    crs=tls_gdf.crs,
)
bins_in_tls = gpd.sjoin(bin_points, tls_gdf, predicate='intersects', how='inner')
# A bin centroid on a shared edge between two touching TLS polygons can match both;
# keep one match per bin.
in_tls_idx = bins_in_tls['bin_idx'].drop_duplicates().values


def write_empty_outputs():
    """Write empty per-sample CSVs so downstream cohort scripts stay consistent."""
    if do_bt:
        pd.DataFrame(columns=['cell_id', 'TLS_zone_T_score', 'TLS_zone_B_score', 'TLS_zone']
                     ).to_csv(os.path.join(csv_dir, f'{sample_id}_BT_zone_classification.csv'), index=False)
        pd.DataFrame(columns=['bin_x', 'bin_y', 'bin_ix', 'bin_iy', 'n_cells',
                              'T_cell_zone_score', 'B_cell_zone_score', 'TLS_zone']
                     ).to_csv(os.path.join(csv_dir, f'{sample_id}_BT_zone_bins.csv'), index=False)
    if do_gc:
        pd.DataFrame(columns=['cell_id', 'GC_zone_score', 'GC_zone', 'discrete_GC']
                     ).to_csv(os.path.join(csv_dir, f'{sample_id}_GC_zone_classification.csv'), index=False)
        pd.DataFrame(columns=['bin_x', 'bin_y', 'bin_ix', 'bin_iy', 'n_cells',
                              'GC_score', 'GC_zone', 'discrete_GC']
                     ).to_csv(os.path.join(csv_dir, f'{sample_id}_GC_zone_bins.csv'), index=False)


if len(in_tls_idx) == 0:
    print(f'No bins fall inside TLS regions for {sample_id}; nothing to score.')
    write_empty_outputs()
    print('Script successfully completed (no TLS bins)')
    raise SystemExit(0)

print(f'{len(in_tls_idx)} of {len(nonempty)} bins intersect TLS regions')


# --- Build pseudo-bulk anndata for in-TLS bins only

bin_adata = ad.AnnData(
    X=bin_counts[in_tls_idx],
    obs=pd.DataFrame({
        'n_cells': bin_cell_count[in_tls_idx],
        'bin_x': bin_x[in_tls_idx],
        'bin_y': bin_y[in_tls_idx],
        'bin_ix': bin_ix[in_tls_idx],
        'bin_iy': bin_iy[in_tls_idx],
    }, index=[f'bin_{i}' for i in in_tls_idx]),
    var=adata.var.copy(),
)
bin_adata.obsm['spatial'] = np.column_stack([bin_adata.obs['bin_x'], bin_adata.obs['bin_y']])

# Normalize per bin to a fixed library size (mirrors the visium HD per-spot scaling), then log1p
sc.pp.normalize_total(bin_adata, target_sum=1e4)
sc.pp.log1p(bin_adata)


# --- Load niche signatures and score each bin

niche_map = {
    'B cell zone':      'B_cell_zone',
    'T cell zone':      'T_cell_zone',
    'Germinal center':  'GC',
}
score_keys = {'B_cell_zone': 'B_cell_zone_score',
              'T_cell_zone': 'T_cell_zone_score',
              'GC':          'GC_score'}

niche_sig = pd.read_csv(os.path.join(project_dir, 'resources', 'gene_lists',
                                     'TLS_niche_signatures', 'niche_signatures.csv'))

panel_set = set(bin_adata.var_names)
genes_used = {}

for niche_name, internal in niche_map.items():
    if internal in ('B_cell_zone', 'T_cell_zone') and not do_bt:
        continue
    if internal == 'GC' and not do_gc:
        continue

    sig = niche_sig[niche_sig['niche'] == niche_name].sort_values('rank')
    ranked_genes = [g for g in sig['gene_symbol'].tolist() if isinstance(g, str)]
    # Top N signature genes that are present on the Xenium panel: restrict to panel
    # genes FIRST, then take the top N by rank. (Selecting the top N first and then
    # dropping off-panel genes would leave fewer than N.)
    ranked_on_panel = [g for g in ranked_genes if g in panel_set]
    present = ranked_on_panel[:top_genes]
    print(f'{niche_name}: using top {len(present)} on-panel signature genes '
          f'({len(ranked_on_panel)}/{len(ranked_genes)} signature genes are on the panel)')

    sc.tl.score_genes(bin_adata, present, random_state=seed, score_name=score_keys[internal])
    genes_used[internal] = present


# --- Classify bins: B / T / mixed

if do_bt:
    T_key = score_keys['T_cell_zone']
    B_key = score_keys['B_cell_zone']

    t_above = bin_adata.obs[T_key] > bt_threshold
    b_above = bin_adata.obs[B_key] > bt_threshold

    bt_class = pd.Series(pd.NA, index=bin_adata.obs_names, dtype='object')
    bt_class[b_above] = 'B_cell_zone'
    bt_class[t_above] = 'T_cell_zone'
    bt_class[t_above & b_above] = 'mixed'
    bin_adata.obs['TLS_zone'] = bt_class
    print('--- B/T zone bins ---')
    print(bin_adata.obs['TLS_zone'].value_counts(dropna=False))


# --- Classify bins: GC, then group into discrete GCs

def merge_close_components(comp_labels, points_ix, points_iy, merge_gap_bins):
    """Union-find merge of components whose nearest bins are within (merge_gap_bins + 1)
    Chebyshev steps of each other. Returns consecutively re-numbered labels (1..K)."""
    n = int(comp_labels.max()) if len(comp_labels) else 0
    if merge_gap_bins <= 0 or n <= 1:
        unique = sorted(set(int(c) for c in comp_labels))
        remap = {old: new + 1 for new, old in enumerate(unique)}
        return np.array([remap[int(c)] for c in comp_labels])

    points = np.column_stack([points_ix, points_iy])
    tree = cKDTree(points)
    pairs = tree.query_pairs(r=merge_gap_bins + 1, p=np.inf)

    parent = list(range(n + 1))
    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    def union(a, b):
        pa, pb = find(a), find(b)
        if pa != pb:
            parent[pa] = pb

    for a, b in pairs:
        ca, cb = int(comp_labels[a]), int(comp_labels[b])
        if ca != cb:
            union(ca, cb)

    merged = np.array([find(int(c)) for c in comp_labels])
    unique = sorted(set(int(m) for m in merged))
    remap = {old: new + 1 for new, old in enumerate(unique)}
    return np.array([remap[int(m)] for m in merged])


if do_gc:
    gc_above = bin_adata.obs[score_keys['GC']] > gc_threshold
    gc_class = pd.Series(pd.NA, index=bin_adata.obs_names, dtype='object')
    gc_class[gc_above] = 'GC'
    bin_adata.obs['GC_zone'] = gc_class

    # Discrete GCs: connected components of GC-positive bins on the bin grid
    bin_adata.obs['discrete_GC'] = pd.NA
    gc_mask = gc_class.notna().values
    if gc_mask.any():
        g_ix = bin_adata.obs['bin_ix'].astype(int).values[gc_mask]
        g_iy = bin_adata.obs['bin_iy'].astype(int).values[gc_mask]

        if cc_connectivity == 8:
            structure = np.ones((3, 3), dtype=int)
        else:
            structure = np.array([[0, 1, 0], [1, 1, 1], [0, 1, 0]], dtype=int)

        ix0, iy0 = g_ix.min(), g_iy.min()
        mask2d = np.zeros((g_iy.max() - iy0 + 1, g_ix.max() - ix0 + 1), dtype=bool)
        mask2d[g_iy - iy0, g_ix - ix0] = True
        labeled, _ = cc_label(mask2d, structure=structure)
        initial_labels = labeled[g_iy - iy0, g_ix - ix0]
        merged_labels = merge_close_components(initial_labels, g_ix, g_iy, merge_gap_bins)

        # Drop components smaller than min_gc_bins, re-number the survivors 1..K
        sizes = np.bincount(merged_labels)
        keep = {lab for lab in range(1, len(sizes)) if sizes[lab] >= min_gc_bins}
        remap = {lab: new + 1 for new, lab in enumerate(sorted(keep))}
        final = np.array([remap.get(int(m), 0) for m in merged_labels])

        gc_component = pd.Series(pd.NA, index=bin_adata.obs_names, dtype='object')
        gc_component.iloc[np.where(gc_mask)[0]] = [c if c > 0 else pd.NA for c in final]
        bin_adata.obs['discrete_GC'] = gc_component

    n_discrete = int(pd.to_numeric(bin_adata.obs['discrete_GC'], errors='coerce').dropna().nunique())
    print('--- GC zone bins ---')
    print(bin_adata.obs['GC_zone'].value_counts(dropna=False))
    print(f'Discrete GCs (conn={cc_connectivity}, merge_gap_bins={merge_gap_bins}, '
          f'min_gc_bins={min_gc_bins}): {n_discrete}')


# --- Map bin classifications back to cells

flat_to_nonempty = -np.ones(n_bins_full, dtype=int)
flat_to_nonempty[nonempty] = np.arange(len(nonempty))
cell_nonempty_idx = flat_to_nonempty[cell_bin_flat]  # >= 0 (every cell sits in a bin)

if do_bt:
    bt_lookup = pd.Series(bin_adata.obs['TLS_zone'].values, index=in_tls_idx)
    t_lookup  = pd.Series(bin_adata.obs[score_keys['T_cell_zone']].values, index=in_tls_idx)
    b_lookup  = pd.Series(bin_adata.obs[score_keys['B_cell_zone']].values, index=in_tls_idx)
    adata.obs['TLS_zone_T_score'] = t_lookup.reindex(cell_nonempty_idx).values
    adata.obs['TLS_zone_B_score'] = b_lookup.reindex(cell_nonempty_idx).values
    adata.obs['TLS_zone'] = pd.Categorical(bt_lookup.reindex(cell_nonempty_idx).values,
                                           categories=['B_cell_zone', 'T_cell_zone', 'mixed'])

if do_gc:
    gc_lookup    = pd.Series(bin_adata.obs['GC_zone'].values, index=in_tls_idx)
    gc_sc_lookup = pd.Series(bin_adata.obs[score_keys['GC']].values, index=in_tls_idx)
    gc_id_lookup = pd.Series(bin_adata.obs['discrete_GC'].values, index=in_tls_idx)
    adata.obs['GC_zone_score'] = gc_sc_lookup.reindex(cell_nonempty_idx).values
    adata.obs['GC_zone'] = pd.Categorical(gc_lookup.reindex(cell_nonempty_idx).values,
                                          categories=['GC'])
    adata.obs['discrete_GC'] = gc_id_lookup.reindex(cell_nonempty_idx).values


# --- Plot helpers (micron space)

def draw_bin_rectangles(ax, xs, ys, colors, grid_um, alpha=0.85, linewidth=0):
    rects = [Rectangle((x - grid_um / 2, y - grid_um / 2), grid_um, grid_um)
             for x, y in zip(xs, ys)]
    pc = PatchCollection(rects, facecolors=colors, edgecolors='none',
                         alpha=alpha, linewidth=linewidth, rasterized=True)
    ax.add_collection(pc)


def setup_spatial_ax(ax, crop_coord):
    ax.set_xlim(crop_coord[0], crop_coord[1])
    ax.set_ylim(crop_coord[2], crop_coord[3])
    ax.invert_yaxis()
    ax.set_aspect('equal')
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_xlabel(''); ax.set_ylabel('')


def grey_cells(ax):
    ax.scatter(xy[:, 0], xy[:, 1], s=2, c='lightgray', marker='.', linewidths=0, alpha=0.6,
               rasterized=True)


vmax_score = 0.8
score_cmap = plt.get_cmap('viridis')
score_norm = mpl.colors.Normalize(vmin=0, vmax=vmax_score)


# --- Plot: B/T zones

if do_bt:
    bt_color_dict = {
        'B_cell_zone': mpl.colors.cnames['cornflowerblue'],
        'T_cell_zone': mpl.colors.cnames['crimson'],
        'mixed':       '#921e95',
    }

    fig, axs = plt.subplots(1, 3, figsize=(9 * 3, 7))
    for ax, (label, key) in zip(axs[:2], [('T cell zone', T_key), ('B cell zone', B_key)]):
        cols = score_cmap(score_norm(bin_adata.obs[key].values.clip(0, vmax_score)))
        draw_bin_rectangles(ax, bin_adata.obs['bin_x'], bin_adata.obs['bin_y'], cols, grid_um, alpha=1.0)
        setup_spatial_ax(ax, crop_coord)
        plt.colorbar(ScalarMappable(norm=score_norm, cmap=score_cmap), ax=ax, shrink=0.6, label='Score')
        ax.set_title(f'{sample_id} ({subtype})\n{label} score (per bin)', fontweight='bold')

    ax = axs[2]
    grey_cells(ax)
    plot_bins = bin_adata.obs.dropna(subset=['TLS_zone']).copy()
    if len(plot_bins):
        cols = [bt_color_dict[z] for z in plot_bins['TLS_zone']]
        draw_bin_rectangles(ax, plot_bins['bin_x'], plot_bins['bin_y'], cols, grid_um, alpha=0.85)
    tls_gdf.plot(ax=ax, facecolor='none', edgecolor='black', linewidth=0.8)
    setup_spatial_ax(ax, crop_coord)
    ax.set_title(f'{sample_id} ({subtype})\nBin classification (threshold={bt_threshold})', fontweight='bold')
    handles = [mlines.Line2D([], [], color=bt_color_dict[z], marker='s', lw=0, markersize=10,
                             label=z.replace('_', ' ')) for z in ['B_cell_zone', 'T_cell_zone', 'mixed']]
    ax.legend(handles=handles, loc='upper left', bbox_to_anchor=(1.02, 1.0), frameon=False, title='TLS zone')
    plt.suptitle(f'{sample_id} ({subtype}) - TLS B/T zone classification per bin '
                 f'(grid={int(grid_um)}um, top_genes={top_genes})', y=1.05, fontweight='bold')


# --- Plot: GC zones

if do_gc:
    gc_color = mpl.colors.cnames['darkorchid']

    fig, axs = plt.subplots(1, 2, figsize=(9 * 2, 7))
    ax = axs[0]
    cols = score_cmap(score_norm(bin_adata.obs[score_keys['GC']].values.clip(0, vmax_score)))
    draw_bin_rectangles(ax, bin_adata.obs['bin_x'], bin_adata.obs['bin_y'], cols, grid_um, alpha=1.0)
    tls_gdf.plot(ax=ax, facecolor='none', edgecolor='black', linewidth=0.6)
    setup_spatial_ax(ax, crop_coord)
    plt.colorbar(ScalarMappable(norm=score_norm, cmap=score_cmap), ax=ax, shrink=0.6, label='Score')
    ax.set_title(f'{sample_id} ({subtype})\nGC score (per bin)', fontweight='bold')

    ax = axs[1]
    grey_cells(ax)
    gc_plot_bins = bin_adata.obs.dropna(subset=['GC_zone']).copy()
    if len(gc_plot_bins):
        draw_bin_rectangles(ax, gc_plot_bins['bin_x'], gc_plot_bins['bin_y'],
                            [gc_color] * len(gc_plot_bins), grid_um, alpha=0.85)
    tls_gdf.plot(ax=ax, facecolor='none', edgecolor='black', linewidth=0.8)
    setup_spatial_ax(ax, crop_coord)
    n_disc = int(pd.to_numeric(bin_adata.obs['discrete_GC'], errors='coerce').dropna().nunique())
    ax.set_title(f'{sample_id} ({subtype})\nGC bins (threshold={gc_threshold}, {n_disc} discrete GCs)',
                 fontweight='bold')
    handles = [mlines.Line2D([], [], color=gc_color, marker='s', lw=0, markersize=10, label='GC')]
    ax.legend(handles=handles, loc='upper left', bbox_to_anchor=(1.02, 1.0), frameon=False, title='GC zone')
    plt.suptitle(f'{sample_id} ({subtype}) - GC zone classification per bin '
                 f'(grid={int(grid_um)}um, top_genes={top_genes})', y=1.05, fontweight='bold')


# --- Save outputs

if do_bt:
    cell_out = adata.obs.loc[adata.obs['TLS_zone'].notna(),
                             ['TLS_zone_T_score', 'TLS_zone_B_score', 'TLS_zone']].copy()
    cell_out.index.name = 'cell_id'
    cell_out.to_csv(os.path.join(csv_dir, f'{sample_id}_BT_zone_classification.csv'))

    bin_out = bin_adata.obs[['bin_x', 'bin_y', 'bin_ix', 'bin_iy', 'n_cells',
                             T_key, B_key, 'TLS_zone']].copy()
    bin_out.to_csv(os.path.join(csv_dir, f'{sample_id}_BT_zone_bins.csv'), index=False)

    for internal in ('T_cell_zone', 'B_cell_zone'):
        pd.DataFrame({'genes_used': genes_used[internal]}).to_csv(
            os.path.join(csv_dir, f'{sample_id}_{internal}_genes_used.csv'), index=False)

if do_gc:
    gc_cell_out = adata.obs.loc[adata.obs['GC_zone'].notna(),
                                ['GC_zone_score', 'GC_zone', 'discrete_GC']].copy()
    gc_cell_out.index.name = 'cell_id'
    gc_cell_out.to_csv(os.path.join(csv_dir, f'{sample_id}_GC_zone_classification.csv'))

    gc_bin_out = bin_adata.obs[['bin_x', 'bin_y', 'bin_ix', 'bin_iy', 'n_cells',
                                score_keys['GC'], 'GC_zone', 'discrete_GC']].copy()
    gc_bin_out.rename(columns={score_keys['GC']: 'GC_score'}, inplace=True)
    gc_bin_out.to_csv(os.path.join(csv_dir, f'{sample_id}_GC_zone_bins.csv'), index=False)

    pd.DataFrame({'genes_used': genes_used['GC']}).to_csv(
        os.path.join(csv_dir, f'{sample_id}_GC_genes_used.csv'), index=False)

pdf_path = os.path.join(plots_dir, f'{sample_id}_zone_classification.pdf')
save_multi_image(pdf_path)
print(f'Saved: {pdf_path}')


# --- Finish

print("Script successfully completed")
