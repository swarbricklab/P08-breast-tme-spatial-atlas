#!/usr/bin/env python3

# --- Aim
# Classify each intra-TLS HEV niche by its distance to the nearest germinal
# center (GC) and label it close (<= 500 um) or far (> 500 um) from a GC.
# Upstream for Extended Data Fig. 10i,j (HEV endothelial cells stratified by GC proximity).
#
# Only HEV niches that sit within a TLS are considered (close_to_TLS == True,
# defined as <= 100 um from a TLS in 05/02_identify_HEV_regions). For every such
# niche we save:
#   - centroid_x, centroid_y     (mean constituent-cell position, in microns)
#   - min_distance_to_GC         (continuous, in um; inf when the sample has no GCs)
#   - close_to_GC                (bool, distance <= near_gc_um; False when inf)
#
# Samples with no GC bins (e.g. the GC-negative samples in the cohort) still
# produce a valid CSV in which every niche has min_distance_to_GC = inf and
# close_to_GC = False. Downstream cohort scripts auto-skip such samples by
# detecting that no niche has a finite GC distance.

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

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
from matplotlib.cm import ScalarMappable

mpl.rcParams['pdf.fonttype'] = 42
sc.set_figure_params(facecolor='white', figsize=(8, 8), fontsize=14)


# --- Argparse arguments

parser = argparse.ArgumentParser(
    description="Per-sample classification of intra-TLS HEV niches by distance to GC",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
)

parser.add_argument("--h5ad_dir",    required=True, help="Path to .h5ad input dir")
parser.add_argument("--sample_id",   required=True, help="Sample (section_id) to process")
parser.add_argument("--project_dir", required=True, help="Project root directory")
parser.add_argument("--res_dir",     required=True, help="Path to results directory")
parser.add_argument("--sample_sheet", default=None,
                    help="Path to sample sheet CSV (falls back to "
                         "config/sample_sheet_xenium.csv)")
parser.add_argument("--gc_threshold", type=float, default=0.8,
                    help="GC score threshold used upstream in 07_TLS_zonation "
                         "(recorded for provenance)")
parser.add_argument("--grid_um", type=float, default=50.0,
                    help="GC bin size in microns (matches 07_TLS_zonation grid)")
parser.add_argument("--near_gc_um", type=float, default=500.0,
                    help="HEV niche is 'close to GC' if min distance <= this (um)")
parser.add_argument("--gc_bins_dir", default=None,
                    help="Directory holding {sample}_GC_zone_bins.csv (default: "
                         "res_dir/05_xenium_mechanisms/07_TLS_zonation/"
                         "01_classify_zones/csv_files)")

args = vars(parser.parse_args())

h5ad_dir       = args["h5ad_dir"]
sample_id      = args["sample_id"]
project_dir    = args["project_dir"]
res_dir        = args["res_dir"]
gc_threshold   = args["gc_threshold"]
grid_um        = args["grid_um"]
near_gc_um     = args["near_gc_um"]
gc_bins_dir    = args["gc_bins_dir"]

threshold_tag = str(gc_threshold).replace('.', '')


# --- Local helper import

sys.path.append(os.path.join(project_dir, 'scripts'))
from helper_functions import save_multi_image


# --- Sample metadata

if args.get("sample_sheet") and str(args["sample_sheet"]) != 'None':
    sample_sheet_path = args["sample_sheet"]
else:
    sample_sheet_path = os.path.join(project_dir, 'config', 'sample_sheet_xenium.csv')

sample_sheet    = pd.read_csv(sample_sheet_path)
tissue_metadata = pd.read_csv(os.path.join(project_dir, 'resources', 'tissue_metadata.csv'))

tissue_id = sample_sheet.loc[sample_sheet['section_id'] == sample_id, 'tissue_id'].values[0]
subtype   = tissue_metadata.loc[tissue_metadata['tissue_id'].astype(str) == str(tissue_id),
                                'clinical_subtype'].values[0]


# --- Load data (xenium adata; .X = log1p counts)

adata = ad.read_h5ad(os.path.join(h5ad_dir, 'xenium', f'{sample_id}.h5ad'))
# All diagnostic plotting is in native micron space (x_centroid/y_centroid), matching the
# niche centroids and GC bins; the crop is the micron bounding box of all cells.
cell_um = adata.obs[['x_centroid', 'y_centroid']].to_numpy()
_pad = 100.0
crop_coord = [cell_um[:, 0].min() - _pad, cell_um[:, 0].max() + _pad,
              cell_um[:, 1].min() - _pad, cell_um[:, 1].max() + _pad]


# --- Load HEV niche assignments and info (05/02)

hev_dir = os.path.join(res_dir, '05_xenium_mechanisms', '02_identify_HEV_regions')
hev_most_recent = hev_dir
print(f'HEV niches loaded from: {hev_most_recent}')

niche_info = pd.read_csv(
    os.path.join(hev_most_recent, 'csv_files', f'{sample_id}_HEV_niches_info.csv'))
cell_to_niche = pd.read_csv(
    os.path.join(hev_most_recent, 'csv_files', f'{sample_id}_HEV_niches_cell_to_niche.csv'),
    index_col=0,
)

intra_tls_ids = set(niche_info.loc[niche_info['close_to_TLS'], 'niche_id'].astype(str))
print(f'{len(intra_tls_ids)} intra-TLS HEV niches '
      f'({len(niche_info) - len(intra_tls_ids)} non-TLS dropped)')


# --- Compute niche centroids (mean of constituent cell positions, in microns)

cell_to_niche = cell_to_niche.dropna(subset=['niche_id'])
cell_to_niche['niche_id'] = cell_to_niche['niche_id'].astype(str)

valid_cells = cell_to_niche.index.intersection(adata.obs_names)
if len(valid_cells) < len(cell_to_niche):
    print(f'WARNING: {len(cell_to_niche) - len(valid_cells)} cell barcodes from '
          f'cell_to_niche not found in adata; skipping those')

cell_idx = pd.Series(np.arange(adata.n_obs), index=adata.obs_names)
# Native micron centroids (x_centroid/y_centroid), matching the GC bins' micron frame;
# obsm['spatial'] is image/pixel space in this project and must not be mixed with the bins.
positions = adata.obs[['x_centroid', 'y_centroid']].to_numpy()[cell_idx[valid_cells].values]
niches_for_cells = cell_to_niche.loc[valid_cells, 'niche_id'].values

centroid_df = (pd.DataFrame({'x': positions[:, 0], 'y': positions[:, 1],
                             'niche_id': niches_for_cells})
                 .groupby('niche_id')
                 .agg(centroid_x=('x', 'mean'),
                      centroid_y=('y', 'mean'),
                      n_cells=('x', 'size')))

# Restrict to intra-TLS HEV niches
centroid_df = centroid_df.loc[centroid_df.index.intersection(intra_tls_ids)]
print(f'Computed centroids for {len(centroid_df)} intra-TLS niches')


# --- Load GC bins for this sample
#
# GC bins come from 07_TLS_zonation/01_classify_zones.py
# ({s}_GC_zone_bins.csv, columns bin_x, bin_y, GC_zone, discrete_GC). Missing or
# empty files are treated as "sample has no GCs", so GC-negative samples remain
# valid and are auto-skipped downstream.

if gc_bins_dir is None:
    gc_bins_dir = os.path.join(res_dir, '05_xenium_mechanisms', '07_TLS_zonation',
                               '01_classify_zones', 'csv_files')
gc_csv = os.path.join(gc_bins_dir, f'{sample_id}_GC_zone_bins.csv')

if os.path.isfile(gc_csv):
    gc_bins = pd.read_csv(gc_csv)
    if 'GC_zone' in gc_bins.columns:
        gc_bins = gc_bins.dropna(subset=['GC_zone'])
else:
    print(f'WARNING: no GC bins file at {gc_csv} — treating sample as GC-negative')
    gc_bins = pd.DataFrame(columns=['bin_x', 'bin_y', 'GC_zone'])
print(f'{len(gc_bins)} GC bins for {sample_id} at threshold {gc_threshold}')


# --- Distance from each niche centroid to the nearest GC bin edge

n_niches = len(centroid_df)
if n_niches == 0:
    print('No intra-TLS niches — writing empty output and exiting.')
    min_distances = np.array([], dtype=float)
elif len(gc_bins) == 0:
    print('No GC bins for this sample — all distances set to inf.')
    min_distances = np.full(n_niches, np.inf)
else:
    niche_xy = centroid_df[['centroid_x', 'centroid_y']].values        # (n_niches, 2)
    gc_xy    = gc_bins[['bin_x', 'bin_y']].values                       # (n_gc, 2)
    dx = np.abs(niche_xy[:, None, 0] - gc_xy[None, :, 0])
    dy = np.abs(niche_xy[:, None, 1] - gc_xy[None, :, 1])
    r  = grid_um / 2
    min_distances = np.sqrt(np.maximum(dx - r, 0) ** 2 + np.maximum(dy - r, 0) ** 2).min(axis=1)


# --- Build output table

distance_series = pd.Series(min_distances, index=centroid_df.index)

output = centroid_df.copy()
output['sample_id']          = sample_id
output['subtype']            = subtype
output['min_distance_to_GC'] = distance_series.values
output['close_to_GC']        = distance_series.values <= near_gc_um  # inf > val -> False
output.index.name            = 'niche_id'

print('\n--- Summary ---')
print(f'Niches close_to_GC (<= {near_gc_um} um): '
      f'{int(output["close_to_GC"].sum())} / {len(output)}')


# --- Diagnostic spatial plots
#   Page 1: intra-TLS HEV niche centroids colored by continuous min distance to GC
#   Page 2: binary close (<= near_gc_um) vs far classification
# GC bins are drawn as dark-orchid squares; all cells in light gray.

BINARY_COLORS = {'close to GC': '#1f77b4', 'far from GC': '#8c6d31'}


def setup_spatial_ax(ax, title):
    # Faint micron-space scatter of all cells for tissue context (replaces the
    # image-space sc.pl.spatial background so the overlay shares the micron frame).
    ax.scatter(cell_um[:, 0], cell_um[:, 1], s=1, c='lightgray',
               linewidths=0, rasterized=True, zorder=1)
    if len(gc_bins) > 0:
        for _, row in gc_bins.iterrows():
            bx, by = row['bin_x'], row['bin_y']
            ax.add_patch(mpatches.Rectangle(
                (bx - grid_um / 2, by - grid_um / 2), grid_um, grid_um,
                facecolor='darkorchid', edgecolor='none', alpha=0.55, zorder=2))
    ax.set_xlim(crop_coord[0], crop_coord[1])
    ax.set_ylim(crop_coord[2], crop_coord[3])
    ax.invert_yaxis()
    ax.set_aspect('equal')
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_xlabel(''); ax.set_ylabel('')
    ax.set_title(title, fontweight='bold')


# --- Page 1: continuous distance

fig, ax = plt.subplots(figsize=(9, 8))
setup_spatial_ax(ax, f'{sample_id} ({subtype})\n'
                 'Intra-TLS HEV niches by min distance to GC (um)')
if len(output) > 0:
    finite = np.isfinite(output['min_distance_to_GC'].values)
    finite_dist = output.loc[finite, 'min_distance_to_GC'].values
    if finite.any():
        vmax = max(near_gc_um * 3, float(np.nanmax(finite_dist)))
        cmap = plt.get_cmap('viridis_r')
        norm = mpl.colors.Normalize(vmin=0, vmax=vmax)
        ax.scatter(output.loc[finite, 'centroid_x'],
                   output.loc[finite, 'centroid_y'],
                   c=finite_dist, cmap=cmap, norm=norm,
                   s=60, edgecolors='black', linewidths=0.5, zorder=3)
        plt.colorbar(ScalarMappable(norm=norm, cmap=cmap), ax=ax,
                     shrink=0.5, label='min distance to GC (um)')
    inf_mask = ~finite
    if inf_mask.any():
        ax.scatter(output.loc[inf_mask, 'centroid_x'],
                   output.loc[inf_mask, 'centroid_y'],
                   c='lightgray', s=60, edgecolors='black', linewidths=0.5,
                   zorder=3)
plt.suptitle(f'{sample_id} ({subtype}) — intra-TLS HEV vs GC distance '
             f'(GC threshold={gc_threshold}, grid={int(grid_um)} um)',
             y=1.05, fontweight='bold')


# --- Page 2: binary close/far classification

fig, ax = plt.subplots(figsize=(9, 8))
setup_spatial_ax(ax, f'{sample_id} ({subtype})\n'
                 f'Intra-TLS HEV niches: close (<={int(near_gc_um)} um) vs far from GC')
if len(output) > 0:
    finite = np.isfinite(output['min_distance_to_GC'].values)
    for is_close, label in [(True, 'close to GC'), (False, 'far from GC')]:
        mask = finite & (output['close_to_GC'].values == is_close)
        if mask.any():
            ax.scatter(output.loc[mask, 'centroid_x'],
                       output.loc[mask, 'centroid_y'],
                       c=BINARY_COLORS[label], s=60, edgecolors='black',
                       linewidths=0.5, zorder=3, label=label)
    inf_mask = ~finite
    if inf_mask.any():
        ax.scatter(output.loc[inf_mask, 'centroid_x'],
                   output.loc[inf_mask, 'centroid_y'],
                   c='lightgray', s=60, edgecolors='black', linewidths=0.5,
                   zorder=3, label='no GC in sample')
    handles = [mlines.Line2D([], [], color=BINARY_COLORS[l], marker='o',
                              markeredgecolor='black', linestyle='', label=l,
                              markersize=8)
               for l in ['close to GC', 'far from GC']]
    ax.legend(handles=handles, loc='center left', bbox_to_anchor=(1.0, 0.5),
              frameon=False, fontsize=10)
plt.suptitle(f'{sample_id} ({subtype}) — intra-TLS HEV binary classification '
             f'(threshold={int(near_gc_um)} um)',
             y=1.05, fontweight='bold')


# --- Save

out_base = os.path.join(res_dir, '05_xenium_mechanisms', '10_HEV_GC_proximity',
                        '01_classify_HEV_by_GC_dist')
csv_dir   = os.path.join(out_base, 'csv_files')
plots_dir = os.path.join(out_base, 'plots')
os.makedirs(csv_dir, exist_ok=True)
os.makedirs(plots_dir, exist_ok=True)

csv_path = os.path.join(csv_dir, f'{sample_id}_intra_TLS_HEV_GC_dist.csv')
output.to_csv(csv_path)

pdf_path = os.path.join(plots_dir, f'{sample_id}_intra_TLS_HEV_GC_dist.pdf')
save_multi_image(pdf_path)

print(f'\nSaved: {csv_path}')
print(f'Saved: {pdf_path}')
print("Script successfully completed")
