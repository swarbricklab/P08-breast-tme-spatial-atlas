#!/usr/bin/env python3

# --- Aim
# Figure S9H: spatial chemokine gradients moving out from the HEVs. For every
# sample we extract density-vs-distance profiles and aggregate them into the
# cohort figure.
#
# The headline S9H figure shows, as a broken-axis profile (mean +/- s.e.m., n=6):
#   LEFT  = distance from the HEV centre within the T cell zone
#   RIGHT = distance from the T-cell-zone edge into the surrounding non-TLS tissue
# for BOTH
#   - CCL19 / CCL21 transcript density (transcripts / mm^2), and
#   - imPVL s1 (PVL Immature s1) / iCAF s2 (CAFs MSC iCAF-like s2) cell density
#     (cells / mm^2).
#
# Per-sample extraction records, for every point (transcript / cell / area-grid
# point), its distance to three anchors:
#   'hev'         : nearest cell of a T-zone-located HEV niche (vessel boundary).
#   'hev_center'  : nearest T-zone-located HEV-niche centroid (radial from core).
#   'tzone_edge'  : nearest T-cell-zone tissue (for the non-TLS segment).
# and tags it with the zone group of its 50um bin (T_cell_zone / B_cell_zone /
# mixed / TLS_unclassified / outside_TLS). HEVs are "T-zone-located" if >= 25% of
# a niche's cells fall in T_cell_zone bins. Distances are tallied in fine 5um bins
# and re-binned to 15um display bins (to 150um per segment) at plot time; per-
# sample bins with < 0.02 mm^2 tissue area or < 15 transcripts are dropped before
# averaging.
#
# --replot_only skips the (slow) SpatialData extraction when the per-sample CSVs
# already exist and just regenerates the cohort figure.

# --- Define environment
# env: spatdat

# --- Load packages

import argparse
import os
import sys
import shutil
import tempfile

import numpy as np
import pandas as pd
import anndata as ad
import spatialdata as sd
from scipy.spatial import cKDTree
from scipy import stats

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import ConnectionPatch, Rectangle, Polygon
from matplotlib.ticker import FuncFormatter

mpl.rcParams['pdf.fonttype'] = 42


# --- Argparse arguments

parser = argparse.ArgumentParser(
    description="S9H: CCL19/CCL21 transcript + imPVL s1/iCAF s2 cell density vs distance from HEV",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
)
# --- runner contract (cohort script)
parser.add_argument("-s", "--samples", nargs="*", type=str, required=True, help="List of samples")
parser.add_argument("--h5ad_dir", default=None, help="defaults to <project_dir>/data/anndata_objects")
parser.add_argument("--project_dir", required=True)
parser.add_argument("--res_dir", required=True)
parser.add_argument("--sample_sheet", default=None,
                    help="defaults to <project_dir>/config/sample_sheet_xenium.csv")
parser.add_argument("--replot_only", action="store_true",
                    help="skip per-sample extraction when the per-sample CSVs already exist")

# --- upstream input dirs
parser.add_argument("--spatdata_dir", default=None,
                    help="defaults to <project_dir>/data/spatialdata_objects/xenium")
parser.add_argument("--zone_dir", default=None,
                    help="TLS-zonation csv_files dir (produced by 07_TLS_zonation); "
                         "defaults to <res_dir>/05_xenium_mechanisms/07_TLS_zonation/"
                         "01_classify_zones/csv_files")

# --- profile / display options kept on the CLI (the analysis parameters are
#     fixed to the manuscript Methods values as module constants, defined just
#     below the argparse block)
parser.add_argument("--profile_genes", nargs="*", default=['CCL21', 'CCL19'],
                    help="transcript genes shown on the density profiles; filename encodes them")
parser.add_argument("--band", choices=['sem', 'ci95'], default='sem',
                    help="shaded band across samples: 'sem' = +/-1 s.e.m.; 'ci95' = t-based 95%% CI")

args = vars(parser.parse_args())
sample_list   = args["samples"]
project_dir   = args["project_dir"]
res_dir       = args["res_dir"]
h5ad_dir      = args["h5ad_dir"]     or os.path.join(project_dir, 'data', 'anndata_objects')
spatdata_dir  = args["spatdata_dir"] or os.path.join(project_dir, 'data', 'spatialdata_objects', 'xenium')
zone_dir      = args["zone_dir"]     or os.path.join(res_dir, '05_xenium_mechanisms', '07_TLS_zonation', '01_classify_zones', 'csv_files')
replot_only   = args["replot_only"]

# --- Analysis parameters (fixed to the values used in the manuscript Methods)
grid_um       = 50.0     # 50 um zonation grid
threshold     = 0.3      # B/T zone score threshold (zonation tag)
dist_bin_um   = 5.0      # fine distance-bin width for tallying
max_dist_um   = 400.0    # max distance profiled from an anchor
area_grid_um  = 5.0      # background area-grid spacing
min_frac      = 0.25     # >=25% of an HEV niche's cells in T-zone bins => T-zone HEV
source_max    = 30.0     # max transcript->source-cell distance for CCL attribution
coexpr_radius = 15.0     # co-expression radius for a different center-gene transcript
reassign_cut  = 10.0     # UNASSIGNED transcripts within this inherit nearest cell label
hev_region    = 20.0     # 'HEV' region distance for regional composition
exclude_other = 20.0     # other-HEV-excluded variant distance
center_genes  = ['CCL19', 'CCL21']                # chemokines that get per-source attribution
marker_genes  = ['TSPAN7', 'MEOX2', 'RAMP3', 'PECAM1', 'NOTCH3', 'ACTA2', 'MCAM', 'SPP1']
fine_bin_um   = 15.0     # display bin width (uniform 15 um)
coarse_bin_um = 15.0
fine_until_um = 0.0      # 0 => uniform display bins
seg_max_um    = 150.0    # x-range per broken-axis segment
ctx_max_um    = 400.0    # x-range for context figures
hev_excl_um   = 0.0      # 0 => plot the no-other-HEV-exclusion variant
min_bin_area_um2 = 20000.0   # drop per-sample density bins < 0.02 mm^2 before averaging
min_bin_cells    = 20.0      # drop per-sample proportion bins < 20 cells
min_bin_tx       = 15.0      # drop per-sample transcript-density bins < 15 transcripts
all_genes     = center_genes + [g for g in marker_genes if g not in center_genes]
threshold_tag = str(threshold).replace('.', '')

# --- Local helper functions
sys.path.append(os.path.join(project_dir, 'scripts'))
from helper_functions import save_multi_image   # noqa: E402

# Tangram / label-transfer subset labels
IMPVL, ICAF = 'PVL Immature s1', 'CAFs MSC iCAF-like s2'
LEC = 'Endothelial Lymphatic LYVE1'   # lymphatic EC, a 4th CCL source for the regional composition
CELLTYPES = {
    'imPVL s1':  IMPVL,                          # perivascular / HEV-associated, CCL19/CCL21+
    'iCAF s2':   ICAF,                           # T-cell-zone FRC, CCL19-only
    'venous EC': 'Endothelial ACKR1',            # venous / HEV endothelium (control)
    'LEC':       'Endothelial Lymphatic LYVE1',  # lymphatic endothelium (control)
}
SOURCE_LABELS = ['imPVL s1', 'iCAF s2', 'other']
COEXPR_LABELS = ['co', 'solo']
ANCHORS = ['hev', 'hev_center', 'tzone_edge']
T_ZONE = 'T_cell_zone'

TX_COLS = ['sample_id', 'subtype', 'anchor', 'hev_excl_um', 'dist_bin_center', 'zone_group', 'gene', 'source',
           'n_transcripts', 'area_um2', 'density_per_mm2']
CT_COLS = ['sample_id', 'subtype', 'anchor', 'hev_excl_um', 'dist_bin_center', 'zone_group', 'cell_type',
           'tangram_subset', 'n_cells', 'n_total_cells', 'area_um2', 'density_per_mm2', 'proportion']
COMP_COLS   = ['sample_id', 'subtype', 'region', 'category', 'source4', 'n_transcripts', 'region_area_um2']
TXTAB_COLS  = ['sample_id', 'subtype', 'region', 'feature_name', 'd_partner', 'd_imPVL', 'd_iCAF',
               'd_lymphEC', 'same_cell_other', 'same_cell_own']


# --- Sample metadata

if args.get("sample_sheet"):
    sample_sheet_path = args["sample_sheet"]
else:
    sample_sheet_path = os.path.join(project_dir, 'config', 'sample_sheet_xenium.csv')
sample_sheet    = pd.read_csv(sample_sheet_path)
tissue_metadata = pd.read_csv(os.path.join(project_dir, 'resources', 'tissue_metadata.csv'))
merged          = sample_sheet.merge(tissue_metadata, on='tissue_id')
subtype_dict    = dict(zip(merged['section_id'], merged['clinical_subtype']))

csv_dir = os.path.join(res_dir, '05_xenium_mechanisms', '09_chemokine_sources',
                       '05_chemokine_gradients', 'csv_files')
os.makedirs(csv_dir, exist_ok=True)


def celltype_csv(sample, level='subset'):
    """Label-transfer predicted cell types (col 'pred_cell_type')."""
    p = os.path.join(res_dir, '02_cell_type_mapping', 'label_transfer', 'xenium',
                     f'{level}_level', f'{sample}_{level}_pred_cell_types.csv')
    return p


# --- Per-sample extraction

def extract_one_sample(sample):
    subtype = subtype_dict.get(sample)
    print(f'\n=== extract {sample} ({subtype}) ===')

    def bin_center(d):
        return (np.floor(np.asarray(d, dtype=float) / dist_bin_um) + 0.5) * dist_bin_um

    # --- Zone bins (50um grid, T/B classification from 07_TLS_zonation)
    zone_csv = os.path.join(zone_dir, f'{sample}_BT_zone_bins.csv')
    bin_df = pd.read_csv(zone_csv)
    xmin = float(bin_df['bin_x'].iloc[0] - (bin_df['bin_ix'].iloc[0] + 0.5) * grid_um)
    ymin = float(bin_df['bin_y'].iloc[0] - (bin_df['bin_iy'].iloc[0] + 0.5) * grid_um)
    bin_zone = bin_df.set_index(['bin_ix', 'bin_iy'])['TLS_zone'].fillna('TLS_unclassified')

    def zone_group_of(x, y):
        ix = np.floor((np.asarray(x, dtype=float) - xmin) / grid_um).astype(int)
        iy = np.floor((np.asarray(y, dtype=float) - ymin) / grid_um).astype(int)
        g = bin_zone.reindex(pd.MultiIndex.from_arrays([ix, iy])).to_numpy()
        return pd.Series(g).fillna('outside_TLS').to_numpy()

    n_tzone_bins = int((bin_df['TLS_zone'] == T_ZONE).sum())
    print(f'T-cell-zone bins: {n_tzone_bins} ({n_tzone_bins * grid_um**2 / 1e6:.3f} mm^2)')

    # --- Cells: centroids, subset label, HEV niche membership
    adata   = ad.read_h5ad(os.path.join(h5ad_dir, 'xenium', f'{sample}.h5ad'), backed='r')
    # micron centroids (obsm['spatial'] is image space in this project; the zone bins
    # from 07_TLS_zonation are built on these native micron coordinates, and the
    # SpatialData transcripts are in the same micron space)
    cell_xy = pd.DataFrame(adata.obs[['x_centroid', 'y_centroid']].to_numpy(), columns=['x', 'y'],
                           index=adata.obs_names.astype(str))
    adata.file.close()

    subset_df = pd.read_csv(celltype_csv(sample), index_col=0)
    subset_df.index = subset_df.index.astype(str)

    hev_root = os.path.join(res_dir, '05_xenium_mechanisms', '02_identify_HEV_regions')
    hev_recent = hev_root
    niche_df = pd.read_csv(os.path.join(hev_recent, 'csv_files',
                                        f'{sample}_HEV_niches_cell_to_niche.csv'), index_col=0)
    niche_df.index = niche_df.index.astype(str)
    print(f'Using HEV niches from {os.path.basename(hev_recent)}')

    cells = cell_xy.join(subset_df['pred_cell_type'].rename('subset')).join(niche_df['niche_id'].rename('niche_id'))

    # --- T-zone-located HEV niches -> anchors (boundary cells + niche centroids)
    hev_cells = cells[cells['niche_id'].notna() & (cells['niche_id'].astype(str) != 'nan')].copy()
    hev_cells['zone_group'] = zone_group_of(hev_cells['x'], hev_cells['y'])
    frac_tzone = hev_cells.groupby('niche_id')['zone_group'].apply(lambda z: (z == T_ZONE).mean())
    tzone_niches = frac_tzone[frac_tzone >= min_frac].index.tolist()
    print(f'HEV niches total={frac_tzone.size}, T-zone-located (>= {min_frac:.0%})={len(tzone_niches)}')

    anchor_cells = hev_cells[hev_cells['niche_id'].isin(tzone_niches)]
    anchor = anchor_cells[['x', 'y']].to_numpy()
    _centroid_df = anchor_cells.groupby('niche_id')[['x', 'y']].mean()
    centroids = _centroid_df.to_numpy()
    centroid_niche_ids = _centroid_df.index.to_numpy()

    _cen = anchor_cells.groupby('niche_id')[['x', 'y']].transform('mean')
    _r   = np.sqrt((anchor_cells['x'] - _cen['x'])**2 + (anchor_cells['y'] - _cen['y'])**2)
    _sz  = anchor_cells.groupby('niche_id').size()
    _per_max = _r.groupby(anchor_cells['niche_id']).max()
    _per_max = _per_max[_per_max.index.isin(_sz[_sz >= 3].index)]
    boundary_radius_um = float(np.median(_per_max)) if len(_per_max) else float('nan')
    niche_radius_median_um = (float(_r.groupby(anchor_cells['niche_id']).median().median())
                              if len(anchor_cells) else float('nan'))

    meta = {
        'sample_id': sample, 'subtype': subtype,
        'grid_um': grid_um, 'threshold': threshold, 'dist_bin_um': dist_bin_um,
        'max_dist_um': max_dist_um, 'area_grid_um': area_grid_um,
        'min_niche_tzone_frac': min_frac, 'source_max_dist': source_max,
        'coexpr_radius_um': coexpr_radius, 'reassign_cutoff': reassign_cut,
        'hev_region_um': hev_region, 'exclude_other_hev_um': exclude_other,
        'n_tzone_bins': n_tzone_bins, 'tzone_area_mm2': n_tzone_bins * grid_um**2 / 1e6,
        'n_hev_niches_total': int(frac_tzone.size), 'n_hev_niches_tzone': len(tzone_niches),
        'n_anchor_cells': int(len(anchor)), 'hev_niche_date': os.path.basename(hev_recent),
        'boundary_radius_um': boundary_radius_um, 'niche_radius_median_um': niche_radius_median_um,
        'center_genes': ','.join(center_genes), 'marker_genes': ','.join(marker_genes),
    }

    if len(anchor) == 0:
        print('No T-zone-located HEV niches; writing empty outputs.')
        pd.DataFrame(columns=TX_COLS).to_csv(os.path.join(csv_dir, f'{sample}_transcript_density_vs_dist.csv'), index=False)
        pd.DataFrame(columns=CT_COLS).to_csv(os.path.join(csv_dir, f'{sample}_celltype_density_vs_dist.csv'), index=False)
        pd.DataFrame(columns=COMP_COLS).to_csv(os.path.join(csv_dir, f'{sample}_coexpr_source_composition.csv'), index=False)
        pd.DataFrame(columns=TXTAB_COLS).to_csv(os.path.join(csv_dir, f'{sample}_coexpr_source_transcripts.csv'), index=False)
        pd.DataFrame([meta]).to_csv(os.path.join(csv_dir, f'{sample}_distance_profile_meta.csv'), index=False)
        return

    hev_tree  = cKDTree(anchor)
    hevc_tree = cKDTree(centroids)
    allhev_tree  = cKDTree(hev_cells[['x', 'y']].to_numpy())
    allhev_niche = hev_cells['niche_id'].to_numpy()

    def near_other_hev(xy):
        if exclude_other <= 0 or len(xy) == 0:
            return np.zeros(len(xy), dtype=bool)
        assigned   = centroid_niche_ids[hevc_tree.query(xy, k=1)[1]]
        d_near, i_near = allhev_tree.query(xy, k=1)
        return (d_near < exclude_other) & (allhev_niche[i_near] != assigned)

    def cell_tree(label):
        pts = cells.loc[cells['subset'] == label, ['x', 'y']].dropna().to_numpy()
        return (cKDTree(pts), len(pts)) if len(pts) else (None, 0)

    imp_tree, n_imp = cell_tree(IMPVL)
    ica_tree, n_ica = cell_tree(ICAF)
    print(f'source cells: imPVL s1={n_imp:,}, iCAF s2={n_ica:,}')

    def attribute_source(xy):
        n = len(xy)
        d_imp = imp_tree.query(xy, k=1)[0] if imp_tree is not None else np.full(n, np.inf)
        d_ica = ica_tree.query(xy, k=1)[0] if ica_tree is not None else np.full(n, np.inf)
        src = np.full(n, 'other', dtype=object)
        imp_win = (d_imp <= d_ica) & (d_imp <= source_max)
        ica_win = (d_ica <  d_imp) & (d_ica <= source_max)
        src[imp_win] = 'imPVL s1'
        src[ica_win] = 'iCAF s2'
        return src

    lec_tree, n_lec = cell_tree(LEC)
    print(f'source cells (composition): +lymphatic EC={n_lec:,}')

    def attribute_source4(xy):
        n = len(xy)
        trees = {'imPVL s1': imp_tree, 'iCAF s2': ica_tree, 'lymphatic EC': lec_tree}
        labs = list(trees)
        D = np.vstack([t.query(xy, k=1)[0] if t is not None else np.full(n, np.inf) for t in trees.values()])
        dmin = D.min(axis=0)
        src = np.array(labs, dtype=object)[D.argmin(axis=0)]
        src[dmin > source_max] = 'other'
        return src

    # --- Background grid over the anchor bbox +/- max_dist; distances to each anchor
    pts_bbox = np.vstack([anchor, centroids])
    ax0, ay0 = pts_bbox.min(axis=0)
    ax1, ay1 = pts_bbox.max(axis=0)
    gx = np.arange(ax0 - max_dist_um, ax1 + max_dist_um, area_grid_um)
    gy = np.arange(ay0 - max_dist_um, ay1 + max_dist_um, area_grid_um)
    GX, GY = np.meshgrid(gx, gy)
    grid_pts  = np.column_stack([GX.ravel(), GY.ravel()])
    grid_zone = zone_group_of(grid_pts[:, 0], grid_pts[:, 1])
    tzone_tree = cKDTree(grid_pts[grid_zone == T_ZONE])

    def dists(xy):
        return {'hev': hev_tree.query(xy, k=1)[0],
                'hev_center': hevc_tree.query(xy, k=1)[0],
                'tzone_edge': tzone_tree.query(xy, k=1)[0]}

    grid = pd.DataFrame({'zone_group': grid_zone, **dists(grid_pts)})
    grid['near_other'] = near_other_hev(grid_pts)
    print(f'Background grid: {len(grid):,} points ({area_grid_um:.0f}um spacing); '
          f'{grid["near_other"].mean():.1%} within {exclude_other:.0f}um of another HEV')

    def assign_region(hev_d, zone):
        hev_d = np.asarray(hev_d); zone = np.asarray(zone)
        return np.where(hev_d <= hev_region, 'HEV',
                        np.where((zone == T_ZONE) & (hev_d > hev_region), 'T_cell_zone', 'none'))

    _grid_region = assign_region(grid['hev'], grid['zone_group'])
    region_area = (pd.Series(_grid_region).value_counts() * area_grid_um**2)
    meta['hev_region_area_um2']   = float(region_area.get('HEV', 0.0))
    meta['tzone_region_area_um2'] = float(region_area.get('T_cell_zone', 0.0))

    # --- Cells: distances + zone group
    cells_all = cells.dropna(subset=['x', 'y']).copy()
    cxy = cells_all[['x', 'y']].to_numpy()
    cells_all['zone_group'] = zone_group_of(cells_all['x'], cells_all['y'])
    for a, d in dists(cxy).items():
        cells_all[a] = d
    cells_all['near_other'] = near_other_hev(cxy)

    def build_tx(dist_col, gdf, tdf):
        g = gdf[gdf[dist_col] < max_dist_um].copy()
        g['c'] = bin_center(g[dist_col])
        area = (g.groupby(['c', 'zone_group']).size() * area_grid_um**2).rename('area_um2').reset_index()
        t = tdf[tdf[dist_col] < max_dist_um].copy()
        t['c'] = bin_center(t[dist_col])
        rows = []
        for gene in all_genes:
            sub = t[t['feature_name'] == gene]
            specs = [('all', None)]
            if gene in center_genes:
                specs += [(s, 'source') for s in SOURCE_LABELS]
                specs += [(s, 'coexpr') for s in COEXPR_LABELS]
            for src, col in specs:
                s = sub if col is None else sub[sub[col] == src]
                n = s.groupby(['c', 'zone_group']).size().rename('n_transcripts').reset_index()
                m = area.merge(n, on=['c', 'zone_group'], how='left')
                m['n_transcripts'] = m['n_transcripts'].fillna(0).astype(int)
                m['gene'] = gene
                m['source'] = src
                rows.append(m)
        tx = pd.concat(rows, ignore_index=True)
        tx['density_per_mm2'] = tx['n_transcripts'] / tx['area_um2'] * 1e6
        return tx.rename(columns={'c': 'dist_bin_center'})

    def build_ct(dist_col, gdf, cdf):
        g = gdf[gdf[dist_col] < max_dist_um].copy()
        g['c'] = bin_center(g[dist_col])
        area = (g.groupby(['c', 'zone_group']).size() * area_grid_um**2).rename('area_um2').reset_index()
        c = cdf[cdf[dist_col] < max_dist_um].copy()
        c['c'] = bin_center(c[dist_col])
        total = c.groupby(['c', 'zone_group']).size().rename('n_total_cells').reset_index()
        rows = []
        for short, label in CELLTYPES.items():
            n = c[c['subset'] == label].groupby(['c', 'zone_group']).size().rename('n_cells').reset_index()
            m = area.merge(n, on=['c', 'zone_group'], how='left').merge(total, on=['c', 'zone_group'], how='left')
            m['n_cells']       = m['n_cells'].fillna(0).astype(int)
            m['n_total_cells'] = m['n_total_cells'].fillna(0).astype(int)
            m['cell_type']      = short
            m['tangram_subset'] = label
            rows.append(m)
        ct = pd.concat(rows, ignore_index=True)
        ct['density_per_mm2'] = ct['n_cells'] / ct['area_um2'] * 1e6
        ct['proportion'] = np.where(ct['n_total_cells'] > 0, ct['n_cells'] / ct['n_total_cells'], np.nan)
        return ct.rename(columns={'c': 'dist_bin_center'})

    # --- Load spatialdata transcripts. Use a persistent .zarr if present;
    #     otherwise unzip the .zip to a scratch tmp dir and clean up afterwards.
    zip_path     = os.path.join(spatdata_dir, f'{sample}.zip')
    zarr_persist = os.path.join(spatdata_dir, f'{sample}.zarr')
    tmp_root = None
    try:
        if os.path.isdir(zarr_persist):
            zarr_path = zarr_persist
            print(f'Reading existing zarr {zarr_path}')
        else:
            tmp_root = tempfile.mkdtemp(prefix=f'chemokine_grad_{sample}_')
            zarr_path = os.path.join(tmp_root, f'{sample}.zarr')
            print(f'Unzipping {zip_path} -> {zarr_path}')
            shutil.unpack_archive(zip_path, zarr_path)
        sdata = sd.read_zarr(zarr_path)
        tx    = sdata.points['transcripts']
        center_tx = (tx[tx['feature_name'].isin(all_genes)][['x', 'y', 'feature_name', 'cell_id']]
                     .compute().reset_index(drop=True))
    finally:
        if tmp_root is not None:
            shutil.rmtree(tmp_root, ignore_errors=True)
    center_tx['feature_name'] = center_tx['feature_name'].astype(str)
    center_tx['cell_id'] = center_tx['cell_id'].astype(str)
    print('Transcripts loaded: ' +
          ', '.join(f'{g}={int((center_tx["feature_name"] == g).sum()):,}' for g in all_genes))

    txy = center_tx[['x', 'y']].to_numpy()
    center_tx['zone_group'] = zone_group_of(center_tx['x'].to_numpy(), center_tx['y'].to_numpy())
    for a, d in dists(txy).items():
        center_tx[a] = d
    is_ccl = center_tx['feature_name'].isin(center_genes).to_numpy()
    center_tx['source'] = 'all'
    center_tx.loc[is_ccl, 'source'] = attribute_source(txy[is_ccl])

    center_tx['coexpr'] = 'all'
    center_tx['d_partner'] = np.inf
    ccl_trees = {g: cKDTree(txy[(center_tx['feature_name'] == g).to_numpy()])
                 for g in center_genes if (center_tx['feature_name'] == g).any()}
    for g in center_genes:
        m = (center_tx['feature_name'] == g).to_numpy()
        if not m.any():
            continue
        d_partner = np.full(int(m.sum()), np.inf)
        for o in center_genes:
            if o != g and o in ccl_trees:
                d_partner = np.minimum(d_partner, ccl_trees[o].query(txy[m], k=1)[0])
        center_tx.loc[m, 'd_partner'] = d_partner
        center_tx.loc[m, 'coexpr'] = np.where(d_partner <= coexpr_radius, 'co', 'solo')

    center_tx['same_cell_other'] = 0
    center_tx['same_cell_own']   = 0
    if len(center_genes) == 2:
        g0, g1 = center_genes
        other_of = {g0: g1, g1: g0}
        cid = center_tx['cell_id'].to_numpy().astype(object)
        eff_cid = cid.copy()
        un = is_ccl & (cid == 'UNASSIGNED')
        if un.any():
            cell_valid = cells.dropna(subset=['x', 'y'])
            ctree = cKDTree(cell_valid[['x', 'y']].to_numpy())
            cbar  = cell_valid.index.to_numpy().astype(str)
            dd, ii = ctree.query(txy[un], k=1)
            near = dd <= reassign_cut
            eff_cid[np.where(un)[0][near]] = cbar[ii[near]]
        assigned = center_tx.loc[is_ccl & (cid != 'UNASSIGNED'), ['feature_name', 'cell_id']]
        per_cell = assigned.groupby(['cell_id', 'feature_name']).size().unstack(fill_value=0)
        for g in center_genes:
            if g not in per_cell.columns:
                per_cell[g] = 0
        sco  = np.zeros(len(center_tx))
        sown = np.zeros(len(center_tx))
        for g in center_genes:
            m = is_ccl & (center_tx['feature_name'].to_numpy() == g)
            sco[m]  = per_cell[other_of[g]].reindex(eff_cid[m]).fillna(0).to_numpy()
            sown[m] = per_cell[g].reindex(eff_cid[m]).fillna(0).to_numpy()
        center_tx['same_cell_other'] = sco.astype(int)
        center_tx['same_cell_own']   = sown.astype(int)

    center_tx['near_other'] = near_other_hev(txy)

    variants = [(0.0, grid, center_tx, cells_all)]
    if exclude_other > 0:
        variants.append((exclude_other,
                         grid[~grid['near_other']], center_tx[~center_tx['near_other']],
                         cells_all[~cells_all['near_other']]))

    tx_frames, ct_frames = [], []
    for excl, gdf, tdf, cdf in variants:
        for a in ANCHORS:
            tx_a = build_tx(a, gdf, tdf); tx_a['anchor'] = a; tx_a['hev_excl_um'] = excl; tx_frames.append(tx_a)
            ct_a = build_ct(a, gdf, cdf); ct_a['anchor'] = a; ct_a['hev_excl_um'] = excl; ct_frames.append(ct_a)
    tx_df = pd.concat(tx_frames, ignore_index=True)
    ct_df = pd.concat(ct_frames, ignore_index=True)

    # --- Regional composition of CCL transcripts (feeds the retired HEV-vs-T-zone barplots)
    ccl = center_tx.loc[is_ccl, ['feature_name', 'coexpr', 'hev', 'zone_group', 'd_partner',
                                 'same_cell_other', 'same_cell_own']].copy()
    xy_ccl = txy[is_ccl]
    ccl['d_imPVL']   = imp_tree.query(xy_ccl, k=1)[0] if imp_tree is not None else np.inf
    ccl['d_iCAF']    = ica_tree.query(xy_ccl, k=1)[0] if ica_tree is not None else np.inf
    ccl['d_lymphEC'] = lec_tree.query(xy_ccl, k=1)[0] if lec_tree is not None else np.inf
    ccl['region'] = assign_region(ccl['hev'], ccl['zone_group'])
    ccl['source4'] = attribute_source4(xy_ccl)
    ccl['category'] = np.where(ccl['coexpr'].to_numpy() == 'co', 'co-expressed',
                               ccl['feature_name'].astype(str) + ' mono')

    in_region = ccl[ccl['region'] != 'none'].copy()
    tx_table = in_region[['region', 'feature_name', 'd_partner', 'd_imPVL', 'd_iCAF',
                          'd_lymphEC', 'same_cell_other', 'same_cell_own']].copy()
    comp = (in_region.groupby(['region', 'category', 'source4']).size()
            .rename('n_transcripts').reset_index())
    comp['region_area_um2'] = comp['region'].map(region_area).astype(float)

    # --- Save per-sample outputs
    tx_df['sample_id'] = sample; tx_df['subtype'] = subtype
    ct_df['sample_id'] = sample; ct_df['subtype'] = subtype
    comp['sample_id'] = sample; comp['subtype'] = subtype
    tx_table['sample_id'] = sample; tx_table['subtype'] = subtype

    tx_df[TX_COLS].to_csv(os.path.join(csv_dir, f'{sample}_transcript_density_vs_dist.csv'), index=False)
    ct_df[CT_COLS].to_csv(os.path.join(csv_dir, f'{sample}_celltype_density_vs_dist.csv'), index=False)
    comp[COMP_COLS].to_csv(os.path.join(csv_dir, f'{sample}_coexpr_source_composition.csv'), index=False)
    tx_table[TXTAB_COLS].to_csv(os.path.join(csv_dir, f'{sample}_coexpr_source_transcripts.csv'), index=False)
    pd.DataFrame([meta]).to_csv(os.path.join(csv_dir, f'{sample}_distance_profile_meta.csv'), index=False)
    print(f'Saved distance profiles for {sample}')


# --- Run extraction (loop over samples), skipped by --replot_only if CSVs exist

for s in sample_list:
    have = os.path.exists(os.path.join(csv_dir, f'{s}_transcript_density_vs_dist.csv'))
    if replot_only and have:
        print(f'--replot_only: reusing existing per-sample CSVs for {s}')
        continue
    extract_one_sample(s)


# =====================================================================================
# --- Aggregate the per-sample profiles + plot the cohort S9H figure
# =====================================================================================

# (display/aggregation parameters are set as module constants above)
band_kind        = args["band"]


def load_concat(suffix):
    frames = []
    for s in sample_list:
        f = os.path.join(csv_dir, f'{s}_{suffix}.csv')
        if os.path.exists(f):
            df = pd.read_csv(f)
            if len(df):
                frames.append(df)
        else:
            print(f'WARNING: missing {f}')
    if not frames:
        raise FileNotFoundError(f'No {suffix} CSVs found for the requested samples')
    return pd.concat(frames, ignore_index=True)

tx_df   = load_concat('transcript_density_vs_dist')
ct_df   = load_concat('celltype_density_vs_dist')
meta_df = load_concat('distance_profile_meta')

if 'hev_excl_um' in tx_df.columns:
    avail = sorted(tx_df['hev_excl_um'].unique())
    if hev_excl_um not in avail:
        print(f'WARNING: hev_excl_um={hev_excl_um} not in {avail}; using {avail[0]}')
        hev_excl_um = avail[0]
    tx_df = tx_df[tx_df['hev_excl_um'] == hev_excl_um].copy()
    ct_df = ct_df[ct_df['hev_excl_um'] == hev_excl_um].copy()
elif hev_excl_um:
    print('WARNING: CSVs lack hev_excl_um; plotting unfiltered data')
    hev_excl_um = 0.0

n_samples = tx_df['sample_id'].nunique()
hev_frac  = float(meta_df['min_niche_tzone_frac'].iloc[0])
grid_res  = float(meta_df['area_grid_um'].iloc[0])
hev_boundary_um = (float(meta_df['boundary_radius_um'].mean())
                   if 'boundary_radius_um' in meta_df.columns else float('nan'))
print(f'Loaded profiles for {n_samples} samples; avg HEV boundary radius = {hev_boundary_um:.1f}um')


# --- Colors

_plot_col_dir = os.path.join(project_dir, 'resources', 'plotting_colors')
_tx_col = pd.read_csv(os.path.join(_plot_col_dir, 'transcript_colors.csv')).set_index('transcript')['color'].to_dict()
CCL_COLORS = {'CCL21': _tx_col.get('CCL21', '#1f77b4'), 'CCL19': _tx_col.get('CCL19', '#d62728')}
PROFILE_GENES = list(args["profile_genes"])
GENE_COLORS = {g: _tx_col.get(g, '#333333') for g in PROFILE_GENES}
GENE_TAG = '_'.join(PROFILE_GENES)

_subset_col = pd.read_csv(os.path.join(_plot_col_dir, 'miniatlas_subset_colors.csv')).set_index('miniatlas_subset')['color'].to_dict()
_short2label = ct_df.drop_duplicates('cell_type').set_index('cell_type')['tangram_subset'].to_dict()
_ct_fallback = {'imPVL s1': '#e377c2', 'iCAF s2': '#2ca02c', 'venous EC': '#8c564b', 'LEC': '#7f7f7f'}
CT_COLORS = {s: _subset_col.get(_short2label.get(s), fb) for s, fb in _ct_fallback.items()}
ZONE_LABEL = {'T_cell_zone': 'T cell zone', 'B_cell_zone': 'B-cell zone', 'mixed': 'mixed',
              'TLS_unclassified': 'TLS (unclassified)', 'outside_TLS': 'outside TLS'}
ZONE_COLORS = {'T_cell_zone': '#4c72b0', 'B_cell_zone': '#dd8452', 'mixed': '#8172b3',
               'TLS_unclassified': '#c7c7c7', 'outside_TLS': '#333333'}
HEV_ANCHOR = 'hev_center'

tx_all = tx_df[tx_df['source'] == 'all'].copy()

SOURCE_COLORS = {'imPVL s1': CT_COLORS['imPVL s1'], 'iCAF s2': CT_COLORS['iCAF s2'], 'other': '#b0b0b0'}
MARKER_CLASS  = {'EC (TSPAN7/MEOX2/RAMP3/PECAM1)': ['TSPAN7', 'MEOX2', 'RAMP3', 'PECAM1'],
                 'PVL (NOTCH3/ACTA2/MCAM)':        ['NOTCH3', 'ACTA2', 'MCAM']}
MARKER_COLORS = {'EC (TSPAN7/MEOX2/RAMP3/PECAM1)': '#c1272d', 'PVL (NOTCH3/ACTA2/MCAM)': '#2b6cb0'}
src_max = float(meta_df['source_max_dist'].iloc[0])


# --- Display bins (uniform by default; fine-near-anchor optional)

def make_edges(fine_bin, fine_until, coarse_bin, max_c):
    if fine_until <= 0 or abs(fine_bin - coarse_bin) < 1e-9:
        edges = np.arange(0, max_c + coarse_bin, coarse_bin)
    else:
        e1 = np.arange(0, fine_until + fine_bin / 2, fine_bin)
        e2 = np.arange(fine_until + coarse_bin, max_c + coarse_bin, coarse_bin)
        edges = np.unique(np.concatenate([e1, e2]))
    return edges[edges <= max_c + coarse_bin / 2]

EDGES   = make_edges(fine_bin_um, fine_until_um, coarse_bin_um, ctx_max_um)
CENTERS = (EDGES[:-1] + EDGES[1:]) / 2.0
EDGES_SEG2 = make_edges(coarse_bin_um, 0, coarse_bin_um, ctx_max_um)

if fine_until_um <= 0 or abs(fine_bin_um - coarse_bin_um) < 1e-9:
    BIN_STR = f'{coarse_bin_um:.0f}um bins'
    BIN_TAG = f'bin_{coarse_bin_um:.0f}um'
else:
    BIN_STR = f'{fine_bin_um:.0f}um<={fine_until_um:.0f}um, {coarse_bin_um:.0f}um>{fine_until_um:.0f}um'
    BIN_TAG = f'bin_{fine_bin_um:.0f}-{coarse_bin_um:.0f}um_fine{fine_until_um:.0f}'
EXCL_STR = f'excl. other-HEV <{hev_excl_um:.0f}um' if hev_excl_um else 'no HEV exclusion'
BIN_TAG += f'_excl{hev_excl_um:.0f}'
BAND_STR = 'line = mean, shaded = 95% CI' if band_kind == 'ci95' else 'line = mean, shaded = +/-1 s.e.m.'
INFO = (f'T-zone HEVs >={hev_frac:.0%} - {BIN_STR} - {grid_res:.0f}um area grid - '
        f'{EXCL_STR} - bins >={min_bin_area_um2/1e6:.2f}mm2 & >={min_bin_tx:.0f} tx - '
        f'n={n_samples} - {BAND_STR}')


# --- Per-sample re-binning, then mean +/- SEM across samples

def rebin(df, value_kind, group_col, by_zone, edges=None):
    edges = EDGES if edges is None else edges
    centers = (edges[:-1] + edges[1:]) / 2.0
    d = df.copy()
    idx = np.digitize(d['dist_bin_center'].to_numpy(), edges) - 1
    keep = (idx >= 0) & (idx < len(centers))
    d = d[keep].copy()
    d['disp'] = centers[idx[keep]]
    keys = ['sample_id', group_col, 'disp'] + (['zone_group'] if by_zone else [])
    if value_kind in ('transcript', 'celldensity'):
        ncol = 'n_transcripts' if value_kind == 'transcript' else 'n_cells'
        g = d.groupby(keys).agg(n=(ncol, 'sum'), area=('area_um2', 'sum')).reset_index()
        g['val'] = np.where(g['area'] > 0, g['n'] / g['area'] * 1e6, np.nan)
        g.loc[g['area'] < min_bin_area_um2, 'val'] = np.nan
        if value_kind == 'transcript' and min_bin_tx > 0:
            g.loc[g['n'] < min_bin_tx, 'val'] = np.nan
    else:
        g = d.groupby(keys).agg(n=('n_cells', 'sum'), tot=('n_total_cells', 'sum')).reset_index()
        g['val'] = np.where(g['tot'] > 0, g['n'] / g['tot'], np.nan)
        g.loc[g['tot'] < min_bin_cells, 'val'] = np.nan
    return g

def across(per_sample, group_col):
    keys = [group_col, 'disp'] + (['zone_group'] if 'zone_group' in per_sample.columns else [])
    return per_sample.groupby(keys)['val'].agg(mean='mean', sem='sem', n='count').reset_index()


def band_half(sem, n):
    sem = np.asarray(sem, dtype=float)
    n   = np.asarray(n, dtype=float)
    if band_kind == 'ci95':
        with np.errstate(invalid='ignore'):
            mult = np.where(n >= 2, stats.t.ppf(0.975, np.maximum(n - 1, 1)), 0.0)
    else:
        mult = 1.0
    return np.nan_to_num(mult * sem, nan=0.0)


def draw(ax, agg, group_col, series, zone=None, xmax=None):
    for s in series:
        sub = agg[agg[group_col] == s['group']]
        if zone is not None and 'zone_group' in sub.columns:
            sub = sub[sub['zone_group'] == zone]
        sub = sub[sub['n'] >= 1].sort_values('disp')
        if xmax is not None:
            sub = sub[sub['disp'] <= xmax]
        if sub.empty:
            continue
        x = sub['disp'].to_numpy(); y = sub['mean'].to_numpy()
        e = band_half(sub['sem'].to_numpy(), sub['n'].to_numpy())
        c = s['color']; a = s.get('alpha', 1.0)
        ax.plot(x, y, color=c, ls=s.get('ls', '-'), lw=2, label=s['label'], alpha=a, zorder=3)
        ax.fill_between(x, y - e, y + e, color=c, alpha=0.15 * a, lw=0, zorder=2)


def _endpoint(agg, group_col, grp, zone, which):
    sub = agg[agg[group_col] == grp]
    if zone is not None and 'zone_group' in sub.columns:
        sub = sub[sub['zone_group'] == zone]
    sub = sub[(sub['n'] >= 1) & (sub['disp'] <= seg_max_um)].sort_values('disp')
    if sub.empty:
        return None
    r = sub.iloc[-1] if which == 'last' else sub.iloc[0]
    return (float(r['disp']), float(r['mean']))


# --- S9H inputs (top row): CCL19 / CCL21 transcript density, HEV -> T cell zone -> non-TLS

tx_hev   = across(rebin(tx_all[tx_all['anchor'] == HEV_ANCHOR],        'transcript', 'gene', by_zone=True),  'gene')
seg2_src = tx_all[(tx_all['anchor'] == 'tzone_edge') & (tx_all['zone_group'] == 'outside_TLS')]
tx_tzone = across(rebin(seg2_src, 'transcript', 'gene', by_zone=False, edges=EDGES_SEG2), 'gene')


# --- S9H inputs (bottom row): imPVL s1 / iCAF s2 cell density, HEV -> T cell zone -> non-TLS

ct_hev   = across(rebin(ct_df[ct_df['anchor'] == HEV_ANCHOR], 'celldensity', 'cell_type', by_zone=True), 'cell_type')
ct2_src  = ct_df[(ct_df['anchor'] == 'tzone_edge') & (ct_df['zone_group'] == 'outside_TLS')]
ct_tzone = across(rebin(ct2_src, 'celldensity', 'cell_type', by_zone=False, edges=EDGES_SEG2), 'cell_type')

# --- Output dir

out_dir = os.path.join(res_dir, '05_xenium_mechanisms', '09_chemokine_sources',
                       '05_chemokine_gradients', 'plots')
os.makedirs(out_dir, exist_ok=True)


# --- Combined stacked S9H figure (CCL19/CCL21 transcripts over imPVL s1/iCAF s2 density)

_zone_col_f = os.path.join(_plot_col_dir, 'TLS_zone_colors.csv')
_zone_col = (pd.read_csv(_zone_col_f).set_index('zone')['color'].to_dict()
             if os.path.exists(_zone_col_f) else {})
REGION_COLORS = {'HEV': _zone_col.get('HEV', '#604a8b'),
                 'T cell zone': _zone_col.get('T_cell_zone', mpl.colors.cnames['crimson']),
                 'non-TLS': _zone_col.get('non_TLS', '#9e9e9e')}
COMBINED_ROWS = [
    ('gene', tx_hev, tx_tzone,
     [{'group': g, 'color': CCL_COLORS[g], 'ls': '-', 'label': g} for g in ['CCL21', 'CCL19']],
     'Transcript density\n(tx / mm2)', True),
    ('cell_type', ct_hev, ct_tzone,
     [{'group': c, 'color': CT_COLORS[c], 'ls': '-', 'label': c} for c in ['imPVL s1', 'iCAF s2']],
     'Cell density\n(cells / mm2)', False),
]
figC, axsC = plt.subplots(2, 2, figsize=(4.2, 2.7), sharex='col', sharey='row',
                          gridspec_kw={'wspace': 0.03, 'hspace': 0.08})
for ri, (gcol, seg1, seg2, series, ylab, sci_y) in enumerate(COMBINED_ROWS):
    axL, axR = axsC[ri, 0], axsC[ri, 1]
    draw(axL, seg1, gcol, series, zone=T_ZONE, xmax=seg_max_um)
    draw(axR, seg2, gcol, series, zone=None,   xmax=seg_max_um)
    axL.set_xlim(0, seg_max_um); axR.set_xlim(0, seg_max_um)
    axL.set_ylabel(ylab, fontweight='bold', fontsize=8.5)
    axL.tick_params(labelsize=8); axR.tick_params(labelsize=8, left=False)
    if sci_y:
        axL.yaxis.set_major_formatter(FuncFormatter(
            lambda v, pos: f'{v / 1e4:.1g}' if v < 1e4 else f'{v / 1e4:g}'))
        axL.text(0.0, 1.02, '1e4', transform=axL.transAxes, ha='left', va='bottom',
                 fontsize=7, color='black')
    _vals = []
    for s in series:
        for _agg, _zone in [(seg1, T_ZONE), (seg2, None)]:
            _sub = _agg[_agg[gcol] == s['group']]
            if _zone is not None and 'zone_group' in _sub.columns:
                _sub = _sub[_sub['zone_group'] == _zone]
            _sub = _sub[(_sub['n'] >= 1) & (_sub['disp'] <= seg_max_um)]
            _vals += _sub['mean'].dropna().tolist()
    if _vals:
        _ymin = float(min(_vals))
        axL.set_ylim(bottom=_ymin)
        _top = axL.get_ylim()[1]
        _tk = [t for t in axL.get_yticks() if _ymin + 1e-6 < t <= _top]
        axL.set_yticks([_ymin] + _tk)
        if not sci_y:
            axL.set_yticklabels([f'{_ymin:.0f}'] + [f'{t:.0f}' for t in _tk])
    else:
        axL.set_ylim(bottom=0)
    for s in series:
        p1 = _endpoint(seg1, gcol, s['group'], T_ZONE, 'last')
        p2 = _endpoint(seg2, gcol, s['group'], None,   'first')
        if p1 and p2:
            figC.add_artist(ConnectionPatch(xyA=p1, coordsA=axL.transData, xyB=p2, coordsB=axR.transData,
                                            ls=':', color=s['color'], lw=1.3, zorder=3))
    if np.isfinite(hev_boundary_um):
        axL.axvline(hev_boundary_um, ls='--', color='0.5', lw=1.0, zorder=1, clip_on=False)
    axR.axvline(0, ls='--', color='0.5', lw=1.0, zorder=1, clip_on=False)
    axL.spines['top'].set_visible(False); axR.spines['top'].set_visible(False)
    axL.spines['right'].set_visible(False); axR.spines['left'].set_visible(False); axR.spines['right'].set_visible(False)
    dk = 0.02
    axL.plot((1 - dk, 1 + dk), (-dk, dk), transform=axL.transAxes, color='k', clip_on=False, lw=1)
    axR.plot((-dk, dk), (-dk, dk), transform=axR.transAxes, color='k', clip_on=False, lw=1)
    _lh, _ll = axR.get_legend_handles_labels()
    if gcol == 'gene':
        _lh, _ll = _lh[::-1], _ll[::-1]
    axR.legend(_lh, _ll, frameon=False, fontsize=7.5, loc='upper right')

BY0, BH = 1.16, 0.08
BY1 = BY0 + BH
if np.isfinite(hev_boundary_um):
    b = hev_boundary_um
    w = min(12.0, 0.4 * b)
    g = 2.5
    axT, axTR = axsC[0, 0], axsC[0, 1]
    trL, trR = axT.get_xaxis_transform(), axTR.get_xaxis_transform()
    axT.add_patch(Polygon([(0, BY0), (b + w - g, BY0), (b - w - g, BY1), (0, BY1)], closed=True,
                          transform=trL, clip_on=False, facecolor=REGION_COLORS['HEV'],
                          edgecolor='white', lw=0.6, zorder=5))
    axT.add_patch(Polygon([(b + w + g, BY0), (seg_max_um, BY0), (seg_max_um, BY1), (b - w + g, BY1)], closed=True,
                          transform=trL, clip_on=False, facecolor=REGION_COLORS['T cell zone'],
                          edgecolor='white', lw=0.6, zorder=5))
    axTR.add_patch(Rectangle((0, BY0), seg_max_um, BH, transform=trR, clip_on=False,
                             facecolor=REGION_COLORS['non-TLS'], edgecolor='white', lw=0.6, zorder=5))
    for ax, tr, xc, lab in [(axT, trL, max(b * 0.4, 9.0), 'HEV'),
                            (axT, trL, (b + w + seg_max_um) / 2, 'T cell zone'),
                            (axTR, trR, seg_max_um / 2, 'non-TLS')]:
        ax.text(xc, BY1 + 0.03, lab, transform=tr, ha='center', va='bottom',
                fontsize=7.5, fontweight='bold', color='black', clip_on=False)

axsC[1, 0].set_xlabel('From HEV center (um)\n(within T cell zone)', fontweight='bold', fontsize=8)
axsC[1, 1].set_xlabel('From T cell zone edge (um)\n(into non-TLS)', fontweight='bold', fontsize=8)
_t = np.arange(0, seg_max_um + 1, 50)
axsC[1, 0].set_xticks(_t[_t < seg_max_um]); axsC[1, 1].set_xticks(_t[_t > 0])
_ip = INFO.split(' - '); _m = (len(_ip) + 1) // 2
_info_2l = ' - '.join(_ip[:_m]) + '\n' + ' - '.join(_ip[_m:])
figC.suptitle('S9H — CCL19/CCL21 transcripts & imPVL s1/iCAF s2 density vs distance from HEV\n' + _info_2l,
              y=1.52, ha='center', fontweight='bold', fontsize=7)
figC.tight_layout()

BAND_SUFFIX = '' if band_kind == 'sem' else f'_{band_kind}'
pdf_path = os.path.join(out_dir, f'S9H_{GENE_TAG}_density_vs_dist_{BIN_TAG}{BAND_SUFFIX}.pdf')
save_multi_image(pdf_path)
plt.close('all')

print(f'\nSaved: {pdf_path}')
print("Script successfully completed")
