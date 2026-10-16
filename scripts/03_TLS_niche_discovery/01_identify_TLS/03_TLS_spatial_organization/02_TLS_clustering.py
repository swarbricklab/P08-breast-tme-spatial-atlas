#!/usr/bin/env python3

# --- Aim
# Quantify TLS spatial clustering per Visium sample with DBSCAN on the pairwise
# edge-to-edge (min spot-to-spot) distance matrix between TLS, using a
# biologically-motivated absolute eps (default 500 µm). Labels each TLS
# 'clustered' (member of a DBSCAN cluster) or 'lone' (noise). Writes one row per
# TLS (clustering_per_TLS.csv, consumed by the cohort overview) and one row per
# sample (clustering_per_sample.csv, aggregate summary).

# --- Define environment
# env: spatial-pilot

# --- Load packages

import argparse
import os

import numpy as np
import pandas as pd
import anndata as ad
from scipy.spatial import ConvexHull
from scipy.spatial.distance import cdist
from sklearn.cluster import DBSCAN

from discrete_tls import define_discrete_tls


# --- Argparse arguments

parser = argparse.ArgumentParser(
    description="Per-sample TLS spatial clustering via DBSCAN (edge-to-edge)",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
)
parser.add_argument("-s", "--samples", nargs="*", type=str, default=None,
                    help="List of section_ids (default: all in the sample sheet)")
parser.add_argument("--project_dir", required=True, help="Project directory (root)")
parser.add_argument("--h5ad_dir", required=True,
                    help="Path to directory for .h5ad files")
parser.add_argument("-r", "--res_dir", help="Path to results directory")
parser.add_argument("--sample_sheet", help="Path to sample sheet CSV file")
parser.add_argument("--TLS_threshold", type=float, default=0.01)
parser.add_argument("--dbscan_eps_um", type=str, default="500",
                    help="Comma-separated DBSCAN eps values in µm (edge-to-edge)")
parser.add_argument("--min_TLS_spots", type=int, default=3,
                    help="Minimum directly-adjacent spots for a discrete TLS "
                         "(Extended Data Fig. 3b stringent definition)")

args = vars(parser.parse_args())

project_dir       = args["project_dir"]
h5ad_dir          = args["h5ad_dir"]
res_dir           = args["res_dir"] or os.path.join(project_dir, "results")
TLS_threshold     = args["TLS_threshold"]
TLS_threshold_adj = str(TLS_threshold).split('.')[-1]
dbscan_eps_list   = [float(e) for e in args["dbscan_eps_um"].split(",")]
min_TLS_spots     = args["min_TLS_spots"]


# --- Sample metadata

if args.get("sample_sheet"):
    sample_sheet_path = args["sample_sheet"]
else:
    sample_sheet_path = os.path.join(project_dir, 'config', 'sample_sheet_visium.csv')

sample_sheet    = pd.read_csv(sample_sheet_path)
sample_lut      = sample_sheet.set_index('section_id')
tissue_metadata = pd.read_csv(os.path.join(project_dir, 'resources', 'tissue_metadata.csv'))
sample_meta = pd.merge(
    sample_sheet, tissue_metadata[['tissue_id', 'clinical_subtype']], on='tissue_id',
)
subtype_dict = dict(zip(sample_meta['section_id'], sample_meta['clinical_subtype']))

sample_list = args["samples"] or list(sample_sheet['section_id'])


# --- Locate the TLS_groups directory

groups_parent = os.path.join(
    res_dir, '03_TLS_niche_discovery', '01_identify_TLS', '02_group_TLS_high',
    f'threshold_{TLS_threshold_adj}', 'csv_files',
)
groups_dir = groups_parent
print(f"Reading TLS_groups from: {groups_dir}")


# --- Compute metrics per sample

GROUP_COL = 'group_singleton_removed'

rows_sample    = []
rows_tls       = []
missing_h5ad   = []
missing_groups = []

for sample in sample_list:
    groups_file = os.path.join(groups_dir, f'{sample}_TLS_groups.csv')
    h5ad_path   = os.path.join(h5ad_dir, 'visium', f'{sample}.h5ad')

    if not os.path.exists(h5ad_path):
        missing_h5ad.append(sample)
        continue
    if not os.path.exists(groups_file):
        missing_groups.append(sample)
        continue

    # Backed read — we only need obsm['spatial'] and uns['spatial']
    adata = ad.read_h5ad(h5ad_path, backed='r')
    spot_diameter_pixel = adata.uns['spatial'][sample]['scalefactors']['spot_diameter_fullres']
    scalefactor_um = spot_diameter_pixel / 55  # pixels per µm
    coords_um = pd.DataFrame(
        adata.obsm['spatial'] / scalefactor_um,
        index=adata.obs_names, columns=['x', 'y'],
    )
    adata.file.close()
    del adata

    # Observation window = convex hull of all in-tissue spots
    hull = ConvexHull(coords_um.values)
    tissue_area_um2 = float(hull.volume)   # scipy: 2D "volume" = polygon area
    tissue_area_mm2 = tissue_area_um2 / 1e6

    groups_df  = pd.read_csv(groups_file, index_col=0)
    # Discrete TLS (Extended Data Fig. 3b): >=min_TLS_spots adjacent spots, single-spot-gap merge
    tls_groups = define_discrete_tls(groups_df, coords_um, group_col=GROUP_COL,
                                     min_spots=min_TLS_spots)
    tls_ids    = list(tls_groups.unique())
    n_tls      = len(tls_ids)

    tls_coord_arr = {
        g: coords_um.loc[tls_groups.index[tls_groups == g]].values
        for g in tls_ids
    }
    tls_sizes = np.array([len(tls_coord_arr[g]) for g in tls_ids], dtype=int)

    # Pairwise edge-to-edge (min spot-to-spot) distance matrix
    D_edge = np.zeros((n_tls, n_tls))
    for i in range(n_tls):
        for j in range(i + 1, n_tls):
            d = float(cdist(tls_coord_arr[tls_ids[i]],
                            tls_coord_arr[tls_ids[j]]).min())
            D_edge[i, j] = d
            D_edge[j, i] = d

    # --- Assemble sample-level row (init DBSCAN metrics as NaN)
    row = {
        'sample_id':           sample,
        'subtype':             subtype_dict.get(sample, np.nan),
        'slide_size':          sample_lut.loc[sample, 'slide_size'],
        'n_TLS':               n_tls,
        'tissue_area_mm2':     tissue_area_mm2,
        'TLS_density_per_mm2': (n_tls / tissue_area_mm2) if tissue_area_mm2 > 0 else np.nan,
    }
    for eps in dbscan_eps_list:
        eps_tag = str(int(eps))
        for s in ('n_lone', 'n_clustered', 'n_clusters', 'frac_lone'):
            row[f'dbscan_eps{eps_tag}_{s}'] = np.nan

    # Per-TLS row stubs (init DBSCAN labels as NaN)
    tls_row_stubs = []
    for i, g in enumerate(tls_ids):
        stub = {
            'sample_id':  sample,
            'subtype':    row['subtype'],
            'TLS_group':  g,
            'size_spots': int(tls_sizes[i]),
        }
        for eps in dbscan_eps_list:
            eps_tag = str(int(eps))
            stub[f'dbscan_label_eps{eps_tag}'] = np.nan
            stub[f'dbscan_class_eps{eps_tag}'] = np.nan
        tls_row_stubs.append(stub)

    # --- DBSCAN — runs at any n_tls >= 1 (single TLS just gets label -1)
    if n_tls >= 1:
        for eps in dbscan_eps_list:
            eps_tag = str(int(eps))
            labels = DBSCAN(eps=eps, min_samples=2,
                            metric='precomputed').fit_predict(D_edge)
            n_lone      = int((labels == -1).sum())
            n_clustered = int((labels != -1).sum())
            n_clusters  = int(labels.max() + 1) if (labels != -1).any() else 0
            frac_lone   = n_lone / n_tls

            for i in range(n_tls):
                two_class = 'lone' if labels[i] == -1 else 'clustered'
                tls_row_stubs[i][f'dbscan_label_eps{eps_tag}'] = int(labels[i])
                tls_row_stubs[i][f'dbscan_class_eps{eps_tag}'] = two_class

            row[f'dbscan_eps{eps_tag}_n_lone']      = n_lone
            row[f'dbscan_eps{eps_tag}_n_clustered'] = n_clustered
            row[f'dbscan_eps{eps_tag}_n_clusters']  = n_clusters
            row[f'dbscan_eps{eps_tag}_frac_lone']   = frac_lone

    rows_sample.append(row)
    rows_tls.extend(tls_row_stubs)

    print(f'  {sample}: n_TLS={n_tls}, area={tissue_area_mm2:.1f} mm²')

if missing_h5ad:
    print(f'WARNING: missing h5ad for: {missing_h5ad}')
if missing_groups:
    print(f'WARNING: missing TLS_groups csv for: {missing_groups}')

sample_df = pd.DataFrame(rows_sample)
tls_df    = pd.DataFrame(rows_tls)
print(f'\nTotal: {len(sample_df)} samples, {len(tls_df)} TLS rows')


# --- Save outputs

out_dir = os.path.join(
    res_dir, '03_TLS_niche_discovery', '01_identify_TLS',
    '03_TLS_spatial_organization', '02_TLS_clustering',
    f'threshold_{TLS_threshold_adj}', 'csv_files',
)
os.makedirs(out_dir, exist_ok=True)
sample_out = os.path.join(out_dir, 'clustering_per_sample.csv')
tls_out    = os.path.join(out_dir, 'clustering_per_TLS.csv')
sample_df.to_csv(sample_out, index=False)
tls_df.to_csv(tls_out, index=False)
print(f'Saved: {sample_out}')
print(f'Saved: {tls_out}')
print("Script successfully completed")
