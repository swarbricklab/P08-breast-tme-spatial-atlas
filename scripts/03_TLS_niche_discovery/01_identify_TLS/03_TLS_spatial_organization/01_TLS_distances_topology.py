#!/usr/bin/env python3

# --- Aim
# For each TLS group in each Visium sample, classify its topology relative to
# cancer-annotated spots (intratumoral / peritumoral / stromal_distant) and
# record per-TLS distances to the nearest other TLS and to the nearest cancer
# spot (euclidean µm + hop distance on the in-tissue spot graph). Writes one
# row per TLS (distances_per_TLS.csv) for the cohort overview figure.

# --- Define environment
# env: spatial-pilot

# --- Load packages

import argparse
import os

import pandas as pd
import numpy as np
import anndata as ad
from scipy.spatial import cKDTree
from scipy.spatial.distance import cdist
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import dijkstra

from discrete_tls import define_discrete_tls


# --- Argparse arguments

parser = argparse.ArgumentParser(
    description="Compute per-TLS topology + distances across the Visium cohort",
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
parser.add_argument("--min_TLS_spots", type=int, default=3,
                    help="Minimum directly-adjacent spots for a discrete TLS "
                         "(Figure S3B stringent definition)")

args = vars(parser.parse_args())

project_dir       = args["project_dir"]
h5ad_dir          = args["h5ad_dir"]
res_dir           = args["res_dir"] or os.path.join(project_dir, "results")
TLS_threshold     = args["TLS_threshold"]
TLS_threshold_adj = str(TLS_threshold).split('.')[-1]
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


# --- Compute distances per sample

CANCER_LABELS = ('invasive', 'in_situ')  # epithelial pathology labels marking cancer
GROUP_COL     = 'group_singleton_removed'

rows = []
missing_h5ad = []
missing_groups = []
missing_anno = []

for sample in sample_list:
    groups_file = os.path.join(groups_dir, f'{sample}_TLS_groups.csv')
    anno_file   = os.path.join(
        project_dir, 'data', 'pathology_annotations', 'visium',
        f'{sample}_pathology_annotations.csv',
    )
    h5ad_path = os.path.join(h5ad_dir, 'visium', f'{sample}.h5ad')

    if not os.path.exists(h5ad_path):
        missing_h5ad.append(sample)
        continue
    if not os.path.exists(groups_file):
        missing_groups.append(sample)
        continue

    # Backed mode: skip loading .X to keep memory low
    adata = ad.read_h5ad(h5ad_path, backed='r')
    spot_diameter_pixel = adata.uns['spatial'][sample]['scalefactors']['spot_diameter_fullres']
    scalefactor_um = spot_diameter_pixel / 55  # pixels per µm
    coords_um = pd.DataFrame(
        adata.obsm['spatial'] / scalefactor_um,
        index=adata.obs_names, columns=['x', 'y'],
    )
    adata.file.close()
    del adata

    groups_df  = pd.read_csv(groups_file, index_col=0)
    # Discrete TLS (Fig S3B): >=min_TLS_spots adjacent spots, single-spot-gap merge
    tls_groups = define_discrete_tls(groups_df, coords_um, group_col=GROUP_COL,
                                     min_spots=min_TLS_spots)
    tls_group_ids = list(tls_groups.unique())

    if os.path.exists(anno_file):
        anno_df = pd.read_csv(anno_file, index_col=0)
        cancer_barcodes = anno_df.index[anno_df['epithelial'].isin(CANCER_LABELS)]
        cancer_barcodes_in_adata = coords_um.index.intersection(cancer_barcodes)
        cancer_coords = coords_um.loc[cancer_barcodes_in_adata].values
    else:
        missing_anno.append(sample)
        cancer_barcodes_in_adata = pd.Index([])
        cancer_coords = np.empty((0, 2))

    # --- Build a tight in-tissue spot graph (radius = 1 spot pitch) for hop
    # distance. Use cKDTree on the µm coords with a small numerical slack.
    coords_arr = coords_um.values
    n_spots = coords_arr.shape[0]
    tree = cKDTree(coords_arr)
    pairs = list(tree.query_pairs(r=110.0))  # 100 µm pitch + 10 µm slack
    if pairs:
        r = [p[0] for p in pairs] + [p[1] for p in pairs]
        c = [p[1] for p in pairs] + [p[0] for p in pairs]
        graph = csr_matrix((np.ones(len(r)), (r, c)), shape=(n_spots, n_spots))
    else:
        graph = csr_matrix((n_spots, n_spots))

    name_to_idx = {n: i for i, n in enumerate(coords_um.index)}
    cancer_idxs = np.array(
        [name_to_idx[n] for n in cancer_barcodes_in_adata if n in name_to_idx],
        dtype=int,
    )
    cancer_idx_set = set(cancer_idxs.tolist())

    # Multi-source BFS from all cancer spots: hops from each in-tissue spot
    # to the nearest cancer spot.
    if len(cancer_idxs):
        d = dijkstra(graph, indices=cancer_idxs, directed=False, unweighted=True)
        hops_to_cancer_per_spot = d.min(axis=0)
    else:
        hops_to_cancer_per_spot = np.full(n_spots, np.inf)

    # Indices per TLS group + their union (for "any other TLS" lookups)
    tls_idx_per_group = {
        g: np.array([name_to_idx[b] for b in tls_groups.index[tls_groups == g]])
        for g in tls_group_ids
    }
    all_tls_idxs = (
        np.concatenate(list(tls_idx_per_group.values()))
        if tls_idx_per_group else np.array([], dtype=int)
    )

    def cancer_quadrant_count(centroid, cancer_idxs_in_ring):
        """Number of quadrants (0-4) around the TLS centroid that contain at
        least one cancer-annotated ring spot. Diagnostic column only; not used
        by the classifier."""
        if not cancer_idxs_in_ring:
            return 0
        pts = coords_arr[list(cancer_idxs_in_ring)]
        dx = pts[:, 0] - centroid[0]
        dy = pts[:, 1] - centroid[1]
        angles = (np.arctan2(dy, dx) + 2 * np.pi) % (2 * np.pi)
        quads = (angles / (np.pi / 2)).astype(int)
        quads = np.clip(quads, 0, 3)
        return len(set(quads.tolist()))

    def classify_tls(tls_idx_set, ring1_set, ring_intra_set, centroid):
        """Classify a TLS as intratumoral / peritumoral / stromal_distant.

        - intratumoral: >=30% of TLS spots are cancer-annotated, OR >50% of
          `ring_intra_set` spots are cancer-annotated. At >50% ring cancer the
          spatial distribution cannot be confined to a single side of the TLS,
          so an explicit directionality check would be redundant.
        - peritumoral: >=1 cancer spot in `ring1_set` (the 1-hop ring).
        - stromal_distant: neither.

        ring1_set drives peritumoral; ring_intra_set drives the % criterion for
        intratumoral. Pass ring1_set == ring_intra_set for the 1-hop variant."""
        if tls_idx_set:
            pct_in_tls = len(tls_idx_set & cancer_idx_set) / len(tls_idx_set)
            if pct_in_tls >= 0.30:
                return 'intratumoral'
        if ring_intra_set:
            cancer_in_ring = ring_intra_set & cancer_idx_set
            pct_intra = len(cancer_in_ring) / len(ring_intra_set)
            if pct_intra > 0.50:
                return 'intratumoral'
        if ring1_set and (ring1_set & cancer_idx_set):
            return 'peritumoral'
        return 'stromal_distant'

    for tls_id, tls_idxs in tls_idx_per_group.items():
        this_coords = coords_um.iloc[tls_idxs].values
        centroid    = this_coords.mean(axis=0)

        # Hops from THIS TLS to every spot (drives both ring identification and
        # the hop-to-other-TLS metric).
        d_from_tls         = dijkstra(graph, indices=tls_idxs, directed=False,
                                      unweighted=True)
        hops_from_this_tls = d_from_tls.min(axis=0)
        ring1_set = set(np.where(
            (hops_from_this_tls > 0) & (hops_from_this_tls <= 1)
        )[0].tolist())
        ring2_set = set(np.where(
            (hops_from_this_tls > 0) & (hops_from_this_tls <= 2)
        )[0].tolist())

        tls_idx_set = set(tls_idxs.tolist())

        topo_1hop = classify_tls(tls_idx_set, ring1_set, ring1_set, centroid)
        topo_2hop = classify_tls(tls_idx_set, ring1_set, ring2_set, centroid)

        # Diagnostic numbers — what the classifier actually saw
        cancer_in_ring1 = ring1_set & cancer_idx_set
        cancer_in_ring2 = ring2_set & cancer_idx_set
        n_quads_ring1 = cancer_quadrant_count(centroid, cancer_in_ring1)
        n_quads_ring2 = cancer_quadrant_count(centroid, cancer_in_ring2)

        # Spot-to-spot µm
        other_idxs = np.setdiff1d(all_tls_idxs, tls_idxs)
        if len(other_idxs):
            other_coords     = coords_um.iloc[other_idxs].values
            dist_to_next_tls = float(cdist(this_coords, other_coords).min())
        else:
            dist_to_next_tls = np.nan
        if len(cancer_coords):
            dist_to_cancer     = float(cdist(this_coords, cancer_coords).min())
            centroid_to_cancer = float(cdist(centroid[None, :], cancer_coords).min())
        else:
            dist_to_cancer     = np.nan
            centroid_to_cancer = np.nan

        # Hop distance
        if len(cancer_idxs):
            h = hops_to_cancer_per_spot[tls_idxs].min()
            hop_to_cancer = float(h) if np.isfinite(h) else np.nan
        else:
            hop_to_cancer = np.nan

        if len(other_idxs):
            h = d_from_tls[:, other_idxs].min()
            hop_to_next_tls = float(h) if np.isfinite(h) else np.nan
        else:
            hop_to_next_tls = np.nan

        # Cancer composition stats (helpful for diagnostic plots)
        n_cancer_in_tls   = len(tls_idx_set & cancer_idx_set)
        n_cancer_ring1    = len(cancer_in_ring1)
        n_cancer_ring2    = len(cancer_in_ring2)
        pct_cancer_in_tls = n_cancer_in_tls / len(tls_idx_set) if tls_idx_set else 0.0
        pct_cancer_ring1  = n_cancer_ring1 / len(ring1_set) if ring1_set else 0.0
        pct_cancer_ring2  = n_cancer_ring2 / len(ring2_set) if ring2_set else 0.0

        rows.append({
            'sample_id': sample,
            'subtype': subtype_dict.get(sample, np.nan),
            'TLS_group': tls_id,
            'size_spots': int(len(tls_idxs)),
            'dist_to_next_TLS_um':   dist_to_next_tls,
            'dist_to_cancer_um':     dist_to_cancer,
            'hops_to_next_TLS':      hop_to_next_tls,
            'hops_to_cancer':        hop_to_cancer,
            'centroid_to_cancer_um': centroid_to_cancer,
            'n_cancer_in_TLS':       n_cancer_in_tls,
            'pct_cancer_in_TLS':     pct_cancer_in_tls,
            'n_cancer_in_ring1':     n_cancer_ring1,
            'n_cancer_in_ring2':     n_cancer_ring2,
            'pct_cancer_in_ring1':   pct_cancer_ring1,
            'pct_cancer_in_ring2':   pct_cancer_ring2,
            'cancer_quadrants_ring1': n_quads_ring1,
            'cancer_quadrants_ring2': n_quads_ring2,
            'tls_topology_1hop':     topo_1hop,
            'tls_topology_2hop':     topo_2hop,
        })

    print(f'  {sample}: {len(tls_group_ids)} TLS, cancer spots: {len(cancer_coords)}')

if missing_h5ad:
    print(f'WARNING: missing h5ad for: {missing_h5ad}')
if missing_groups:
    print(f'WARNING: missing TLS_groups csv for: {missing_groups}')
if missing_anno:
    print(f'WARNING: missing pathology annotation for: {missing_anno}')

dist_df = pd.DataFrame(rows)
print(f'\nTotal: {len(dist_df)} TLS rows across {dist_df["sample_id"].nunique()} samples')


# --- Save outputs

out_dir = os.path.join(
    res_dir, '03_TLS_niche_discovery', '01_identify_TLS',
    '03_TLS_spatial_organization', '01_TLS_distances_topology',
    f'threshold_{TLS_threshold_adj}', 'csv_files',
)
os.makedirs(out_dir, exist_ok=True)
out_path = os.path.join(out_dir, 'distances_per_TLS.csv')
dist_df.to_csv(out_path, index=False)
print(f'Saved: {out_path}')
print("Script successfully completed")
