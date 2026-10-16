#!/usr/bin/env python3

# --- Aim
# DEG analysis (input to the Extended Data Fig. 10i volcano + Extended Data Fig. 10j dotplot): for ACKR1+ endothelial
# cells in intra-TLS HEV niches, compare GC-proximal (close, <= 500 um from a GC)
# vs GC-distal (far) cells, pooled across the GC-positive samples. The unit is the
# cell (Extended Data Fig. 10i reports cell counts, e.g. ~795 close / ~4021 far).
#
# For every gene detected in >= 10 % of the pooled cells, close vs far log-transformed
# expression is compared with a two-sided Mann-Whitney U test, Benjamini-Hochberg
# corrected; the effect size is the mean log-transformed difference (close - far). The
# per-gene table is written to DEG_close_vs_far.csv and rendered by
# 03_volcano_dotplot.py.
#
# Only GC-positive samples contribute (auto-skip: a sample is dropped when its
# per-sample HEV-GC distance CSV has no niche with a finite GC distance). No
# sample IDs are hard-coded.

# --- Define environment
# env: xen-hev

# --- Load packages

import argparse
import os

import anndata as ad
import pandas as pd
import numpy as np
import scipy.sparse as sp
from scipy import stats
from statsmodels.stats.multitest import multipletests


# --- Argparse arguments

parser = argparse.ArgumentParser(
    description="DEG close vs far from GC (cell-level, ACKR1+ intra-TLS HEVs) -> CSV",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
)
parser.add_argument("-s", "--samples", nargs='+', required=True,
                    help="Cohort samples; GC-negative ones auto-skipped")
parser.add_argument("--h5ad_dir", required=True)
parser.add_argument("--project_dir", required=True)
parser.add_argument("--res_dir", required=True)
parser.add_argument("--sample_sheet", default=None,
                    help="Accepted for runner compatibility (unused here)")
parser.add_argument("--gc_threshold", type=float, default=0.8,
                    help="GC score threshold used upstream (selects per-sample CSVs)")
parser.add_argument("--near_gc_um", type=float, default=500.0,
                    help="Cell is 'close to GC' if its niche distance <= this (um)")
parser.add_argument("--ackr1_label", default='Endothelial ACKR1')
parser.add_argument("--annotation_level", default='subset')
parser.add_argument("--min_pct_detected", type=float, default=0.10,
                    help="Test a gene only if expressed (>0) in >= this fraction "
                         "of the pooled close+far cells")

args = vars(parser.parse_args())
samples          = args["samples"]
h5ad_dir         = args["h5ad_dir"]
project_dir      = args["project_dir"]
res_dir          = args["res_dir"]
gc_threshold     = args["gc_threshold"]
near_gc_um       = args["near_gc_um"]
ackr1_label      = args["ackr1_label"]
annotation_level = args["annotation_level"]
min_pct_detected = args["min_pct_detected"]


# --- Load per-sample HEV-GC distance CSVs and drop samples without GCs

dist_csv_dir = os.path.join(res_dir, '05_xenium_mechanisms', '10_HEV_GC_proximity',
                            '01_classify_HEV_by_GC_dist', 'csv_files')
print(f'HEV-GC distance CSVs from: {dist_csv_dir}')

dist_dfs = []
for s in samples:
    f = os.path.join(dist_csv_dir, f'{s}_intra_TLS_HEV_GC_dist.csv')
    if not os.path.isfile(f):
        print(f'WARNING: missing {f} — skipping {s}')
        continue
    dist_dfs.append(pd.read_csv(f))
if not dist_dfs:
    raise FileNotFoundError('No per-sample HEV-GC distance CSVs found — run '
                            '01_classify_HEV_by_GC_dist.py first.')
dist_df = pd.concat(dist_dfs, ignore_index=True)
dist_df['niche_id'] = dist_df['niche_id'].astype(str)

samples_with_gc = (dist_df.groupby('sample_id')['min_distance_to_GC']
                          .apply(lambda d: np.isfinite(d).any()))
samples_to_keep = samples_with_gc[samples_with_gc].index.tolist()
print(f'GC-positive samples kept: {samples_to_keep}')
dist_df = dist_df[dist_df['sample_id'].isin(samples_to_keep)].copy()


# --- Build per-cell expression table (ACKR1+ intra-TLS HEV endothelial cells)

hev_dir = os.path.join(res_dir, '05_xenium_mechanisms', '02_identify_HEV_regions')
hev_recent = hev_dir
print(f'HEV niches loaded from: {hev_recent}')

label_dir = os.path.join(res_dir, '02_cell_type_mapping',
                         'label_transfer', 'xenium', f'{annotation_level}_level')

all_expr = []
all_meta = []
gene_names_ref = None

for sample in samples_to_keep:
    adata = ad.read_h5ad(os.path.join(h5ad_dir, 'xenium', f'{sample}.h5ad'))

    cell_to_niche = pd.read_csv(
        os.path.join(hev_recent, 'csv_files',
                     f'{sample}_HEV_niches_cell_to_niche.csv'),
        index_col=0,
    ).dropna(subset=['niche_id'])
    cell_to_niche['niche_id'] = cell_to_niche['niche_id'].astype(str)

    sample_dist = dist_df[dist_df['sample_id'] == sample].set_index('niche_id')
    cell_to_niche = cell_to_niche[cell_to_niche['niche_id'].isin(sample_dist.index)]

    annot_path = os.path.join(label_dir,
                              f'{sample}_{annotation_level}_pred_cell_types.csv')
    ann = pd.read_csv(annot_path, index_col=0)['pred_cell_type'].astype(str)
    ackr1_cells = ann[ann == ackr1_label].index

    candidate = cell_to_niche.index.intersection(adata.obs_names)
    valid = candidate.intersection(ackr1_cells)
    if len(valid) == 0:
        print(f'  {sample}: 0 ACKR1+ HEV cells — skip')
        continue
    print(f'  {sample}: {len(valid)} ACKR1+ HEV cells '
          f'({cell_to_niche.loc[valid, "niche_id"].nunique()} niches)')

    sub = adata[valid]
    if gene_names_ref is None:
        gene_names_ref = list(sub.var_names)
    else:
        assert list(sub.var_names) == gene_names_ref, f'gene order mismatch in {sample}'

    niches = cell_to_niche.loc[valid, 'niche_id'].values
    X = sub.X.toarray() if sp.issparse(sub.X) else np.asarray(sub.X)
    cell_index = [f'{sample}__{c}' for c in valid]

    expr_df = pd.DataFrame(X, columns=gene_names_ref, index=cell_index)
    meta = pd.DataFrame({
        'sample_id':          sample,
        'subtype':            sample_dist['subtype'].iloc[0],
        'niche_id':           niches,
        'min_distance_to_GC': sample_dist['min_distance_to_GC'].reindex(niches).values,
    }, index=cell_index)
    meta['close_to_GC'] = meta['min_distance_to_GC'] <= near_gc_um

    all_expr.append(expr_df)
    all_meta.append(meta)
    del adata, sub

expr_df = pd.concat(all_expr, axis=0)
meta_df = pd.concat(all_meta, axis=0)
assert (expr_df.index == meta_df.index).all()
print(f'\nTotal ACKR1+ HEV cells: {len(expr_df)}; genes in panel: {expr_df.shape[1]}')


# --- DEG: Mann-Whitney close vs far per gene, BH-corrected

close_mask = meta_df['close_to_GC'].astype(bool).values
far_mask   = ~close_mask
n_close = int(close_mask.sum())
n_far   = int(far_mask.sum())
print(f'Close cells: {n_close}  Far cells: {n_far}')
if n_close < 3 or n_far < 3:
    raise RuntimeError(f'Need >=3 cells per group; got {n_close}/{n_far}')

X_close = expr_df.values[close_mask]
X_far   = expr_df.values[far_mask]

n_total = X_close.shape[0] + X_far.shape[0]
min_units_req = int(np.ceil(min_pct_detected * n_total))
detected = ((X_close > 0).sum(axis=0) + (X_far > 0).sum(axis=0)) >= min_units_req
genes_kept = [g for g, k in zip(expr_df.columns, detected) if k]
print(f'Genes kept (detected in >={min_pct_detected:.0%} of {n_total} cells, '
      f'i.e. >={min_units_req} cells): {len(genes_kept)} of {expr_df.shape[1]}')

X_close_k = X_close[:, detected]
X_far_k   = X_far[:, detected]

mean_close = X_close_k.mean(axis=0)
mean_far   = X_far_k.mean(axis=0)
mean_diff  = mean_close - mean_far          # log1p space ~ log fold change

pct_close = (X_close_k > 0).mean(axis=0)
pct_far   = (X_far_k   > 0).mean(axis=0)

pvals = np.full(len(genes_kept), np.nan)
for i in range(len(genes_kept)):
    try:
        pvals[i] = stats.mannwhitneyu(X_close_k[:, i], X_far_k[:, i],
                                      alternative='two-sided').pvalue
    except ValueError:
        pvals[i] = np.nan

valid_p = np.isfinite(pvals)
padj = np.full_like(pvals, np.nan)
padj[valid_p] = multipletests(pvals[valid_p], method='fdr_bh')[1]

deg = pd.DataFrame({
    'gene':       genes_kept,
    'n_close':    n_close,
    'n_far':      n_far,
    'mean_close': mean_close,
    'mean_far':   mean_far,
    'mean_diff':  mean_diff,
    'pct_close':  pct_close,
    'pct_far':    pct_far,
    'pval':       pvals,
    'padj':       padj,
}).sort_values('padj').reset_index(drop=True)


# --- Save DEG CSV (input to 03_volcano_dotplot.py)

out_base = os.path.join(res_dir, '05_xenium_mechanisms', '10_HEV_GC_proximity',
                        '02_DEG_close_vs_far')
csv_dir = os.path.join(out_base, 'csv_files')
os.makedirs(csv_dir, exist_ok=True)

deg_csv = os.path.join(csv_dir, 'DEG_close_vs_far.csv')
deg.to_csv(deg_csv, index=False)
print(f'\nSaved: {deg_csv}')


# --- Console summary

sig = (deg['padj'] < 0.05) & (deg['mean_diff'].abs() > 0.05)
print(f'\n{int((sig & (deg["mean_diff"] > 0)).sum())} up (close), '
      f'{int((sig & (deg["mean_diff"] < 0)).sum())} down (far), '
      f'{int((~sig).sum())} ns of {len(deg)} genes')
print('\nTop 10 up in close (by padj):')
print(deg[sig & (deg['mean_diff'] > 0)].head(10)[
    ['gene', 'mean_close', 'mean_far', 'mean_diff', 'pct_close', 'pct_far', 'padj']
].to_string(index=False))
print('\nScript successfully completed')
