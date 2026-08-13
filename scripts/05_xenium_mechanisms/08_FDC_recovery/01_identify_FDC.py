#!/usr/bin/env python3

# --- Aim
# Recover follicular dendritic cells (FDCs) that were misannotated as Tfh by the label-transfer
# step (both express CXCL13 strongly). For each sample, re-cluster the subset-level Tfh cells,
# run DEG, and auto-flag any leiden cluster whose top DEGs include CR1/CR2 (the most specific FDC
# markers on the panel) as FDC. Saves the per-sample FDC/Tfh barcode CSV (the sole FDC input
# consumed downstream), the full DEG table, and a multi-page diagnostic PDF of UMAPs, dotplots,
# violins and a spatial plot. Supports figure S9J (FDC recovery).

# --- Define environment
# env: spatial-pilot

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

mpl.rcParams['pdf.fonttype'] = 42
sc.set_figure_params(facecolor='white', fontsize=10)


# --- Argparse arguments

parser = argparse.ArgumentParser(
    description="Recover FDC cells misannotated as Tfh by per-sample Tfh subclustering",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
)
parser.add_argument("--h5ad_dir",    help="Path to anndata directory")
parser.add_argument("--sample_id",   help="Sample (section_id) to process")
parser.add_argument("--project_dir", help="Project root directory")
parser.add_argument("--res_dir",     help="Path to results directory")
parser.add_argument("--sample_sheet", default=None, help="Path to sample sheet CSV file")

# Methods parameters (hard-coded defaults; exposed for discoverability, never required)
parser.add_argument("--leiden_res",   type=float, default=1.0,  help="Leiden clustering resolution")
parser.add_argument("--n_pcs",        type=int,   default=30,   help="Number of PCA components")
parser.add_argument("--n_neighbors",  type=int,   default=15,   help="kNN neighbours for the graph")
parser.add_argument("--scale_max",    type=float, default=10.0, help="Clip value for sc.pp.scale")
parser.add_argument("--cd3_min_rank", type=int,   default=1000, help="Candidate FDC cluster must have CD3E and CD3G ranked > this (depleted)")
parser.add_argument("--max_cr_sum",   type=int,   default=250,  help="Max combined CR1+CR2 rank for FDC call")
parser.add_argument("--min_cr2_pct",  type=float, default=5.0,  help="Min pct of cells in candidate cluster expressing CR2")
parser.add_argument("--seed",         type=int,   default=0,    help="Random seed")

args = vars(parser.parse_args())

h5ad_dir     = args["h5ad_dir"]
sample_id    = args["sample_id"]
project_dir  = args["project_dir"]
res_dir      = args["res_dir"]
leiden_res   = args["leiden_res"]
n_pcs        = args["n_pcs"]
n_neighbors  = args["n_neighbors"]
scale_max    = args["scale_max"]
cd3_min_rank = args["cd3_min_rank"]
max_cr_sum   = args["max_cr_sum"]
min_cr2_pct  = args["min_cr2_pct"]
seed         = args["seed"]


# --- Local helper functions

sys.path.append(os.path.join(project_dir, 'scripts'))
from helper_functions import auto_crop_xenium, save_multi_image


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
# adata.X is log-normalised in the P08 Xenium objects.

adata = ad.read_h5ad(os.path.join(h5ad_dir, 'xenium', f'{sample_id}.h5ad'))

crop_coord = auto_crop_xenium(adata)


# --- Load subset-level cell-type annotations
# Primary path, falling back to the regenerated label-transfer results if absent.

level = 'subset'
annotation_file = os.path.join(res_dir, '02_cell_type_mapping', 'label_transfer', 'xenium',
                               f'{level}_level', f'{sample_id}_{level}_pred_cell_types.csv')
if not os.path.exists(annotation_file):
    annotation_file = os.path.join(res_dir, '02_cell_type_mapping_regenerated', 'label_transfer',
                                   'xenium', f'{level}_level', f'{sample_id}_{level}_pred_cell_types.csv')
print(f'Using cell-type annotations: {annotation_file}')
annotations = pd.read_csv(annotation_file, index_col=0)

annotation_key = 'pred_cell_type'
adata.obs[annotation_key] = annotations['pred_cell_type']


# --- Markers

TFH_LABEL = 'T_cells_c3_CD4+_Tfh_CXCL13'

# CR1 / CR2 are the only highly FDC-specific genes on the panel; used for auto-pick.
fdc_specific   = ['CR1', 'CR2']
fdc_markers    = ['CR1', 'CR2', 'FCER2', 'VCAM1', 'CXCL13']
t_cell_markers = ['CD3E', 'CD3G', 'CD4', 'PDCD1', 'ICOS', 'BCL6', 'CXCR5', 'IL21', 'FOXP3']

fdc_specific   = [g for g in fdc_specific   if g in adata.var_names]
fdc_markers    = [g for g in fdc_markers    if g in adata.var_names]
t_cell_markers = [g for g in t_cell_markers if g in adata.var_names]


# --- Output directories

out_root = os.path.join(res_dir, '05_xenium_mechanisms', '08_FDC_recovery', 'per_sample')
csv_dir  = os.path.join(out_root, 'csv_files')
plot_dir = os.path.join(out_root, 'plots')
os.makedirs(csv_dir,  exist_ok=True)
os.makedirs(plot_dir, exist_ok=True)


# --- Subset to Tfh cells

tfh_mask = adata.obs[annotation_key] == TFH_LABEL
n_tfh    = int(tfh_mask.sum())
print(f'\n{sample_id}: {n_tfh} cells annotated as {TFH_LABEL}')

if n_tfh < 20:
    print(f'Too few Tfh cells (<20) — skipping clustering for {sample_id}')
    pd.DataFrame(columns=['cell_id', 'group']).to_csv(
        os.path.join(csv_dir, f'{sample_id}_FDC_Tfh_barcodes.csv'), index=False)
    print("Script successfully completed")
    raise SystemExit

tfh = adata[tfh_mask, :].copy()


# --- Re-process the Tfh subset

# adata.X is already log-normalised — scale + PCA on the full panel (no HVG: targeted panel is small).
n_comps = min(n_pcs, tfh.n_obs - 1, tfh.n_vars - 1)
sc.pp.scale(tfh, max_value=scale_max)
sc.pp.pca(tfh, n_comps=n_comps, random_state=seed)
sc.pp.neighbors(tfh, n_neighbors=min(n_neighbors, tfh.n_obs - 1), random_state=seed)
sc.tl.umap(tfh, random_state=seed)
sc.tl.leiden(tfh, resolution=leiden_res, random_state=seed)

clusters = list(tfh.obs['leiden'].cat.categories)
print(f'Leiden clusters: {clusters}')
print(tfh.obs['leiden'].value_counts().sort_index().to_string())


# --- DEG (Wilcoxon, one-vs-rest per leiden cluster)

sc.tl.rank_genes_groups(tfh, 'leiden', method='wilcoxon')

deg_records = []
for cluster in clusters:
    df = sc.get.rank_genes_groups_df(tfh, group=cluster)
    df['cluster'] = cluster
    df['rank'] = range(1, len(df) + 1)
    deg_records.append(df)
deg_df = pd.concat(deg_records, ignore_index=True)


# --- Auto-identify FDC cluster
# Strategy: among clusters where CD3E and CD3G are deeply ranked (i.e. depleted, ruling out
# true Tfh), pick the cluster with the lowest combined CR1+CR2 rank. Call FDC only if BOTH
# (a) that combined rank is below max_cr_sum, and (b) the candidate cluster has more than
# min_cr2_pct % of its cells actually expressing CR2 (avoids picking up FRC-like clusters
# that share CR2 enrichment relative to peers but have very few CR2+ cells in absolute terms).

def _rank(df_c, gene):
    return int(df_c.loc[gene, 'rank']) if gene in df_c.index else 10**9

# Compute pct expressing per cluster from the unscaled log-norm matrix (adata.X)
cr2_pct = {}
cr1_pct = {}
if 'CR2' in adata.var_names:
    x = adata[tfh.obs_names, 'CR2'].X
    cr2_expr = np.asarray(x.toarray()).flatten() if hasattr(x, 'toarray') else np.asarray(x).flatten()
    cr2_pct = (pd.Series(cr2_expr > 0, index=tfh.obs_names)
                 .groupby(tfh.obs['leiden'], observed=True).mean() * 100).to_dict()
if 'CR1' in adata.var_names:
    x = adata[tfh.obs_names, 'CR1'].X
    cr1_expr = np.asarray(x.toarray()).flatten() if hasattr(x, 'toarray') else np.asarray(x).flatten()
    cr1_pct = (pd.Series(cr1_expr > 0, index=tfh.obs_names)
                 .groupby(tfh.obs['leiden'], observed=True).mean() * 100).to_dict()

candidate_scores = {}
for cluster in clusters:
    df_c = deg_df[deg_df['cluster'] == cluster].set_index('names')
    cd3e_r = _rank(df_c, 'CD3E')
    cd3g_r = _rank(df_c, 'CD3G')
    cr1_r  = _rank(df_c, 'CR1')
    cr2_r  = _rank(df_c, 'CR2')
    combined = cr1_r + cr2_r
    cr2_p  = cr2_pct.get(cluster, 0.0)
    cr1_p  = cr1_pct.get(cluster, 0.0)
    print(f'  cluster {cluster}: CR1 rank={cr1_r} (pct={cr1_p:.1f}%), CR2 rank={cr2_r} (pct={cr2_p:.1f}%), '
          f'sum={combined}, CD3E rank={cd3e_r}, CD3G rank={cd3g_r}')
    if cd3e_r > cd3_min_rank and cd3g_r > cd3_min_rank:
        candidate_scores[cluster] = (combined, cr1_r, cr2_r, cr1_p, cr2_p, cd3e_r, cd3g_r)

fdc_clusters = []
if candidate_scores:
    # Pick ALL CD3-depleted clusters that pass the CR1+CR2 rank-sum and CR2-pct
    # thresholds. Selecting all passing clusters (rather than only the single best)
    # keeps samples that harbour multiple discrete FDC clusters correctly counted.
    sorted_clusters = sorted(candidate_scores.items(), key=lambda kv: kv[1][0])
    rejected = []
    for cluster, scores in sorted_clusters:
        combined, cr1_r, cr2_r, cr1_p, cr2_p, cd3e_r, cd3g_r = scores
        if combined > max_cr_sum:
            rejected.append((cluster, f'CR1+CR2 sum={combined} > {max_cr_sum}'))
            continue
        if cr2_p < min_cr2_pct:
            rejected.append((cluster, f'CR2 pct={cr2_p:.1f}% < {min_cr2_pct}%'))
            continue
        fdc_clusters.append(cluster)
        print(f'  -> FDC cluster {cluster}: CR1+CR2 sum={combined}, CR2 pct={cr2_p:.1f}%, '
              f'CR1 pct={cr1_p:.1f}%, CD3E={cd3e_r}, CD3G={cd3g_r}')
    if not fdc_clusters:
        for cluster, reason in rejected:
            print(f'  -> Rejected cluster {cluster}: {reason}')
        print('  -> No FDC call')
    else:
        for cluster, reason in rejected:
            print(f'  -> Rejected cluster {cluster}: {reason}')
else:
    print(f'  -> No cluster has both CD3E rank > {cd3_min_rank} and CD3G rank > {cd3_min_rank}; no FDC call')

tfh.obs['is_FDC'] = (tfh.obs['leiden'].isin(fdc_clusters)
                                       .map({True: 'FDC', False: 'Tfh'})
                                       .astype('category'))


# --- Save barcodes and DEG table

# Combined FDC + Tfh barcode CSV (for Xenium Explorer overlay)
barcodes_df = pd.DataFrame({
    'cell_id': tfh.obs.index,
    'group':   tfh.obs['is_FDC'].astype(str).values,   # 'FDC' or 'Tfh'
})
barcodes_df.to_csv(os.path.join(csv_dir, f'{sample_id}_FDC_Tfh_barcodes.csv'), index=False)
n_fdc_out = int((barcodes_df['group'] == 'FDC').sum())
n_tfh_out = int((barcodes_df['group'] == 'Tfh').sum())
print(f'\nSaved barcode CSV: {n_fdc_out} FDC, {n_tfh_out} Tfh ({len(barcodes_df)} total)')

deg_df.to_csv(os.path.join(csv_dir, f'{sample_id}_Tfh_subcluster_DEGs.csv'), index=False)


# --- Plots

ctx = f'{sample_id} ({subtype}) - {n_tfh} Tfh cells, leiden res={leiden_res}, FDC clusters={fdc_clusters or "none"}'

# Plot 1: UMAP coloured by leiden
sc.pl.umap(tfh, color='leiden', legend_loc='on data', show=False, title=f'{sample_id} ({subtype})\nTfh leiden')
plt.suptitle(ctx, y=1.05, fontweight='bold')

# Plot 2: UMAP — FDC markers
if fdc_markers:
    sc.pl.umap(tfh, color=fdc_markers, ncols=3, color_map='Reds', show=False)
    plt.suptitle(f'{ctx}\nFDC markers', y=1.02, fontweight='bold')

# Plot 3: UMAP — T cell markers
if t_cell_markers:
    sc.pl.umap(tfh, color=t_cell_markers, ncols=3, color_map='Blues', show=False)
    plt.suptitle(f'{ctx}\nT cell markers', y=1.02, fontweight='bold')

# Plot 4: Dotplot of all markers per leiden
sc.pl.dotplot(tfh, var_names=fdc_markers + t_cell_markers, groupby='leiden',
              standard_scale='var', show=False)
plt.suptitle(f'{ctx}\nMarker expression per leiden cluster', y=1.05, fontweight='bold')

# Plot 5: After FDC call — violins grouped by FDC vs Tfh (clusters collapsed)
if fdc_clusters:
    sc.pl.stacked_violin(tfh, var_names=fdc_markers + t_cell_markers, groupby='is_FDC', show=False)
    plt.suptitle(f'{ctx}\nFDC vs remaining Tfh (violins)', y=1.05, fontweight='bold')

# Plot 7b: Same comparison as a dotplot (shows pct expressing + mean expression)
if fdc_clusters:
    sc.pl.dotplot(tfh, var_names=fdc_markers + t_cell_markers, groupby='is_FDC',
                  standard_scale='var', show=False)
    plt.suptitle(f'{ctx}\nFDC vs remaining Tfh (dotplot)', y=1.05, fontweight='bold')

# Plot 7: Spatial — H&E (left) + FDC/Tfh cells on the tissue, no image (right)
adata.obs['Tfh_subtype'] = 'other'
adata.obs.loc[tfh.obs.index, 'Tfh_subtype'] = tfh.obs['is_FDC'].astype(str).values
adata.obs['Tfh_subtype'] = adata.obs['Tfh_subtype'].astype('category')

palette = {'Tfh': '#1f77b4', 'FDC': '#d62728', 'other': 'lightgrey'}
groups_to_show = [g for g in ['FDC', 'Tfh'] if g in adata.obs['Tfh_subtype'].cat.categories]

fig, axs = plt.subplots(1, 2, figsize=(9, 4))
# Left: aligned H&E only
sc.pl.spatial(adata, color=None, crop_coord=crop_coord, ax=axs[0], show=False)
axs[0].set_title(f'{sample_id} ({subtype}) — H&E', fontweight='bold')
# Right: the recovered FDC & Tfh cells on the same tissue coordinates, without the H&E
# background (img_key=None) and with larger spots so the sparse cells are visible.
sub = adata[adata.obs['Tfh_subtype'].isin(groups_to_show)].copy()
sc.pl.spatial(sub, color='Tfh_subtype', palette=palette, groups=groups_to_show,
              img_key=None, spot_size=150, crop_coord=crop_coord, ax=axs[1], show=False,
              na_in_legend=False, title='FDC vs Tfh on tissue')
fig.subplots_adjust(top=0.80, wspace=0.2)
fig.suptitle(ctx, y=1.02, fontsize=9, fontweight='bold')


# --- Save outputs

pdf_path = os.path.join(plot_dir, f'{sample_id}_FDC_analysis.pdf')
save_multi_image(pdf_path)
plt.close('all')

print(f'\nSaved: {pdf_path}')
print("Script successfully completed")
