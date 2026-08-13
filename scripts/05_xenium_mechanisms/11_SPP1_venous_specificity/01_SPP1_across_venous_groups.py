#!/usr/bin/env python3

# --- Aim
# Compare SPP1 log-norm expression across four ACKR1+ venous EC groups and
# produce Figure S10C (SPP1 + venous/HEV-marker per-gene-max dot plot with
# per-group cell-count bars) plus the supporting violin/box panels.
#
# Stratification follows the HEV-microenvironment approach: HEV-vs-non-HEV is
# determined by HEV-niche membership and the niche-level `close_to_TLS` flag
# (not by per-cell in_TLS), so the (HEV / non-HEV) cell counts are consistent
# across analyses. Non-HEV venous cells (those outside any HEV niche) are then
# further split by per-cell in_TLS to address the SPP1-in-TLS question that's
# specific to this analysis.
#
#   1. Non-HEV venous (stroma)  ACKR1+, NOT in HEV niche, in_TLS=False
#   2. HEV-like (stroma)        ACKR1+ in HEV niche, niche close_to_TLS=False
#   3. Non-HEV venous (in TLS)  the novel question: ACKR1+, NOT in HEV niche,
#                               in_TLS=True
#   4. HEV venous (in TLS)      ACKR1+ in HEV niche, niche close_to_TLS=True
#
# Every ACKR1+ cell falls into exactly one of the four groups (groups 2+4
# together = all ACKR1+ in a niche; groups 1+3 together = all ACKR1+ outside
# any niche).
#
# Hypothesis: HEV venous (in TLS) cells have the highest SPP1 expression
# (canonical observation). Non-HEV venous (in TLS) cells may be intermediate
# if the TLS microenvironment alone is enough to drive SPP1; otherwise their
# expression should resemble the outside-TLS baseline.
#
# Output (one PDF): cohort dot plot of SPP1 + HEV markers (ACKR1, TSPAN7, RAMP3,
# MEOX2), per-gene-max scaled, across the four groups, with per-group cell-count
# bars  [S10C].


# --- Define environment
# env: xen-hev

# --- Load packages

import argparse
import os
import sys
from datetime import datetime

import anndata as ad
import pandas as pd
import numpy as np
import scipy.sparse as sp

import matplotlib as mpl
import matplotlib.pyplot as plt

mpl.rcParams['pdf.fonttype'] = 42
TICK_FS, SUP_FS, ANNO_FS = 8, 14, 7
PANEL_TITLE_FS = 11


# --- Argparse arguments

parser = argparse.ArgumentParser(
    description="SPP1 expression across ACKR1+ venous EC groups (TLS / HEV "
                "splits) — Figure S10C",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
)
parser.add_argument("-s", "--samples", nargs='+', required=True,
                    help="Xenium sample IDs (section_id) to include")
parser.add_argument("--h5ad_dir", required=True,
                    help="Path to anndata_objects directory")
parser.add_argument("--project_dir", required=True, help="Project directory (root)")
parser.add_argument("--res_dir", required=True, help="Path to results directory")
parser.add_argument("--sample_sheet", default=None,
                    help="Path to sample sheet CSV file")
parser.add_argument("--spp1_gene", default='SPP1',
                    help="Probe name for SPP1 in the panel")
parser.add_argument("--ackr1_label", default='Endothelial ACKR1',
                    help="Cell-type label (subset level) for venous ECs")
parser.add_argument("--annotation_level", default='subset',
                    help="Annotation level used to identify ACKR1+ cells")
parser.add_argument("--dotplot_genes", nargs='+',
                    default=['SPP1', 'ACKR1', 'TSPAN7', 'RAMP3', 'MEOX2'],
                    help="Genes shown on the per-gene-max dotplot")

args = vars(parser.parse_args())
samples           = args["samples"]
h5ad_dir          = args["h5ad_dir"]
project_dir       = args["project_dir"]
res_dir           = args["res_dir"]
spp1_gene         = args["spp1_gene"]
ackr1_label       = args["ackr1_label"]
annotation_level  = args["annotation_level"]
dotplot_genes     = args["dotplot_genes"]


# --- Local imports

sys.path.append(os.path.join(project_dir, 'scripts'))
from helper_functions import save_multi_image

sys.path.append(os.path.join(project_dir, 'scripts', '05_xenium_mechanisms'))
import plotting_functions


# --- Group / colour conventions

# Order from "no expected SPP1" → "highest expected SPP1"
group_order = [
    'Non-HEV venous (stroma)',
    'HEV-like (stroma)',
    'Non-HEV venous (in TLS)',
    'HEV venous (in TLS)',
]

# Pull HEV-related colors from existing palettes for consistency; pick distinct
# greys / accents for the two non-HEV groups so they don't clash with the
# blue/brown HEV colors.
def _load_group_colors():
    color_dir = os.path.join(project_dir, 'resources', 'plotting_colors')
    tls_hev = pd.read_csv(os.path.join(color_dir, 'TLS_HEV_colors.csv'),
                          index_col=0)['color'].to_dict()
    return {
        'Non-HEV venous (stroma)': '#bdbdbd',
        'HEV-like (stroma)':            tls_hev.get('non_TLS_HEV', '#8c6d31'),
        'Non-HEV venous (in TLS)':      '#e6ab02',
        'HEV venous (in TLS)':          tls_hev.get('near_TLS_HEV', '#1f77b4'),
    }


group_colors = _load_group_colors()


# --- Helpers

def find_most_recent(dir_path):
    valid = []
    for d in os.listdir(dir_path):
        if not os.path.isdir(os.path.join(dir_path, d)):
            continue
        try:
            datetime.strptime(d, '%Y%m%d')
            valid.append(d)
        except ValueError:
            pass
    if not valid:
        raise FileNotFoundError(f'No YYYYMMDD subdirectory in {dir_path}')
    return os.path.join(dir_path,
                        max(valid, key=lambda x: datetime.strptime(x, '%Y%m%d')))


# --- Sample sheet & subtype mapping

if args.get("sample_sheet"):
    sample_sheet_path = args["sample_sheet"]
else:
    sample_sheet_path = os.path.join(project_dir, 'config',
                                     'sample_sheet_xenium.csv')

sample_sheet = pd.read_csv(sample_sheet_path)
tissue_metadata = pd.read_csv(
    os.path.join(project_dir, 'resources', 'tissue_metadata.csv'))

section_to_subtype = (tissue_metadata.set_index('tissue_id')['clinical_subtype']
                                      .to_dict())
sample_sheet['subtype'] = sample_sheet['tissue_id'].map(section_to_subtype)
subtype_dict = sample_sheet.set_index('section_id')['subtype'].to_dict()


# --- Locate inputs (TLS annotations + HEV niche outputs)

tls_csv_dir = os.path.join(res_dir, '05_xenium_mechanisms',
                           '01_identify_TLS_regions', 'csv_files')

hev_root = os.path.join(res_dir, '05_xenium_mechanisms',
                        '02_identify_HEV_regions')
hev_recent = find_most_recent(hev_root)
print(f'HEV niche CSVs from: {hev_recent}')


# --- Build cell-level table per sample

per_sample_dfs = []
per_sample_adata = []
samples_no_spp1 = []

for sample in samples:

    print(f'\n=== {sample} ===')

    # Load anndata (.X = log-norm)
    adata = ad.read_h5ad(os.path.join(h5ad_dir, 'xenium', f'{sample}.h5ad'))

    # Subset-level annotations
    annot_path = os.path.join(
        res_dir, '02_cell_type_mapping', 'label_transfer', 'xenium',
        f'{annotation_level}_level',
        f'{sample}_{annotation_level}_pred_cell_types.csv')
    if not os.path.exists(annot_path):
        # Fallback: regenerated label-transfer outputs (materialised for all samples)
        annot_path = os.path.join(
            res_dir, '02_cell_type_mapping_regenerated', 'label_transfer',
            'xenium', f'{annotation_level}_level',
            f'{sample}_{annotation_level}_pred_cell_types.csv')
    ann = pd.read_csv(annot_path, index_col=0)['pred_cell_type'].astype(str)
    ann = ann.reindex(adata.obs_names)

    # TLS per-cell annotation
    tls_df = pd.read_csv(
        os.path.join(tls_csv_dir, f'{sample}_cell_TLS_annotations.csv'),
        index_col=0,
    )
    in_tls = (tls_df['in_TLS'].reindex(adata.obs_names)
                              .fillna(False).astype(bool))

    # HEV niche assignments + close_to_TLS flag
    cell_to_niche = pd.read_csv(
        os.path.join(hev_recent, 'csv_files',
                     f'{sample}_HEV_niches_cell_to_niche.csv'),
        index_col=0,
    )
    cell_to_niche['niche_id'] = cell_to_niche['niche_id'].astype('object')
    niche_info = pd.read_csv(
        os.path.join(hev_recent, 'csv_files',
                     f'{sample}_HEV_niches_info.csv'),
    )
    niche_info['niche_id'] = niche_info['niche_id'].astype(str)
    close_to_tls = niche_info.set_index('niche_id')['close_to_TLS'].to_dict()

    niche_per_cell = cell_to_niche['niche_id'].reindex(adata.obs_names)
    in_niche = niche_per_cell.notna()
    niche_close = (niche_per_cell[in_niche].astype(str)
                                            .map(close_to_tls)
                                            .reindex(adata.obs_names)
                                            .fillna(False).astype(bool))

    # Assign venous group (only ACKR1+ cells get a non-null group). Groups 2+4
    # use the niche-level close_to_TLS flag — same approach as the HEV
    # microenvironment analysis — so the HEV / non-HEV totals stay consistent.
    # Per-cell in_TLS is only used to split the non-HEV side (groups 1 vs 3).
    ackr1_mask = ann.eq(ackr1_label).reindex(adata.obs_names, fill_value=False)
    venous_group = pd.Series(np.nan, index=adata.obs_names, dtype=object)

    # 1. Non-HEV venous (stroma): ACKR1+, no niche, in_TLS=False
    m1 = ackr1_mask & (~in_niche) & (~in_tls)
    venous_group.loc[m1] = 'Non-HEV venous (stroma)'
    # 2. HEV-like (stroma): ACKR1+ in niche AND niche far from TLS
    m2 = ackr1_mask & in_niche & (~niche_close)
    venous_group.loc[m2] = 'HEV-like (stroma)'
    # 3. Non-HEV venous (in TLS): ACKR1+, no niche, in_TLS=True
    m3 = ackr1_mask & (~in_niche) & in_tls
    venous_group.loc[m3] = 'Non-HEV venous (in TLS)'
    # 4. HEV venous (in TLS): ACKR1+ in niche AND niche close to TLS
    m4 = ackr1_mask & in_niche & niche_close
    venous_group.loc[m4] = 'HEV venous (in TLS)'

    # SPP1 log-norm expression
    if spp1_gene not in adata.var_names:
        print(f'  WARNING: {spp1_gene} not in panel for {sample}; skipping')
        samples_no_spp1.append(sample)
        continue
    gi = np.where(adata.var_names == spp1_gene)[0][0]
    v = adata.X[:, gi]
    spp1_expr = (v.toarray().ravel() if sp.issparse(v)
                 else np.asarray(v).ravel()).astype(float)

    # Keep only ACKR1+ cells that fall into one of the four groups
    keep_mask = venous_group.isin(group_order).values
    cell_ids = adata.obs_names[keep_mask]
    coords = adata.obs[['x_centroid', 'y_centroid']].to_numpy()

    sample_df = pd.DataFrame({
        'cell_id': cell_ids,
        'sample_id': sample,
        'subtype': subtype_dict.get(sample, '?'),
        'group': venous_group.values[keep_mask],
        'x': coords[keep_mask, 0],
        'y': coords[keep_mask, 1],
        'in_TLS': in_tls.values[keep_mask],
        'in_HEV_niche': in_niche.values[keep_mask],
        'SPP1_lognorm': spp1_expr[keep_mask],
    })
    per_sample_dfs.append(sample_df)
    print(f'  ACKR1+ ECs in groups: {len(sample_df)} ('
          + ', '.join(f'{g}={int((sample_df["group"] == g).sum())}'
                       for g in group_order) + ')')

    # Sub-anndata of the kept cells for the cohort dotplot
    sub = adata[keep_mask].copy()
    sub.obs['group'] = pd.Categorical(venous_group.values[keep_mask],
                                       categories=group_order, ordered=True)
    sub.obs['sample_id'] = sample
    per_sample_adata.append(sub)

    del adata

if not per_sample_dfs:
    raise SystemExit('No samples with SPP1 — nothing to plot')

all_cells = pd.concat(per_sample_dfs, ignore_index=True)
print(f'\nSamples included:  {sorted(all_cells["sample_id"].unique())}')
print(f'Samples dropped:   {samples_no_spp1}')
print(f'Total ACKR1+ ECs across cohort: {len(all_cells)}')


# --- Plot config

# --- Cohort dotplot (SPP1 + HEV markers), per-gene-max scaled [S10C]

all_adata = ad.concat(per_sample_adata, merge='same')
all_adata.obs['group'] = pd.Categorical(
    all_adata.obs['group'], categories=group_order, ordered=True)

genes_in_panel = [g for g in dotplot_genes if g in all_adata.var_names]
missing = [g for g in dotplot_genes if g not in all_adata.var_names]
if missing:
    print(f'WARNING: dotplot genes missing from panel: {missing}')
print(f'Dotplot genes: {genes_in_panel}')

fig, ax = plt.subplots(figsize=(5, 5))
plotting_functions.per_gene_norm_dotplot_with_legends(
    all_adata, var_names=genes_in_panel, groupby='group', ax=ax,
    cmap='Reds', vmin=0, vmax=1,
    standard_scale='var_max', size_normalize='max',
    max_dot_size=100, space_per_dot=0.04, width_ratio=0.8,
    count_palette=group_colors, count_title='Cells',
    count_ax_title='ACKR1+ ECs',
    cbar_title='Mean expression\n(scaled to max)',
    size_title='% expressing\n(scaled to max)',
    transpose=True,
)
plt.suptitle(
    f'S10C — {spp1_gene} & HEV markers in ACKR1+ venous EC groups '
    f'(mean expr & dot size are per-gene fraction of max)',
    fontweight='bold', y=1.1, va='bottom')


# --- Save

out_base = os.path.join(res_dir, '05_xenium_mechanisms',
                        '11_SPP1_venous_specificity')
plots_dir = os.path.join(out_base, 'plots')
os.makedirs(plots_dir, exist_ok=True)

pdf_path = os.path.join(plots_dir, 'S10C_SPP1_across_venous_groups.pdf')
save_multi_image(pdf_path)

print(f'\nSaved: {pdf_path}')
print('Script successfully completed')
