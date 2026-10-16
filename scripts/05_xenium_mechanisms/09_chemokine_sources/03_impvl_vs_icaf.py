#!/usr/bin/env python3

# --- Aim
# Extended Data Fig. 9b: inside TLS there are multiple CCL19+ stromal populations that are distinct lineages —
# iCAF s2 (CAFs MSC iCAF-like s2) carries fibroblast markers, imPVL s1 (PVL Immature s1)
# carries mural/perivascular markers. Pooling cells across all samples, draw the per-gene-max
# scaled dot plot contrasting TLS-resident imPVL s1 vs iCAF s2:
#   Chemokine       CCL19, CCL21
#   Fibroblast/CAF  FAP, LUM, PDGFRA, PDPN
#   Mural/PVL       KCNJ8, MCAM, NOTCH3, RGS5
# Colour = mean expression scaled to per-gene max; dot size = fraction of per-gene max; a
# per-group cell-count bar, colourbar and size legend are rendered to the right. A per-sample
# cell-count table (focus subset x TLS/non-TLS x sample) is written alongside.

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
import geopandas as gpd

import matplotlib as mpl
import matplotlib.pyplot as plt

mpl.rcParams['pdf.fonttype'] = 42
sc.set_figure_params(facecolor='white', fontsize=10)


# --- Argparse arguments
# Cohort: dispatches through submit_multi_sample_job_nci.sh -> run_multi_sample_job.sh, which
# passes --h5ad_dir --project_dir --res_dir --sample_sheet --samples.

parser = argparse.ArgumentParser(
    description="Extended Data Fig. 9b imPVL s1 vs iCAF s2 per-gene-max dot plot inside TLS (cohort)",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
)
parser.add_argument("--h5ad_dir", required=True, help="Path to anndata objects dir")
parser.add_argument("--project_dir", required=True, help="Project root directory")
parser.add_argument("--res_dir", required=True, help="Path to results directory")
parser.add_argument("--samples", required=True, nargs='+', help="List of sample IDs")
parser.add_argument("--sample_sheet", default=None, help="Path to sample sheet CSV (else default)")

args = vars(parser.parse_args())

h5ad_dir    = args["h5ad_dir"]
project_dir = args["project_dir"]
res_dir     = args["res_dir"]
samples     = args["samples"]

if args.get("sample_sheet"):
    sample_sheet_path = args["sample_sheet"]
else:
    sample_sheet_path = os.path.join(project_dir, 'config', 'sample_sheet_xenium.csv')


# --- Local modules

sys.path.append(os.path.join(project_dir, 'scripts', '05_xenium_mechanisms'))
import plotting_functions

sys.path.append(os.path.join(project_dir, 'scripts'))
from helper_functions import save_multi_image


# --- The two focus populations (full subset strings -> short label)

FOCUS = {
    'CAFs MSC iCAF-like s2': 'iCAF s2',    # T-zone FRC candidate (CCL19)
    'PVL Immature s1':       'imPVL s1',   # HEV-associated perivascular (CCL19+CCL21)
}
GROUP_ORDER = ['iCAF s2 (TLS)', 'imPVL s1 (TLS)',
               'iCAF s2 (non-TLS)', 'imPVL s1 (non-TLS)']
TLS_COMPARISON = ['iCAF s2 (TLS)', 'imPVL s1 (TLS)']


# --- Marker panel, grouped by functional category (all on the Xenium Prime 5K panel)

GENE_GROUPS_TARGETED = {
    'Chemokine':        ['CCL19', 'CCL21'],
    'Fibroblast/\nCAF': ['PDGFRA', 'PDPN', 'LUM', 'FAP'],
    'Mural/\nPVL':      ['NOTCH3', 'MCAM', 'RGS5', 'KCNJ8'],
}


# --- Per-gene-max dotplot styling

SPACE_PER_DOT_Y = 0.055   # figure-fraction per dot along the group (y) axis
SPACE_PER_DOT_X = 0.028   # figure-fraction per dot along the gene (x) axis
LABEL_FS        = 9
DOT_MAX         = 150
EDGE_BUFFER     = 0.75    # aligns the count barh rows exactly with the dot rows


# --- Sample metadata

sample_sheet    = pd.read_csv(sample_sheet_path)
tissue_metadata = pd.read_csv(os.path.join(project_dir, 'resources', 'tissue_metadata.csv'))
merged          = sample_sheet.merge(tissue_metadata, on='tissue_id')
subtype_dict    = dict(zip(merged['section_id'], merged['clinical_subtype']))


# --- Cell-type classifications: derive CAF + PVL lineage context (per project convention)

classif = pd.read_csv(os.path.join(project_dir, 'resources', 'cell_types', 'miniatlas',
                                   'miniatlas_subset_classifications.csv'))
classif = classif.drop_duplicates('miniatlas_subset')
stromal_context = classif[classif['miniatlas_major'].isin(['CAFs', 'PVL'])]
stromal_context = stromal_context.sort_values(['miniatlas_major', 'miniatlas_subset'])
context_subsets = stromal_context['miniatlas_subset'].tolist()
subset_to_short = dict(zip(stromal_context['miniatlas_subset'],
                           stromal_context['miniatlas_subset_short']))


# --- Colors (consistent with the rest of the project)

subset_colors = pd.read_csv(os.path.join(project_dir, 'resources', 'plotting_colors',
                                         'miniatlas_subset_colors.csv'))
color_lut = dict(zip(subset_colors['miniatlas_subset'], subset_colors['color']))
FOCUS_COLORS = {short: color_lut.get(full, 'grey') for full, short in FOCUS.items()}
SHORT_COLORS = {subset_to_short[full]: color_lut.get(full, 'grey') for full in context_subsets}


def label_color(label):
    for short, c in FOCUS_COLORS.items():
        if label.startswith(short):
            return c
    return SHORT_COLORS.get(label, 'grey')


# --- Cell-type annotation loader (label transfer)

def load_subset_labels(sample):
    path = os.path.join(res_dir, '02_cell_type_mapping', 'label_transfer',
                        'xenium', 'subset_level', f'{sample}_subset_pred_cell_types.csv')
    df = pd.read_csv(path, index_col=0)
    df.index = df.index.astype(str)
    return df['pred_cell_type']


# --- Output directories

out_root = os.path.join(res_dir, '05_xenium_mechanisms', '09_chemokine_sources',
                        '03_impvl_vs_icaf')
csv_dir  = os.path.join(out_root, 'csv_files')
plot_dir = os.path.join(out_root, 'plots')
os.makedirs(csv_dir,  exist_ok=True)
os.makedirs(plot_dir, exist_ok=True)


# --- Per-gene-max multi-axis dotplot
# Colour = mean expr scaled to per-gene max; dot size = fraction of per-gene max. Count bars
# (horizontal), a vertical colourbar and a size legend are rendered to the right.
def fancy_dotplot(adata_sub, gene_categories, groupby, title,
                  cbar_title='Mean expression\n(scaled to max)', size_normalize='max',
                  sort_genes=True):
    gene_categories = {c: gs for c, gs in gene_categories.items() if gs}
    # Genes within each functional group are plotted alphabetically.
    if sort_genes:
        gene_categories = {c: sorted(gs) for c, gs in gene_categories.items()}
    genes_flat  = [g for gs in gene_categories.values() for g in gs]
    group_order = adata_sub.obs[groupby].cat.categories.tolist()
    n_genes = len(genes_flat)
    n_cats  = len(gene_categories)
    n_grp   = len(group_order)

    # Fit the dot grid + count bars + colourbar + size legend inside x in [0.06, 0.97]
    # and y in [0.10, ~0.90]. The legend block to the right consumes a fixed budget.
    legend_budget = 0.012 + 0.05 + 0.012 + 0.05 + 0.07          # bar gap+cbar gap/w+size gap/w
    spx = min(SPACE_PER_DOT_X,
              (0.97 - 0.06 - (n_cats - 1) * 0.01 - legend_budget) / (n_genes + 3))
    spy = min(SPACE_PER_DOT_Y, 0.70 / n_grp)
    figw = float(np.clip(0.22 / spx, 6, 14))
    row_in = 0.17 if n_grp > 2 else 0.22
    figh = float(np.clip(row_in / spy, 2.6, 9))

    group_colors = {g: label_color(g) for g in group_order}

    fig, ax = plt.subplots(figsize=(figw, figh))
    ax.set_position([0.06, 0.12, 0.4, 0.4])   # helper reads x0/y0 as the grid origin
    pgn = plotting_functions.per_gene_norm_dotplot_with_legends(
        adata_sub, var_names=genes_flat, groupby=groupby, ax=ax,
        cmap='Reds', vmin=0, vmax=1, standard_scale='var_max',
        max_dot_size=DOT_MAX, space_per_dot=spy, width_ratio=spx / spy,
        count_palette=group_colors, count_title='',
        label_fontsize=LABEL_FS - 1,
        size_normalize=size_normalize, color_mode='standard_scale',
        transpose=True, gene_categories=gene_categories,
        category_title_rotation=0,
    )
    for sub_ax in pgn['dotplot_axes']:
        for t in sub_ax.get_xticklabels():
            t.set_fontsize(LABEL_FS - 1)
            t.set_rotation(90)
        sub_ax.grid(False)
    for t in pgn['dotplot_ax'].get_yticklabels():
        t.set_fontsize(LABEL_FS)

    # Drop the helper's auxiliary axes and re-render count bars / cbar / size legend.
    fig.delaxes(pgn['count_ax'])
    fig.delaxes(pgn['cbar_ax'])
    fig.delaxes(pgn['size_ax'])

    dotplot_bbox = pgn['dotplot_bbox']

    def _fmt_count(n):
        n = int(n)
        return f'{n / 1000:.1f}k' if n >= 1000 else f'{n}'

    counts    = adata_sub.obs[groupby].value_counts().reindex(group_order)
    count_max = int(counts.max())
    xmax_global = count_max * 1.25

    bar_pad = 0.012
    bar_w   = 3 * spx
    bar_ax  = fig.add_axes([dotplot_bbox.x1 + bar_pad, dotplot_bbox.y0,
                            bar_w, dotplot_bbox.height])
    bar_ys     = list(range(n_grp))
    bar_colors = [group_colors[g] for g in group_order]
    bar_ax.barh(bar_ys, counts.values, color=bar_colors, edgecolor='black',
                height=0.7, linewidth=0.5)
    bar_ax.set_xlim(0, xmax_global)
    bar_ax.set_ylim(n_grp - 1 + EDGE_BUFFER, -EDGE_BUFFER)
    bar_ax.set_yticks([])
    for v, yp in zip(counts.values, bar_ys):
        bar_ax.text(v + xmax_global * 0.03, yp, _fmt_count(v), va='center',
                    ha='left', fontsize=LABEL_FS - 1)
    for sp_ in ('top', 'right', 'left'):
        bar_ax.spines[sp_].set_visible(False)
    bar_ax.spines['bottom'].set_visible(True)
    bar_ax.spines['bottom'].set_edgecolor('black')
    bar_ax.tick_params(axis='x', length=3, width=0.8, labelsize=LABEL_FS - 1)
    bar_ax.set_xticks([0, count_max])
    bar_ax.set_xticklabels(['0', _fmt_count(count_max)])
    bar_ax.set_title('n cells', fontsize=LABEL_FS, pad=6)
    bar_ax.grid(False)

    # Vertical colourbar spanning the full dot-grid height, right of the count bars
    fig.canvas.draw()
    inv = fig.transFigure.inverted()
    bar_right = inv.transform(bar_ax.get_tightbbox(fig.canvas.get_renderer()))[1][0]
    if n_grp > 2:
        cbar_h  = 3 * spy
        cbar_y0 = dotplot_bbox.y0 + (dotplot_bbox.height - cbar_h) / 2
    else:
        cbar_h  = dotplot_bbox.height
        cbar_y0 = dotplot_bbox.y0
    cax = fig.add_axes([bar_right + 0.05, cbar_y0, 0.012, cbar_h])
    cbar = mpl.colorbar.ColorbarBase(
        cax, cmap=plt.get_cmap('Reds'),
        norm=mpl.colors.Normalize(vmin=0, vmax=1), orientation='vertical',
    )
    cbar.set_label(cbar_title, fontsize=LABEL_FS, labelpad=3)
    cax.yaxis.set_label_position('left')
    cax.tick_params(labelsize=LABEL_FS - 1)
    cax.set_yticks([0, 0.5, 1.0])
    cax.grid(False)

    # Vertical size legend, right of the colourbar
    cbar_pos = cax.get_position()
    size_h = cbar_h * 0.85
    size_y = cbar_y0 + (cbar_h - size_h) / 2
    size_ax = fig.add_axes([cbar_pos.x1 + 0.05, size_y, 0.07, size_h])
    size_ax.set_xlim(0, 1); size_ax.set_ylim(0, 1)
    size_ax.axis('off')
    if size_normalize == 'max':
        size_title = '% expressing\n(scaled to max)'
        ref_sizes  = [0.25, 0.5, 1.0]            # dot area as fraction of DOT_MAX
        ref_labels = ['25%', '50%', '100%']
    else:
        genes_flat = [g for gs in gene_categories.values() for g in gs
                      if g in adata_sub.var_names]
        sub = adata_sub[:, genes_flat]
        X = sub.X.toarray() if hasattr(sub.X, 'toarray') else np.asarray(sub.X)
        gvals = adata_sub.obs[groupby].values
        max_frac = 1e-6
        for grp in group_order:
            m = (gvals == grp)
            if m.sum():
                max_frac = max(max_frac, float((X[m] > 0).mean(axis=0).max()))
        ref_fracs  = [max_frac * 0.25, max_frac * 0.5, max_frac]
        ref_sizes  = [(f / max_frac) ** 1.5 for f in ref_fracs]
        ref_labels = [f'{int(round(f * 100))}%' for f in ref_fracs]
        size_title = 'Fraction of\ncells expressing'
    size_ax.set_title(size_title, fontsize=LABEL_FS, ha='left', x=0, pad=4)
    y_positions = np.linspace(0.88, 0.12, len(ref_sizes))
    for sval, lbl, y in zip(ref_sizes, ref_labels, y_positions):
        size_ax.scatter(0.15, y, s=sval * DOT_MAX, c='lightgrey',
                        edgecolors='black', linewidths=0.5,
                        transform=size_ax.transAxes, clip_on=False)
        size_ax.text(0.45, y, lbl, fontsize=LABEL_FS,
                     va='center', ha='left', transform=size_ax.transAxes)

    title_x = (dotplot_bbox.x0 + dotplot_bbox.x1) / 2
    plt.suptitle(title, fontweight='bold', fontsize=LABEL_FS + 1,
                 x=title_x, y=pgn['suptitle_y'], va='bottom')


# --- Load samples: keep only CAF + PVL cells (all genes), tag subset + TLS membership
# in-TLS is a point-in-polygon against the TLS region polygons, which are built in
# x_centroid/y_centroid (micron) space by 05/01_identify_TLS_regions.

focus_subsets = list(FOCUS.keys())

pooled = []
count_rows = []
for sample in samples:
    print(f'loading {sample}')
    adata = ad.read_h5ad(os.path.join(h5ad_dir, 'xenium', f'{sample}.h5ad'))
    adata.obs_names = adata.obs_names.astype(str)

    # subset label from the label transfer
    adata.obs['subset'] = load_subset_labels(sample).reindex(adata.obs_names).values

    # keep only the stromal context (CAF + PVL) to bound memory
    keep = adata.obs['subset'].isin(context_subsets).values
    adata = adata[keep].copy()
    if adata.n_obs == 0:
        print(f'  no CAF/PVL cells, skipping')
        continue

    # in-TLS via point-in-polygon against the TLS region polygons (x_centroid/y_centroid)
    tls_geojson = os.path.join(res_dir, '05_xenium_mechanisms', '01_identify_TLS_regions',
                               'csv_files', f'{sample}_TLS_regions.geojson')
    if os.path.exists(tls_geojson):
        tls_gdf = gpd.read_file(tls_geojson)
        pts = gpd.GeoDataFrame(
            geometry=gpd.points_from_xy(adata.obs['x_centroid'].values,
                                        adata.obs['y_centroid'].values),
            index=adata.obs_names, crs=tls_gdf.crs,
        )
        joined = gpd.sjoin(pts, tls_gdf[['tls_id', 'geometry']], how='left', predicate='within')
        joined = joined[~joined.index.duplicated(keep='first')]
        adata.obs['in_TLS'] = joined['tls_id'].reindex(adata.obs_names).notna().values
    else:
        print(f'  WARNING no TLS geojson for {sample}, marking all cells non-TLS')
        adata.obs['in_TLS'] = False

    adata.obs['sample']  = sample
    adata.obs['subtype'] = subtype_dict.get(sample)
    adata.obs['short']   = adata.obs['subset'].map(subset_to_short)

    pooled.append(adata)

    for sub in focus_subsets:
        m = adata.obs['subset'] == sub
        count_rows.append({
            'sample':  sample,
            'subtype': subtype_dict.get(sample),
            'subset':  FOCUS[sub],
            'n_TLS':     int((m & adata.obs['in_TLS']).sum()),
            'n_non_TLS': int((m & ~adata.obs['in_TLS']).sum()),
        })

pool = ad.concat(pooled, join='outer', merge='same')
pool.obs_names_make_unique()
print(f'\nPooled {pool.n_obs:,} CAF+PVL cells across {len(pooled)} samples')

counts_df = pd.DataFrame(count_rows)
counts_df.to_csv(os.path.join(csv_dir, 'focus_cell_counts_per_sample.csv'), index=False)
print('\n' + counts_df.to_string(index=False))


# --- Build the focus group column: {short} ({TLS|non-TLS})

pool.obs['tls_tag'] = np.where(pool.obs['in_TLS'], 'TLS', 'non-TLS')
is_focus = pool.obs['subset'].isin(focus_subsets)
pool.obs['group'] = np.where(
    is_focus,
    pool.obs['short'].astype(str) + ' (' + pool.obs['tls_tag'].astype(str) + ')',
    np.nan,
)
pool.obs['group'] = pd.Categorical(pool.obs['group'], categories=GROUP_ORDER)

group_totals = pool.obs['group'].value_counts().reindex(GROUP_ORDER)
print('\nFocus group sizes:\n' + group_totals.to_string())


# --- Restrict the marker panel to genes actually present

present = pool.var_names
def present_groups(groups):
    out = {cat: [g for g in gs if g in present] for cat, gs in groups.items()}
    return {c: gs for c, gs in out.items() if gs}

gene_groups_targeted = present_groups(GENE_GROUPS_TARGETED)
all_defined = {g for gs in GENE_GROUPS_TARGETED.values() for g in gs}
missing = sorted(g for g in all_defined if g not in present)
if missing:
    print(f'\nMarkers absent from panel (dropped): {missing}')


# --- Extended Data Fig. 9b: main comparison, imPVL s1 vs iCAF s2 inside TLS (per-gene max)

tls_pool = pool[pool.obs['group'].isin(TLS_COMPARISON)].copy()
tls_pool.obs['group'] = tls_pool.obs['group'].cat.remove_unused_categories()

fancy_dotplot(tls_pool, gene_groups_targeted, 'group',
              'imPVL s1 vs iCAF s2 inside TLS\n(per-gene max)')


# --- Save outputs

for fnum in plt.get_fignums():
    for ax in plt.figure(fnum).axes:
        ax.grid(False)

pdf_path = os.path.join(plot_dir, 'ED9B_impvl_s1_vs_icaf_s2_TLS.pdf')
save_multi_image(pdf_path)
plt.close('all')

print(f'\nSaved: {pdf_path}')
print(f'Saved: {os.path.join(csv_dir, "focus_cell_counts_per_sample.csv")}')
print("Script successfully completed")
