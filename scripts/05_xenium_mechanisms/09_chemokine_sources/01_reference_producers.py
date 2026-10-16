#!/usr/bin/env python3

# --- Aim
# Ground the niche-chemokine source attribution in the scRNA-seq reference (miniatlas).
# This single cohort script merges three source analyses:
#   (a) CCL19 / CCL21 producer subsets  (miniatlas subsets with >= X% of cells expressing),
#   (b) CXCL13 producer subsets         (same metric),
#   (c) the Extended Data Fig. 9a dot plot of the TLS-niche chemokines / receptors / HEV markers across every
#       reference subset, with candidate producers (> 15% of cells expressing any of
#       CCL19 / CXCL13 / CCL21) highlighted.
# For every subset we compute the % of cells expressing (raw counts > 0) and the mean
# log-norm expression. Producer lists are emitted at BOTH thresholds: over_10 (>= 10%, used
# by the downstream Extended Data Fig. 9i reference-producer correlation) and over_15 (>= 15%, the headline).
# FDC are NOT in the reference and are always added downstream as a CXCL13 producer regardless
# of this list.

# --- Define environment
# env: xen-hev

# --- Load packages

import argparse
import os
import sys

import numpy as np
import pandas as pd
import scipy.sparse as sp

import anndata as ad
import scanpy as sc

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
from matplotlib.colors import to_rgb, to_hex

mpl.rcParams['pdf.fonttype'] = 42
sc.set_figure_params(facecolor="white", figsize=(8, 8), fontsize=10)
sc.settings.verbosity = 3


# --- Argparse arguments
# --samples / --h5ad_dir / --sample_sheet are accepted so this dispatches unchanged through
# the cohort runner (submit_multi_sample_job_nci.sh -> run_multi_sample_job.sh), but the
# figure and producer lists are built purely from the miniatlas reference.

parser = argparse.ArgumentParser(
    description="Reference (miniatlas) CCL19/CCL21 & CXCL13 producers + Extended Data Fig. 9a niche-gene dot plot",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
)
parser.add_argument("--h5ad_dir", default=None, help="Path to anndata objects dir (holds sc_ref/)")
parser.add_argument("--project_dir", required=True, help="Project directory (root)")
parser.add_argument("--res_dir", required=True, help="Path to results directory")
parser.add_argument("--sample_sheet", default=None, help="(inert) accepted for the runner contract")
parser.add_argument("-s", "--samples", nargs="*", default=None,
                    help="(inert) accepted for the runner contract")
parser.add_argument("--min_pct_exprs", nargs="*", type=float, default=[10.0, 15.0],
                    help="thresholds (%%): a subset is a producer if >= this %% of its cells "
                         "express; one output subdir over_<pct> per value")
parser.add_argument("--bold_pct", type=float, default=15.0,
                    help="Extended Data Fig. 9a: bold a subset if > this %% of its cells express any of "
                         "CCL19/CXCL13/CCL21 (headline candidate-producer threshold)")

args = vars(parser.parse_args())

h5ad_dir      = args["h5ad_dir"]
project_dir   = args["project_dir"]
res_dir       = args["res_dir"]
MIN_PCT_EXPRS = list(args["min_pct_exprs"])
BOLD_PCT      = args["bold_pct"]

if h5ad_dir is None:
    h5ad_dir = os.path.join(project_dir, 'data', 'anndata_objects')

LINEAGE_ORDER = ['Epithelial', 'Immune', 'Stromal']
SUBSET_COL    = 'celltype_subset'      # miniatlas obs column holding the subset label
CYCLING_COL   = 'miniatlas_cycling'    # classifications column (major class w/ Cycling split out)
BOLD_GENES    = ['CCL19', 'CXCL13', 'CCL21']   # Extended Data Fig. 9a: producers of any of these get a bold label


# --- Output directories (plots / csv_files split)

out_base = os.path.join(res_dir, '05_xenium_mechanisms', '09_chemokine_sources',
                        '01_reference_producers')
csv_root  = os.path.join(out_base, 'csv_files')
plot_root = os.path.join(out_base, 'plots')
os.makedirs(csv_root,  exist_ok=True)
os.makedirs(plot_root, exist_ok=True)


# --- Local plotting helper

sys.path.append(os.path.join(project_dir, 'scripts'))
from helper_functions import save_multi_image


# --- Load & prepare the miniatlas reference
# miniatlas_GEO.h5ad already stores raw counts in layers['counts'] and log-norm
# (normalize_total 1e4 + log1p) in .X, with ensembl-id var_names + var['original_name'] gene
# symbols. % expressing needs raw counts -> we use layers['counts'].

miniatlas_file = os.path.join(h5ad_dir, 'sc_ref', 'miniatlas_GEO.h5ad')
adata = ad.read_h5ad(miniatlas_file)
print(f'Loaded miniatlas: {adata.shape[0]:,} cells x {adata.shape[1]:,} genes')

adata.var_names = adata.var['original_name'].astype(str)
adata.var_names_make_unique()
if 'counts' not in adata.layers:
    raise KeyError("miniatlas_GEO.h5ad has no layers['counts'] (raw counts) — % expressing "
                   "needs raw counts; regenerate the reference with raw counts stored.")
adata.layers['log_norm'] = adata.X.copy()   # .X is already log-norm

if SUBSET_COL not in adata.obs:
    raise KeyError(f"{SUBSET_COL} not in obs; available: {list(adata.obs.columns)}")


# --- Classifications: subset ordering (lineage -> major -> subset) + color luts

color_dir = os.path.join(project_dir, 'resources', 'plotting_colors')

classif = pd.read_csv(os.path.join(project_dir, 'resources', 'cell_types', 'miniatlas',
                                   'miniatlas_subset_classifications.csv'))
classif['lineage'] = pd.Categorical(classif['lineage'], categories=LINEAGE_ORDER, ordered=True)
classif = classif.sort_values(['lineage', 'miniatlas_major', 'miniatlas_subset'])
classif_dedup = classif[~classif['miniatlas_subset'].duplicated(keep='first')].copy()
classif_dedup.index = classif_dedup['miniatlas_subset']

lineage_lut = pd.read_csv(os.path.join(color_dir, 'miniatlas_lineage_colors.csv')).set_index(
    'lineage')['color'].to_dict()
miniatlas_major_lut = pd.read_csv(os.path.join(color_dir, 'miniatlas_major_colors.csv'))
miniatlas_major_lut = dict(zip(miniatlas_major_lut['miniatlas_major'], miniatlas_major_lut['color']))
subset_lineage = classif.set_index('miniatlas_subset')['lineage'].astype(str).to_dict()

# Cell types actually present in the reference, kept in the lineage-sorted order
present = set(adata.obs[SUBSET_COL].dropna().unique())
subset_order = [s for s in classif_dedup['miniatlas_subset'] if s in present]
adata.obs[SUBSET_COL] = pd.Categorical(adata.obs[SUBSET_COL], categories=subset_order, ordered=True)


# --- Per-subset expression: % expressing (counts > 0) + mean log-norm, for any gene list

def expression_table(genes):
    """Per-subset % expressing (raw counts > 0) and mean log-norm for `genes`."""
    missing = [g for g in genes if g not in adata.var_names]
    if missing:
        raise KeyError(f"genes missing from reference: {missing}")
    gi     = [adata.var_names.get_loc(g) for g in genes]
    counts = adata.layers['counts'][:, gi]
    logn   = adata.layers['log_norm'][:, gi]
    counts = counts.toarray() if sp.issparse(counts) else np.asarray(counts)
    logn   = logn.toarray()   if sp.issparse(logn)   else np.asarray(logn)

    tbl = pd.DataFrame({'subset': adata.obs[SUBSET_COL].to_numpy()})
    for i, g in enumerate(genes):
        tbl[f'{g}_pos']  = (counts[:, i] > 0).astype(float)
        tbl[f'{g}_logn'] = logn[:, i]
    grp = tbl.groupby('subset', observed=True)

    res = pd.DataFrame({'subset': subset_order})
    res['lineage'] = [subset_lineage.get(s, 'NA') for s in subset_order]
    res['n_cells'] = grp.size().reindex(subset_order).to_numpy()
    for g in genes:
        res[f'{g}_pct']  = (grp[f'{g}_pos'].mean() * 100).reindex(subset_order).to_numpy()
        res[f'{g}_mean'] = grp[f'{g}_logn'].mean().reindex(subset_order).to_numpy()
    return res


# --- Producer emit: per-threshold producer lists + a simple scatter dot plot

def emit_producers(genes, name, expr_csv, prod_csv, pdf_name, legend_anchor):
    """For each threshold in MIN_PCT_EXPRS, write producer lists (over_<pct>/) + a dot plot.
    `name` -> subdir under csv_files/ & plots/; file names are kept stable so the
    downstream attribution / correlation scripts can read them unchanged."""
    res  = expression_table(genes)
    n    = len(subset_order)
    SMAX = 190.0                                     # dot area at 100% expressing
    vmax = float(np.nanmax(res[[f'{g}_mean' for g in genes]].to_numpy()))
    ypos = {s: (n - 1 - i) for i, s in enumerate(subset_order)}   # first subset at the top

    for min_pct_expr in MIN_PCT_EXPRS:
        rt = res.copy()
        for g in genes:
            rt[f'{g}_producer'] = rt[f'{g}_pct'] >= min_pct_expr
        for g in genes:
            prods = rt.loc[rt[f'{g}_producer'], 'subset'].tolist()
            print(f'[{name} over_{int(min_pct_expr)}] {g} producers '
                  f'(>= {min_pct_expr:.0f}%, n={len(prods)}): {prods}')

        csv_dir  = os.path.join(csv_root,  name, f'over_{int(min_pct_expr)}')
        plot_dir = os.path.join(plot_root, name, f'over_{int(min_pct_expr)}')
        os.makedirs(csv_dir,  exist_ok=True)
        os.makedirs(plot_dir, exist_ok=True)
        rt.to_csv(os.path.join(csv_dir, expr_csv), index=False)
        rt.to_csv(os.path.join(csv_dir, prod_csv), index=False)   # read by downstream attribution

        fig, ax = plt.subplots(figsize=(3.0 + 0.4 * len(genes), max(3.0, n * 0.17)))
        sc_h = None
        for xi, g in enumerate(genes):
            sizes = rt[f'{g}_pct'].to_numpy() / 100.0 * SMAX + 3
            ys    = np.array([ypos[s] for s in rt['subset']])
            sc_h  = ax.scatter(np.full(n, xi), ys, s=sizes, c=rt[f'{g}_mean'].to_numpy(),
                               cmap='viridis', vmin=0, vmax=vmax, edgecolors='none', zorder=2)
            prod = rt[f'{g}_producer'].to_numpy()
            if prod.any():
                ax.scatter(np.full(prod.sum(), xi), ys[prod], s=sizes[prod],
                           facecolors='none', edgecolors='#d1352b', linewidths=1.0, zorder=3)
        ax.set_xlim(-0.5, len(genes) - 0.5)
        ax.set_ylim(-0.7, n - 0.3)
        ax.set_xticks(range(len(genes)))
        ax.set_xticklabels(genes, fontsize=9, fontweight='bold')
        ax.set_yticks([ypos[s] for s in subset_order])
        ax.set_yticklabels(subset_order, fontsize=6)
        prod_any = rt.set_index('subset')[[f'{g}_producer' for g in genes]].any(axis=1).to_dict()
        for tick, s in zip(ax.get_yticklabels(), subset_order):
            if prod_any.get(s):
                tick.set_fontweight('bold'); tick.set_color('#d1352b')
        ax.tick_params(length=0)
        for sp_ in ax.spines.values():
            sp_.set_visible(False)
        ax.set_axisbelow(True)
        for s in subset_order:   # lineage color strip left of the labels
            ax.scatter(-0.45, ypos[s], marker='s', s=14,
                       color=lineage_lut.get(subset_lineage.get(s), '#cccccc'),
                       transform=ax.transData, clip_on=False, zorder=4)
        for p in [5, int(min_pct_expr), 50, 100]:
            ax.scatter([], [], s=p / 100.0 * SMAX + 3, color='0.5', edgecolors='none', label=f'{p}%')
        leg = ax.legend(title='% expressing', loc='upper left', bbox_to_anchor=(legend_anchor, 1.0),
                        frameon=False, labelspacing=1.1, fontsize=7, title_fontsize=7.5, borderpad=0.8)
        leg._legend_box.align = 'left'
        cbar = fig.colorbar(sc_h, ax=ax, fraction=0.025, pad=0.02, aspect=30)
        cbar.set_label('mean log-norm', fontsize=7)
        cbar.ax.tick_params(labelsize=6)
        fig.suptitle(f"miniatlas {'/'.join(genes)} producers\n"
                     f"(>= {min_pct_expr:.0f}% of cells expressing = producer, ringed)",
                     y=1.04, fontweight='bold', fontsize=9)
        save_multi_image(os.path.join(plot_dir, pdf_name))
        plt.close('all')
        print(f'  saved {csv_dir} / {plot_dir}')


emit_producers(['CCL19', 'CCL21'], 'ccl_producers',
               'ccl_producer_expression.csv', 'ccl_producers.csv',
               'ccl_producers.pdf', legend_anchor=1.25)
emit_producers(['CXCL13'], 'cxcl13_producers',
               'cxcl13_producer_expression.csv', 'cxcl13_producers.csv',
               'cxcl13_producers.pdf', legend_anchor=1.6)


# ===========================================================================================
# --- Extended Data Fig. 9a: TLS-niche chemokines / receptors / HEV markers across the reference
# Faithful port of the reference per-niche split dotplots. Per-gene min-max-scaled mean
# expression (color), fraction of cells expressing (dot size), swap_axes so genes are rows
# and cell types are columns, a lineage / major-class row-color matrix, and L-R brackets.
# Subset labels are bold if > BOLD_PCT% of cells express ANY of CCL19 / CXCL13 / CCL21.
# ===========================================================================================

use_data = 'log_norm'
level = 'subset'
sorted_cell_types = subset_order

# Major-class categories present (2nd legend), lineage-ordered, Cycling last
major_present = list(dict.fromkeys(classif_dedup.loc[sorted_cell_types, CYCLING_COL]))
if 'Cycling' in major_present:
    major_present.remove('Cycling')
    major_present.append('Cycling')

# Row-color matrix (row 0 = lineage, row 1 = major class), aligned to sorted_cell_types
row_colors_df = pd.DataFrame(
    {
        'lineage': classif_dedup['lineage'].map(lineage_lut),
        'major':   classif_dedup[CYCLING_COL].map(miniatlas_major_lut),
    },
    index=classif_dedup.index,
).loc[sorted_cell_types]

row_colors_df_T = row_colors_df.T
color_matrix = np.empty(row_colors_df_T.shape + (4,), dtype=float)  # RGBA
for row_idx, row in enumerate(row_colors_df_T.columns):
    for col_idx, col in enumerate(row_colors_df_T.index):
        color_matrix[col_idx, row_idx] = mcolors.to_rgba(row_colors_df_T.loc[col, row])

subset_short_lut = dict(zip(classif['miniatlas_subset'], classif['miniatlas_subset_short']))


# --- Genes per niche

PLOT3_GENES = {
    'T cell zone': ['CCL19', 'CCR7', 'TCF7'],
    'B cell zone': ['CXCL13', 'CXCR5'],
    'HEVs':        ['CCL21', 'TSPAN7'],
}
PLOT4_GENES = {
    'HEVs':       ['TSPAN7', 'SELE', 'SELP'],
    'Chemokines': ['CCL19', 'CCL21'],
}


def present_only(genes_dict):
    """Drop genes absent from the reference; drop niches left empty."""
    out = {}
    for niche, genes in genes_dict.items():
        keep = [g for g in genes if g in adata.var_names]
        dropped = [g for g in genes if g not in adata.var_names]
        if dropped:
            print(f'WARNING: {niche}: genes absent from reference (dropped): {dropped}')
        if keep:
            out[niche] = keep
    return out


PLOT3_GENES = present_only(PLOT3_GENES)
PLOT4_GENES = present_only(PLOT4_GENES)


# --- Per-cell-type expression tables (all niche genes): % expressing, mean log-norm, scaled mean
# scaled_df reproduces scanpy's standard_scale='var'; we recolor the dots with it because this
# repo's scanpy mis-renders that scaling under swap_axes (colors per column instead of per gene).

all_niche_genes = list(dict.fromkeys(
    g for d in (PLOT3_GENES, PLOT4_GENES) for gs in d.values() for g in gs))
gi = [adata.var_names.get_loc(g) for g in all_niche_genes]
counts = adata.layers['counts'][:, gi]
logn   = adata.layers['log_norm'][:, gi]
counts = counts.toarray() if sp.issparse(counts) else np.asarray(counts)
logn   = logn.toarray()   if sp.issparse(logn)   else np.asarray(logn)

sub_arr = adata.obs[SUBSET_COL].to_numpy()
mean_df = pd.DataFrame(index=sorted_cell_types, columns=all_niche_genes, dtype=float)
pct_df  = pd.DataFrame(index=sorted_cell_types, columns=all_niche_genes, dtype=float)
for ct in sorted_cell_types:
    m = (sub_arr == ct)
    mean_df.loc[ct] = logn[m].mean(axis=0)
    pct_df.loc[ct]  = (counts[m] > 0).mean(axis=0) * 100.0

denom = (mean_df.max() - mean_df.min()).replace(0, np.nan)
scaled_df = ((mean_df - mean_df.min()) / denom).fillna(0.0)

# Cell types to bold: > BOLD_PCT % expressing ANY of the bold genes
bold_genes_present = [g for g in BOLD_GENES if g in pct_df.columns]
missing_bold = [g for g in BOLD_GENES if g not in pct_df.columns]
if missing_bold:
    print(f'WARNING: bold genes absent from reference (ignored for bolding): {missing_bold}')
bold_mask = (pct_df[bold_genes_present] > BOLD_PCT).any(axis=1)
bold_set = set(bold_mask[bold_mask].index)
print(f'Extended Data Fig. 9a: cell types bolded (> {BOLD_PCT:.0f}% expressing any of {bold_genes_present}, '
      f'n={len(bold_set)}):')
for s in sorted_cell_types:
    if s in bold_set:
        print(f'  {s}: ' + ', '.join(f'{g}={pct_df.loc[s, g]:.1f}%' for g in bold_genes_present))


# --- Functions for plotting

def lighten_color(color, amount=0.5):
    r, g, b = to_rgb(color)
    return to_hex([1 - (1 - c) * (1 - amount) for c in (r, g, b)])


def fix_dot_colors(main_ax, genes):
    """Recolor a scanpy dotplot's dots per gene (min-max scaled mean), overwriting the
    per-column coloring that this repo's scanpy applies under swap_axes. Dot positions and
    sizes are left untouched; each dot is matched to its (cell type, gene) by nearest tick."""
    xticks = np.asarray(main_ax.get_xticks())
    yticks = np.asarray(main_ax.get_yticks())
    xlabs  = [t.get_text() for t in main_ax.get_xticklabels()]
    ylabs  = [t.get_text() for t in main_ax.get_yticklabels()]
    coll = max(main_ax.collections, key=lambda c: len(c.get_offsets()))
    offs = np.asarray(coll.get_offsets())
    vals = np.empty(len(offs))
    for k, (ox, oy) in enumerate(offs):
        ct = xlabs[int(np.argmin(np.abs(xticks - ox)))]
        gn = ylabs[int(np.argmin(np.abs(yticks - oy)))]
        vals[k] = scaled_df.loc[ct, gn] if (ct in scaled_df.index and gn in scaled_df.columns) else 0.0
    coll.set_array(vals)
    coll.set_clim(0, 1)
    coll.set_cmap(cmap)


def join_lr_multi(axs, lr_pairs, offset=0.3, offset_start=0.2):
    """Draw brackets for multiple ligand-receptor pairs, all using the same furthest x_left."""
    if not isinstance(axs, (list, tuple, np.ndarray)):
        axs = [axs]
    label_to_ax, label_to_idx, label_to_bbox, ax_inv = {}, {}, {}, {}
    for ax in axs:
        inv = ax.transData.inverted()
        ax_inv[ax] = inv
        for idx, label in enumerate(ax.get_yticklabels()):
            label_to_ax[label.get_text()] = ax
            label_to_idx[label.get_text()] = idx
            label_to_bbox[label.get_text()] = label.get_window_extent()
    x_lefts, label_objs = [], []
    for lr1, lr2 in lr_pairs:
        ax_lig = label_to_ax[lr1]
        bbox_lig = label_to_bbox[lr1]
        x_start = ax_inv[ax_lig].transform((bbox_lig.x0, bbox_lig.y0))[0] - offset
        x_lefts.append(x_start - offset_start)
        label_objs.append((lr1, lr2))
    x_left_global = min(x_lefts)
    for lr1, lr2 in label_objs:
        ax_lig, ax_rec = label_to_ax[lr1], label_to_ax[lr2]
        idx_lig, idx_rec = label_to_idx[lr1], label_to_idx[lr2]
        bbox_lig, bbox_rec = label_to_bbox[lr1], label_to_bbox[lr2]
        x_start = ax_inv[ax_lig].transform((bbox_lig.x0, bbox_lig.y0))[0] - offset
        x_right = ax_inv[ax_rec].transform((bbox_rec.x0, bbox_rec.y0))[0] - offset
        y_lig_data, y_rec_data = idx_lig + 0.5, idx_rec + 0.5
        ax_lig.plot([x_start, x_left_global, x_left_global, x_right],
                    [y_lig_data, y_lig_data, y_rec_data, y_rec_data],
                    color='black', linewidth=2, zorder=10, clip_on=False, solid_joinstyle='miter')
    return x_left_global


def auto_resize_dotplot(ax, space_per_dot=0.1):
    pos = ax.get_position()
    xticks = len(ax.get_xticks())
    yticks = len(ax.get_yticks())
    ax.set_position([pos.x0, pos.y0, xticks * space_per_dot, yticks * space_per_dot])


def reposition_dotplot_axs(axs, y_offset=0.1):
    ax_pos_prev = None
    for n, ax in enumerate(axs):
        pos = ax.get_position()
        ypos = pos.y0 if n == 0 else ax_pos_prev.y0 - pos.height - y_offset
        ax.set_position([pos.x0, ypos, pos.width, pos.height])
        ax_pos_prev = ax.get_position()
        if n < len(axs) - 1:
            ax.set_xticks([])


def position_legend(ax):
    fig = ax.figure
    fig.canvas.draw()
    bbox = ax.get_tightbbox(fig.canvas.get_renderer())
    bbox_fig = fig.transFigure.inverted().transform(bbox)
    return bbox_fig[0][1], bbox_fig[1][1], bbox_fig[0][0], bbox_fig[1][0]


def add_color_legend(fig, labels, color_dict, title='', x=0, y=0, offset=0.01,
                     spacing_per_label=0.05, loc='center left', fontsize=10):
    handles = [
        mpatches.Patch(facecolor=color_dict[label], label=label, edgecolor='black', linewidth=0.8)
        for label in labels if label in color_dict
    ]
    legend = fig.legend(
        handles=handles, title=title, loc=loc, bbox_to_anchor=(x + offset, y),
        borderaxespad=0, frameon=False, handlelength=0.8,
        bbox_transform=fig.transFigure, fontsize=fontsize,
    )
    legend.get_title().set_fontweight('bold')
    legend.get_title().set_fontsize(fontsize)
    legend._legend_box.align = "left"
    return legend


def bold_producer_labels(heatmap_ax):
    """Rename x-tick labels to short names; bold those in bold_set (the producers)."""
    orig_texts = [lbl.get_text() for lbl in heatmap_ax.get_xticklabels()]
    heatmap_ax.set_xticklabels([subset_short_lut.get(t, t) for t in orig_texts])
    for tick, orig in zip(heatmap_ax.get_xticklabels(), orig_texts):
        if orig in bold_set:
            tick.set_fontweight('bold')


# --- Shared niche colors (rectangles on the left of the dot grids)

color_dict_score = {
    'B cell zone': mpl.colors.cnames["cornflowerblue"],
    'T cell zone': mpl.colors.cnames["crimson"],
    'HEVs':        mpl.colors.cnames["indigo"],
}
for key in color_dict_score:
    color_dict_score[key] = lighten_color(color_dict_score[key], amount=0.3)
color_dict_score['Chemokines'] = lighten_color(mpl.colors.cnames["cornflowerblue"], amount=0.3)


# --- Shared dot size / color scale

dot_max, vmin, vmax, cmap = 1, 0, 1, 'Reds'
space_per_dot = 0.04
sc.set_figure_params(fontsize=11)
plt.rcParams['axes.grid'] = False


def make_niche_dotplot(top_genes_dict, lr_pairs, x_left_global=None):
    """Reproduce the per-niche split dotplot: one swap_axes dotplot per niche, stacked, with
    the row-color matrix, ligand-receptor brackets, niche rectangles, legends and bold
    producer labels. Returns x_left_global for reuse (Plot 4 reuses Plot 3's)."""
    n_keys = len(top_genes_dict)
    fig, axs = plt.subplots(ncols=1, nrows=n_keys, figsize=(5, 5))
    axs = np.atleast_1d(axs)

    heatmap_axs = []
    size_ax = cbar_ax = None
    for i, (zone, genes) in enumerate(top_genes_dict.items()):
        dp = sc.pl.dotplot(
            adata, var_names=genes, groupby=SUBSET_COL, layer=use_data, ax=axs[i],
            dot_max=dot_max, cmap=cmap, vmin=vmin, vmax=vmax, show=False,
            swap_axes=True, standard_scale='var',
        )
        fix_dot_colors(dp['mainplot_ax'], genes)   # correct the per-column mis-coloring
        auto_resize_dotplot(dp['mainplot_ax'], space_per_dot=space_per_dot)
        heatmap_axs.append(dp['mainplot_ax'])
        buffer = 0.25
        dp['mainplot_ax'].set_ylim(len(dp['mainplot_ax'].get_yticks()) + buffer, -buffer)
        if i < n_keys - 1:
            for name, dp_ax in dp.items():
                if name != 'mainplot_ax':
                    dp_ax.set_visible(False)
        else:
            size_ax = dp['size_legend_ax']
            cbar_ax = dp['color_legend_ax']

    reposition_dotplot_axs(heatmap_axs, y_offset=0.01)

    # Colorbar + size legend to the right of the last dot grid
    ax_pos = heatmap_axs[-1].get_position()
    cbar_height = space_per_dot * 0.8
    cbar_width = space_per_dot * 5
    cbar_ax.set_position([ax_pos.x1 + 0.03, ax_pos.y0, cbar_width, cbar_height])
    cbar_pos = cbar_ax.get_position()
    size_ax.set_position([cbar_pos.x0, cbar_pos.y1 + 0.12, cbar_width, cbar_height * 3])

    # Ligand-receptor brackets (only pairs whose both genes are plotted)
    plotted = {g for gs in top_genes_dict.values() for g in gs}
    lr_use = [(a, b) for (a, b) in lr_pairs if a in plotted and b in plotted]
    if lr_use:
        x_left_global = join_lr_multi(heatmap_axs, lr_use, offset=0.25)

    # Row-color matrix above the first dot grid
    row_colors_ax = plt.axes([0.1, 0.1, 0.2, 0.2])
    new_heatmap_pos = heatmap_axs[0].get_position()
    presence_ax_height = space_per_dot * len(row_colors_df_T.index)
    row_colors_ax.set_position([new_heatmap_pos.x0, new_heatmap_pos.y1 + 0.01,
                                new_heatmap_pos.width, presence_ax_height])
    row_colors_ax.imshow(color_matrix, aspect="auto")
    row_colors_ax.set_xticks([])
    row_colors_ax.set_xticks(np.arange(-0.5, len(row_colors_df_T.columns), 1), minor=True)
    row_colors_ax.set_yticks(np.arange(-0.5, len(row_colors_df_T.index), 1), minor=True)
    row_colors_ax.grid(which="minor", color="white", linestyle='-', linewidth=1)
    row_colors_ax.tick_params(which="minor", bottom=False, left=False)
    row_colors_ax.set_yticks(np.arange(len(row_colors_df_T.index)))
    row_colors_ax.set_yticklabels(['Lineage', 'Major class'], fontsize=9)
    row_colors_ax.grid(False)

    # Niche rectangles + vertical labels down the left of each dot grid
    for n, ax in enumerate(heatmap_axs):
        label = list(top_genes_dict.keys())[n]
        face_color = color_dict_score[label]
        x_fig, _ = fig.transFigure.inverted().transform(ax.transData.transform((x_left_global, 0)))
        pos = ax.get_position()
        y0, y1 = pos.y0, pos.y1
        height = abs(y1 - y0)
        width = 0.02
        x1 = x_fig - 0.02
        fig.patches.append(mpl.patches.Rectangle(
            (x1 - width, y0), width, height, transform=fig.transFigure,
            facecolor=face_color, edgecolor='black', linewidth=0, alpha=1))
        label_txt = label.replace(' zone', '\nzone') if 'zone' in label else label
        fig.text(x1 - width - 0.01, (y0 + y1) / 2, label_txt,
                 ha='right', va='center', color='black')

    # Lineage + major-class legends, right of the size legend
    _, top_edge, left_edge, right_edge = position_legend(size_ax)
    legend1 = add_color_legend(fig, labels=lineage_lut.keys(), color_dict=lineage_lut,
                               title='Lineage', x=left_edge, y=top_edge, offset=0,
                               spacing_per_label=0.05, loc='lower left', fontsize=9)
    fig.canvas.draw()
    legend_bbox_fig = fig.transFigure.inverted().transform(
        legend1.get_window_extent(fig.canvas.get_renderer()))
    add_color_legend(fig, labels=major_present, color_dict=miniatlas_major_lut,
                     title='Major class', x=right_edge, y=legend_bbox_fig[1][1], offset=0.02,
                     spacing_per_label=0.05, loc='upper left', fontsize=9)

    # Short + bold cell-type labels on the last (only labeled) dot grid
    bold_producer_labels(heatmap_axs[-1])

    row_colors_ax_pos = row_colors_ax.get_position()
    plt.suptitle(f'Extended Data Fig. 9a — miniatlas {level} cell types (use_data={use_data})\n'
                 f'bold = > {BOLD_PCT:.0f}% expressing any of ' + '/'.join(bold_genes_present),
                 x=(row_colors_ax_pos.x0 + row_colors_ax_pos.x1) / 2,
                 y=row_colors_ax_pos.y1 + 0.1, fontweight='bold')
    return x_left_global


# --- Build the Extended Data Fig. 9a niche dotplot

make_niche_dotplot(PLOT3_GENES, [('CCL19', 'CCR7'), ('CXCL13', 'CXCR5')])


# --- Save Extended Data Fig. 9a outputs

niche_csv_dir  = os.path.join(csv_root,  'niche_genes')
niche_plot_dir = os.path.join(plot_root, 'niche_genes')
os.makedirs(niche_csv_dir,  exist_ok=True)
os.makedirs(niche_plot_dir, exist_ok=True)

expr = pct_df.add_suffix('_pct').join(mean_df.add_suffix('_mean')).join(scaled_df.add_suffix('_scaled'))
expr['bold_producer'] = bold_mask.reindex(sorted_cell_types).values
expr.to_csv(os.path.join(niche_csv_dir, 'niche_gene_pct_expressing.csv'))

for fnum in plt.get_fignums():
    for ax in plt.figure(fnum).axes:
        ax.grid(False)

save_multi_image(os.path.join(niche_plot_dir, 'ED9A_miniatlas_niche_genes_dotplots.pdf'))
plt.close('all')

print(f'Saved Extended Data Fig. 9a: {niche_plot_dir}')
print("Script successfully completed")
