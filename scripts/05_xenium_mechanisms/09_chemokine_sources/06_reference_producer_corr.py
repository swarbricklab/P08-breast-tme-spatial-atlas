#!/usr/bin/env python3

# --- Aim
# Extended Data Fig. 9i: within-TLS spatial correlation of reference-producer cell abundance vs CCL19 / CCL21
# transcript density, per 50 um bin. The candidate cell types are the scRNA-seq reference
# CCL19 / CCL21 producers (subsets with >= 10 % of cells expressing, the `over_10` list from
# 01_reference_producers). For every 50 um bin we count each producer subset's cells and sum
# the per-cell CCL19 / CCL21 transcript counts, then correlate producer-cell density against
# chemokine density. Two coefficients per (sample, cell type, gene): raw Spearman, and a
# partial Spearman that controls for the total number of cells per bin (so a subset must track
# a chemokine beyond the generic "dense bins have more of everything" effect). Per-sample
# coefficients are tested against 0 with a two-sided Wilcoxon signed-rank (n = 6). The headline
# panel (Extended Data Fig. 9i) is the T-cell-zone partial Spearman; the all-TLS and B-cell-zone strata are shown
# alongside. Reference producers of a given gene get a bold red x-tick label and a red outline
# on that gene's % expressing annotation cell. Read-only: per-cell transcript counts, cell
# centroids, TLS zone bins, cell-type labels and the reference producer list.

# --- Define environment
# env: xen-hev

# --- Load packages

import argparse
import os
import sys

import numpy as np
import pandas as pd
import anndata as ad
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
from scipy.stats import rankdata, spearmanr, wilcoxon

mpl.rcParams['pdf.fonttype'] = 42
plt.rcParams['axes.grid'] = False

TICK_FS, SUP_FS, ANNO_FS = 8, 12, 7
PRODUCER_C = '#d1352b'   # highlight colour for reference producers (matches 01_reference_producers dot plot)


# --- Argparse arguments
# Cohort: dispatches through submit_multi_sample_job_nci.sh -> run_multi_sample_job.sh, which
# passes --h5ad_dir --project_dir --res_dir --sample_sheet --samples. The Methods parameters
# (50 um grid, over_10 producer set, min 10 bins per stratum) are hard-coded defaults.

parser = argparse.ArgumentParser(
    description="Extended Data Fig. 9i: reference-producer abundance vs CCL19/CCL21 density correlation per 50um bin (cohort)",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("-s", "--samples", nargs='+', required=True)
parser.add_argument("--h5ad_dir", default=None, help="defaults to <project_dir>/data/anndata_objects")
parser.add_argument("--project_dir", required=True)
parser.add_argument("--res_dir", required=True)
parser.add_argument("--sample_sheet", default=None, help="Path to sample sheet CSV (else default)")
parser.add_argument("--grid_um", type=float, default=50.0)
parser.add_argument("--min_bins", type=int, default=10, help="skip a (sample, stratum) with fewer bins")
parser.add_argument("--producer_tag", default="over_10",
                    help="reference-producer threshold subdir to read producer lists from + tag outputs")
parser.add_argument("--genes", nargs='+', default=['CCL19', 'CCL21'],
                    help="chemokines to include (one panel each); the candidate set is the union of "
                         "these genes' reference producers")

args = vars(parser.parse_args())
samples      = args["samples"]
project_dir  = args["project_dir"]
res_dir      = args["res_dir"]
h5ad_dir     = args["h5ad_dir"] or os.path.join(project_dir, 'data', 'anndata_objects')
grid_um      = args["grid_um"]
min_bins     = args["min_bins"]
producer_tag = args["producer_tag"]
GENES        = list(args["genes"])

if args.get("sample_sheet"):
    sample_sheet_path = args["sample_sheet"]
else:
    sample_sheet_path = os.path.join(project_dir, 'config', 'sample_sheet_xenium.csv')

subtype_order = ['Luminal A', 'Luminal B', 'Luminal B (HER2+)', 'HER2+', 'TNBC']
LINEAGE_ORDER = ['Epithelial', 'Immune', 'Stromal']
genes_tag = '' if set(GENES) == {'CCL19', 'CCL21'} else '_' + '_'.join(GENES)


# --- Local plotting helper

sys.path.append(os.path.join(project_dir, 'scripts'))
from helper_functions import save_multi_image   # noqa: E402


# --- Reference cell-type classifications: short display name + lineage + major type per subset
# x-tick labels and every internal cell key use the classifications' `miniatlas_subset_short`.

_classif = pd.read_csv(os.path.join(project_dir, 'resources', 'cell_types', 'miniatlas',
                                    'miniatlas_subset_classifications.csv'))
_classif = _classif.drop_duplicates(subset='miniatlas_subset', keep='first')
_sub2short = dict(zip(_classif['miniatlas_subset'], _classif['miniatlas_subset_short']))
_sub2lin   = dict(zip(_classif['miniatlas_subset'], _classif['lineage']))
_sub2maj   = dict(zip(_classif['miniatlas_subset'], _classif['miniatlas_major']))
disp = lambda sub: _sub2short.get(sub, sub)
_lin_lut = dict(pd.read_csv(os.path.join(project_dir, 'resources', 'plotting_colors',
                                         'miniatlas_lineage_colors.csv'))[['lineage', 'color']].values)
_maj_lut = dict(pd.read_csv(os.path.join(project_dir, 'resources', 'plotting_colors',
                                         'miniatlas_major_colors.csv'))[['miniatlas_major', 'color']].values)


# --- Reference-grounded producer sets (per gene) from 01_reference_producers

_prod_csv = os.path.join(res_dir, '05_xenium_mechanisms', '09_chemokine_sources',
                         '01_reference_producers', 'csv_files', 'ccl_producers',
                         producer_tag, 'ccl_producers.csv')
_prod_tbl = pd.read_csv(_prod_csv).drop_duplicates(subset='subset')
PRODUCER_SETS = {g: _prod_tbl.loc[_prod_tbl[f'{g}_producer'], 'subset'].tolist() for g in GENES}
_union = sorted(set().union(*(set(PRODUCER_SETS[g]) for g in GENES)))   # producers of ANY included gene


# order cell types by lineage (Epithelial, Immune, Stromal) then major type then subset — an
# interpretive grouping (per global convention: derive the sort from the classifications file)
def _lineage_sort_key(full):
    lin = _sub2lin.get(full, 'NA')
    li = LINEAGE_ORDER.index(lin) if lin in LINEAGE_ORDER else len(LINEAGE_ORDER)
    return (li, str(_sub2maj.get(full, 'zzz')), full)
_cells_full = sorted(_union, key=_lineage_sort_key)

SUBSETS       = {full: disp(full) for full in _cells_full}     # full reference name -> short display
CELL_ORDER    = [disp(full) for full in _cells_full]
SHORT_TO_FULL = {short: full for full, short in SUBSETS.items()}
# per-(display cell, gene) reference-producer flag
PRODUCER_FLAG = {(disp(full), g): (full in PRODUCER_SETS[g]) for full in _cells_full for g in GENES}
PROD_LABEL = f"≥{producer_tag.replace('over_', '')}%"

# per-cell lineage / major-type colours (for the annotation strips + legends)
CELL_LINEAGE = {ct: _sub2lin.get(SHORT_TO_FULL[ct], 'NA') for ct in CELL_ORDER}
CELL_MAJOR   = {ct: _sub2maj.get(SHORT_TO_FULL[ct], 'NA') for ct in CELL_ORDER}
LINEAGE_C    = {ct: _lin_lut.get(CELL_LINEAGE[ct], '#cccccc') for ct in CELL_ORDER}
MAJOR_C      = {ct: _maj_lut.get(CELL_MAJOR[ct], '#cccccc') for ct in CELL_ORDER}

# per-producer % of cells expressing each chemokine in the reference (annotation row + shared cbar)
_pct_tbl = _prod_tbl.set_index('subset')
PCT_EXPR = {disp(full): {g: (float(_pct_tbl.loc[full, f'{g}_pct']) if full in _pct_tbl.index else np.nan)
                         for g in GENES} for full in _cells_full}
PCT_VMAX = float(np.nanmax([PCT_EXPR[ct][g] for ct in CELL_ORDER for g in GENES]))

print(f'Reference producers ({producer_tag}):')
for g in GENES:
    print(f"  {g}: {[disp(s) for s in PRODUCER_SETS[g]]}")
print(f'Cell types (x-axis, lineage-ordered): {CELL_ORDER}')

# standard transcript colours from the project palette
_tx_lut = dict(pd.read_csv(os.path.join(project_dir, 'resources', 'plotting_colors',
                                        'transcript_colors.csv'))[['transcript', 'color']].values)
GENE_C = {g: _tx_lut.get(g, '#888888') for g in GENES}

BIN_TAG = f'{int(grid_um)}µm bins · in TLS'


def suptitle(fig, text, y=1.05):
    fig.suptitle(text, x=0.5, y=y, ha='center', va='bottom', fontweight='bold', fontsize=SUP_FS)


def star(p):
    if not np.isfinite(p): return 'ns'
    if p < 1e-3: return '***'
    if p < 1e-2: return '**'
    if p < 5e-2: return '*'
    return 'ns'


def boxplot_panel(ax, values_per_x, colors, samples_per_x, box_w):
    """Colour-by-x boxes at integer positions (width box_w) + per-sample jittered dots and
    connecting lines. Caller owns xlim / xticks / strips (positions are 0..n-1)."""
    n = len(values_per_x)
    bp = ax.boxplot([np.asarray(v) for v in values_per_x],
                    positions=list(range(n)), widths=box_w,
                    patch_artist=True, showfliers=False,
                    medianprops=dict(color='black', linewidth=1.2), zorder=2)
    for patch, color in zip(bp['boxes'], colors):
        rgb = mpl.colors.to_rgb(color)
        patch.set_facecolor((*rgb, 0.55)); patch.set_edgecolor(color); patch.set_linewidth(0.8)
    for elem in (*bp['whiskers'], *bp['caps']):
        elem.set_color('black'); elem.set_linewidth(0.8)
    rng = np.random.default_rng(seed=42)
    jit = 0.14 * box_w
    all_samples = sorted({s for sl in samples_per_x for s in sl})
    jitter = {(s, xi): float(rng.uniform(-jit, jit)) for xi in range(n) for s in all_samples}
    for xi, (sl, vl) in enumerate(zip(samples_per_x, values_per_x)):
        xs = [xi + jitter[(s, xi)] for s in sl]
        ax.scatter(xs, vl, c='black', s=12, linewidths=0, zorder=3)
    ax.tick_params(axis='y', labelsize=TICK_FS); ax.grid(False)
    for sp_ in ('top', 'right'):
        ax.spines[sp_].set_visible(False)


def partial_spearman(x, y, z):
    """Spearman partial correlation of x,y controlling for z (rank residuals)."""
    x, y, z = np.asarray(x, float), np.asarray(y, float), np.asarray(z, float)
    if len(x) < 4:
        return np.nan
    rx, ry, rz = rankdata(x), rankdata(y), rankdata(z)
    A = np.vstack([rz, np.ones_like(rz)]).T
    def resid(a):
        coef, *_ = np.linalg.lstsq(A, a, rcond=None)
        return a - A @ coef
    ex, ey = resid(rx), resid(ry)
    if ex.std() == 0 or ey.std() == 0:
        return np.nan
    return float(np.corrcoef(ex, ey)[0, 1])


# --- Metadata

sample_sheet    = pd.read_csv(sample_sheet_path)
tissue_metadata = pd.read_csv(os.path.join(project_dir, 'resources', 'tissue_metadata.csv'))
merged          = sample_sheet.merge(tissue_metadata, on='tissue_id')
subtype_dict    = dict(zip(merged['section_id'], merged['clinical_subtype']))

zone_dir = os.path.join(res_dir, '05_xenium_mechanisms', '07_TLS_zonation',
                        '01_classify_zones', 'csv_files')
tx_dir   = os.path.join(res_dir, '05_xenium_mechanisms', '09_chemokine_sources',
                        '02_extract_transcripts', 'csv_files')


def celltype_csv(sample, level='subset'):
    """Label-transfer predicted cell types (col 'pred_cell_type')."""
    p = os.path.join(res_dir, '02_cell_type_mapping', 'label_transfer', 'xenium',
                     f'{level}_level', f'{sample}_{level}_pred_cell_types.csv')
    return p


# --- Per-sample: build per-bin table, compute correlations

corr_rows = []      # per (sample, stratum, cell_type, gene): raw + partial r
bin_tables = {}     # (sample) -> per-bin df (for scatter examples)
for s in samples:
    binf = pd.read_csv(os.path.join(zone_dir, f'{s}_BT_zone_bins.csv'))
    if len(binf) == 0:
        print(f'{s}: no TLS zone bins — skipping')
        continue
    # recover the 50um grid origin from the bin-centre coords (bin_x = xmin + (bin_ix+0.5)*grid)
    xmin = float(binf['bin_x'].iloc[0] - (binf['bin_ix'].iloc[0] + 0.5) * grid_um)
    ymin = float(binf['bin_y'].iloc[0] - (binf['bin_iy'].iloc[0] + 0.5) * grid_um)
    bin_zone = binf.set_index(['bin_ix', 'bin_iy'])['TLS_zone']
    bin_zone.index.names = ['ix', 'iy']   # match the cell bin indices

    # cells -> bin (micron centroids; obsm['spatial'] is image space in this project)
    adata = ad.read_h5ad(os.path.join(h5ad_dir, 'xenium', f'{s}.h5ad'), backed='r')
    cxy = adata.obs[['x_centroid', 'y_centroid']].copy()
    cxy.columns = ['x', 'y']
    cxy.index = cxy.index.astype(str)
    adata.file.close()

    sub = pd.read_csv(celltype_csv(s), index_col=0)
    sub.index = sub.index.astype(str)
    cxy['subset'] = sub['pred_cell_type']
    cxy['ix'] = np.floor((cxy['x'] - xmin) / grid_um).astype(int)
    cxy['iy'] = np.floor((cxy['y'] - ymin) / grid_um).astype(int)
    cxy['short'] = cxy['subset'].map(SUBSETS)
    n_total = cxy.groupby(['ix', 'iy']).size().rename('n_total_cells')
    cell_counts = (cxy.dropna(subset=['short']).groupby(['ix', 'iy', 'short']).size().unstack(fill_value=0)
                   .reindex(columns=CELL_ORDER, fill_value=0))

    # per-cell chemokine transcript counts -> bin (integer counts from 02_extract_transcripts)
    for gene in GENES:
        tx_path = os.path.join(tx_dir, f'{s}_{gene}_per_cell.csv')
        if os.path.isfile(tx_path):
            gser = pd.read_csv(tx_path).set_index('cell_id')[f'n_{gene}_transcripts']
            gser.index = gser.index.astype(str)
            cxy[gene] = cxy.index.map(gser).fillna(0).astype(float)
        else:
            print(f'  WARNING: missing {tx_path} — {gene} density set to 0 for {s}')
            cxy[gene] = 0.0
    tx_counts = cxy.groupby(['ix', 'iy'])[GENES].sum()

    # restrict to valid TLS bins (those in the zone-bin table)
    tab = pd.DataFrame(index=bin_zone.index)
    tab['zone'] = bin_zone
    tab = tab.join(cell_counts).join(tx_counts).join(n_total)
    tab = tab.fillna(0)

    bin_tables[s] = tab
    for stratum, mask in [('all TLS', tab['zone'].notna()),
                          ('T cell zone', tab['zone'] == 'T_cell_zone'),
                          ('B cell zone', tab['zone'] == 'B_cell_zone')]:
        t = tab[mask]
        if len(t) < min_bins:
            continue
        for ct in CELL_ORDER:
            for gene in GENES:
                raw = spearmanr(t[ct], t[gene]).correlation if t[ct].std() and t[gene].std() else np.nan
                par = partial_spearman(t[ct], t[gene], t['n_total_cells'])
                corr_rows.append({'sample_id': s, 'subtype': subtype_dict.get(s), 'stratum': stratum,
                                  'cell_type': ct, 'gene': gene, 'is_producer': PRODUCER_FLAG[(ct, gene)],
                                  'n_bins': len(t), 'spearman': raw, 'partial_spearman': par})
    print(f'{s}: done ({len(tab)} TLS bins)')

cdf = pd.DataFrame(corr_rows)


# --- Correlation boxplots: one panel per gene (shared y), x = candidate producers ordered by
# lineage then major type. Boxes in the gene's transcript colour, per-sample dots + lines,
# one-sample Wilcoxon (r != 0) star per box. Above each panel a 3-row colour grid annotates each
# producer — Lineage (top), Major type, then % of cells expressing that chemokine in the reference
# (bottom, shared vmin/vmax across genes, own colorbar). Reference producers of the panel's gene
# get a bold-red x-tick label. All text one size (BODY_FS) except gene + main titles.

BODY_FS  = 8                         # every non-title label uses this size
TITLE_FS = 11                        # gene panel titles + main title (the only larger text)
BOX_W    = 0.82                      # box / annotation-cell width in x-category units (fat boxes)
X_PAD    = BOX_W / 2 + 0.35          # buffer between the outer boxes and the panel edge (x-units)
PER_CAT  = 0.30                      # inches of panel width per producer
PCT_CMAP = 'BuPu'
_pct_norm = mpl.colors.Normalize(vmin=0, vmax=PCT_VMAX)


def corr_boxplot_figure(stratum, metric, title):
    sub = cdf[cdf['stratum'] == stratum]
    if sub.empty:
        return
    n = len(CELL_ORDER)
    vmn, vmx = sub[metric].min(), sub[metric].max()
    ymin = min(0, vmn if np.isfinite(vmn) else 0)
    ymax = max(0, vmx if np.isfinite(vmx) else 0)
    pad = 0.16 * ((ymax - ymin) or 1)
    cmap = plt.get_cmap(PCT_CMAP)

    # --- geometry (inches). Square annotation cells with equal H/V gaps: one x-category unit maps
    # to `xu` inches; a cell is BOX_W units wide/tall, the gap is (1-BOX_W) units in both directions.
    panel_w = n * PER_CAT
    xrange  = (n - 1) + 2 * X_PAD
    xu      = panel_w / xrange
    h_y, gap_y = 0.40 * BOX_W, 0.5 * (1 - BOX_W)   # short rects; small vertical gap
    pitch      = h_y + gap_y
    YTOP       = 3 * h_y + 2 * gap_y             # strip ax data ylim so 1 y-unit == xu inches
    strip_H    = YTOP * xu
    MAIN_H, SGAP, TGAP, NOTE_H = 1.35, 0.06, 0.16, 0.74
    BOT, ML, PGAP, RIGHT = 1.5, 1.15, 0.08, 2.2
    n_g = len(GENES)
    panels_w = n_g * panel_w + (n_g - 1) * PGAP            # total width spanned by the gene panels
    fig_w = ML + panels_w + RIGHT
    fig_h = BOT + MAIN_H + SGAP + strip_H + TGAP + NOTE_H
    strip_bottom = BOT + MAIN_H + SGAP

    def _pos(x, y, w, h):
        return [x / fig_w, y / fig_h, w / fig_w, h / fig_h]

    fig = plt.figure(figsize=(fig_w, fig_h))
    for ci, gene in enumerate(GENES):
        pl = ML + ci * (panel_w + PGAP)
        ax_main  = fig.add_axes(_pos(pl, BOT, panel_w, MAIN_H))
        ax_strip = fig.add_axes(_pos(pl, strip_bottom, panel_w, strip_H))

        vals, samps = [], []
        for ct in CELL_ORDER:
            sd = sub[(sub['cell_type'] == ct) & (sub['gene'] == gene)].dropna(subset=[metric])
            vals.append(sd[metric].to_numpy()); samps.append(sd['sample_id'].tolist())
        boxplot_panel(ax_main, vals, [GENE_C[gene]] * n, samps, BOX_W)
        ax_main.axhline(0, color='0.6', lw=0.6, zorder=0)
        ax_main.set_ylim(ymin - pad, ymax + pad)
        ax_main.set_xlim(-X_PAD, (n - 1) + X_PAD)
        ax_main.set_xticks(range(n))
        ax_main.set_xticklabels(CELL_ORDER, rotation=90, fontsize=BODY_FS)
        ax_main.tick_params(axis='y', labelsize=BODY_FS)
        for tick, ct in zip(ax_main.get_xticklabels(), CELL_ORDER):      # bold-red producers of this gene
            if PRODUCER_FLAG[(ct, gene)]:
                tick.set_fontweight('bold'); tick.set_color(PRODUCER_C)
        for xi, v in enumerate(vals):                                    # Wilcoxon r != 0 stars
            vv = v[np.isfinite(v)]
            p = np.nan
            if len(vv) >= 5 and np.any(vv != 0):
                try:
                    p = wilcoxon(vv).pvalue
                except ValueError:
                    p = np.nan
            ax_main.text(xi, ymax + pad * 0.35, star(p), ha='center', va='bottom', fontsize=BODY_FS)
        if ci == 0:
            ax_main.set_ylabel(metric.replace('_', ' ').capitalize() + '\n(per bin)',
                               fontsize=BODY_FS, fontweight='bold')
        else:
            ax_main.tick_params(labelleft=False)                          # shared y -> drop right labels

        # 3-row annotation grid (short rects): Lineage (top), Major type, % expressing (bottom)
        for xi, ct in enumerate(CELL_ORDER):
            x0 = xi - BOX_W / 2
            ax_strip.add_patch(mpatches.Rectangle((x0, 2 * pitch), BOX_W, h_y,
                                                  facecolor=LINEAGE_C[ct], edgecolor='black', linewidth=0.5))
            ax_strip.add_patch(mpatches.Rectangle((x0, pitch), BOX_W, h_y,
                                                  facecolor=MAJOR_C[ct], edgecolor='black', linewidth=0.5))
            pv = PCT_EXPR[ct][gene]; _pe = PRODUCER_FLAG[(ct, gene)]   # producer -> red outline (like x-ticks)
            ax_strip.add_patch(mpatches.Rectangle((x0, 0.0), BOX_W, h_y,
                                                  facecolor=(cmap(_pct_norm(pv)) if np.isfinite(pv) else '#dddddd'),
                                                  edgecolor=(PRODUCER_C if _pe else 'black'),
                                                  linewidth=(1.1 if _pe else 0.5), zorder=(3 if _pe else 2)))
        ax_strip.set_xlim(-X_PAD, (n - 1) + X_PAD)
        ax_strip.set_ylim(0, YTOP); ax_strip.set_yticks([])
        ax_strip.tick_params(axis='x', labelbottom=False, bottom=False)
        ax_strip.grid(False)
        for sp_ in ax_strip.spines.values():
            sp_.set_visible(False)
        ax_strip.set_title(gene, fontsize=TITLE_FS, fontweight='bold', pad=4)
        if ci == 0:
            for y_c, lab in [(2 * pitch + h_y / 2, 'Lineage'), (pitch + h_y / 2, 'Major type'),
                             (h_y / 2, '% expressing (ref)')]:
                ax_strip.text(-0.02, y_c, lab, transform=ax_strip.get_yaxis_transform(),
                              ha='right', va='center', fontsize=BODY_FS)

    # --- right column: legends (Lineage, Major type, x-tick colour key) stacked flush, then cbar
    rx = (ML + panels_w + 0.12) / fig_w
    lin_present = [l for l in LINEAGE_ORDER if l in set(CELL_LINEAGE.values())]
    maj_present = list(dict.fromkeys(CELL_MAJOR[ct] for ct in CELL_ORDER))
    lin_h = [mpatches.Patch(facecolor=_lin_lut.get(l, '#ccc'), edgecolor='black', linewidth=0.5, label=l)
             for l in lin_present]
    maj_h = [mpatches.Patch(facecolor=_maj_lut.get(m, '#ccc'), edgecolor='black', linewidth=0.5, label=m)
             for m in maj_present]
    top_y = (strip_bottom + strip_H) / fig_h
    _lkw = dict(loc='upper left', frameon=False, fontsize=BODY_FS,
                title_fontproperties={'weight': 'bold', 'size': BODY_FS}, borderpad=0.1, labelspacing=0.25)
    leg1 = fig.legend(handles=lin_h, title='Lineage', bbox_to_anchor=(rx, top_y),
                      handlelength=1.0, handleheight=1.0, **_lkw)
    leg1._legend_box.align = 'left'; fig.add_artist(leg1)
    fig.canvas.draw()
    bb1 = leg1.get_window_extent().transformed(fig.transFigure.inverted())
    leg2 = fig.legend(handles=maj_h, title='Major type', bbox_to_anchor=(rx, bb1.y0 - 0.005),
                      handlelength=1.0, handleheight=1.0, **_lkw)
    leg2._legend_box.align = 'left'; fig.add_artist(leg2)
    fig.canvas.draw()
    bb2 = leg2.get_window_extent().transformed(fig.transFigure.inverted())
    pct = producer_tag.replace('over_', '')
    expr_h = [Line2D([], [], linestyle='none', marker='none', label=f'≥{pct}%'),
              Line2D([], [], linestyle='none', marker='none', label=f'<{pct}%')]
    leg3 = fig.legend(handles=expr_h, title='Producer call (x-axis)', bbox_to_anchor=(rx, bb2.y0 - 0.005),
                      handlelength=0, handletextpad=0, **_lkw)
    leg3._legend_box.align = 'left'; fig.add_artist(leg3)
    _et = leg3.get_texts()
    _et[0].set_color(PRODUCER_C); _et[0].set_fontweight('bold')   # >= threshold -> bold red x-ticks
    _et[1].set_color('black')                                     # < threshold  -> plain black x-ticks

    # short horizontal colorbar for the % expressing row, just below the Producer Call legend
    fig.canvas.draw()
    bb3 = leg3.get_window_extent().transformed(fig.transFigure.inverted())
    sm = mpl.cm.ScalarMappable(norm=_pct_norm, cmap=cmap); sm.set_array([])
    cax = fig.add_axes([rx, bb3.y0 - 0.07, 0.6 / fig_w, 0.03])
    cbar = fig.colorbar(sm, cax=cax, orientation='horizontal')
    cbar.set_label('% expressing (ref)', fontsize=BODY_FS)
    cbar.ax.tick_params(labelsize=BODY_FS)

    # shared x-axis label centred under the panels
    center_x = (ML + panels_w / 2) / fig_w
    fig.text(center_x, 0.12 / fig_h, 'Candidate producers', ha='center', va='bottom',
             fontsize=BODY_FS, fontweight='bold')

    # title + note in the reserved top band (top-anchored)
    note = (f'{BIN_TAG} · reference {PROD_LABEL} · bold red = predicted producer of this gene, '
            f'plain = predicted non-producer (expect ns) · * = Wilcoxon r≠0 (n={sub["sample_id"].nunique()})')
    fig.text(center_x, (fig_h - 0.14) / fig_h, f'{title} — {stratum}', ha='center', va='top',
             fontweight='bold', fontsize=TITLE_FS)
    fig.text(center_x, (fig_h - 0.34) / fig_h, note, ha='center', va='top', fontsize=BODY_FS)


# Extended Data Fig. 9i: reference-producer vs chemokine correlation in the T cell zone
corr_boxplot_figure('T cell zone', 'partial_spearman',
                    'Extended Data Fig. 9i — reference-producer vs chemokine correlation (density-controlled)')


# --- Report + save

print('\nCohort mean partial Spearman (T cell zone):')
_tz = cdf[cdf['stratum'] == 'T cell zone']
if not _tz.empty:
    piv = (_tz.groupby(['cell_type', 'gene'])['partial_spearman']
           .mean().unstack().reindex(CELL_ORDER)[GENES])
    print(piv.round(2).to_string())

out_base = os.path.join(res_dir, '05_xenium_mechanisms', '09_chemokine_sources',
                        '06_reference_producer_corr')
os.makedirs(os.path.join(out_base, 'plots'), exist_ok=True)
os.makedirs(os.path.join(out_base, 'csv_files'), exist_ok=True)
cdf.to_csv(os.path.join(out_base, 'csv_files', f'reference_producer_corr{genes_tag}.csv'), index=False)
pdf_path = os.path.join(out_base, 'plots', f'ED9I_reference_producer_corr{genes_tag}.pdf')
save_multi_image(pdf_path)
print(f'\nSaved: {pdf_path}')
print('Script successfully completed')
