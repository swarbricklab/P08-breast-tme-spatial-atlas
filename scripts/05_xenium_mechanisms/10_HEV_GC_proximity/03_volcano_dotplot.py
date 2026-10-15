#!/usr/bin/env python3

# --- Aim
# S10I + S10J: render the close-vs-far DEG table (DEG_close_vs_far.csv from
# 02_DEG_close_vs_far.py) as:
#   - S10I: a volcano of the mean log-transformed difference (close - far) vs -log10 BH-P,
#           coloured by significance (purple = higher close to GC, brown = lower,
#           grey = ns), with canonical HEV markers outlined + bold-labelled and all
#           significant genes labelled (adjustText avoids overlaps).
#   - S10J: a transposed dot plot of the significant genes, split into up-/down-close
#           blocks, with close (<=500 um) vs far (>500 um) rows, dot size = % of cells
#           expressing and colour = per-gene-max-scaled mean expression.
# Both are built purely from the DEG CSV (no expression reload).

# --- Define environment
# env: shapely

# --- Load packages

import argparse
import os

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from adjustText import adjust_text

mpl.rcParams['pdf.fonttype'] = 42

# Colour scheme: darkorchid = higher close to GC (up), brown = higher far from GC (down).
CLOSE_COLOR = '#9932CC'   # darkorchid
FAR_COLOR   = '#8c6d31'   # brown


# --- Argparse arguments (Methods-fixed thresholds are hard-coded defaults)

parser = argparse.ArgumentParser(
    description="S10I volcano + S10J dot plot from the close-vs-far DEG table.",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
)
parser.add_argument("-s", "--samples", nargs="*", help="Accepted for runner compatibility (unused).")
parser.add_argument("--h5ad_dir", default=None, help="Accepted for runner compatibility (unused).")
parser.add_argument("--project_dir", required=True)
parser.add_argument("--res_dir", required=True)
parser.add_argument("--sample_sheet", default=None, help="Accepted for runner compatibility (unused).")
parser.add_argument("--near_gc_um", type=float, default=500.0)
parser.add_argument("--lfc_threshold", type=float, default=0.05,
                    help="Below |mean diff| this, dots are grey (ns).")
parser.add_argument("--padj_threshold", type=float, default=0.05,
                    help="Above BH-adjusted p this, dots are grey (ns).")
parser.add_argument("--hev_markers", nargs='+',
                    default=['TSPAN7', 'RAMP3', 'MEOX2', 'SPNS2', 'IL33', 'ACKR1',
                             'SELE', 'SELP', 'TMEM173', 'LTBR', 'SPP1'],
                    help="Bolded labels + outlined dots on the volcano.")
parser.add_argument("--min_label_distance", type=float, default=0.15,
                    help="Min L2 distance (mixed data units) between annotated points "
                         "before they're allowed to share a label.")

args = vars(parser.parse_args())
project_dir        = args["project_dir"]
res_dir            = args["res_dir"]
near_gc_um         = args["near_gc_um"]
lfc_threshold      = args["lfc_threshold"]
padj_threshold     = args["padj_threshold"]
hev_markers        = args["hev_markers"]
min_label_distance = args["min_label_distance"]

unit_label = 'cells'   # DEG is cell-level (S10I reports cell counts)


# --- Paths

out_base  = os.path.join(res_dir, '05_xenium_mechanisms', '10_HEV_GC_proximity')
deg_csv   = os.path.join(out_base, '02_DEG_close_vs_far', 'csv_files', 'DEG_close_vs_far.csv')
plots_dir = os.path.join(out_base, '03_volcano_dotplot', 'plots')
os.makedirs(plots_dir, exist_ok=True)

print(f'Reading DEG CSV: {deg_csv}')
deg = pd.read_csv(deg_csv)
deg['neg_log10_padj'] = -np.log10(np.clip(deg['padj'].values, 1e-300, None))
n_close = int(deg['n_close'].iloc[0])
n_far   = int(deg['n_far'].iloc[0])


# --- Volcano helpers (adjustText labelling + directional arrows)

def adjust_text_func(ax, df, x_val, y_val, gene_col='gene',
                     min_distance=0.15, label_size=7, color='#0f3c5a', lw=0.4,
                     bold_genes=None, priority_genes=None):
    """Annotate genes in df. priority_genes are ALWAYS labelled; the rest are
    filtered by min_distance to avoid collisions, then adjustText repositions
    everything so labels don't overlap."""
    priority_set = set(priority_genes) if priority_genes else set()
    texts = []
    annotated_coords = []

    def _label(row, ignore_distance=False):
        label = row[gene_col]
        x, y = row[x_val], row[y_val]
        if not ignore_distance:
            for x0, y0 in annotated_coords:
                if np.hypot(x - x0, y - y0) < min_distance:
                    return
        fw = 'bold' if (bold_genes and label in bold_genes) else 'normal'
        text = ax.text(x, y, label, fontsize=label_size, ha='right', va='bottom',
                       color='black', fontweight=fw)
        texts.append(text)
        annotated_coords.append((x, y))

    for _, row in df[df[gene_col].isin(priority_set)].iterrows():
        _label(row, ignore_distance=True)
    for _, row in df[~df[gene_col].isin(priority_set)].iterrows():
        _label(row, ignore_distance=False)

    for t in texts:
        t.set_zorder(200)
    if texts:
        adjust_text(
            texts, ax=ax,
            arrowprops=dict(arrowstyle='-', color=color, lw=lw),
            expand_points=(1.2, 1.4), expand_text=(1.2, 1.4),
            force_text=(0.5, 0.5),
        )


def draw_directional_arrows(ax, left_label='Far from GC', right_label='Close to GC',
                            label_fontsize=10, y_offset=1.03):
    arrowprops = dict(arrowstyle='-|>', color='black', lw=1,
                      shrinkA=0, shrinkB=0, joinstyle='miter')
    ax.annotate('', xy=(1, y_offset), xytext=(0.7, y_offset),
                xycoords='axes fraction', textcoords='axes fraction',
                arrowprops=arrowprops, clip_on=False)
    ax.annotate('', xy=(0, y_offset), xytext=(0.3, y_offset),
                xycoords='axes fraction', textcoords='axes fraction',
                arrowprops=arrowprops, clip_on=False)
    midpoint = 0.3 / 2
    ax.text(1 - midpoint, y_offset + 0.01, right_label,
            fontsize=label_fontsize, ha='center', va='bottom',
            transform=ax.transAxes, clip_on=False)
    ax.text(0 + midpoint, y_offset + 0.01, left_label,
            fontsize=label_fontsize, ha='center', va='bottom',
            transform=ax.transAxes, clip_on=False)


def draw_volcano(deg, annotate_df, title, figsize=(4, 4),
                 label_size=7, min_distance=0.15,
                 highlight_outline=None, priority_genes=None):
    fig, ax = plt.subplots(figsize=figsize)

    sig_mask = (deg['padj'] < padj_threshold) & (deg['mean_diff'].abs() > lfc_threshold)
    sig_up   = sig_mask & (deg['mean_diff'] > 0)
    sig_down = sig_mask & (deg['mean_diff'] < 0)
    ns_mask  = ~sig_mask

    DOT_Z = 2
    ax.scatter(deg.loc[ns_mask, 'mean_diff'], deg.loc[ns_mask, 'neg_log10_padj'],
               s=6, c='#bbbbbb', alpha=0.6, linewidths=0, zorder=DOT_Z,
               label=f'ns ({int(ns_mask.sum())})')
    ax.scatter(deg.loc[sig_up, 'mean_diff'], deg.loc[sig_up, 'neg_log10_padj'],
               s=14, c=CLOSE_COLOR, alpha=0.85, linewidths=0, zorder=DOT_Z,
               label=f'↑ close (n={int(sig_up.sum())})')
    ax.scatter(deg.loc[sig_down, 'mean_diff'], deg.loc[sig_down, 'neg_log10_padj'],
               s=14, c=FAR_COLOR, alpha=0.85, linewidths=0, zorder=DOT_Z,
               label=f'↓ close (n={int(sig_down.sum())})')

    LINE_Z = 0
    ax.axhline(y=-np.log10(padj_threshold), color='gray', linestyle='--',
               linewidth=0.6, alpha=0.7, zorder=LINE_Z)
    ax.axvline(x=lfc_threshold,  color='gray', linestyle='--',
               linewidth=0.6, alpha=0.7, zorder=LINE_Z)
    ax.axvline(x=-lfc_threshold, color='gray', linestyle='--',
               linewidth=0.6, alpha=0.7, zorder=LINE_Z)
    ax.axvline(x=0, color='darkgray', linewidth=0.6, linestyle='-', zorder=LINE_Z)
    ax.axhline(y=0, color='darkgray', linewidth=0.6, linestyle='-', zorder=LINE_Z)

    xlim = ax.get_xlim()
    xmax = max(abs(xlim[0]), abs(xlim[1]))
    ax.set_xlim(-xmax, xmax)

    if highlight_outline:
        sel = annotate_df[annotate_df['gene'].isin(highlight_outline)]
        if len(sel) > 0:
            ax.scatter(sel['mean_diff'], sel['neg_log10_padj'],
                       s=22, facecolors='none', edgecolors='#0f3c5a',
                       linewidths=0.9, zorder=4,
                       label=f'HEV markers (n={len(sel)})')

    adjust_text_func(ax, annotate_df, 'mean_diff', 'neg_log10_padj',
                     gene_col='gene', min_distance=min_distance,
                     label_size=label_size, bold_genes=set(hev_markers),
                     priority_genes=priority_genes)

    ax.set_xlabel('Mean expression difference (close − far)')
    ax.set_ylabel('−log10(BH-adjusted p)')
    draw_directional_arrows(ax, left_label='Far from GC', right_label='Close to GC')

    xlim = ax.get_xlim()
    ax.set_xlim(xlim[0] * 1.05, xlim[1] * 1.05)
    ylim = ax.get_ylim()
    ax.set_ylim(ylim[0], ylim[1] * 1.05)
    ax.grid(False)
    ax.set_position([0.15, 0.15, 0.70, 0.70])

    xtick_fs = mpl.rcParams['xtick.labelsize']
    if isinstance(xtick_fs, str):
        xtick_fs = mpl.font_manager.FontProperties(size=xtick_fs).get_size_in_points()
    ax.legend(loc='center left', bbox_to_anchor=(1.02, 0.5),
              fontsize=xtick_fs, markerscale=1.8, frameon=False, borderaxespad=0)

    plt.suptitle(title, y=1.10, fontweight='bold', fontsize=8)
    return fig


def draw_deg_dotplot_transposed(deg, gene_list, title, bold_genes=None, near_gc_um=500,
                                space_per_dot=0.04, width_ratio=0.8, max_dot_size=100,
                                label_fontsize=8, figsize=(5, 5),
                                cbar_height_factor=0.4, cbar_width_factor=3.4,
                                cbar_to_size_gap=0.085, cbar_title_pad=10,
                                bar_fill_alpha=0.55):
    """Transposed dot plot: genes on the y-axis (split into up-/down-close blocks,
    most-significant at top), close/far on the x-axis, vertical count bars on top,
    vertical colourbar + size legend on the right. Built from the DEG CSV."""
    df = deg[deg['gene'].isin(gene_list)].copy()
    if df.empty:
        print('No genes for transposed DEG dotplot — skipping.')
        return None
    df['direction'] = np.where(df['mean_diff'] > 0, 'up', 'down')
    gmax = df[['mean_close', 'mean_far']].max(axis=1).clip(lower=1e-9)
    df['close_norm'] = df['mean_close'] / gmax
    df['far_norm']   = df['mean_far']   / gmax
    df['signed_z'] = (np.sign(df['mean_diff'])
                      * -np.log10(np.clip(df['pval'].values, 1e-300, None)))
    up_df = df[df['direction'] == 'up'  ].sort_values('signed_z', ascending=False)
    dn_df = df[df['direction'] == 'down'].sort_values('signed_z', ascending=True)
    df = pd.concat([up_df, dn_df], ignore_index=False)

    n_up = int((df['direction'] == 'up').sum())
    n_dn = len(df) - n_up
    bold_set = set(bold_genes) if bold_genes else set()

    fig = plt.figure(figsize=figsize)

    ax_x0 = 0.42
    ax_y0 = 0.30
    sub_gap = space_per_dot * 0.5
    dot_w = 2 * space_per_dot   # 2 categories (close, far) on x

    cmap = plt.cm.Reds
    data_min = float(min(df['close_norm'].min(), df['far_norm'].min()))
    data_max = float(max(df['close_norm'].max(), df['far_norm'].max()))
    if data_min == data_max:
        data_min = max(data_min - 1e-3, 0.0)
        data_max = min(data_max + 1e-3, 1.0)
    norm = mpl.colors.Normalize(vmin=data_min, vmax=data_max)

    cursor_y = ax_y0
    dn_ax = up_ax = None
    if n_dn > 0:
        dn_ax = fig.add_axes([ax_x0, cursor_y, dot_w, n_dn * space_per_dot * width_ratio])
        cursor_y += n_dn * space_per_dot * width_ratio + (sub_gap if n_up > 0 else 0)
    if n_up > 0:
        up_ax = fig.add_axes([ax_x0, cursor_y, dot_w, n_up * space_per_dot * width_ratio])
    sub_axes = []
    if up_ax is not None: sub_axes.append((up_ax, up_df, '↑ close'))
    if dn_ax is not None: sub_axes.append((dn_ax, dn_df, '↓ close'))
    bottom_ax = dn_ax if dn_ax is not None else up_ax

    for (sub_ax, sub_df, _sub_label) in sub_axes:
        n_g = len(sub_df)
        for y_idx, (_, row) in enumerate(sub_df.iterrows()):
            sub_ax.scatter(0, y_idx, s=row['pct_close'] * max_dot_size,
                           c=[cmap(norm(row['close_norm']))],
                           edgecolor='black', linewidth=0.4, zorder=3)
            sub_ax.scatter(1, y_idx, s=row['pct_far'] * max_dot_size,
                           c=[cmap(norm(row['far_norm']))],
                           edgecolor='black', linewidth=0.4, zorder=3)
        sub_ax.set_yticks(range(n_g))
        sub_ax.set_yticklabels(sub_df['gene'].tolist(), fontsize=label_fontsize)
        for tl in sub_ax.get_yticklabels():
            if tl.get_text() in bold_set:
                tl.set_fontweight('bold')
        sub_ax.set_ylim(-0.5, n_g - 0.5)
        sub_ax.invert_yaxis()                     # most-significant gene at top
        sub_ax.set_xlim(-0.6, 1.6)
        if sub_ax is bottom_ax:
            sub_ax.set_xticks([0, 1])
            sub_ax.set_xticklabels([f'<{int(near_gc_um)} µm', f'>{int(near_gc_um)} µm'],
                                   rotation=90, fontsize=label_fontsize)
            sub_ax.set_xlabel('HEV-GC\ndistance', fontsize=label_fontsize)
        else:
            sub_ax.set_xticks([])
        for sp in ('top', 'right', 'bottom', 'left'):
            sub_ax.spines[sp].set_visible(True)
        sub_ax.tick_params(axis='x', labelsize=label_fontsize)

    fig.canvas.draw()
    ax_pos = mpl.transforms.Bbox.union([sa[0].get_position() for sa in sub_axes])

    # Count bars — vertical, flush on top of the dotplot.
    count_h = 2.0 * space_per_dot
    count_ax = fig.add_axes([ax_pos.x0, ax_pos.y1, ax_pos.width, count_h])
    counts_arr = [n_close, n_far]                 # x=0 close, x=1 far
    edge_colors = [CLOSE_COLOR, FAR_COLOR]
    fill_colors = [mpl.colors.to_rgba(c, alpha=bar_fill_alpha) for c in edge_colors]
    for x_i, (val, fill, edge) in enumerate(zip(counts_arr, fill_colors, edge_colors)):
        count_ax.bar(x_i, val, width=0.7, color=fill, edgecolor=edge, linewidth=0.9)
    annot_gap = max(counts_arr) * 0.05
    for i, v in enumerate(counts_arr):
        lbl = f'{v/1000:.1f}k' if v >= 1000 else f'{int(v)}'
        count_ax.text(i, v + annot_gap, lbl, ha='center', va='bottom',
                      fontsize=label_fontsize, color='black')
    count_ax.set_xlim(-0.6, 1.6)
    count_ax.set_ylim(0, max(counts_arr) * 1.35)
    count_ax.set_xticks([]); count_ax.set_yticks([])
    for sp in count_ax.spines.values():
        sp.set_visible(False)
    count_ax.set_xlabel(f'n {unit_label}', fontsize=label_fontsize)
    count_ax.xaxis.set_label_position('top')

    # Vertical colourbar (right).
    cbar_pad = 0.03
    cbar_x0 = ax_pos.x1 + cbar_pad
    cbar_width  = space_per_dot * cbar_height_factor
    cbar_height = space_per_dot * cbar_width_factor
    cbar_ax = fig.add_axes([cbar_x0, ax_pos.y0, cbar_width, cbar_height])
    sm = mpl.cm.ScalarMappable(cmap=cmap, norm=norm); sm.set_array([])
    cb = fig.colorbar(sm, cax=cbar_ax, orientation='vertical')
    cb.set_ticks([data_min, data_max])
    cb.set_ticklabels([f'{data_min:.2f}', f'{data_max:.2f}'])
    cbar_ax.tick_params(labelsize=label_fontsize)
    cbar_ax.set_title('Mean expression\n(scaled to max)', loc='left',
                      fontsize=label_fontsize, pad=cbar_title_pad, color='black')

    # Size legend — vertical stack above the cbar.
    size_ax_y0 = ax_pos.y0 + cbar_height + cbar_to_size_gap
    size_ax = fig.add_axes([cbar_x0, size_ax_y0, 0.14, space_per_dot * cbar_width_factor])
    size_ax.axis('off'); size_ax.set_xlim(0, 1); size_ax.set_ylim(0, 1)
    for y_pos, pct in zip([0.2, 0.5, 0.8], [0.25, 0.5, 1.0]):
        size_ax.scatter(0.12, y_pos, s=pct * max_dot_size, c='lightgray',
                        edgecolor='black', linewidth=0.4, zorder=3)
        size_ax.text(0.32, y_pos, f'{int(pct*100)}%', ha='left', va='center',
                     fontsize=label_fontsize - 1, color='black')
    size_ax.set_title('% expressing', loc='left', fontsize=label_fontsize, pad=2, color='black')

    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    tops = []
    for art in [count_ax, size_ax] + [sa[0] for sa in sub_axes]:
        try:
            bb = art.get_tightbbox(renderer)
            tops.append(fig.transFigure.inverted().transform(bb)[1][1])
        except Exception:
            pass
    suptitle_y = (max(tops) if tops else 0.95) + 0.008
    plt.suptitle(title, x=(ax_pos.x0 + ax_pos.x1) / 2, y=suptitle_y, fontweight='bold',
                 fontsize=label_fontsize + 1, va='bottom', color='black')
    return fig


# --- Plot

plt.rcParams['axes.axisbelow'] = True
hev_set = set(hev_markers)

# S10I — volcano: annotate all sig genes + HEV markers (HEV markers always shown,
# even if non-sig). Both HEV markers and every sig gene are priority so no dot is
# dropped by min_distance; adjustText repositions overlaps.
sig_mask = (deg['padj'] < padj_threshold) & (deg['mean_diff'].abs() > lfc_threshold)
sig_genes = set(deg.loc[sig_mask, 'gene'])
annotate = deg[sig_mask | deg['gene'].isin(hev_set)]
print(f'S10I annotate set: {len(annotate)} genes ({int(sig_mask.sum())} sig + HEV markers)')
fig_v = draw_volcano(
    deg, annotate,
    title=(f'S10I — HEVs close vs far from GC (≤{int(near_gc_um)} µm = close)\n'
           f'{n_close} close vs {n_far} far {unit_label}'),
    figsize=(4, 4), label_size=7, min_distance=min_label_distance,
    highlight_outline=hev_set, priority_genes=hev_set | sig_genes,
)
volcano_pdf = os.path.join(plots_dir, 'S10I_DEG_close_vs_far_volcano.pdf')
fig_v.savefig(volcano_pdf, bbox_inches='tight', dpi=300)
plt.close(fig_v)
print(f'Saved: {volcano_pdf}')

# S10J — transposed dot plot of the significant genes (up/down-close blocks).
sig_gene_list = deg.loc[sig_mask, 'gene'].tolist()
print(f'S10J dotplot: {len(sig_gene_list)} sig genes')
fig_d = draw_deg_dotplot_transposed(
    deg, sig_gene_list,
    title=(f'S10J — ACKR1+ HEV ECs, close vs far from GC\n'
           f'significant genes (BH-P<{padj_threshold}, |diff|>{lfc_threshold})'),
    bold_genes=hev_set, near_gc_um=near_gc_um,
)
if fig_d is not None:
    dotplot_pdf = os.path.join(plots_dir, 'S10J_DEG_close_vs_far_dotplot.pdf')
    fig_d.savefig(dotplot_pdf, bbox_inches='tight', dpi=300)
    plt.close(fig_d)
    print(f'Saved: {dotplot_pdf}')

print('Script successfully completed')
