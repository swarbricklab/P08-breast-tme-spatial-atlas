#!/usr/bin/env python3

# --- Aim
# Cohort-level summary of the TLS zonation classifications from 01_classify_zones.py.
# Produces (one multi-page PDF):
#   * Figure S9C  - spatial map of B/T cell zone bins per sample
#   * Figure S9D  - per-bin T:B cell-ratio distribution, split by zone
#   * Figure S10G (support) - spatial overview of GC bins per sample
#   * Figure S10H - number of discrete GCs per sample, coloured by clinical subtype


# --- Define environment
# env: xen-hev

# --- Load packages

import argparse
import os
import math

import anndata as ad
import pandas as pd
import numpy as np

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle

import geopandas as gpd

mpl.rcParams['pdf.fonttype'] = 42


# --- Argparse arguments

parser = argparse.ArgumentParser(description="Cohort-level summary of TLS zonation (B/T + GC)",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

# Runner contract (submit_multi_sample_job_nci.sh -> run_multi_sample_job.sh)
parser.add_argument("--h5ad_dir", help="Path to directory for .h5ad input files")
parser.add_argument("--project_dir", help="Path to the root of the project directory")
parser.add_argument("--res_dir", help="Path to results directory")
parser.add_argument("--sample_sheet", help="Path to sample sheet CSV file", default=None)
parser.add_argument("-s", "--samples", nargs="*", type=str, help="List of samples (section_id)")

# Analysis parameters (Methods-fixed defaults; exposed but never required)
parser.add_argument("--grid_um", type=float, default=50.0, help="Bin size used in step 01 (microns)")
parser.add_argument("--annotation_level", default='major',
                    help="Cell-type level used for the T:B per-bin ratio")

args = vars(parser.parse_args())

h5ad_dir         = args["h5ad_dir"]
project_dir      = args["project_dir"]
res_dir          = args["res_dir"]
sample_list      = args["samples"]
grid_um          = args["grid_um"]
annotation_level = args["annotation_level"]


# --- Font sizes
TICK_FS  = 8
TITLE_FS = 12
SUP_FS   = 14
ANNO_FS  = 7


# --- Sample metadata

if args.get("sample_sheet"):
    sample_sheet_path = args["sample_sheet"]
else:
    sample_sheet_path = os.path.join(project_dir, 'config', 'sample_sheet_xenium.csv')

sample_sheet    = pd.read_csv(sample_sheet_path)
tissue_metadata = pd.read_csv(os.path.join(project_dir, 'resources', 'tissue_metadata.csv'))

merged = sample_sheet.merge(tissue_metadata, on='tissue_id')
subtype_dict = dict(zip(merged['section_id'], merged['clinical_subtype']))

# Clinical-subtype ordering (global convention)
subtype_order = ['Luminal A', 'Luminal B', 'Luminal B (HER2+)', 'HER2+', 'TNBC']
present_subtypes = sorted({subtype_dict[s] for s in sample_list if s in subtype_dict},
                          key=lambda x: subtype_order.index(x) if x in subtype_order else float('inf'))
samples_by_subtype = {st: sorted([s for s in sample_list if subtype_dict.get(s) == st])
                      for st in present_subtypes}
sorted_samples = [s for st in present_subtypes for s in samples_by_subtype[st]]


# --- Subtype colours (from SCSubtype_colors.csv)

SUBTYPE_COLOR_SYNONYMS = {
    'LumA':         'Luminal A',
    'LumB':         'Luminal B',
    'LumB (HER2+)': 'Luminal B (HER2+)',
    'LumB (ER+)':   'Luminal B (ER+)',
}
color_lut = pd.read_csv(os.path.join(project_dir, 'resources', 'plotting_colors',
                                     'SCSubtype_colors.csv'), index_col='subtype')


def subtype_color(subtype):
    key = SUBTYPE_COLOR_SYNONYMS.get(subtype, subtype)
    if key in color_lut.index:
        return color_lut.loc[key, 'color']
    return 'dimgray'


# --- Zone palettes

bt_color_dict = {
    'B_cell_zone':  mpl.colors.cnames['cornflowerblue'],
    'T_cell_zone':  mpl.colors.cnames['crimson'],
    'mixed':        '#921e95',
    'unclassified': mpl.colors.cnames['lightgray'],
}
gc_color = mpl.colors.cnames['darkorchid']
classified_zones = ['B_cell_zone', 'T_cell_zone', 'mixed']


def pretty_zone(z):
    return z.replace('_', ' ')


# --- Cell-type helper (major-level T/B call)

def load_major_labels(sample):
    """Return a Series (index = cell_id) of major-level cell types."""
    path = os.path.join(res_dir, '02_cell_type_mapping', 'label_transfer',
                        'xenium', f'{annotation_level}_level',
                        f'{sample}_{annotation_level}_pred_cell_types.csv')
    df = pd.read_csv(path, index_col=0)
    return df['pred_cell_type']


# --- Load per-sample step-01 outputs

zone_dir = os.path.join(res_dir, '05_xenium_mechanisms', '07_TLS_zonation', '01_classify_zones', 'csv_files')
tls_regions_dir = os.path.join(res_dir, '05_xenium_mechanisms', '01_identify_TLS_regions', 'csv_files')

bt_bin_dfs   = {}
bt_cell_dfs  = {}
gc_bin_dfs   = {}
tls_gdf_dict = {}
in_tls_dict  = {}

for sample in sorted_samples:
    bt_bin  = pd.read_csv(os.path.join(zone_dir, f'{sample}_BT_zone_bins.csv'))
    bt_cell = pd.read_csv(os.path.join(zone_dir, f'{sample}_BT_zone_classification.csv'))
    if not bt_cell.empty:
        bt_cell.rename(columns={bt_cell.columns[0]: 'cell_id'}, inplace=True)
        major = load_major_labels(sample)
        bt_cell['major'] = bt_cell['cell_id'].map(major)
    gc_bin  = pd.read_csv(os.path.join(zone_dir, f'{sample}_GC_zone_bins.csv'))

    bt_bin_dfs[sample]  = bt_bin
    bt_cell_dfs[sample] = bt_cell
    gc_bin_dfs[sample]  = gc_bin
    tls_gdf_dict[sample] = gpd.read_file(os.path.join(tls_regions_dir, f'{sample}_TLS_regions.geojson'))
    in_tls_dict[sample]  = pd.read_csv(os.path.join(tls_regions_dir,
                                                    f'{sample}_cell_TLS_annotations.csv'), index_col=0)

print(f'Loaded zonation outputs for {len(sorted_samples)} samples')


# --- Micron crop + grey backdrop helpers (x_centroid/y_centroid frame)

def square_crop(x, y, border=0.1):
    x0, x1 = float(x.min()), float(x.max())
    y0, y1 = float(y.min()), float(y.max())
    xr, yr = x1 - x0, y1 - y0
    cx, cy = x0 + xr / 2, y0 + yr / 2
    half = max(xr, yr) / 2
    b = (border / 2) * max(xr, yr)
    return (cx - half - b, cx + half + b, cy - half - b, cy + half + b)


def spatial_grid(zone_col, palette, title, legend_handles, suptitle):
    """3-column per-sample spatial grid: grey cells + classified bins + TLS outlines."""
    ncols = 3
    nrows = math.ceil(len(sorted_samples) / ncols)
    panel = 4.5
    fig, axs = plt.subplots(nrows, ncols, figsize=(panel * ncols, panel * nrows), squeeze=False)
    axs_flat = axs.ravel()

    for n, sample in enumerate(sorted_samples):
        ax = axs_flat[n]
        adata = ad.read_h5ad(os.path.join(h5ad_dir, 'xenium', f'{sample}.h5ad'))
        xy = adata.obs[['x_centroid', 'y_centroid']].values.astype(float)
        crop = square_crop(xy[:, 0], xy[:, 1])
        ax.scatter(xy[:, 0], xy[:, 1], s=1, c='lightgray', marker='.', linewidths=0, alpha=0.5,
                   rasterized=True)

        bin_df = gc_bin_dfs[sample] if zone_col == 'GC_zone' else bt_bin_dfs[sample]
        bin_df = bin_df.dropna(subset=[zone_col])
        if len(bin_df):
            rects = [Rectangle((x - grid_um / 2, y - grid_um / 2), grid_um, grid_um)
                     for x, y in zip(bin_df['bin_x'], bin_df['bin_y'])]
            if isinstance(palette, dict):
                colors = [palette[z] for z in bin_df[zone_col]]
            else:
                colors = [palette] * len(bin_df)
            ax.add_collection(PatchCollection(rects, facecolors=colors, edgecolors='none',
                                              alpha=0.85, rasterized=True))
        tls_gdf_dict[sample].plot(ax=ax, facecolor='none', edgecolor='black', linewidth=0.6)

        ax.set_xlim(crop[0], crop[1]); ax.set_ylim(crop[2], crop[3]); ax.invert_yaxis()
        ax.set_aspect('equal'); ax.set_xticks([]); ax.set_yticks([])
        ax.set_title(f'{sample}\n({subtype_dict.get(sample, "?")})', fontweight='bold', fontsize=TITLE_FS)
        del adata

    for k in range(len(sorted_samples), len(axs_flat)):
        axs_flat[k].axis('off')

    fig.legend(handles=legend_handles, loc='lower center', ncol=len(legend_handles),
               frameon=False, fontsize=11)
    plt.suptitle(suptitle, y=1.02, fontweight='bold', fontsize=SUP_FS)
    plt.subplots_adjust(wspace=0.05, hspace=0.25, bottom=0.08)


# --- Figure S9C: spatial map of B/T cell zone bins

bt_handles = [mlines.Line2D([], [], color=bt_color_dict[z], marker='s', lw=0, markersize=12,
                            label=pretty_zone(z)) for z in classified_zones]
spatial_grid('TLS_zone', bt_color_dict, 'B/T zones', bt_handles,
             f'B/T cell zone classification across cohort (grid={int(grid_um)}um)')


# --- Figure S9D: per-bin T:B cell-ratio distribution, split by zone
# Cells in the same bin share (T_score, B_score), so that pair identifies a bin.

zones_for_ratio = ['B_cell_zone', 'T_cell_zone']
ratio_rows = []
for sample in sorted_samples:
    cdf = bt_cell_dfs[sample]
    if cdf.empty:
        continue
    cdf = cdf.dropna(subset=['TLS_zone'])
    for (_, _), bin_cells in cdf.groupby(['TLS_zone_T_score', 'TLS_zone_B_score']):
        zone = bin_cells['TLS_zone'].iloc[0]
        if zone not in zones_for_ratio:
            continue
        n_t = int((bin_cells['major'] == 'T-cells').sum())
        n_b = int((bin_cells['major'] == 'B-cells').sum())
        if (n_t + n_b) == 0:
            continue
        ratio_rows.append({'sample_id': sample, 'zone': zone, 'n_T': n_t, 'n_B': n_b,
                           'T_prop_of_TB': n_t / (n_t + n_b)})
ratio_df = pd.DataFrame(ratio_rows)

BODY_ALPHA = 0.55
fig, ax = plt.subplots(figsize=(4.0, 4.0))
ax.set_position([0.15, 0.30, 0.45, 0.30])
edges = np.linspace(0, 1, 26)
centers = (edges[:-1] + edges[1:]) / 2
bin_w = (edges[1] - edges[0]) * 0.85
ax.axvline(0.5, color='lightgray', lw=1, zorder=0)
if not ratio_df.empty:
    for zone in zones_for_ratio:
        vals = ratio_df.loc[ratio_df['zone'] == zone, 'T_prop_of_TB'].values
        rgb = mpl.colors.to_rgb(bt_color_dict[zone])
        counts, _ = np.histogram(vals, bins=edges)
        props = counts / counts.sum() if counts.sum() else counts
        ax.bar(centers, props, width=bin_w, color=(*rgb, BODY_ALPHA),
               edgecolor=(*rgb, 1.0), linewidth=1.0, label=pretty_zone(zone), zorder=2)
ax.set_xlim(0, 1)
ax.set_xticks([0, 0.25, 0.5, 0.75, 1.0])
ax.set_xticklabels(['0/100', '25/75', '50/50', '75/25', '100/0'], fontsize=TICK_FS)
ax.set_xlabel('T : B cells per bin', fontsize=TICK_FS, fontweight='bold')
ax.set_ylabel('Proportion of bins', fontsize=TICK_FS)
ax.tick_params(labelsize=TICK_FS)
ax.grid(False)
for sp in ['top', 'right']:
    ax.spines[sp].set_visible(False)
ax.legend(loc='center left', bbox_to_anchor=(1.02, 0.5), frameon=False, fontsize=TICK_FS,
          handlelength=0.8, handleheight=0.8)
ax.set_title('Per-bin T:B cell ratio distribution by zone', fontsize=SUP_FS, fontweight='bold', pad=6)


# --- Figure S10G (support): spatial overview of GC bins

gc_handles = [mlines.Line2D([], [], color=gc_color, marker='s', lw=0, markersize=12, label='GC')]
spatial_grid('GC_zone', gc_color, 'GC zones', gc_handles,
             f'GC zone classification across cohort (grid={int(grid_um)}um)')


# --- Figure S10H: number of discrete GCs per sample, coloured by clinical subtype

display_name_dict = {s: s for s in sample_list}


def make_barplot_figure_transposed(n_samples):
    """The discrete-GC barplot rotated 90°: horizontal bars, samples on the y-axis, a
    thin subtype-colour strip on the RIGHT, legend further right. The bar AXES are an
    exact rotation of the original ax: new ax width = original ax height, new ax height
    = original ax width."""
    # Original bar-ax dimensions (from make_barplot_figure), in inches
    per_sample_in    = 0.30
    orig_avail       = 0.85 - 0.30
    orig_bar_h_frac  = orig_avail - orig_avail * 0.05 - 0.005
    orig_ax_h_in     = orig_bar_h_frac * 2.9        # original ax HEIGHT (count direction)
    orig_ax_w_in     = n_samples * per_sample_in    # original ax WIDTH  (sample direction)

    ax_w_in = orig_ax_h_in       # rotate: new width  = original height
    ax_h_in = orig_ax_w_in       #         new height = original width

    left_in    = 0.75            # room for the sample y-labels
    strip_gap  = 0.06
    strip_w_in = 0.10
    legend_in  = 1.35            # room for the subtype legend on the right
    top_in     = 0.42            # title
    bottom_in  = 0.15

    fig_w_in = left_in + ax_w_in + strip_gap + strip_w_in + legend_in
    fig_h_in = top_in + ax_h_in + bottom_in
    fig = plt.figure(figsize=(fig_w_in, fig_h_in))

    bar_ax = fig.add_axes([left_in / fig_w_in, bottom_in / fig_h_in,
                           ax_w_in / fig_w_in, ax_h_in / fig_h_in])
    strip_x_in = left_in + ax_w_in + strip_gap
    strip_ax = fig.add_axes([strip_x_in / fig_w_in, bottom_in / fig_h_in,
                             strip_w_in / fig_w_in, ax_h_in / fig_h_in], sharey=bar_ax)
    return fig, bar_ax, strip_ax


def subtype_ordered_barplot_transposed(values_by_sample, bar_ax, strip_ax, title,
                                       annotate=True, value_fmt=None):
    """Original barplot rotated 90°: samples on the y-axis (kept), horizontal bars,
    count axis on the bottom, subtype-colour rectangles as a vertical strip on the
    RIGHT, subtype legend to the right of that."""
    values = np.array([values_by_sample[s] for s in sorted_samples], dtype=float)
    colors = [subtype_color(subtype_dict.get(s, '?')) for s in sorted_samples]
    y = np.arange(len(sorted_samples))

    bar_ax.barh(y, values, color=colors, alpha=0.85, edgecolor='black', linewidth=0.4)
    bar_ax.set_yticks(y)
    bar_ax.set_yticklabels([display_name_dict.get(s, s) for s in sorted_samples],
                           fontsize=TICK_FS)
    bar_ax.set_ylim(-0.5, len(y) - 0.5); bar_ax.invert_yaxis()   # first sample at top
    bar_ax.tick_params(axis='y', labelsize=TICK_FS)
    bar_ax.set_xticks([])                 # x-axis removed; counts shown as annotations
    if annotate and len(values):
        xmax = float(values.max()) if values.max() > 0 else 1.0
        for yi, v in enumerate(values):
            label = value_fmt(v) if value_fmt is not None else f'{int(v)}'
            bar_ax.text(v + xmax * 0.02, yi, label, ha='left', va='center', fontsize=ANNO_FS)
        bar_ax.set_xlim(0, xmax * 1.20)
    bar_ax.grid(False)
    for sp_ in ('top', 'right', 'bottom'):   # remove x-axis; keep left (y-axis)
        bar_ax.spines[sp_].set_visible(False)
    bar_ax.set_title(title, fontweight='bold', fontsize=TITLE_FS, pad=6)

    # Subtype-colour strip on the RIGHT (one stacked rectangle per sample)
    for yi, col in enumerate(colors):
        strip_ax.add_patch(mpatches.Rectangle((0, yi - 0.4), 1, 0.8,
                                              facecolor=col, edgecolor='black',
                                              linewidth=0.3, alpha=0.85))
    strip_ax.set_xlim(0, 1)
    strip_ax.set_xticks([])
    strip_ax.tick_params(axis='y', left=False, labelleft=False)
    strip_ax.grid(False)
    for spine in strip_ax.spines.values():
        spine.set_visible(False)

    legend_handles = [mpatches.Patch(facecolor=subtype_color(st), edgecolor='black',
                                     label=st, alpha=0.85)
                      for st in present_subtypes]
    bar_ax.legend(handles=legend_handles, title='Subtype', loc='upper left',
                  bbox_to_anchor=(1.30, 1.0), frameon=False, fontsize=TICK_FS,
                  title_fontproperties={'weight': 'bold', 'size': TICK_FS},
                  alignment='left', handlelength=0.8, handleheight=0.8)


def n_discrete_gcs(gc_bin_df):
    if 'discrete_GC' not in gc_bin_df.columns or gc_bin_df.empty:
        return 0
    return int(pd.to_numeric(gc_bin_df['discrete_GC'], errors='coerce').dropna().nunique())


gc_counts = {s: n_discrete_gcs(gc_bin_dfs[s]) for s in sorted_samples}
print('Discrete GC counts per sample:')
print(pd.Series(gc_counts))

fig, bar_ax, strip_ax = make_barplot_figure_transposed(len(sorted_samples))
subtype_ordered_barplot_transposed(gc_counts, bar_ax, strip_ax, title='GCs per sample')


# --- Save outputs

out_dir   = os.path.join(res_dir, '05_xenium_mechanisms', '07_TLS_zonation', '02_zonation_overview')
plots_dir = os.path.join(out_dir, 'plots')
csv_dir   = os.path.join(out_dir, 'csv_files')
os.makedirs(plots_dir, exist_ok=True)
os.makedirs(csv_dir, exist_ok=True)

pd.Series(gc_counts, name='n_discrete_GCs').to_csv(
    os.path.join(csv_dir, 'all_samples_n_discrete_GCs.csv'), index_label='sample_id')
if not ratio_df.empty:
    ratio_df.to_csv(os.path.join(csv_dir, 'per_bin_TB_ratio.csv'), index=False)

pdf_path = os.path.join(plots_dir, 'cohort_zonation_overview.pdf')
pp = PdfPages(pdf_path)
for fig_num in plt.get_fignums():
    plt.figure(fig_num).savefig(pp, format='pdf', bbox_inches='tight', dpi=300)
pp.close()
print(f'Saved: {pdf_path}')


# --- Finish

print('Script successfully completed')
