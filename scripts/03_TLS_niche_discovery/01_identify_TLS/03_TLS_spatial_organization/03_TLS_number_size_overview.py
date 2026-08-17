#!/usr/bin/env python3

# --- Aim
# Figure-quality cohort overview of per-sample TLS characteristics across the
# Visium cohort (Supplementary Figure S3B): a per-sample sorted bar of TLS number
# with stacked panels for TLS size, spatial clustering (DBSCAN clustered vs lone)
# and topology (intratumoral / peritumoral / stromal_distant). Reads the
# per-sample TLS-model spot counts and the per-TLS distances/clustering tables
# produced by 01 and 02.

# --- Define environment
# env: spatial-pilot

# --- Load packages

import argparse
import os
import sys
from datetime import datetime

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.patches as mpatches
import seaborn as sns

mpl.rcParams['pdf.fonttype'] = 42

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..'))
from helper_functions import save_multi_image


# --- Legend stacking helpers
# Figure-level legends anchored just right of the widest plotted axes.

def position_legend(axs):
    """Return (right_edge, top_edge) of the tight bbox of axs in figure coords."""
    if not isinstance(axs, (list, tuple, np.ndarray)):
        axs = [axs]
    fig = axs[0].figure
    fig.canvas.draw()
    right_edges, top_edges = [], []
    inv = fig.transFigure.inverted()
    for ax in axs:
        bbox = ax.get_tightbbox(fig.canvas.get_renderer())
        bbox_fig = inv.transform(bbox)
        right_edges.append(bbox_fig[1][0])
        top_edges.append(bbox_fig[1][1])
    return max(right_edges), max(top_edges)


def add_color_legend(fig, labels, color_dict, title='', x=0, y=0,
                     loc='upper left', fontsize=10, edgecolor='black',
                     linewidth=0.6, alpha=1.0):
    """Add a figure-level legend anchored at (x, y) in figure coords."""
    handles = [
        mpatches.Patch(facecolor=color_dict[label], label=label,
                       edgecolor=edgecolor, linewidth=linewidth, alpha=alpha)
        for label in labels if label in color_dict
    ]
    legend = fig.legend(
        handles=handles, title=title, loc=loc,
        bbox_to_anchor=(x, y), bbox_transform=fig.transFigure,
        borderaxespad=0, borderpad=0, frameon=False,
        handlelength=0.8, handleheight=0.8,
        fontsize=fontsize, title_fontsize=fontsize,
    )
    legend.get_title().set_fontweight('bold')
    legend._legend_box.align = "left"
    return legend


def legend_bottom_y(legend, fig):
    """Bottom y of a legend in figure coords (call after fig.canvas.draw())."""
    return fig.transFigure.inverted().transform(
        legend.get_window_extent(fig.canvas.get_renderer())
    )[0][1]


# --- Argparse

parser = argparse.ArgumentParser(
    description="Cohort overview of TLS number, size, clustering and topology",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
)
parser.add_argument("--project_dir", required=True, help="Project directory (root)")
parser.add_argument("-r", "--res_dir", help="Path to results directory")
parser.add_argument("--sample_sheet", help="Path to sample sheet CSV file")
parser.add_argument("--TLS_threshold", type=float, default=0.01)

args = vars(parser.parse_args())

project_dir = args["project_dir"]
res_dir = args["res_dir"] or os.path.join(project_dir, "results")
TLS_threshold = args["TLS_threshold"]
TLS_threshold_adj = str(TLS_threshold).split('.')[-1]


# --- Sample sheet + tissue metadata + subtype colors

if args.get("sample_sheet"):
    sample_sheet_path = args["sample_sheet"]
else:
    sample_sheet_path = os.path.join(project_dir, 'config', 'sample_sheet_visium.csv')

sample_sheet = pd.read_csv(sample_sheet_path)
tissue_metadata = pd.read_csv(os.path.join(project_dir, 'resources', 'tissue_metadata.csv'))
sample_meta = pd.merge(
    sample_sheet, tissue_metadata[['tissue_id', 'clinical_subtype']], on='tissue_id',
)
subtype_dict = dict(zip(sample_meta['section_id'], sample_meta['clinical_subtype']))
donor_dict = dict(zip(sample_sheet['section_id'], sample_sheet['donor_id']))

color_lut = pd.read_csv(os.path.join(project_dir, 'resources', 'plotting_colors',
                                     'SCSubtype_colors.csv')).set_index('subtype')
subtype_color_dict = color_lut['color'].to_dict()

subtype_order = ['Luminal A', 'Luminal B', 'Luminal B (HER2+)', 'HER2+', 'TNBC']


# --- Locate inputs
# Number / size / clustering / topology are all reported on the discrete-TLS set
# (Figure S3B definition) written by 01 and 02. Total spot counts per section come
# from the per-sample TLS-model CSVs.

def find_most_recent(parent):
    candidates = []
    for d in os.listdir(parent):
        full = os.path.join(parent, d)
        if not os.path.isdir(full):
            continue
        try:
            datetime.strptime(d, '%Y%m%d')
            candidates.append(d)
        except ValueError:
            continue
    if not candidates:
        raise FileNotFoundError(f"No dated subdirectory found in {parent}")
    return os.path.join(parent, max(candidates, key=lambda x: datetime.strptime(x, '%Y%m%d')))


tls_model_parent = os.path.join(
    res_dir, '03_TLS_niche_discovery', '01_identify_TLS', '01_run_TLS_model',
)
tls_model_dir = find_most_recent(tls_model_parent)
print(f"Per-sample TLS_model from: {tls_model_dir}")

# Per-TLS distances/topology (produced by 01_TLS_distances_topology.py); optional.
distances_path = os.path.join(
    res_dir, '03_TLS_niche_discovery', '01_identify_TLS',
    '03_TLS_spatial_organization', '01_TLS_distances_topology',
    f'threshold_{TLS_threshold_adj}', 'csv_files', 'distances_per_TLS.csv',
)
if os.path.exists(distances_path):
    distances_df = pd.read_csv(distances_path)
    print(f"Loaded distances for {distances_df['sample_id'].nunique()} samples.")
else:
    distances_df = pd.DataFrame(
        columns=['sample_id', 'subtype', 'TLS_group', 'size_spots',
                 'dist_to_next_TLS_um', 'dist_to_cancer_um', 'tls_topology_2hop'],
    )
    print("No distances_per_TLS.csv found — topology panels will be empty.")

# Per-TLS spatial-clustering labels (produced by 02_TLS_clustering.py); optional.
clustering_path = os.path.join(
    res_dir, '03_TLS_niche_discovery', '01_identify_TLS',
    '03_TLS_spatial_organization', '02_TLS_clustering',
    f'threshold_{TLS_threshold_adj}', 'csv_files', 'clustering_per_TLS.csv',
)
if os.path.exists(clustering_path):
    clustering_df = pd.read_csv(clustering_path)
    print(f"Loaded clustering labels for {clustering_df['sample_id'].nunique()} samples.")
else:
    clustering_df = pd.DataFrame(
        columns=['sample_id', 'subtype', 'TLS_group', 'size_spots',
                 'dbscan_label_eps500', 'dbscan_class_eps500'],
    )
    print("No clustering_per_TLS.csv found — clustering panel will be empty.")


# --- Build per-sample summary + per-TLS sizes from the discrete-TLS tables
# One row per discrete TLS with its size (spots). Prefer the clustering table
# (carries size + DBSCAN labels); fall back to the distances table.

_disc_cols = ['sample_id', 'subtype', 'TLS_group', 'size_spots']
if not clustering_df.empty and 'size_spots' in clustering_df.columns:
    discrete_tls_tbl = clustering_df[_disc_cols].copy()
elif not distances_df.empty and 'size_spots' in distances_df.columns:
    discrete_tls_tbl = distances_df[_disc_cols].copy()
else:
    discrete_tls_tbl = pd.DataFrame(columns=_disc_cols)

slide_size_dict = dict(zip(sample_sheet['section_id'], sample_sheet['slide_size']))
area_small = 6.5 * 6.5  # mm^2 for 6.5 mm capture areas
area_large = 11 * 11

summary_rows = []
sizes_rows = []

for section_id, subtype in subtype_dict.items():
    tls_file = os.path.join(tls_model_dir, f'{section_id}_TLS_model.csv')
    if not os.path.exists(tls_file):
        continue

    total_spots = len(pd.read_csv(tls_file, index_col=0))
    slide_size = slide_size_dict.get(section_id)
    area = area_large if slide_size == '11 mm' else area_small

    sub = discrete_tls_tbl[discrete_tls_tbl['sample_id'] == section_id]
    n_aggregates = len(sub)
    median_size = float(sub['size_spots'].median()) if n_aggregates else np.nan

    summary_rows.append({
        'sample_id': section_id,
        'subtype': subtype,
        'donor_id': donor_dict.get(section_id),
        'slide_size': slide_size,
        'area': area,
        'total_spots': total_spots,
        'nunique': n_aggregates,
        'normalized_n': n_aggregates / area * 100,
        'n_normalized_total_spots': n_aggregates / total_spots * 5000 if total_spots else 0,
        'median_group_size': median_size,
    })

    for _, r in sub.iterrows():
        sizes_rows.append({
            'sample_id': section_id,
            'subtype': subtype,
            'group_singleton_removed': r['TLS_group'],
            'count': int(r['size_spots']),
        })

summary_df = (pd.DataFrame(summary_rows) if summary_rows
              else pd.DataFrame(columns=['sample_id', 'subtype', 'donor_id',
                                         'slide_size', 'area', 'total_spots',
                                         'nunique', 'normalized_n',
                                         'n_normalized_total_spots',
                                         'median_group_size']))
sizes_df = (pd.DataFrame(sizes_rows) if sizes_rows
            else pd.DataFrame(columns=['sample_id', 'subtype',
                                       'group_singleton_removed', 'count']))

# Drop rows whose section_id isn't in the current sample sheet (extra safety)
summary_df = summary_df.dropna(subset=['subtype']).copy()
sizes_df = sizes_df.dropna(subset=['subtype']).copy()

print(f"Built summary for {len(summary_df)} samples, {len(sizes_df)} discrete TLS.")

# Always plot the canonical 5-subtype order; keep empty slots for subtypes
# with no observations.
palette = {s: subtype_color_dict[s] for s in subtype_order}


# --- Styling

sns.set_style('ticks')
mpl.rcParams.update({'font.size': 10, 'axes.titleweight': 'bold'})

bar_alpha = 0.7
box_face_alpha = 0.4   # softer fill so the dots behind it stay visible
err_kws = {'color': 'black', 'linewidth': 0.8}
HORIZ_YLABEL = dict(rotation=0, ha='right', va='center', labelpad=10)


def _soften_bars(ax):
    """Pastel-bar + dark-edge look, applied to every patch on the axes."""
    for bar in ax.patches:
        face = bar.get_facecolor()
        if len(face) != 4:
            continue
        r, g, b, _ = face
        bar.set_facecolor((r, g, b, bar_alpha))
        bar.set_edgecolor('black')


# --- Per-sample sorted overview (one bar per section, coloured by subtype)
# Stack from bottom up:
#   bar (count) → size box → TLS clustering stack → TLS topology stack
#   → subtype strip → capture-area strip + title.

capture_area_colors = {'6.5 mm': '#d9d9d9', '11 mm': '#525252'}
topo_categories = ['intratumoral', 'peritumoral', 'stromal_distant']
topo_colors = {
    'intratumoral':    '#b30000',
    'peritumoral':     '#fc8d59',
    'stromal_distant': '#74add1',
}
# Pretty labels for the legend (data stays lowercase_with_underscores).
topo_display_labels = {
    'intratumoral':    'Intratumoral',
    'peritumoral':     'Peritumoral',
    'stromal_distant': 'Stromal/distant',
}


TOPO_COL = 'tls_topology_2hop'  # 2-hop ring: paired with the stringent >=30% / >50% thresholds


def _draw_topology_stack(ax, dist_df, sample_order, ylabel):
    """Stacked proportion bar per sample: intratumoral / peritumoral / stromal."""
    n = len(sample_order)
    bottoms = np.zeros(n)
    for cat in topo_categories:
        heights = np.zeros(n)
        for i, s in enumerate(sample_order):
            sub = dist_df[dist_df['sample_id'] == s]
            if len(sub):
                heights[i] = (sub[TOPO_COL] == cat).sum() / len(sub) * 100
        ax.bar(range(n), heights, bottom=bottoms,
               color=topo_colors[cat], edgecolor='black', linewidth=0.3,
               label=cat)
        bottoms += heights
    ax.set_ylim(0, 100)
    ax.set_ylabel(ylabel, **HORIZ_YLABEL)
    sns.despine(ax=ax)


def _draw_box_per_sample(ax, source_df, value_col, sample_order, bar_colors,
                         ylabel, log=False):
    """Boxplot of `value_col` per sample, with strip overlay of raw values."""
    box_data = [
        source_df.loc[source_df['sample_id'] == s, value_col].dropna().values
        for s in sample_order
    ]
    bp = ax.boxplot(
        box_data, positions=range(len(sample_order)), widths=0.7,
        patch_artist=True, showfliers=False,
        medianprops={'color': 'black', 'linewidth': 1.2},
    )
    # Recolour each box (fill semi-transparent, outline in the subtype colour)
    # so identity stays readable despite the faded fill. Whiskers stay black.
    for i, color in enumerate(bar_colors):
        r, g, b, _ = mpl.colors.to_rgba(color)
        bp['boxes'][i].set_facecolor((r, g, b, box_face_alpha))
        bp['boxes'][i].set_edgecolor(color)
        bp['boxes'][i].set_linewidth(1.0)

    rng = np.random.default_rng(0)
    for i, vals in enumerate(box_data):
        if len(vals) == 0:
            continue
        jitter = (rng.random(len(vals)) - 0.5) * 0.3
        # zorder=1: dots draw behind the box (default ~2). The semi-transparent
        # box fill keeps them visible while preserving the outline.
        ax.scatter(np.full(len(vals), i) + jitter, vals,
                   s=8, color='black', linewidth=0, zorder=1)

    if log:
        ax.set_yscale('log')
    ax.set_ylabel(ylabel, **HORIZ_YLABEL)
    ax.tick_params(labelbottom=False)
    sns.despine(ax=ax)


# --- TLS spatial-clustering composition (from 02_TLS_clustering.py)
# DBSCAN at eps=500 µm on edge-to-edge (min spot-to-spot) distances. Each
# distinct cluster is drawn as its own shade of blue (darkest = largest); every
# lone TLS is amber. Legend collapses to two entries: clustered + lone.

CLUSTER_EPS        = 500
_cluster_label_col = f'dbscan_label_eps{CLUSTER_EPS}'
_cluster_class_col = f'dbscan_class_eps{CLUSTER_EPS}'

_blues = mpl.colormaps['Blues']
CLUSTER_SHADES       = [mpl.colors.to_hex(_blues(v)) for v in np.linspace(0.85, 0.35, 8)]
CLUSTER_LEGEND_COLOR = CLUSTER_SHADES[1]   # representative blue for the legend swatch
LONE_COLOR           = '#fdae61'           # amber


def cluster_stack(sample_id):
    """Ordered (heights, colors) for the fraction-normalised clustering stack.
    Bottom -> top: cluster 1 (largest) .. cluster k as blue shades, then one
    amber segment per lone TLS. Every segment is a single spatial group, so its
    height is that group's TLS count (a cluster's members, or 1 for a lone TLS);
    the amber block is thus segmented into individual TLS by black lines, exactly
    as the blue clustered region is segmented into distinct clusters."""
    sub = clustering_df[clustering_df['sample_id'] == sample_id]
    if sub.empty or _cluster_class_col not in sub.columns:
        return [], []
    heights, colors = [], []
    clustered = sub[sub[_cluster_class_col] == 'clustered']
    if len(clustered):
        sizes = (clustered[_cluster_label_col].astype(int)
                 .value_counts().sort_values(ascending=False))
        for rank, count in enumerate(sizes.values):
            heights.append(int(count))
            colors.append(CLUSTER_SHADES[rank % len(CLUSTER_SHADES)])
    n_lone = int((sub[_cluster_class_col] == 'lone').sum())
    # One unit-height band per lone TLS (a lone TLS is a group of size 1) so the
    # amber block is divided by black lines into individual TLS, mirroring the
    # per-cluster segmentation of the blue clustered region.
    for _ in range(n_lone):
        heights.append(1)
        colors.append(LONE_COLOR)
    total = sum(heights)
    if total:
        heights = [h / total for h in heights]
    return heights, colors


def _draw_cluster_stack(ax, sample_order, ylabel):
    """Per-sample fraction-normalised stacked bar of clustered (blue shades,
    one per distinct cluster) vs lone (amber) TLS. The number of distinct
    clusters reads straight off the number of blue shades in each bar."""
    for i, s in enumerate(sample_order):
        heights, colors = cluster_stack(s)
        bottom = 0.0
        for h, c in zip(heights, colors):
            ax.bar(i, h, bottom=bottom, color=c, edgecolor='black',
                   linewidth=0.3, width=0.8)
            bottom += h
    ax.set_ylim(0, 1.05)
    ax.set_yticks([0, 0.5, 1.0])
    ax.set_yticklabels(['0', '0.5', '1'])
    ax.set_ylabel(ylabel, **HORIZ_YLABEL)
    ax.tick_params(labelbottom=False)
    sns.despine(ax=ax)


def per_sample_barplot(df, plot_key, ylabel, title, round_annot=True,
                       label_map=None):
    df = df.sort_values(by=plot_key, ascending=False).copy()
    sample_order = df['sample_id'].tolist()
    display_labels = [
        label_map.get(s, s) if label_map else s for s in sample_order
    ]
    # Compact replicate suffix on the x-tick labels: '_Rep1' -> '_R1'.
    display_labels = [lbl.replace('_Rep', '_R') for lbl in display_labels]
    n = len(sample_order)
    bar_colors = [palette[s] for s in df['subtype']]
    slide_colors_per_sample = [capture_area_colors.get(sz, '#cccccc')
                               for sz in df['slide_size']]

    _space_per_bar = 0.04  # figure-fraction units per sample
    _ax_w          = _space_per_bar * n
    _bar_h         = 0.20
    _box_h         = 0.20
    _strip_h       = 0.04
    _gap           = 0.03
    _gap_strip     = 0.01
    _x0            = 0.10
    _y0_bar        = 0.30

    # Layout (top → bottom of the figure): cap strip → subt strip →
    # count bar → size box → clustering stack → topology stack. xticklabels go
    # on the topology stack (bottommost data panel).
    _y0_topo       = _y0_bar
    _y0_cluster    = _y0_topo       + _box_h   + _gap
    _y0_size_box   = _y0_cluster    + _box_h   + _gap
    _y0_count_bar  = _y0_size_box   + _box_h   + _gap
    _y0_subt_strip = _y0_count_bar  + _bar_h   + _gap_strip
    _y0_cap_strip  = _y0_subt_strip + _strip_h + _gap_strip

    fig = plt.figure(figsize=(4, 3))
    ax_topo        = fig.add_axes([_x0, _y0_topo,        _ax_w, _box_h])
    ax_cluster     = fig.add_axes([_x0, _y0_cluster,     _ax_w, _box_h],   sharex=ax_topo)
    ax_size        = fig.add_axes([_x0, _y0_size_box,    _ax_w, _box_h],   sharex=ax_topo)
    ax_bar         = fig.add_axes([_x0, _y0_count_bar,   _ax_w, _bar_h],   sharex=ax_topo)
    ax_subt_strip  = fig.add_axes([_x0, _y0_subt_strip,  _ax_w, _strip_h], sharex=ax_topo)
    ax_cap_strip   = fig.add_axes([_x0, _y0_cap_strip,   _ax_w, _strip_h], sharex=ax_topo)

    def _draw_strip(ax, colors_per_sample, edge='black'):
        for i, color in enumerate(colors_per_sample):
            ax.add_patch(mpatches.Rectangle((i - 0.4, 0), 0.8, 1,
                                            facecolor=color, edgecolor=edge,
                                            linewidth=0.3))
        ax.set_ylim(0, 1)
        ax.set_yticks([])
        ax.tick_params(axis='x', labelbottom=False, bottom=False)
        for spine in ax.spines.values():
            spine.set_visible(False)

    _draw_strip(ax_cap_strip, slide_colors_per_sample)
    _draw_strip(ax_subt_strip, bar_colors)
    n_samples = len(df)
    n_tls     = int(df['nunique'].sum())
    ax_cap_strip.set_title(
        f'{title}\nn = {n_samples} sections, {n_tls} TLS',
        fontweight='bold', pad=4,
    )

    # TLS topology classification (2-hop ring); see 01_TLS_distances_topology.py.
    _draw_topology_stack(ax_topo, distances_df, sample_order,
                         ylabel='TLS topology\n(% of TLS)')
    # TLS spatial-clustering composition (DBSCAN eps=500 µm): clustered (blue
    # shades, one per distinct cluster) vs lone (amber). A directly interpretable
    # measure of whether TLS occur in little clusters or spread out as lone TLS.
    _draw_cluster_stack(ax_cluster, sample_order,
                        ylabel='TLS clustering\n(clustered / lone)')
    _draw_box_per_sample(ax_size, sizes_df, 'count',
                         sample_order, bar_colors,
                         ylabel='TLS size\n(# spots, log)', log=True)

    # --- Count bar (sits at the top of the data panels)
    ax_bar.bar(
        range(n), df[plot_key].values, color=bar_colors,
        edgecolor='black', linewidth=0.6, alpha=bar_alpha,
    )
    for i, v in enumerate(df[plot_key].values):
        label = f'{v:.1f}' if round_annot else f'{int(v)}'
        ax_bar.text(i, v, label, ha='center', va='bottom', fontsize=9)

    ax_bar.set_xlim(-0.6, n - 0.4)
    ax_bar.set_ylim(0, df[plot_key].max() * 1.25)
    ax_bar.set_xlabel('')
    ax_bar.tick_params(axis='x', labelbottom=False, bottom=False)

    ax_bar.set_ylabel(ylabel, **HORIZ_YLABEL)
    sns.despine(ax=ax_bar)

    # xticklabels go on the topology ax (bottommost data panel)
    ax_topo.set_xticks(range(n))
    ax_topo.set_xticklabels(display_labels, rotation=90)
    ax_topo.tick_params(axis='x', labelbottom=True, bottom=True)

    # Stack the legends: capture-area, then subtype, then clustering, then
    # topology, all anchored to a common left edge just right of the widest ax.
    right_edge, top_edge = position_legend(
        [ax_bar, ax_size, ax_cluster, ax_topo,
         ax_subt_strip, ax_cap_strip],
    )
    leg_x = right_edge + 0.02
    leg_pad = 0.02

    capture_legend = add_color_legend(
        fig, labels=list(capture_area_colors.keys()),
        color_dict=capture_area_colors, title='Capture area',
        x=leg_x, y=top_edge, loc='upper left',
    )
    fig.canvas.draw()

    subtype_legend = add_color_legend(
        fig, labels=subtype_order, color_dict=palette, title='Subtype',
        x=leg_x, y=legend_bottom_y(capture_legend, fig) - leg_pad,
        loc='upper left', alpha=bar_alpha,
    )
    fig.canvas.draw()

    # TLS clustering: two entries only (shades in the bar = distinct clusters).
    clustered_label = f'Clustered (<{CLUSTER_EPS} µm)'
    clustering_legend = add_color_legend(
        fig, labels=[clustered_label, 'Lone'],
        color_dict={clustered_label: CLUSTER_LEGEND_COLOR, 'Lone': LONE_COLOR},
        title='TLS clustering',
        x=leg_x, y=legend_bottom_y(subtype_legend, fig) - leg_pad,
        loc='upper left',
    )
    fig.canvas.draw()

    topo_legend_labels = [topo_display_labels[c] for c in topo_categories]
    topo_legend_colors = {topo_display_labels[c]: topo_colors[c]
                          for c in topo_categories}
    add_color_legend(
        fig, labels=topo_legend_labels, color_dict=topo_legend_colors,
        title='TLS topology',
        x=leg_x, y=legend_bottom_y(clustering_legend, fig) - leg_pad,
        loc='upper left',
    )

    return fig


per_sample_barplot(summary_df, 'nunique', 'Number of TLS',
                   'TLS per sample (sorted)', round_annot=False)


# --- Save outputs

out_dir = os.path.join(
    res_dir, '03_TLS_niche_discovery', '01_identify_TLS',
    '03_TLS_spatial_organization', '03_TLS_number_size_overview',
    f'threshold_{TLS_threshold_adj}',
)
out_plots = os.path.join(out_dir, 'plots')
out_csv = os.path.join(out_dir, 'csv_files')
os.makedirs(out_plots, exist_ok=True)
os.makedirs(out_csv, exist_ok=True)

save_multi_image(os.path.join(out_plots, 'TLS_number_size_overview.pdf'))
plt.close('all')

# Save the tables we plotted (subset to relevant columns)
summary_out = summary_df[[
    'sample_id', 'subtype', 'donor_id', 'slide_size', 'area', 'total_spots',
    'nunique', 'n_normalized_total_spots', 'median_group_size',
]].rename(columns={'nunique': 'n_aggregates'})
summary_out.to_csv(os.path.join(out_csv, 'summary_per_sample.csv'), index=False)

sizes_df[['sample_id', 'subtype', 'group_singleton_removed', 'count']].rename(
    columns={'group_singleton_removed': 'aggregate_id', 'count': 'size_spots'},
).to_csv(os.path.join(out_csv, 'sizes_per_aggregate.csv'), index=False)

print("Script successfully completed")
