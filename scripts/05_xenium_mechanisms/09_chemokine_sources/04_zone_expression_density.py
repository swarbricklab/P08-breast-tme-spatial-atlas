#!/usr/bin/env python3

# --- Aim
# Extended Data Fig. 9e: cohort library-normalised (CP10K) expression density of the niche chemokines CXCL13 &
# CCL19 across three zone categories, per mm² of occupied 50 µm bins:
#   - B_cell_zone   (TLS, B zone)
#   - T_cell_zone   (TLS, T zone)
#   - non_TLS       (everything outside the classified TLS bins; per-sample subsampled to
#                    roughly match the B+T bin count as a baseline)
# One panel per gene, per-sample dots, two-sided paired Wilcoxon signed-rank (n=6) for
# B-T, B-non_TLS, T-non_TLS. Density = sum(filtered_counts) / bin area; the `filtered_counts`
# layer is library-normalised to 1e4/cell (CP10K), so this is CP10K signal / mm².

# --- Define environment
# env: xen-hev

# --- Load packages

import argparse
import os

import numpy as np
import pandas as pd
import anndata as ad
import scipy.sparse as sp
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
from scipy import stats

mpl.rcParams['pdf.fonttype'] = 42
TICK_FS, SUP_FS, ANNO_FS = 8, 12, 7


# --- Argparse arguments
# Cohort: dispatches through submit_multi_sample_job_nci.sh -> run_multi_sample_job.sh, which
# passes --h5ad_dir --project_dir --res_dir --sample_sheet --samples. The Methods parameters
# (50 µm grid, B/T threshold 0.3, CP10K, avg_BT baseline) are hard-coded defaults.

parser = argparse.ArgumentParser(
    description="Extended Data Fig. 9e CXCL13/CCL19 lib-norm (CP10K) expression density per zone (cohort)",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("-s", "--samples", nargs='+', required=True)
parser.add_argument("--h5ad_dir", required=True)
parser.add_argument("--project_dir", required=True)
parser.add_argument("--res_dir", required=True)
parser.add_argument("--sample_sheet", default=None, help="Path to sample sheet CSV (else default)")
parser.add_argument("--grid_um", type=float, default=50.0)
parser.add_argument("--threshold", type=float, default=0.3)
parser.add_argument("--genes", nargs='+', default=['CXCL13', 'CCL19'])
parser.add_argument("--seed", type=int, default=42)
parser.add_argument("--subsample_strategy", default='avg_BT',
                    choices=['avg_BT', 'sum_BT', 'max_BT', 'min_BT'],
                    help="Per-sample non_TLS subsample size: avg/sum/max/min of n_B and n_T")
args = vars(parser.parse_args())
samples            = args["samples"]
h5ad_dir           = args["h5ad_dir"]
project_dir        = args["project_dir"]
res_dir            = args["res_dir"]
grid_um            = args["grid_um"]
threshold          = args["threshold"]
genes              = args["genes"]
seed               = args["seed"]
subsample_strategy = args["subsample_strategy"]

if args.get("sample_sheet"):
    sample_sheet_path = args["sample_sheet"]
else:
    sample_sheet_path = os.path.join(project_dir, 'config', 'sample_sheet_xenium.csv')

# Upstream per-bin zone labels (B/T/mixed) from 07_TLS_zonation/01_classify_zones.py.
# NB the 50 µm bin grid here MUST match the one 07/01 used (same coordinate source =
# x_centroid/y_centroid, same origin = per-sample min, same grid_um) so bin_ix/bin_iy align.
bin_csv_dir = os.path.join(res_dir, '05_xenium_mechanisms', '07_TLS_zonation',
                           '01_classify_zones', 'csv_files')
bin_area_mm2 = (grid_um ** 2) / 1e6

zones_of_i = ['B_cell_zone', 'T_cell_zone', 'non_TLS']
zone_palette = {
    'B_cell_zone': mpl.colors.cnames['cornflowerblue'],
    'T_cell_zone': mpl.colors.cnames['crimson'],
    'non_TLS':     'dimgray',
}


# --- Subtype lookup for ordering

sample_sheet    = pd.read_csv(sample_sheet_path)
tissue_metadata = pd.read_csv(os.path.join(project_dir, 'resources', 'tissue_metadata.csv'))
merged          = sample_sheet.merge(tissue_metadata, on='tissue_id')
subtype_dict    = dict(zip(merged['section_id'], merged['clinical_subtype']))


# --- Per-sample binning: lib-norm (CP10K) chemokine sums per 50 µm bin, tagged by zone

def load_sample(sample):
    bin_csv = os.path.join(bin_csv_dir, f'{sample}_BT_zone_bins.csv')
    if not os.path.isfile(bin_csv):
        print(f'  WARNING: missing {bin_csv} — skipping {sample}')
        return None
    bin_zone_df = pd.read_csv(bin_csv)[['bin_ix', 'bin_iy', 'TLS_zone']].rename(
        columns={'TLS_zone': 'zone'})

    adata = ad.read_h5ad(os.path.join(h5ad_dir, 'xenium', f'{sample}.h5ad'))
    xy = adata.obs[['x_centroid', 'y_centroid']].to_numpy()
    xmin, ymin = xy.min(axis=0)
    xmax, ymax = xy.max(axis=0)
    nx = int(np.ceil((xmax - xmin) / grid_um))
    ny = int(np.ceil((ymax - ymin) / grid_um))
    ix = np.minimum(((xy[:, 0] - xmin) / grid_um).astype(int), nx - 1)
    iy = np.minimum(((xy[:, 1] - ymin) / grid_um).astype(int), ny - 1)

    cell_df = pd.DataFrame({'bin_ix': ix, 'bin_iy': iy,
                            'cell_id': adata.obs.index.values})
    for g in genes:
        if g in adata.var_names:
            gi = np.where(adata.var_names == g)[0][0]
            v = adata.layers['filtered_counts'][:, gi]
            v = v.toarray().ravel() if sp.issparse(v) else np.asarray(v).ravel()
            cell_df[g] = v.astype(float)
        else:
            cell_df[g] = 0.0

    agg = (cell_df.groupby(['bin_ix', 'bin_iy'], observed=True)
                 .agg(n_cells=('bin_ix', 'size'),
                      **{g: (g, 'sum') for g in genes})
                 .reset_index())

    agg = agg.merge(bin_zone_df, on=['bin_ix', 'bin_iy'], how='left')
    agg['zone'] = agg['zone'].fillna('non_TLS')

    agg['sample_id'] = sample
    agg['subtype']   = subtype_dict.get(sample)
    del adata
    return agg


print(f'Loading {len(samples)} samples...')
parts = [load_sample(s) for s in samples]
all_bins = pd.concat([p for p in parts if p is not None], ignore_index=True)
print(f'Total bins across cohort: {len(all_bins):,}')
print(all_bins.groupby(['sample_id', 'zone']).size().unstack(fill_value=0).to_string())


# --- Per-sample subsample of non_TLS bins (baseline matched in bin count)

rng_master = np.random.default_rng(seed)
sampled_parts = []
print(f'\nSubsampling non_TLS bins per sample (strategy={subsample_strategy}):')
for s, sub in all_bins.groupby('sample_id', observed=True):
    n_B = int((sub['zone'] == 'B_cell_zone').sum())
    n_T = int((sub['zone'] == 'T_cell_zone').sum())
    n_other = int((sub['zone'] == 'non_TLS').sum())
    if subsample_strategy == 'avg_BT':
        target = int(round((n_B + n_T) / 2))
    elif subsample_strategy == 'sum_BT':
        target = n_B + n_T
    elif subsample_strategy == 'max_BT':
        target = max(n_B, n_T)
    else:
        target = min(n_B, n_T)
    target = min(target, n_other)
    non_TLS_bins = sub[sub['zone'] == 'non_TLS']
    if target > 0 and len(non_TLS_bins) > 0:
        sampled = non_TLS_bins.sample(
            n=target,
            random_state=int(rng_master.integers(0, 2 ** 31 - 1)))
    else:
        sampled = non_TLS_bins.iloc[:0]
    sampled_parts.append(pd.concat([sub[sub['zone'] != 'non_TLS'], sampled],
                                   ignore_index=True))
    print(f'  {s:<18s} B={n_B:>5d}  T={n_T:>5d}  '
          f'non_TLS(pool)={n_other:>5d}  non_TLS(sampled)={len(sampled):>5d}')

bins = pd.concat(sampled_parts, ignore_index=True)
bins = bins[bins['zone'].isin(zones_of_i)].copy()


# --- Per-sample × zone density (lib-norm CP10K signal / mm²)

density = (bins.groupby(['sample_id', 'subtype', 'zone'], observed=True)
                .agg(n_bins=('n_cells', 'size'),
                     **{g: (g, 'sum') for g in genes})
                .reset_index())
density['area_mm2'] = density['n_bins'] * bin_area_mm2
for g in genes:
    density[f'{g}_per_mm2'] = density[g] / density['area_mm2']


# --- Stats helpers

def stat_label(p):
    if p < 1e-4: return '****'
    if p < 1e-3: return '***'
    if p < 1e-2: return '**'
    if p < 5e-2: return '*'
    return 'ns'


def paired_w(vA, vB, sA, sB):
    """Two-sided paired Wilcoxon signed-rank across the samples common to both zones."""
    common = sorted(set(sA) & set(sB))
    if len(common) < 3:
        return 1.0, 'ns'
    vAa = np.array([vA[list(sA).index(s)] for s in common])
    vBa = np.array([vB[list(sB).index(s)] for s in common])
    try:
        p = float(stats.wilcoxon(vAa, vBa, zero_method='wilcox').pvalue)
    except ValueError:
        p = 1.0
    return p, stat_label(p)


# --- Plot helpers

def boxplot_panel(ax, values_per_x, colors, x_labels, samples_per_x=None):
    """Boxplot with overlaid jittered dots. If `samples_per_x` is provided, also draws gray
    paired-sample connecting lines across x-positions using deterministic jitter."""
    bp = ax.boxplot([np.asarray(v) for v in values_per_x],
                    positions=list(range(len(x_labels))), widths=0.55,
                    patch_artist=True, showfliers=False,
                    medianprops=dict(color='black', linewidth=1.2),
                    zorder=2)
    for patch, color in zip(bp['boxes'], colors):
        rgb = mpl.colors.to_rgb(color)
        patch.set_facecolor((*rgb, 0.55))
        patch.set_edgecolor(color)
        patch.set_linewidth(0.8)
    for elem in (*bp['whiskers'], *bp['caps']):
        elem.set_color('black')
        elem.set_linewidth(0.8)

    if samples_per_x is not None:
        rng = np.random.default_rng(seed=42)
        all_samples = sorted({s for sl in samples_per_x for s in sl})
        jitter = {(s, xi): float(rng.uniform(-0.12, 0.12))
                  for xi in range(len(x_labels)) for s in all_samples}
        per_x_lookup = [dict(zip(sl, vl)) for sl, vl in
                        zip(samples_per_x, values_per_x)]
        for s in all_samples:
            vals = [d.get(s, np.nan) for d in per_x_lookup]
            if all(np.isfinite(v) for v in vals):
                xs = [xi + jitter[(s, xi)] for xi in range(len(x_labels))]
                ax.plot(xs, vals, color='gray', alpha=0.4, linewidth=0.6,
                        zorder=1)
        for xi, (sl, vl) in enumerate(zip(samples_per_x, values_per_x)):
            xs = [xi + jitter[(s, xi)] for s in sl]
            ax.scatter(xs, vl, c='black', s=14, linewidths=0, zorder=3)
    else:
        rng = np.random.default_rng(seed=42)
        for xi, vals in enumerate(values_per_x):
            jit = rng.uniform(-0.12, 0.12, size=len(vals))
            ax.scatter(xi + jit, vals, c='black', s=14, linewidths=0, zorder=3)

    ax.set_xticks(range(len(x_labels)))
    ax.set_xticklabels(x_labels, rotation=90, fontsize=TICK_FS)
    ax.set_xlim(-0.5, len(x_labels) - 0.5)
    ax.tick_params(axis='y', labelsize=TICK_FS)
    ax.grid(False)
    for sp_ in ['top', 'right']:
        ax.spines[sp_].set_visible(False)


def add_stat_bracket(ax, x1, x2, y_top, label, line_h_frac=0.025, lw=0.6,
                     fontsize=None):
    yr = ax.get_ylim()
    line_h = (yr[1] - yr[0]) * line_h_frac
    ax.plot([x1, x1, x2, x2], [y_top, y_top + line_h, y_top + line_h, y_top],
            color='black', linewidth=lw, clip_on=False)
    fig = ax.figure
    axis_h_in = ax.get_position().height * fig.get_figheight()
    data_per_inch = (yr[1] - yr[0]) / axis_h_in if axis_h_in > 0 else 0
    glyph_drop = ((fontsize or 10) * 0.30 / 72.0) * data_per_inch
    ax.text((x1 + x2) / 2, y_top + line_h - glyph_drop, label,
            ha='center', va='baseline', fontsize=fontsize, clip_on=False)


# --- Multi-panel density figure builder (one panel per gene)

x_labels = ['B cell zone', 'T cell zone', 'non-TLS']


def make_density_figure(panel_specs, ylabel, suptitle):
    """panel_specs: list of (title, metric_col) tuples."""
    n = len(panel_specs)
    panel_w_in   = 0.9
    wspace_in    = 0.35
    margin_left  = 0.85
    margin_right = 0.55
    fig_w_in     = margin_left + n * panel_w_in + (n - 1) * wspace_in + margin_right
    fig_h_in     = 2.9

    fig = plt.figure(figsize=(fig_w_in, fig_h_in))
    panel_w = panel_w_in / fig_w_in
    wspace  = wspace_in  / fig_w_in
    panel_h, panel_y = 0.46, 0.26
    p_left = margin_left / fig_w_in
    axes = []

    for i, (title, metric) in enumerate(panel_specs):
        ax = fig.add_axes([p_left, panel_y, panel_w, panel_h])
        vals  = [density.loc[density['zone'] == z, metric].values for z in zones_of_i]
        samps = [density.loc[density['zone'] == z, 'sample_id'].tolist() for z in zones_of_i]

        boxplot_panel(ax, vals, [zone_palette[z] for z in zones_of_i], x_labels,
                      samples_per_x=samps)
        ax.set_title(title, fontsize=11, fontweight='bold', pad=14)
        ax.set_xlabel('Tissue region', fontsize=TICK_FS, fontweight='bold')
        if i == 0:
            ax.set_ylabel(ylabel, fontsize=TICK_FS, fontweight='bold')

        sci_fmt = ScalarFormatter(useMathText=True)
        sci_fmt.set_scientific(True)
        sci_fmt.set_powerlimits((-2, 3))
        ax.yaxis.set_major_formatter(sci_fmt)
        ax.yaxis.get_offset_text().set_fontsize(TICK_FS)

        ymax_data = max(np.max(v) if len(v) else 0 for v in vals)
        pair_pos = [(0, 1), (0, 2), (1, 2)]
        bracket_info = []
        for (a, b) in pair_pos:
            p_val, lab = paired_w(vals[a], vals[b], samps[a], samps[b])
            print(f'  {title}: {zones_of_i[a]} vs {zones_of_i[b]} -> p={p_val:.3g} ({lab})')
            if lab != 'ns':
                bracket_info.append((a, b, lab))
        n_brackets = len(bracket_info)
        if n_brackets > 0:
            step  = 0.11    # vertical spacing between stacked brackets
            first = 0.10    # gap between top whisker and first bracket
            ax.set_ylim(0, ymax_data * (1.05 + first + step * n_brackets))
            for k, (a, b, lab) in enumerate(bracket_info):
                y_top = ymax_data * (1.0 + first + step * k)
                add_stat_bracket(ax, a, b, y_top, lab,
                                 line_h_frac=0.020, fontsize=TICK_FS + 2)
        else:
            ax.set_ylim(bottom=0)
        axes.append(ax)
        p_left += panel_w + wspace

    fig.canvas.draw()
    title_x = (axes[0].get_position().x0 + axes[-1].get_position().x1) / 2
    top_extent = max(fig.transFigure.inverted().transform(
        ax.get_tightbbox(fig.canvas.get_renderer()))[1][1] for ax in axes)
    fig.suptitle(suptitle, x=title_x, y=top_extent + 0.04,
                 ha='center', va='bottom', fontsize=SUP_FS, fontweight='bold')
    return fig


print()
fig_transcripts = make_density_figure(
    panel_specs=[(g, f'{g}_per_mm2') for g in genes],
    ylabel='Norm. expression / mm²',
    suptitle=(f'Extended Data Fig. 9e — lib-norm (CP10K) density per zone — non-TLS baseline subsampled '
              f'per sample ({subsample_strategy})'),
)


# --- Save

out_dir = os.path.join(res_dir, '05_xenium_mechanisms', '09_chemokine_sources',
                       '04_zone_expression_density')
plot_dir = os.path.join(out_dir, 'plots')
csv_dir  = os.path.join(out_dir, 'csv_files')
os.makedirs(plot_dir, exist_ok=True)
os.makedirs(csv_dir,  exist_ok=True)

suffix = f'{subsample_strategy}_seed{seed}'
pdf_path = os.path.join(plot_dir, f'ED9E_cohort_expression_density_libnorm_{suffix}.pdf')
fig_transcripts.savefig(pdf_path, bbox_inches='tight', dpi=300)
plt.close('all')

density.to_csv(os.path.join(csv_dir, f'cohort_expression_density_{suffix}.csv'), index=False)
print(f'\nSaved -> {pdf_path}')
print('Script successfully completed')
