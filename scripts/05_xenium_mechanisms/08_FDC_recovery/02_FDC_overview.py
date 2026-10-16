#!/usr/bin/env python3

# --- Aim
# Cohort-level overview of the FDC-vs-Tfh recovery across all samples (Extended Data Fig. 9j). Reads the
# per-sample FDC/Tfh barcode CSVs written by 01_identify_FDC.py and produces a multi-page PDF:
#   1. 6-panel spatial figure of FDC (amber) and Tfh (blue) cells, background cells in grey.
#   2. Bar plots of FDC count and pct-of-original-Tfh reannotated per sample.
#   3. Pooled FDC vs Tfh violin and standard dotplot for the marker panel.
#   4. Extended Data Fig. 9j: the per-gene-normalised FDC-vs-Tfh subcluster marker dot plot (FDC CR1/CR2/FCER2,
#      the CXCL13-CXCR5 axis, pan-T CD3E/CD3G/CD4 and Tfh PDCD1/ICOS/BCL6/IL21/TIGIT), dot size
#      scaled to the per-gene max so lowly-detected genes stay visually comparable.
# Samples with no recovered FDC subcluster are handled gracefully (shown as absent / skipped).

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
from matplotlib.lines import Line2D
import seaborn as sns

mpl.rcParams['pdf.fonttype'] = 42
sc.set_figure_params(facecolor='white', fontsize=10)


# --- Argparse arguments

parser = argparse.ArgumentParser(
    description="Cohort overview of FDC vs Tfh recovery (Extended Data Fig. 9j)",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
)
parser.add_argument("--h5ad_dir",    help="Path to anndata directory")
parser.add_argument("--project_dir", help="Project root directory")
parser.add_argument("--res_dir",     help="Path to results directory")
parser.add_argument("--sample_sheet", default=None, help="Path to sample sheet CSV file")
parser.add_argument("--samples",     required=True, nargs='+', help="List of sample IDs")

args = vars(parser.parse_args())

h5ad_dir    = args["h5ad_dir"]
project_dir = args["project_dir"]
res_dir     = args["res_dir"]
samples     = args["samples"]


# --- Local helper functions and shared plotting utilities

sys.path.append(os.path.join(project_dir, 'scripts'))
from helper_functions import auto_crop_xenium, save_multi_image

sys.path.append(os.path.join(project_dir, 'scripts', '05_xenium_mechanisms'))
import plotting_functions


# --- Sample metadata

if args.get("sample_sheet"):
    sample_sheet_path = args["sample_sheet"]
else:
    sample_sheet_path = os.path.join(project_dir, 'config', 'sample_sheet_xenium.csv')

sample_sheet    = pd.read_csv(sample_sheet_path)
tissue_metadata = pd.read_csv(os.path.join(project_dir, 'resources', 'tissue_metadata.csv'))
section_to_tissue = sample_sheet.set_index('section_id')['tissue_id'].to_dict()
tissue_to_subtype = (tissue_metadata.assign(tissue_id=tissue_metadata['tissue_id'].astype(str))
                                     .set_index('tissue_id')['clinical_subtype'].to_dict())

def label_for(sample):
    subtype = tissue_to_subtype.get(str(section_to_tissue.get(sample, '')), '')
    return f'{sample}\n({subtype})' if subtype else sample


# --- Output directories

out_root = os.path.join(res_dir, '05_xenium_mechanisms', '08_FDC_recovery', 'all_samples')
csv_dir  = os.path.join(out_root, 'csv_files')
plot_dir = os.path.join(out_root, 'plots')
os.makedirs(csv_dir,  exist_ok=True)
os.makedirs(plot_dir, exist_ok=True)


# --- Load samples + FDC/Tfh barcodes

fdc_csv_dir = os.path.join(res_dir, '05_xenium_mechanisms', '08_FDC_recovery', 'per_sample', 'csv_files')

adatas = {}
counts = []
for sample in samples:
    print(f'loading {sample}')
    adata = ad.read_h5ad(os.path.join(h5ad_dir, 'xenium', f'{sample}.h5ad'))
    adata.obs['Tfh_group'] = 'other'
    bc_file = os.path.join(fdc_csv_dir, f'{sample}_FDC_Tfh_barcodes.csv')
    if os.path.exists(bc_file):
        bc = pd.read_csv(bc_file)
        fdc_ids = bc.loc[bc['group'] == 'FDC', 'cell_id']
        tfh_ids = bc.loc[bc['group'] == 'Tfh', 'cell_id']
        adata.obs.loc[adata.obs_names.isin(fdc_ids), 'Tfh_group'] = 'FDC'
        adata.obs.loc[adata.obs_names.isin(tfh_ids), 'Tfh_group'] = 'Tfh'
    else:
        print(f'  no FDC/Tfh barcode CSV for {sample} — all cells marked "other"')
    counts.append({
        'sample':            sample,
        'n_fdc':             int((adata.obs['Tfh_group'] == 'FDC').sum()),
        'n_tfh':             int((adata.obs['Tfh_group'] == 'Tfh').sum()),
        'n_total_orig_tfh':  int((adata.obs['Tfh_group'].isin(['FDC', 'Tfh'])).sum()),
        'total_cells':       adata.n_obs,
    })
    adatas[sample] = adata

counts_df = pd.DataFrame(counts)
counts_df['pct_reannotated_FDC'] = (counts_df['n_fdc'] / counts_df['n_total_orig_tfh'].replace(0, np.nan) * 100).fillna(0)
counts_df.to_csv(os.path.join(csv_dir, 'FDC_reannotation_counts.csv'), index=False)
print('\n' + counts_df.to_string(index=False))


# --- Markers grouped by functional category (Extended Data Fig. 9j panel)
GENE_GROUPS = {
    'FDC':    ['CR1', 'CR2', 'FCER2'],
    'CXCL13': ['CXCL13', 'CXCR5'],
    'T cell': ['CD3E', 'CD3G', 'CD4'],
    'Tfh':    ['PDCD1', 'ICOS', 'BCL6', 'IL21', 'TIGIT'],
}
panel_var_names = next(iter(adatas.values())).var_names
all_markers = [g for gs in GENE_GROUPS.values() for g in gs if g in panel_var_names]


# --- Plot 1: Cohort spatial (3 cols x 2 rows, shared legend)

FDC_COLOR = '#f5a623'   # warm amber-orange
TFH_COLOR = '#888EC6'   # T-cells major colour
BG_COLOR  = 'lightgrey'

S_BG  = 0.15
S_TFH = 1.2
S_FDC = 4.0

ncols, nrows = 3, 2
fig, axes = plt.subplots(nrows, ncols, figsize=(4.2 * ncols, 4.2 * nrows))
axes_flat = axes.flatten()

for i, sample in enumerate(samples):
    ax = axes_flat[i]
    adata = adatas[sample]
    coords = adata.obsm['spatial']
    grp = adata.obs['Tfh_group'].astype(str).values
    crop = auto_crop_xenium(adata)   # (x_min, x_max, y_min, y_max)
    bg_m  = grp == 'other'
    tfh_m = grp == 'Tfh'
    fdc_m = grp == 'FDC'
    ax.scatter(coords[bg_m, 0],  coords[bg_m, 1],  s=S_BG,  c=BG_COLOR,  alpha=0.4, rasterized=True, linewidths=0)
    ax.scatter(coords[tfh_m, 0], coords[tfh_m, 1], s=S_TFH, c=TFH_COLOR, rasterized=True, linewidths=0)
    ax.scatter(coords[fdc_m, 0], coords[fdc_m, 1], s=S_FDC, c=FDC_COLOR, rasterized=True, linewidths=0)
    ax.set_xlim(crop[0], crop[1])
    ax.set_ylim(crop[3], crop[2])   # invert y to match Xenium image orientation
    ax.set_aspect('equal')
    ax.set_title(label_for(sample), fontweight='bold', fontsize=14)
    ax.set_xticks([])
    ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)

for j in range(len(samples), len(axes_flat)):
    axes_flat[j].axis('off')

legend_handles = [
    Line2D([0], [0], marker='o', linestyle='', color=FDC_COLOR, markersize=10, label='FDC'),
    Line2D([0], [0], marker='o', linestyle='', color=TFH_COLOR, markersize=8,  label='Tfh'),
    Line2D([0], [0], marker='o', linestyle='', color=BG_COLOR,  markersize=6,  label='Other cells'),
]
fig.legend(handles=legend_handles, loc='lower center', ncol=3,
           bbox_to_anchor=(0.5, -0.01), frameon=False, fontsize=12)
plt.suptitle('FDC vs Tfh spatial distribution', y=1.1, fontweight='bold', fontsize=16)
plt.subplots_adjust(wspace=0.05, hspace=0.15)


# --- Plot 2: bar plots of FDC count and pct reannotated per sample

counts_plot = counts_df.copy()
n = len(counts_plot)
bar_w = 0.6
# Width: ~0.6" per bar + side padding — keeps the axes narrow as the cohort is small.
fig_w = 0.7 * n + 1.5
fig, axes = plt.subplots(1, 2, figsize=(fig_w * 2, 3.5))
xs = np.arange(n)

ax = axes[0]
ax.bar(xs, counts_plot['n_tfh'], width=bar_w, color=TFH_COLOR, label='Tfh', linewidth=0)
ax.bar(xs, counts_plot['n_fdc'], width=bar_w, bottom=counts_plot['n_tfh'], color=FDC_COLOR, label='FDC', linewidth=0)
ax.set_ylabel('Cells originally annotated as Tfh')
ax.set_title('Reannotation counts', fontweight='bold')
ax.set_xticks(xs)
ax.set_xticklabels(counts_plot['sample'], rotation=45, ha='right')
ax.set_xlim(-0.5, n - 0.5)
ax.grid(False)
ax.legend(frameon=False)
for spine in ('top', 'right'):
    ax.spines[spine].set_visible(False)

ax = axes[1]
ax.bar(xs, counts_plot['pct_reannotated_FDC'], width=bar_w, color=FDC_COLOR, linewidth=0)
ax.set_ylabel('% of original Tfh reannotated as FDC')
ax.set_title('Reannotation rate', fontweight='bold')
ax.set_xticks(xs)
ax.set_xticklabels(counts_plot['sample'], rotation=45, ha='right')
ax.set_xlim(-0.5, n - 0.5)
ax.grid(False)
for spine in ('top', 'right'):
    ax.spines[spine].set_visible(False)

plt.suptitle('FDC reannotation per sample', y=1.1, fontweight='bold')
plt.tight_layout()


# --- Plot 3/4/5: pooled FDC vs Tfh gene expression
# Pool only FDC + Tfh cells across all samples into a single AnnData on the marker subset.

pooled_obs = []
pooled_x = []
for sample in samples:
    adata = adatas[sample]
    mask = adata.obs['Tfh_group'].isin(['FDC', 'Tfh']).values
    if not mask.any():
        continue
    sub = adata[mask, all_markers]
    pooled_obs.append(pd.DataFrame({
        'sample': sample,
        'Tfh_group': adata.obs.loc[mask, 'Tfh_group'].astype(str).values,
    }, index=adata.obs_names[mask]))
    pooled_x.append(sub.X.toarray() if hasattr(sub.X, 'toarray') else np.asarray(sub.X))

pooled_obs = pd.concat(pooled_obs)
pooled_x   = np.vstack(pooled_x)
pool = ad.AnnData(X=pooled_x, obs=pooled_obs, var=pd.DataFrame(index=all_markers))

# Only include groups that are actually present (samples may lack a recovered FDC subcluster).
present_groups = [g for g in ['FDC', 'Tfh'] if (pool.obs['Tfh_group'] == g).any()]
pool.obs['Tfh_group'] = pd.Categorical(pool.obs['Tfh_group'], categories=present_groups)
print(f'\nPooled: {pool.n_obs} cells (FDC={int((pool.obs.Tfh_group=="FDC").sum())}, '
      f'Tfh={int((pool.obs.Tfh_group=="Tfh").sum())})')

# Plot 3: stacked violin (cohort)
sc.pl.stacked_violin(pool, var_names=all_markers, groupby='Tfh_group', show=False)
plt.suptitle('Cohort: FDC vs Tfh (stacked violins)', y=1.1, fontweight='bold')

# Plot 4: standard dotplot (cohort)
sc.pl.dotplot(pool, var_names=all_markers, groupby='Tfh_group',
              standard_scale='var', show=False)
plt.suptitle('Cohort: FDC vs Tfh (standard dotplot)', y=1.1, fontweight='bold')


# Plot 5 (Extended Data Fig. 9j): per-gene-scaled dotplot — gene groups split into horizontally-titled
# sub-axes, horizontal n-cells bars, and a vertical cbar + size legend stacked to the right.
gene_categories_dict = {name: [g for g in gs if g in pool.var_names]
                         for name, gs in GENE_GROUPS.items()}
gene_categories_dict = {k: v for k, v in gene_categories_dict.items() if v}
genes_flat = [g for gs in gene_categories_dict.values() for g in gs]

count_palette = {'FDC': FDC_COLOR, 'Tfh': TFH_COLOR}

# Layout constants (X bumped up slightly so the gene panel reads comfortably
# with only 2 categories on Y).
SPACE_PER_DOT_Y = 0.055
SPACE_PER_DOT_X = 0.028
LABEL_FS        = 9
DOT_MAX         = 150

fig, ax = plt.subplots(figsize=(7, 4))
pgn = plotting_functions.per_gene_norm_dotplot_with_legends(
    pool, var_names=genes_flat, groupby='Tfh_group', ax=ax,
    cmap='Reds', vmin=0, vmax=1, standard_scale='var_max',
    max_dot_size=DOT_MAX, space_per_dot=SPACE_PER_DOT_Y,
    width_ratio=SPACE_PER_DOT_X / SPACE_PER_DOT_Y,
    count_palette=count_palette, count_title='',
    label_fontsize=LABEL_FS - 1,
    size_normalize='max', color_mode='standard_scale',
    transpose=True, gene_categories=gene_categories_dict,
    category_title_rotation=0,
)
for sub_ax in pgn['dotplot_axes']:
    for t in sub_ax.get_xticklabels():
        t.set_fontsize(LABEL_FS - 1)
        t.set_rotation(90)
    # CXCL13 title runs into its neighbours' top tick labels — drop it.
    if sub_ax.get_title() == 'CXCL13':
        sub_ax.set_title('')
for t in pgn['dotplot_ax'].get_yticklabels():
    t.set_fontsize(LABEL_FS)

# --- Drop helper aux axes; re-render the count bars, cbar, and size legend directly.
fig.delaxes(pgn['count_ax'])
fig.delaxes(pgn['cbar_ax'])
fig.delaxes(pgn['size_ax'])

dotplot_bbox = pgn['dotplot_bbox']

# n cells bar (horizontal barh)
EDGE_BUFFER = 0.7
def _fmt_count(nval):
    nval = int(nval)
    return f'{nval/1000:.1f}k' if nval >= 1000 else f'{nval}'

group_order = present_groups
counts_bars = pool.obs['Tfh_group'].value_counts().reindex(group_order)
count_max = int(counts_bars.max())
xmax_global = count_max * 1.25

bar_pad = 0.012
bar_w   = 3 * SPACE_PER_DOT_X
bar_ax = fig.add_axes([dotplot_bbox.x1 + bar_pad, dotplot_bbox.y0,
                       bar_w, dotplot_bbox.height])
bar_ys = list(range(len(group_order)))
bar_colors = [count_palette[g] for g in group_order]
bar_ax.barh(bar_ys, counts_bars.values, color=bar_colors, edgecolor='black',
            height=0.7, linewidth=0.5)
bar_ax.set_xlim(0, xmax_global)
bar_ax.set_ylim(len(group_order) - 1 + EDGE_BUFFER, -EDGE_BUFFER)
bar_ax.set_yticks([])
for v, yp in zip(counts_bars.values, bar_ys):
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

# Vertical cbar to the right of the bar ax, spanning full dotplot height
fig.canvas.draw()
inv = fig.transFigure.inverted()
bar_right = inv.transform(bar_ax.get_tightbbox(fig.canvas.get_renderer()))[1][0]

cbar_y0 = dotplot_bbox.y0
cbar_h  = dotplot_bbox.y1 - dotplot_bbox.y0

cbar_pad = 0.05
cbar_w   = 0.012
cax = fig.add_axes([bar_right + cbar_pad, cbar_y0, cbar_w, cbar_h])
cbar = mpl.colorbar.ColorbarBase(
    cax, cmap=plt.get_cmap('Reds'),
    norm=mpl.colors.Normalize(vmin=0, vmax=1),
    orientation='vertical',
)
cbar.set_label('Mean expression\n(scaled)', fontsize=LABEL_FS, labelpad=3)
cax.yaxis.set_label_position('left')
cax.tick_params(labelsize=LABEL_FS - 1)
cax.set_yticks([0, 0.5, 1.0])
cax.grid(False)

# Vertical size legend right of cbar — 3 markers since the FDC plot is short.
cbar_pos = cax.get_position()
size_h  = cbar_h * 0.85
size_y  = cbar_y0 + (cbar_h - size_h) / 2
size_x  = cbar_pos.x1 + 0.05
size_w  = 0.07
size_ax = fig.add_axes([size_x, size_y, size_w, size_h])
size_ax.set_xlim(0, 1); size_ax.set_ylim(0, 1)
size_ax.axis('off')
size_ax.set_title('% of\nper-gene max', fontsize=LABEL_FS,
                  ha='left', x=0, pad=4)

frac_examples = [0.25, 0.5, 1.0]
y_positions = np.linspace(0.88, 0.12, len(frac_examples))
for frac, y in zip(frac_examples, y_positions):
    size_ax.scatter(0.15, y, s=frac * DOT_MAX, c='lightgrey',
                    edgecolors='black', linewidths=0.5,
                    transform=size_ax.transAxes, clip_on=False)
    size_ax.text(0.45, y, f'{int(frac * 100)}%', fontsize=LABEL_FS,
                 va='center', ha='left', transform=size_ax.transAxes)

# Suptitle centred on the dot grid; fontsize matches the sub-ax category
# titles (label_fontsize + 2 inside the helper = LABEL_FS + 1).
title_x = (dotplot_bbox.x0 + dotplot_bbox.x1) / 2
plt.suptitle('Marker expression\nFDC vs Tfh', fontweight='bold',
             fontsize=LABEL_FS + 1,
             x=title_x, y=pgn['suptitle_y'] - 0.03, va='bottom')


# --- Save outputs

pdf_path = os.path.join(plot_dir, 'FDC_cohort_summary.pdf')
save_multi_image(pdf_path)
plt.close('all')

print(f'\nSaved: {pdf_path}')
print(f'Saved: {os.path.join(csv_dir, "FDC_reannotation_counts.csv")}')
print("Script successfully completed")
