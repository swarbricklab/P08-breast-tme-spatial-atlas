#!/usr/bin/env python3

# --- Aim
# Plot mean CD40 and CD40LG expression across the cell types that comprise at least 1% of
# cells within the identified Xenium TLS regions, min-max scaled 0-1 per gene (S7D).

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

import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.patches as mpatches

mpl.rcParams['pdf.fonttype'] = 42


# --- Argparse arguments

parser = argparse.ArgumentParser(description="CD40/CD40LG expression in Xenium TLS regions (S7D)",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument("-s", "--samples",
    nargs="*",
    type=str,
    help="List of Xenium sample IDs (section_id)")
parser.add_argument("-r", "--res_dir", help="Path to results directory")
parser.add_argument("--h5ad_dir", help="Path to directory with .h5ad objects")
parser.add_argument("--project_dir", help="Project directory (root)")
parser.add_argument("--sample_sheet", help="Path to sample sheet CSV file", default=None)
parser.add_argument("--sc_ref", help="Single-cell reference used for cell-type classifications", default='miniatlas')
parser.add_argument("--annotation_levels", nargs="*", default=['subset'],
                    help="Annotation levels to plot (default: subset only)")
parser.add_argument("--genes", nargs="*", default=['CD40LG', 'CD40'],
                    help="Genes to display on the heatmap")
parser.add_argument("--min_cells_pct", type=float, default=1.0,
                    help="Minimum percentage of in-TLS cells a cell type must comprise to be plotted")

args = vars(parser.parse_args())

project_dir = args["project_dir"]
h5ad_dir = args["h5ad_dir"]
res_dir = args["res_dir"]
sample_list = args["samples"]
sc_ref = args["sc_ref"]
levels = args["annotation_levels"]
genes = args["genes"]
min_cells_pct = args["min_cells_pct"]


# --- Local imports

sys.path.append(os.path.join(project_dir, 'scripts'))
from helper_functions import save_multi_image


# --- Load the cell-type classifications

classif_dir = os.path.join(project_dir, 'resources', 'cell_types', sc_ref)
classif_file = f'{sc_ref}_subset_classifications.csv'
classif_df = pd.read_csv(os.path.join(classif_dir, classif_file))


# --- Load plotting colours

color_dir = os.path.join(project_dir, 'resources', 'plotting_colors')
color_dict = pd.read_csv(os.path.join(color_dir, f'{sc_ref}_major_colors.csv'),
                         index_col=0).to_dict()['color']

for level in levels:
    for cell_type in classif_df[f'{sc_ref}_{level}'].unique():
        if pd.isna(cell_type) or cell_type in color_dict:
            continue
        major_type = classif_df[classif_df[f'{sc_ref}_{level}'] == cell_type][f'{sc_ref}_major'].values
        major_type_in_dict = [m for m in major_type if m in color_dict]
        color_dict[cell_type] = color_dict[major_type_in_dict[0]] if major_type_in_dict else 'lightgrey'

sub_short_dict = classif_df.set_index(f'{sc_ref}_subset')[[f'{sc_ref}_subset_short']].to_dict()


# --- Cell-type ordering (lineage, then major)

lineage_order = ['Epithelial', 'Immune', 'Stromal']

major_cell_types = classif_df[f'{sc_ref}_major'].dropna().unique().tolist()


def get_lineage(cell_type, level='major'):
    if cell_type in classif_df[f'{sc_ref}_{level}'].values:
        return classif_df[classif_df[f'{sc_ref}_{level}'] == cell_type]['lineage'].values[0]
    return 'other'


def get_major(cell_type, level, upper_level='major'):
    if cell_type in classif_df[f'{sc_ref}_{level}'].values:
        return classif_df[classif_df[f'{sc_ref}_{level}'] == cell_type][f'{sc_ref}_{upper_level}'].values[0]
    return 'other'


major_cell_types_sorted = sorted(
    major_cell_types,
    key=lambda x: (lineage_order.index(get_lineage(x)) if get_lineage(x) in lineage_order else len(lineage_order), x))


# --- Load anndata objects, cell-type labels and TLS annotations

def load_label_transfer(sample, level):
    annot_path = os.path.join(res_dir, '02_cell_type_mapping', 'label_transfer', 'xenium',
                              f'{level}_level', f'{sample}_{level}_pred_cell_types.csv')
    if not os.path.exists(annot_path):
        annot_path = os.path.join(res_dir, '02_cell_type_mapping_regenerated', 'label_transfer',
                                  'xenium', f'{level}_level', f'{sample}_{level}_pred_cell_types.csv')
    return pd.read_csv(annot_path, index_col=0)['pred_cell_type']

tls_csv_dir = os.path.join(res_dir, '05_xenium_mechanisms', '01_identify_TLS_regions', 'csv_files')

adata_dict = {}

for sample in sample_list:
    adata = ad.read_h5ad(os.path.join(h5ad_dir, 'xenium', f'{sample}.h5ad'))

    tls_df = pd.read_csv(os.path.join(tls_csv_dir, f'{sample}_cell_TLS_annotations.csv'), index_col=0)
    in_tls = tls_df['in_TLS'].astype(str).map({'True': True, 'False': False})
    adata.obs['in_TLS'] = in_tls.reindex(adata.obs_names).fillna(False).astype(bool)

    for level in levels:
        annotations = load_label_transfer(sample, level)
        adata.obs[level] = annotations.reindex(adata.obs_names)

    # Keep only cells located inside TLS regions
    adata_dict[sample] = adata[adata.obs['in_TLS']].copy()

all_adata = ad.concat(list(adata_dict.values()), merge='same')
all_adata.obs_names_make_unique()


# --- Helper: resize a matrixplot axis based on the number of genes / cell types

def auto_resize_heatmap(ax, space_per_dot=0.04, width_ratio=0.7):
    pos = ax.get_position()
    new_width = len(ax.get_xticks()) * space_per_dot * width_ratio
    new_height = len(ax.get_yticks()) * space_per_dot
    ax.set_position([pos.x0, pos.y0, new_width, new_height])


# --- Plot (S7D): CD40 / CD40LG heatmap, min-max scaled per gene

# Restrict to genes present in the panel
genes_present = [g for g in genes if g in all_adata.var_names]
missing = [g for g in genes if g not in all_adata.var_names]
if missing:
    print(f'Warning: genes missing from the panel and skipped: {missing}')

vmin, vmax = 0, 1

for level in levels:

    # Restrict to cell types comprising at least min_cells_pct% of the in-TLS cells
    counts = all_adata.obs[level].value_counts()
    min_cells = min_cells_pct / 100 * all_adata.n_obs
    keep_types = counts[counts >= min_cells].index.tolist()
    adata_sub = all_adata[all_adata.obs[level].isin(keep_types)].copy()

    # Order cell types by descending number of cells (drop any unused categories)
    level_counts = adata_sub.obs[level].value_counts()
    cells_sorted = level_counts[level_counts > 0].index.tolist()
    adata_sub.obs[level] = pd.Categorical(adata_sub.obs[level], categories=cells_sorted, ordered=True)

    print(f'Annotation level {level}: plotting {adata_sub.n_obs} cells across '
          f'{len(cells_sorted)} cell types')

    fig, ax = plt.subplots(figsize=(5, 5))
    dp = sc.pl.matrixplot(
        adata_sub,
        var_names=genes_present,
        groupby=level,
        ax=ax,
        cmap='viridis',
        vmin=vmin,
        vmax=vmax,
        show=False,
        swap_axes=True,
        standard_scale='var',
    )
    heatmap_ax = dp['mainplot_ax']
    cbar_ax = dp['color_legend_ax']
    auto_resize_heatmap(heatmap_ax, space_per_dot=0.04, width_ratio=0.7)
    heatmap_ax.set_ylim(len(heatmap_ax.get_yticks()), 0)

    # Position the colourbar to the right of the heatmap
    ax_pos = heatmap_ax.get_position()
    cbar_ax.set_position([ax_pos.x1 + 0.03, ax_pos.y0, 0.04 * 3.5, 0.04 * 0.7])
    cbar_ax.set_title(cbar_ax.get_title().replace('Mean', 'Scaled'), fontsize=8)

    # Cell-count barplot on top, coloured by major class
    cell_counts = adata_sub.obs[level].value_counts().reindex(cells_sorted)
    color_list = [color_dict.get(get_major(cell, level), 'lightgrey') for cell in cells_sorted]

    bar_ax = fig.add_axes([ax_pos.x0, ax_pos.y1 + 0.02, ax_pos.width,
                           ax_pos.height / max(len(heatmap_ax.get_yticks()), 1) * 3])
    x = np.arange(len(cells_sorted))
    bar_ax.bar(x, cell_counts.values, color=color_list, linewidth=0.5, edgecolor='black')
    bar_ax.set_xlim(heatmap_ax.get_xlim()[0] - 0.5, heatmap_ax.get_xlim()[1] - 0.5)
    bar_ax.set_xticks([])
    bar_ax.set_yticks([])
    for spine in ['top', 'right', 'left']:
        bar_ax.spines[spine].set_visible(False)
    bar_ax.set_ylabel('# cells\n(x10³)', fontsize=8)
    for rect, count in zip(bar_ax.patches, cell_counts.values):
        height = rect.get_height()
        if height > 0:
            label = f'{count / 1000:.1f}' if count >= 1000 else f'{count}'
            bar_ax.text(rect.get_x() + rect.get_width() / 2, height + bar_ax.get_ylim()[1] * 0.05,
                        label, ha='center', va='bottom', fontsize=7, rotation=90)
    bar_ax.set_ylim(0, bar_ax.get_ylim()[1] * 1.3)
    bar_ax.set_title('Gene expression in TLS regions', fontweight='bold', pad=10)

    # Replace subset-level tick labels with short names
    if level == 'subset':
        new_labels = [sub_short_dict[f'{sc_ref}_subset_short'].get(label.get_text(), label.get_text())
                      for label in heatmap_ax.get_xticklabels()]
        heatmap_ax.set_xticklabels(new_labels)

    # Colour legend for the major classes present
    major_present = [m for m in major_cell_types_sorted
                     if m in {get_major(c, level) for c in cells_sorted}]
    handles = [mpatches.Patch(facecolor=color_dict[m], label=m, edgecolor='black', linewidth=0.8)
               for m in major_present if m in color_dict]
    # Place the legend above the colourbar (left edges aligned), growing upward, so it
    # never overlaps the colourbar or its title.
    fig.canvas.draw()
    cbar_bbox = fig.transFigure.inverted().transform(
        cbar_ax.get_tightbbox(fig.canvas.get_renderer()))
    cbar_left, cbar_top = cbar_bbox[0][0], cbar_bbox[1][1]
    leg = fig.legend(handles=handles, loc='lower left',
                     bbox_to_anchor=(cbar_left, cbar_top + 0.03),
                     bbox_transform=fig.transFigure,
                     frameon=False, handlelength=0.8, fontsize=8)

    # Place the suptitle above the topmost element on the figure (heatmap, cell-count
    # bar, colourbar or legend), centred over the heatmap. va='bottom' anchors the text
    # bottom at that edge so the title grows upward and never overlaps the plot or legend.
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    inv = fig.transFigure.inverted()
    tops = [inv.transform(a.get_tightbbox(renderer))[1][1]
            for a in (heatmap_ax, bar_ax, cbar_ax)]
    tops.append(inv.transform(leg.get_window_extent(renderer))[1][1])
    content_top = max(tops)
    bar_ax_pos = bar_ax.get_position()
    plt.suptitle(f'All samples - Gene expression in TLS regions '
                 f'(min_cells={min_cells_pct}%, annotation_key={level})',
                 x=(bar_ax_pos.x0 + bar_ax_pos.x1) / 2, y=content_top + 0.02,
                 va='bottom', fontweight='bold')


# --- Save the output

out_dir = os.path.join(res_dir, '05_xenium_mechanisms', '01_identify_TLS_regions',
                       'all_samples', 'gene_expression')
plots_dir = os.path.join(out_dir, 'plots')
os.makedirs(plots_dir, exist_ok=True)

genes_str = '_'.join(genes_present)
save_multi_image(os.path.join(plots_dir, f'All_samples_TLS_gene_expression_{genes_str}.pdf'))


# --- Finish

print("Script successfully completed")
