#!/usr/bin/env python3

# --- Aim
# Plot Figure 5C: cell type-specific gene expression within TLS
# Shows expression patterns of key chemokines and receptors in different cell populations inside TLS structures

# --- Define environment
# env: xen-hev

# --- Load packages

import argparse
import os
import anndata as ad
import scanpy as sc

import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.patches as mpatches
from matplotlib.colors import Normalize

mpl.rcParams['pdf.fonttype'] = 42


# --- Argparse arguments

parser = argparse.ArgumentParser(description="Import Visium unfiltered data & plot",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument("-s", "--samples",
    nargs="*",  # 0 or more values expected => creates a list
    type=str,
    help="List of samples")
parser.add_argument("-r", "--res_dir", help="Path to results directory")
parser.add_argument("--h5ad_dir", help="Path to directory for .h5ad output files")
parser.add_argument("--project_dir", help="Project directory (root)")
parser.add_argument("--sample_sheet", help="Path to sample sheet CSV file", default=None)

args = vars(parser.parse_args())

# Arguments
project_dir = args["project_dir"]
h5ad_dir = args["h5ad_dir"]
res_dir = args["res_dir"]
sample_list = args["samples"]


# --- Load the anndata objects

adata_dict = {}

for sample in sample_list:
    adata = ad.read_h5ad(os.path.join(h5ad_dir, 'xenium', f'{sample}.h5ad'))
    adata_dict[sample] = adata


# --- Load the tangram annotations

# Decide which annotation levels to use
level = 'subset'

for sample in sample_list:
    annotation_file = f'{sample}_{level}_pred_cell_types.csv'
    annotation_dir = os.path.join(res_dir, '02_cell_type_mapping', 'label_transfer', 'xenium')
    annotations = pd.read_csv(os.path.join(annotation_dir, f'{level}_level', annotation_file), index_col=0)
    # Add to adata
    adata_dict[sample].obs[level] = annotations['pred_cell_type']


# --- Load plotting colors

color_dir = os.path.join(project_dir, 'resources', 'plotting_colors')

# Use the major level for colors
color_level = 'major'
color_file = f'miniatlas_{color_level}_colors.csv'

color_dict = pd.read_csv(os.path.join(color_dir, color_file), index_col=0).to_dict()['color']


# --- Load the miniatlas classifications

# These are needed to check the alignment between major and minor levels
miniatlas_dir = os.path.join(project_dir, 'resources', 'cell_types', 'miniatlas')
miniatlas_file = 'miniatlas_subset_classifications.csv'

miniatlas_df = pd.read_csv(os.path.join(miniatlas_dir, miniatlas_file),)

# Add the other levels to the color_dict
for cell_type in miniatlas_df[f'miniatlas_{level}'].unique():
    if cell_type not in color_dict.keys():
        # Check that the value is not nan
        if pd.isna(cell_type):
            continue
        print(f'Adding color for {cell_type}')
    # Get the corresponding major cell type
    major_type = miniatlas_df[miniatlas_df[f'miniatlas_{level}'] == cell_type][f'miniatlas_major'].values
    # Find which value is in the color_dict
    major_type_in_dict = [m for m in major_type if m in color_dict.keys()]
    if len(major_type_in_dict) > 0:
        color_dict[cell_type] = color_dict[major_type_in_dict[0]]
    else:
        color_dict[cell_type] = 'lightgrey'

# Order miniatlas_major by lineage
lineage_order = ['Epithelial', 'Immune', 'Stromal']

major_cell_types = miniatlas_df[f'miniatlas_major'].dropna().unique().tolist()

# Sort based on lineage first (recorded in the miniatlas df)
def get_lineage(cell_type, level='major'):
    if cell_type in miniatlas_df[f'miniatlas_{level}'].values:
        lineage = miniatlas_df[miniatlas_df[f'miniatlas_{level}'] == cell_type]['lineage'].values[0]
    else:
        lineage = 'other'
    return lineage

major_cell_types_sorted = sorted(major_cell_types, key=lambda x: (lineage_order.index(get_lineage(x)) if get_lineage(x) in lineage_order else len(lineage_order), x))

# Also define a dict with shorter names for the subset level cell types
miniatlas_sub_dict = miniatlas_df.set_index('miniatlas_subset')[['miniatlas_subset_short']].to_dict()


# --- Load the TLS annotations

tls_dir = os.path.join(res_dir, '05_xenium_mechanisms', '01_identify_TLS_regions', 'csv_files')

tls_df = pd.DataFrame()

for sample in sample_list:
    tls_file = f'{sample}_cell_TLS_annotations.csv'
    tmp_df = pd.read_csv(os.path.join(tls_dir, tls_file), index_col=0)
    # Add to adata
    adata_dict[sample].obs['in_TLS'] = tmp_df['in_TLS']

# Subset all adata objects to only the cells in the TLS annotation df
for sample in sample_list:
    adata_sub = adata_dict[sample][adata_dict[sample].obs['in_TLS']]
    adata_dict[sample] = adata_sub

# Extract a single adata object with all samples
adata_list = list(adata_dict.values())
all_adata = ad.concat(adata_list, merge='same')

all_adata.obs_names_make_unique()


# --- Some functions for plotting

# Automatic resizing function
def auto_resize_dotplot(ax, space_per_dot=0.1, width_ratio=1):
    """
    Automatically resize the dotplot axes based on the number of genes.
    """
    # Get the current position of the axes
    pos = ax.get_position()
    # Resize based on xticks & yticks
    xticks = len(ax.get_xticks())
    yticks = len(ax.get_yticks())
    # Calculate new width based on number of genes
    new_width = xticks * space_per_dot * width_ratio
    new_height = yticks * space_per_dot
    # Set the new position
    ax.set_position([pos.x0, pos.y0, new_width, new_height])

# Function for repositioning the dotplot axes
def reposition_dotplot_axs(axs, y_offset=0.1):
    for n, ax in enumerate(axs):
        # Get the current position of the axes
        pos = ax.get_position()
        if n == 0:
            ypos = pos.y0  # Keep the first plot at the original position
        else:
            ypos = ax_pos_prev.y0 - pos.height - y_offset  # Offset from the previous plot
        # Set the new position with an offset
        ax.set_position([pos.x0, ypos, pos.width, pos.height])
        ax_pos_prev = ax.get_position()  # Store the position for the next iteration
        # Remove the xticks for all but the last plot
        if n < len(axs) - 1:
            ax.set_xticks([])

# Add legends
def position_legend(ax):
    """
    Given a single axis, return the bottom, top, left, and right edges
    in figure coordinates, for robust legend placement.
    Returns (bottom_edge, top_edge, left_edge, right_edge)
    """
    fig = ax.figure
    fig.canvas.draw()
    bbox = ax.get_tightbbox(fig.canvas.get_renderer())
    inv_trans = fig.transFigure.inverted()
    bbox_fig = inv_trans.transform(bbox)
    # bbox_fig: [[x0, y0], [x1, y1]]
    left_edge = bbox_fig[0][0]
    right_edge = bbox_fig[1][0]
    bottom_edge = bbox_fig[0][1]
    top_edge = bbox_fig[1][1]
    return bottom_edge, top_edge, left_edge, right_edge

# Custom dotplot function with per-gene size normalization
def custom_dotplot_per_gene_norm(adata, var_names, groupby, ax, cmap='Reds', vmin=0, vmax=1, 
                                  standard_scale='var', max_dot_size=200, edgecolor='black', linewidth=0.5):
    """
    Create a dotplot with per-gene size normalization (0% to max% for each gene).
    Colors show standardized expression (0-1), sizes show fraction of cells expressing.
    """
    # Get categories
    categories = adata.obs[groupby].cat.categories.tolist()
    
    # Calculate mean expression and fraction for each gene-category combination
    gene_data = []
    for gene in var_names:
        for cat in categories:
            mask = adata.obs[groupby] == cat
            gene_values = adata[mask, gene].X.toarray().flatten() if hasattr(adata.X, 'toarray') else adata[mask, gene].X.flatten()
            
            # Calculate fraction of cells expressing (>0)
            frac = (gene_values > 0).sum() / len(gene_values) if len(gene_values) > 0 else 0
            
            # Calculate mean expression
            mean_expr = gene_values.mean() if len(gene_values) > 0 else 0
            
            gene_data.append({
                'gene': gene,
                'category': cat,
                'mean_expr': mean_expr,
                'fraction': frac
            })
    
    df = pd.DataFrame(gene_data)
    
    # Standardize expression per gene if requested
    if standard_scale == 'var':
        for gene in var_names:
            gene_mask = df['gene'] == gene
            gene_vals = df.loc[gene_mask, 'mean_expr']
            if gene_vals.max() > gene_vals.min():
                df.loc[gene_mask, 'mean_expr_scaled'] = (gene_vals - gene_vals.min()) / (gene_vals.max() - gene_vals.min())
            else:
                df.loc[gene_mask, 'mean_expr_scaled'] = 0
    else:
        df['mean_expr_scaled'] = df['mean_expr']
    
    # Clip to vmin, vmax
    df['mean_expr_scaled'] = df['mean_expr_scaled'].clip(vmin, vmax)
    
    # Normalize sizes per gene (min to max_fraction for each gene)
    for gene in var_names:
        gene_mask = df['gene'] == gene
        gene_fracs = df.loc[gene_mask, 'fraction']
        max_frac = gene_fracs.max()
        min_frac = gene_fracs.min()
        if max_frac > min_frac:
            # Scale from min to max, then map to [0.001, 1.0] so min is barely visible
            normalized = (gene_fracs - min_frac) / (max_frac - min_frac)
            df.loc[gene_mask, 'size_norm'] = 0.001 + normalized * 1
        elif max_frac > 0:
            # All fractions are the same (and non-zero)
            df.loc[gene_mask, 'size_norm'] = 0.6
        else:
            df.loc[gene_mask, 'size_norm'] = 0
    
    # Create color normalization
    norm = Normalize(vmin=vmin, vmax=vmax)
    cmap_obj = plt.get_cmap(cmap)
    
    # Plot
    for idx, row in df.iterrows():
        x_pos = categories.index(row['category'])
        y_pos = var_names.index(row['gene'])
        
        color = cmap_obj(norm(row['mean_expr_scaled']))
        size = row['size_norm'] * max_dot_size
        
        ax.scatter(x_pos, y_pos, s=size, c=[color], edgecolor=edgecolor, linewidth=linewidth, zorder=10)
    
    # Set ticks and labels
    ax.set_xticks(range(len(categories)))
    ax.set_xticklabels(categories, rotation=90, ha='center', fontsize=8)
    ax.set_yticks(range(len(var_names)))
    ax.set_yticklabels(var_names, fontsize=8)
    
    # Set limits with proper buffer
    ax.set_xlim(-0.5, len(categories) - 0.5)
    ax.set_ylim(len(var_names) - 0.5, -0.5)  # Inverted: max on top, min on bottom
    
    # No grid
    ax.grid(False)
    
    # Return max fraction per gene for legend
    max_fracs_per_gene = df.groupby('gene')['fraction'].max()
    overall_max_frac = max_fracs_per_gene.max()
    
    return ax, overall_max_frac

# Function to add legend
def add_color_legend(fig, labels, color_dict, title='', x=0, y=0, offset=0.01, spacing_per_label=0.05, loc='center left', fontsize=10):
    # Create legend handles (colored boxes with labels)
    handles = [
        mpatches.Patch(facecolor=color_dict[label], label=label, edgecolor='black', linewidth=0.8)
        for label in labels if label in color_dict
    ]
    # Estimate height needed
    est_height = spacing_per_label * len(handles)
    bbox_anchor = (x + offset, y)
    # Add the legend
    legend = fig.legend(
        handles=handles,
        title=title,
        loc=loc,
        bbox_to_anchor=bbox_anchor,
        borderaxespad=0,
        frameon=False,
        handlelength=0.8,
        bbox_transform=fig.transFigure,  # <-- This is the key fix!
        fontsize=fontsize,
    )
    legend.get_title().set_fontweight('bold')
    legend.get_title().set_fontsize(fontsize)  # Set the title fontsize
    legend._legend_box.align = "left"
    # Return updated vertical position for next stacked legend
    return legend

def join_lr_multi(axs, lr_pairs, offset=0.3, offset_start=0.2):
    """
    Draws brackets connecting ligands and receptors, extending the bracket
    further left if a ligand has multiple receptors or a receptor has
    multiple ligands.

    This version correctly handles labels that appear on multiple axes,
    counting and offsetting them independently per axis.

    Args:
        axs (list | Axes): The matplotlib Axes object(s).
        lr_pairs (list[tuple]): A list of (ligand, receptor) string pairs.
        offset (float): Horizontal distance from the label's left edge to the end of the bracket.
        offset_start (float): Horizontal distance from the longest label to the main vertical line.
    """
    if not isinstance(axs, (list, tuple, np.ndarray)):
        axs = [axs]
    # 1. Build a lookup that can handle labels on MULTIPLE axes.
    # label_to_info['L1'] -> [{'ax': ax1, 'idx': 2, ...}, {'ax': ax2, 'idx': 4, ...}]
    label_to_info = {}
    all_labels = set(lr for pair in lr_pairs for lr in pair)
    if not all_labels:
        return 0
    renderer = axs[0].get_figure().canvas.get_renderer()
    for ax in axs:
        inv = ax.transData.inverted()
        yticks = ax.get_yticks()
        for idx, label in enumerate(ax.get_yticklabels()):
            text = label.get_text()
            if text in all_labels:
                if text not in label_to_info:
                    label_to_info[text] = []
                # Use actual tick position instead of enumeration index
                y_pos = yticks[idx] if idx < len(yticks) else idx
                info = {
                    'ax': ax,
                    'idx': y_pos,
                    'bbox': label.get_window_extent(renderer),
                    'inv': inv
                }
                label_to_info[text].append(info)
    # 2. Resolve pairs by finding common axes.
    # A pair (L1, R1) is only valid if L1 and R1 are found on the same axis.
    resolved_pairs = []
    for lr1, lr2 in lr_pairs:
        info_list1 = label_to_info.get(lr1, [])
        info_list2 = label_to_info.get(lr2, [])
        # Find all valid (L1, R1) instances on a common axis
        for info1 in info_list1:
            for info2 in info_list2:
                if info1['ax'] == info2['ax']:
                    resolved_pairs.append({
                        'lr1': lr1,
                        'lr2': lr2,
                        'ax': info1['ax'],
                        'info1': info1, # All info for lr1 (ax, idx, bbox, inv)
                        'info2': info2  # All info for lr2 (ax, idx, bbox, inv)
                    })
    if not resolved_pairs:
        return 0 # No valid pairs found
    # 3. Count occurrences PER AXIS, using only resolved pairs.
    ligand_counts_per_ax = {ax: {} for ax in axs}
    receptor_counts_per_ax = {ax: {} for ax in axs}
    for pair in resolved_pairs:
        ax, lr1, lr2 = pair['ax'], pair['lr1'], pair['lr2']
        ligand_counts_per_ax[ax].setdefault(lr1, 0)
        receptor_counts_per_ax[ax].setdefault(lr2, 0)
        ligand_counts_per_ax[ax][lr1] += 1
        receptor_counts_per_ax[ax][lr2] += 1
    # 4. Determine the main vertical alignment line (Global)
    all_x0_coords = []
    for pair in resolved_pairs:
        all_x0_coords.append(pair['info1']['inv'].transform((pair['info1']['bbox'].x0, 0))[0])
        all_x0_coords.append(pair['info2']['inv'].transform((pair['info2']['bbox'].x0, 0))[0])
    if not all_x0_coords:
        return 0
    min_x0 = min(all_x0_coords)
    x_left_align = min_x0 - offset_start - offset
    # 5. Prepare for drawing with PER-AXIS seen counters.
    ligand_seen_per_ax = {ax: {lig: 0 for lig in ax_counts} for ax, ax_counts in ligand_counts_per_ax.items()}
    receptor_seen_per_ax = {ax: {rec: 0 for rec in ax_counts} for ax, ax_counts in receptor_counts_per_ax.items()}
    drawn_verticals_per_ax = {ax: [] for ax in axs}
    indent_step = offset * 2
    # 6. Main drawing loop (using resolved pairs).
    for pair in resolved_pairs:
        # Get all resolved data for this specific pair instance
        ax = pair['ax']
        lr1, lr2 = pair['lr1'], pair['lr2']
        info1, info2 = pair['info1'], pair['info2']
        idx_lig, idx_rec = info1['idx'], info2['idx']
        bbox_lig, bbox_rec = info1['bbox'], info2['bbox']
        inv_lig, inv_rec = info1['inv'], info2['inv'] # inv_lig will equal inv_rec
        # Calculate horizontal end-points
        x_lig_end = inv_lig.transform((bbox_lig.x0, 0))[0] - offset
        x_rec_end = inv_rec.transform((bbox_rec.x0, 0))[0] - offset
        # Get PER-AXIS counts for vertical spread
        n_total_lig = ligand_counts_per_ax[ax].get(lr1, 0)
        n_total_rec = receptor_counts_per_ax[ax].get(lr2, 0)
        n_seen_lig = ligand_seen_per_ax[ax].get(lr1, 0)
        n_seen_rec = receptor_seen_per_ax[ax].get(lr2, 0)
        # Calculate vertical positions (using per-axis counts)
        y_spread_lig, y_spread_rec = 0.3, 0.3
        y_lig_data = idx_lig + ((n_seen_lig - (n_total_lig - 1) / 2) * y_spread_lig if n_total_lig > 1 else 0)
        y_rec_data = idx_rec
        if n_total_rec > 1:
            order = (n_total_rec - 1) - n_seen_rec
            y_rec_data += (order - (n_total_rec - 1) / 2) * y_spread_rec
        # Iteratively find a non-overlapping x-position (per-axis check)
        y_span = tuple(sorted((y_lig_data, y_rec_data)))
        current_x = x_left_align
        if n_total_lig > 1 or n_total_rec > 1:
            current_x -= indent_step
        drawn_verticals = drawn_verticals_per_ax[ax] # Check on common axis
        while True:
            has_overlap = any(
                max(y_span[0], s['span'][0]) < min(y_span[1], s['span'][1])
                for s in drawn_verticals if s['x'] == current_x
            )
            if not has_overlap:
                x_left_draw = current_x
                break
            else:
                current_x -= indent_step
        # Draw the bracket on the common axis
        x_points = [x_lig_end, x_left_draw, x_left_draw, x_rec_end]
        y_points = [y_lig_data, y_lig_data, y_rec_data, y_rec_data]
        ax.plot(x_points, y_points, color='black', linewidth=2, zorder=10, clip_on=False, solid_joinstyle='miter')
        # Record the segment that was just drawn for this axis
        drawn_verticals.append({'x': x_left_draw, 'span': y_span})
        # Update seen counters PER AXIS
        ligand_seen_per_ax[ax][lr1] += 1
        receptor_seen_per_ax[ax][lr2] += 1
    # 7. Collapse all drawn verticals to find the furthest point
    all_drawn_verticals = [d for verticals in drawn_verticals_per_ax.values() for d in verticals]
    if all_drawn_verticals:
        furthest_x = min(d['x'] for d in all_drawn_verticals)
    else:
        furthest_x = x_left_align 
    return furthest_x

# Find the major cell type
def get_major(cell_type, level, upper_level='major'):
    if cell_type in miniatlas_df[f'miniatlas_{level}'].values:
        major = miniatlas_df[miniatlas_df[f'miniatlas_{level}'] == cell_type][f'miniatlas_{upper_level}'].values[0]
    else:
        major = 'other'
    return major


# --- Limit to cells that form at least 1% of all cells in TLS regions
# To avoid having too many cell types with very few cells in TLSs that make the plots messy

# Legend positioning cutoffs
# If total genes < genes_cutoff OR total cell types < cell_types_cutoff, legend will be placed to the right of colorbar
genes_cutoff = 6
cell_types_cutoff = 10

annotation_key = level

min_cells_pct = 1
min_cells = min_cells_pct/100 * all_adata.shape[0]
all_adata_sub = all_adata[all_adata.obs[annotation_key].value_counts()[all_adata.obs[annotation_key]].values >= min_cells, :]

print(f'Plotting {all_adata_sub.shape[0]} cells and {all_adata_sub.obs[annotation_key].nunique()} cell types')

# Filter major_sorted for corresponding cell types still present
cells_present = all_adata_sub.obs[annotation_key].unique().tolist()

major_present = [list(set([get_major(c, level) for c in cells_present]))][0]

# Filter major_sorted for corresponding cell types still present
major_sorted = [m for m in major_cell_types_sorted if m in major_present]

# Sort cells by total cell number
cells_sorted = all_adata_sub.obs[annotation_key].value_counts().index.tolist()

all_adata_sub.obs[annotation_key] = pd.Categorical(
    all_adata_sub.obs[annotation_key],
    categories=cells_sorted,
    ordered=True
)

# Extract colors for the sorted cells
color_list = [color_dict[cell] if cell in color_dict else 'lightgrey' for cell in cells_sorted]


# --- Prepare for plotting

# Define color scale
cmap = 'Reds'
vmin = 0
vmax = 1

top_genes_dict = {
    'T cell\nzone': ['CCL19', 'CCR7', 'TCF7'],
    'B cell\nzone': ['CXCL13', 'CXCR5'],
    'HEVs': ['CCL21', 'TSPAN7'],
}

lr_pairs = [('CCL19', 'CCR7'), ('CCL21', 'CCR7'), ('CXCL13', 'CXCR5'), ('CXCL12', 'CXCR4'), ('SPP1', 'CD44'), ('TNFSF14', 'LTBR'),
            ('TNFSF14', 'TNFRSF14'), ('LTA', 'TNFRSF1A'), ('LTA', 'TNFRSF1B'), ('CD40LG', 'CD40'), ('CXCL9', 'CXCR3'), ('IL7', 'IL7R'),
            ('CXCL14', 'CXCR4')]


# Check if genes are in the object
for zone, genes in top_genes_dict.items():
    genes_present = [gene for gene in genes if gene in all_adata_sub.var_names]
    if len(genes_present) < len(genes):
        missing_genes = set(genes) - set(genes_present)
        print(f'Warning: The following genes are missing from the data and will be skipped in the {zone} plot: {missing_genes}')
    top_genes_dict[zone] = genes_present

# Also identify missing genes in lr_pairs
for lig, rec in lr_pairs:
    if lig not in all_adata_sub.var_names:
        print(f'Warning: Ligand {lig} is missing from the data and will be skipped in the ligand-receptor brackets.')
    if rec not in all_adata_sub.var_names:
        print(f'Warning: Receptor {rec} is missing from the data and will be skipped in the ligand-receptor brackets.')
        lr_pairs.remove((lig, rec))


# Also define some cell types of interest
cells_of_interest = []

n_keys = len(top_genes_dict)


# --- Plot 1: Dotplot with per-gene size normalization (0% to max% per gene)

# Plotting
fig, axs = plt.subplots(ncols=1, nrows=n_keys, figsize=(5, 5))

if n_keys == 1:
    axs = [axs]

heatmap_axs = []
max_fracs = []

max_dot_size = 100

for i, (zone, genes) in enumerate(top_genes_dict.items()):
    ax = axs[i]
    # Use custom dotplot with per-gene size normalization
    custom_ax, max_frac = custom_dotplot_per_gene_norm(
        all_adata_sub,
        var_names=genes,
        groupby=annotation_key,
        ax=ax,
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        standard_scale='var',
        max_dot_size=max_dot_size,
        edgecolor='black',
        linewidth=0.5
    )
    max_fracs.append(max_frac)
    
    # Adjust the size of the axes
    auto_resize_dotplot(ax, space_per_dot=0.04, width_ratio=0.8)
    heatmap_axs.append(ax)
    # Adjust the ylims with proper buffer
    n_genes = len(genes)
    buffer = 0.75
    ax.set_ylim(n_genes - 1 + buffer, -buffer)
    ax.set_xlim(-0.75, len(ax.get_xticks()) - 0.25)

# Reposition the dotplot axes
reposition_dotplot_axs(heatmap_axs, y_offset=0.01)

# Force canvas draw after repositioning axes
fig.canvas.draw()

# Create manual horizontal colorbar
ax_pos = heatmap_axs[-1].get_position()
space_per_dot = 0.04
cbar_height = space_per_dot * 0.6
cbar_width = space_per_dot * 3.5

cbar_ax = fig.add_axes([ax_pos.x1 + 0.03, ax_pos.y0, cbar_width, cbar_height])
norm = Normalize(vmin=vmin, vmax=vmax)
cbar = plt.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), cax=cbar_ax, orientation='horizontal')
cbar.ax.set_title('Mean expression\nin group', fontsize=8)
cbar.ax.tick_params(labelsize=8)

# Create horizontal size legend above colorbar
size_ax = fig.add_axes([ax_pos.x1 + 0.03, ax_pos.y0 + cbar_height + 0.09, cbar_width, space_per_dot*1])
size_ax.set_xlim(0, 1)
size_ax.set_ylim(0, 1)
size_ax.axis('off')

# Add horizontal size legend showing min to max
overall_max_frac = max(max_fracs)
size_values = [0.001, 0.5, 1.0]
x_positions = [0.15, 0.5, 0.85]

for size_val, x_pos in zip(size_values, x_positions):
    size_ax.scatter(x_pos, 0.5, s=size_val * max_dot_size, c='lightgray', edgecolor='black', linewidth=0.5, zorder=10)

size_ax.set_ylim(0, 1)

# Add text labels
size_ax.set_title('Fraction of cells\nexpressing', fontsize=8)
size_ax.text(0.15, 0.05, 'min', ha='center', va='top', fontsize=8)
size_ax.text(0.85, 0.05, 'max', ha='center', va='top', fontsize=8)

# Join the ligand-receptor pairs

# Force canvas draw to update tick label positions before drawing brackets
fig.canvas.draw()

x_left_global = join_lr_multi(heatmap_axs, lr_pairs, offset=0.25)

# Draw rectangles from top to bottom of each ax
for n, ax in enumerate(heatmap_axs):
    # Add vertical text
    label = list(top_genes_dict.keys())[n]
    face_color = 'lightskyblue'
    # Convert x_left_global from data to figure coordinates
    x_fig, _ = fig.transFigure.inverted().transform(ax.transData.transform((x_left_global, 0)))
    pos = ax.get_position()
    y0 = pos.y0
    y1 = pos.y1
    height = abs(y1 - y0)
    width = 0.02  # Width in figure coordinates
    x1 = x_fig - 0.02  # Center the rectangle on x_left_global
    highlight_rect = mpl.patches.Rectangle(
        (x1-width, y0), width, height,
        transform=fig.transFigure,
        facecolor=face_color,
        edgecolor='black',
        linewidth=0,
        alpha=1
    )
    fig.patches.append(highlight_rect)
    # Calculate midpoint y
    y_text = (y0 + y1) / 2
    x_text = x1-width - 0.01
    if 'zone' in label:
        label = label.replace(' zone', '\nzone')
    fig.text(
        x_text, y_text,
        label,
        ha='right', va='center',
        color='black',
    )

# Add barplot showing total cell numbers on top of top heatmap ax
top_ax = heatmap_axs[0]
bar_ax = plt.axes([0.1, 0.1, 0.2, 0.2])
top_ax_pos = top_ax.get_position()
bar_ax.set_position([top_ax_pos.x0, top_ax_pos.y1+0.01, top_ax_pos.width, top_ax_pos.height/len(top_ax.get_yticks())*3])

# Get cell counts
cell_counts = all_adata_sub.obs[annotation_key].value_counts().reindex(cells_sorted)

# Plot bars
x = np.arange(len(cells_sorted))
bar_ax.bar(x, cell_counts[cells_sorted], color=color_list, linewidth=0.5, edgecolor='black')

# Set xlim same as heatmap
bar_ax.set_xlim(top_ax.get_xlim()[0], top_ax.get_xlim()[1])

# Set spines inivisible
bar_ax.spines['top'].set_visible(False)
bar_ax.spines['right'].set_visible(False)
bar_ax.spines['left'].set_visible(False)

# Add ylabel
bar_ax.set_ylabel('# cells\n(x10³)', fontsize=8)

# Remove xticks & yticks
bar_ax.set_xticks([])
bar_ax.set_yticks([])

# Annotate the bars with cell counts
for rect, count in zip(bar_ax.patches, cell_counts[cells_sorted]):
    height = rect.get_height()
    if height > 0:
        # Rework the number format if >1000
        if count >= 1000:
            # Round to 1 decimal
            count = f'{count/1000:.1f}'
        bar_ax.text(
            rect.get_x() + rect.get_width() / 2, height + bar_ax.get_ylim()[1]*0.05,
            f'{count}',
            ha='center', va='bottom',
            fontsize=7, rotation=90
        )

# Add a color legend
# Determine position based on number of genes and cell types
total_genes = sum(len(genes) for genes in top_genes_dict.values())
n_cell_types = len(cells_sorted)

if total_genes < genes_cutoff or n_cell_types < cell_types_cutoff:
    # Position to the right of cbar_ax, aligned with bottom
    bottom_edge, top_edge, left_edge, right_edge = position_legend(cbar_ax)
    legend_x = right_edge
    legend_y = bottom_edge
    legend_offset = 0.01
else:
    # Keep current position (below size_ax)
    bottom_edge, top_edge, left_edge, right_edge = position_legend(size_ax)
    legend_x = left_edge
    legend_y = top_edge
    legend_offset = 0

add_color_legend(
    fig,
    labels=major_sorted,
    color_dict=color_dict,
    title='',
    x=legend_x,
    y=legend_y,
    offset=legend_offset,
    spacing_per_label=0.025,
    loc='lower left',
    fontsize=8
)

# Add a title
bar_ax.set_title('Gene expression in TLS regions', fontweight='bold')

# Add some arrows to the barplot ax to indicate cells of interest
for label in cells_of_interest:
    if label in cells_sorted:
        idx = cells_sorted.index(label)
        xpos = idx  # bar index is position
        ypos = cell_counts[label] + bar_ax.get_ylim()[1]*0.15
        bar_ax.annotate(
            "", xy=(xpos, ypos), xycoords=("data", "data"),
            xytext=(xpos, ypos + bar_ax.get_ylim()[1]*0.2), textcoords=("data", "data"),
            arrowprops=dict(arrowstyle="-|>", color="black", lw=1.5, joinstyle='miter')
        )

# Adjust the ylim to make space for the annotations
bar_ax.set_ylim(0, bar_ax.get_ylim()[1]*1.3)

# If subset level, replace xticklabels of the last heatmap ax by the short names
if 'subset' in level:
    last_heatmap_ax = heatmap_axs[-1]
    xticklabels = last_heatmap_ax.get_xticklabels()
    new_labels = []
    for label in xticklabels:
        cell_type = label.get_text()
        if cell_type in miniatlas_sub_dict['miniatlas_subset_short']:
            new_label = miniatlas_sub_dict['miniatlas_subset_short'][cell_type]
        else:
            new_label = cell_type
        new_labels.append(new_label)
    last_heatmap_ax.set_xticklabels(new_labels)
    # Extract adjusted labels for cells of interest
    cells_of_interest_adjusted = [miniatlas_sub_dict['miniatlas_subset_short'][cell] if cell in miniatlas_sub_dict['miniatlas_subset_short'] else cell for cell in cells_of_interest]
else:
    cells_of_interest_adjusted = cells_of_interest

# Position the title
bar_ax_pos = bar_ax.get_position()
title_x = (bar_ax_pos.x0 + bar_ax_pos.x1) / 2
title_y = bar_ax_pos.y1 + 0.2

plt.suptitle(f'All samples - Gene expression in TLS regions (min_cells={min_cells_pct}%, annotation_key={annotation_key})', x=title_x, y=title_y, fontweight='bold')

# --- Save the output

out_dir = os.path.join(res_dir, 'figures')
os.makedirs(out_dir, exist_ok=True)

filename = 'Fig_5C.pdf'
plt.savefig(os.path.join(out_dir, filename), bbox_inches='tight')

# --- Finish

print("Script successfully completed")