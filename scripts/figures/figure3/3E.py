#!/usr/bin/env python3

# --- Aim
# PVL subtype DEG analysis in scRNA-seq reference - matrixplot with gene category markers (Figure 3E)
# Shows top genes per PVL subtype with presence matrix marking contractile, complement, chemokine, cell cycle, and hypoxia response genes

# --- Define environment
# env: vis-hd

# --- Load packages
import argparse
import os

import anndata as ad
import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
import seaborn as sns
import matplotlib as mpl

sc.set_figure_params(facecolor="white", figsize=(8, 8), fontsize=10)
sc.settings.verbosity = 3

mpl.rcParams['pdf.fonttype'] = 42


# --- Argparse arguments
parser = argparse.ArgumentParser(description="PVL subtype DEG analysis for Figure 3E",
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
res_dir = args["res_dir"] or os.path.join(project_dir, 'results')  # default: <project_dir>/results

print(f'Project directory: {project_dir}')


# --- Load data

miniatlas_dir = os.path.join(project_dir, 'data', 'anndata_objects', 'sc_ref')
miniatlas_file = 'miniatlas_GEO.h5ad'

miniatlas = ad.read_h5ad(os.path.join(miniatlas_dir, miniatlas_file))


# --- Subset miniatlas object to PVLs

decon_level = 'subset'

miniatlas = miniatlas[miniatlas.obs['celltype_major'] == 'PVL'].copy()
miniatlas.obs[f'celltype_{decon_level}'].value_counts()

PVL_types = miniatlas.obs[f'celltype_{decon_level}'].unique().tolist()

# Also make a dict with shortened names to update the yticklabels
classifications = pd.read_csv(os.path.join(project_dir, 'resources', 'cell_types', 'miniatlas', 'miniatlas_subset_classifications.csv'))
subset_short_lut = dict(zip(classifications['miniatlas_subset'], classifications['miniatlas_subset_short']))


# --- Set var_names as gene symbols for interpretability

gene_symbol_name = 'original_name'

miniatlas.var['ensembl_id'] = miniatlas.var_names.copy()
miniatlas.var_names = miniatlas.var[gene_symbol_name].astype(str)


# --- Renormalize after subsetting

if miniatlas.raw is not None:
    miniatlas_renorm = miniatlas.raw.to_adata().copy()
else:   # the P08 miniatlas_GEO.h5ad keeps the raw counts in layers['counts'] instead of .raw
    miniatlas_renorm = sc.AnnData(X=miniatlas.layers['counts'].copy(), obs=miniatlas.obs.copy(), var=miniatlas.var.copy())
miniatlas_renorm.var_names = miniatlas_renorm.var[gene_symbol_name].astype(str)

miniatlas_renorm.layers['norm'] = sc.pp.normalize_total(miniatlas_renorm, target_sum=1e4, inplace=False)['X']
miniatlas_renorm.layers['log_norm'] = sc.pp.log1p(miniatlas_renorm.layers['norm'], copy=True)
miniatlas_renorm.layers['scaled'] = sc.pp.scale(miniatlas_renorm.layers['log_norm'], max_value=10, copy=True)

# Make sure var_names are unique
miniatlas_renorm.var_names_make_unique()


# --- Find DEGs

use_layer = 'log_norm'

rank_genes_key = 'rank_genes_PVL'
groupby = f'celltype_{decon_level}'

sc.tl.rank_genes_groups(miniatlas_renorm, groupby, use_raw=False, method='wilcoxon', key_added=rank_genes_key,
                        pts=True, layer=use_layer)


# --- Extract the DEG table

dfs = []

result = miniatlas_renorm.uns[rank_genes_key]
groups = result['names'].dtype.names
reference = result['params']['reference']

# Iterate through the groups
for group in groups:
    # Create a DataFrame for the current group
    df = pd.DataFrame({
        'group': group,
        'reference': reference,
        'genes': result['names'][group],
        'pvals': result['pvals'][group],
        'logfoldchanges': result['logfoldchanges'][group],
        'pts': result['pts'][group],
        'pvals_adj': result['pvals_adj'][group],
        'z-score': result['scores'][group]
    })
    
    # Calculate pct_expressing and pct_expressing_rest
    gene_list = df['genes'].tolist()
    adata_sub = miniatlas_renorm[miniatlas_renorm.obs[groupby] == group]
    adata_rest = miniatlas_renorm[miniatlas_renorm.obs[groupby] != group]
    
    # Convert to DataFrame and calculate percentage expressing
    df_sub = adata_sub.to_df(layer=use_layer)[gene_list]
    pct_expressing_list = (df_sub > 0).mean(axis=0).tolist()
    
    # Same for rest
    df_rest = adata_rest.to_df(layer=use_layer)[gene_list]
    pct_expressing_rest_list = (df_rest > 0).mean(axis=0).tolist()
    
    # Add to dataframe
    df['pct_expressing'] = pct_expressing_list
    df['pct_expressing_rest'] = pct_expressing_rest_list
    
    # Append the DataFrame to the list
    dfs.append(df)

# Concatenate the DataFrames from all groups into one DataFrame
final_df = pd.concat(dfs, ignore_index=True)

# Extract top 10 genes per group
top = 10
top_genes_per_group = final_df.groupby('group').head(top)

# Make a dict
top_genes_dict = {}
for group in top_genes_per_group['group'].unique():
    top_genes = top_genes_per_group[top_genes_per_group['group'] == group]['genes'].tolist()
    top_genes_dict[group] = top_genes


# --- Create presence matrix for gene categories

# Extract a flat list of all genes
top_genes = list(top_genes_dict.values())
top_genes = [gene for sublist in top_genes for gene in sublist]

presence_df = pd.DataFrame(index=None, columns=top_genes)
presence_df = presence_df.astype(float)

# Add contractile genes
presence_df.loc['Contractile', :] = 0
contractile_genes = ['ACTA2', 'TAGLN', 'MYH11', 'CNN1', 'MYL9', 'MYL12B', 'MYL12A', 'MYLK2', 'MYLK3']
for gene in top_genes:
    if gene in contractile_genes:
        presence_df.loc['Contractile', gene] = 1

# Add complement genes
presence_df.loc['Complement', :] = 0
for gene in top_genes:
    if gene.startswith('C1') or gene.startswith('C3'):
        presence_df.loc['Complement', gene] = 1

# Add chemokines
presence_df.loc['Chemokine', :] = 0
for gene in top_genes:
    if gene.startswith('CCL') or gene.startswith('CXCL'):
        presence_df.loc['Chemokine', gene] = 1

# Add cell cycle genes
presence_df.loc['Cell cycle', :] = 0
cell_cycle_genes = ['MKI67', 'TOP2A', 'PCNA', 'CCNB1', 'CCNB2', 'CDK1', 'CDK2', 'CDC20', 'CDC25C', 'AURKB',
                   'STMN1', 'PCLAF', 'CKS2']
for gene in top_genes:
    if gene in cell_cycle_genes:
        presence_df.loc['Cell cycle', gene] = 1

# Add hypoxia response genes
presence_df.loc['Hypoxia response', :] = 0
hypoxia_genes = ['HIGD1B', 'RGS5', 'NDUFA4L2']
for gene in top_genes:
    if gene in hypoxia_genes:
        presence_df.loc['Hypoxia response', gene] = 1


# --- Plot matrixplot with presence markers

# Edit the cluster names on top
name_dict = {
    'PVL Differentiated s3': 'PVL s3\n(Contractile)',
    'PVL Immature s1': 'PVL s1\n(Inflammatory)',
    'PVL_Immature s2': 'PVL s2\n(Hypoxia response)',
}

sc.pl.rank_genes_groups_matrixplot(miniatlas_renorm, n_genes=top, groupby=groupby, key=rank_genes_key, use_raw=False, layer=use_layer,
                                   standard_scale='var', cmap='viridis')

# Adjust the size of the axes
axs = plt.gcf().get_axes()

heatmap_ax = axs[0]

# Set the desired ratio
heatmap_ax.set_aspect(1.75)

# Adjust the other axes
new_heatmap_pos = heatmap_ax.get_position()

bar_ax = axs[1]
bar_ax_pos = bar_ax.get_position()
bar_ax.set_position([new_heatmap_pos.x1, new_heatmap_pos.y0, bar_ax_pos.width*0.5, new_heatmap_pos.height])

# Move the colorbar
cbar_ax = axs[-1]
bar_ax_pos = bar_ax.get_position()
cbar_ax_pos = cbar_ax.get_position()
cbar_ax.set_position([bar_ax_pos.x1+bar_ax_pos.width*0.3, bar_ax_pos.y0, cbar_ax_pos.width*0.8, cbar_ax_pos.height])

# Adjust the top ax
top_ax = axs[2]
top_ax_pos = top_ax.get_position()
top_ax.set_position([new_heatmap_pos.x0, top_ax_pos.y0, new_heatmap_pos.width, top_ax_pos.height])

# Update cluster labels on top
for label in top_ax.texts:
    text = label.get_text()
    if text in name_dict:
        label.set_text(name_dict[text])
    label.set_rotation(0)
    label.set_ha('center')

# Add an ax for the presence matrix markers
presence_ax = plt.axes([0.1, 0.1, 0.2, 0.2])

# Get xticklabels as a list of labels
xticklabels = [label.get_text() for label in heatmap_ax.get_xticklabels()]

# Drop duplicate columns if any
presence_df = presence_df.loc[:, ~presence_df.columns.duplicated()]
presence_df = pd.concat([presence_df[[col]] for col in xticklabels], axis=1)
presence_df.columns = xticklabels

# Plot the presence matrix
sns.heatmap(presence_df, ax=presence_ax, cmap='Blues', vmin=-0.1, vmax=1.1,
            cbar=False, xticklabels=True, yticklabels=True, lw=1, linecolor='white')

# Move the presence ax to the right position
presence_ax_height = new_heatmap_pos.height/len(PVL_types)*len(presence_df.index)/2
offset = new_heatmap_pos.height/len(PVL_types) * 0.2

presence_ax.set_position([new_heatmap_pos.x0, new_heatmap_pos.y1 + offset, new_heatmap_pos.width, presence_ax_height])

# Disable the grid
presence_ax.grid(False)

# Set spines visible
for spine in presence_ax.spines.values():
    spine.set_visible(True)
    spine.set_linewidth(1.5)

# Remove xticks
presence_ax.set_xticks([])

# Rotate yticklabels
presence_ax.set_yticklabels(presence_ax.get_yticklabels(), rotation=0, fontsize=8)

# Move top ax above presence ax
presence_ax_pos = presence_ax.get_position()
top_ax.set_position([new_heatmap_pos.x0, presence_ax_pos.y1, new_heatmap_pos.width, top_ax_pos.height])

# Add lines separating the clusters
for i, cluster in enumerate(PVL_types):
    x_pos = i * top
    presence_ax.axvline(x_pos, color='lightsteelblue', linewidth=1, linestyle='--')

# Update the yticklabels of the heatmap
new_yticklabels = []

for label in heatmap_ax.get_yticklabels():
    text = label.get_text()
    if text in subset_short_lut:
        new_yticklabels.append(subset_short_lut[text])
    else:
        new_yticklabels.append(text)

heatmap_ax.set_yticklabels(new_yticklabels, rotation=0, fontsize=8)

# Calculate position for plot title
title_x = (new_heatmap_pos.x0 + new_heatmap_pos.x1) / 2
title_y = top_ax.get_position().y1 + 0.3

plt.suptitle(f'Miniatlas PVLs - top {top} genes per PVL subtype', y=title_y, x=title_x, fontweight='bold')


# --- Save the output

out_dir = os.path.join(res_dir, 'figures')
os.makedirs(out_dir, exist_ok=True)

filename = 'Fig_3E.pdf'
plt.savefig(os.path.join(out_dir, filename), bbox_inches='tight', dpi=500)

print(f"Saved plot to: {os.path.join(out_dir, filename)}")
print("Script successfully completed")
