#!/usr/bin/env python3

# --- Aim
# Plot Figure 4C: STING signature scores in venous ECs stratified by environment (canonical venous, stromal HEV-like & TLS-associated HEV)

# --- Define environment
# env: xen-hev

# --- Load packages
import argparse
import os
from datetime import datetime
import sys

import anndata as ad
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

import scanpy as sc

import pandas as pd
import numpy as np

import seaborn as sns

sc.set_figure_params(figsize=(8, 8), fontsize=20)
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

# Import function script as a module
sys.path.append(os.path.join(project_dir, 'scripts', '05_xenium_mechanisms',))
import plotting_functions


# --- Extract some info from the sample sheet

# Use sample_sheet argument if provided, else default
if args.get("sample_sheet"):
    sample_sheet_path = args["sample_sheet"]
else:
    sample_dir = os.path.join(project_dir, 'config')
    sample_sheet_file = 'sample_sheet_xenium.csv'
    sample_sheet_path = os.path.join(sample_dir, sample_sheet_file)

sample_sheet = pd.read_csv(sample_sheet_path)

# Also load tissue metadata
metadata_dir = os.path.join(project_dir, 'resources',)
tissue_metadata_file = 'tissue_metadata.csv'
tissue_metadata = pd.read_csv(os.path.join(metadata_dir, tissue_metadata_file))

# Make a section_id to subtype mapping
section_to_subtype = tissue_metadata.set_index('tissue_id')['clinical_subtype'].to_dict()
# Map section_id to subtype
sample_sheet['subtype'] = sample_sheet['tissue_id'].map(section_to_subtype)

subtype_dict = sample_sheet.set_index('section_id')['subtype'].to_dict()


# --- Load the anndata objects

adata_dict = {}

for sample in sample_list:
    adata = ad.read_h5ad(os.path.join(h5ad_dir, 'xenium', f'{sample}.h5ad'))
    adata_dict[sample] = adata


# --- Load the HEV annotations

# Find most recent results
def find_most_recent(dir):
    # Find most recent results
    subdirs = [d for d in os.listdir(dir) if os.path.isdir(os.path.join(dir, d))]
    # Filter subdirectories to keep only those with date format YYYYmmdd
    date_format = '%Y%m%d'
    valid_subdirs = []
    for subdir in subdirs:
        try:
            datetime.strptime(subdir, date_format)
            valid_subdirs.append(subdir)
        except ValueError:
            pass
    # Sort the valid subdirectories by date and get the most recent one
    if valid_subdirs:
        most_recent_dir = max(valid_subdirs, key=lambda x: datetime.strptime(x, date_format))
        most_recent_dir_path = os.path.join(dir, most_recent_dir)
        print("Most recent results:", most_recent_dir_path)
        return most_recent_dir_path
    else:
        print("No valid directories found.")

distance_dir = os.path.join(res_dir, '05_xenium_mechanisms', '02_identify_HEV_regions')
most_recent_distance = find_most_recent(distance_dir)

for sample in sample_list:
    # Also load HEV niche barcode and add to adata
    filename = f'{sample}_HEV_niches_cell_to_niche.csv'
    niche_barcode_df = pd.read_csv(os.path.join(most_recent_distance, 'csv_files', filename), index_col=0)
    niche_barcode_df.dropna(subset=['niche_id'], inplace=True)
    # Prepend sample name to niche_id
    niche_barcode_df['niche_id'] = sample + '_' + niche_barcode_df['niche_id'].astype(str)
    adata_dict[sample].obs['niche_id'] = niche_barcode_df['niche_id']
    # Also load the distance calculations
    distance_file = f'{sample}_HEV_niches_info.csv'
    distance_df = pd.read_csv(os.path.join(most_recent_distance, 'csv_files', distance_file))
    # Prepend sample name to niche_id
    distance_df['niche_id'] = sample + '_' + distance_df['niche_id'].astype(str)
    # Make a mapping dict of distance to TLS
    distance_dict = distance_df.set_index('niche_id')['close_to_TLS'].to_dict()
    # Map to adata
    adata_dict[sample].obs['close_to_TLS'] = adata_dict[sample].obs['niche_id'].map(distance_dict)


# --- Also load cell annotations

# Decide which annotation levels to use
level = 'subset'

for sample in sample_list:
    annotation_file = f'{sample}_{level}_pred_cell_types.csv'
    annotation_dir = os.path.join(project_dir, 'results', '02_cell_type_mapping', 'label_transfer', 'xenium')
    annotations = pd.read_csv(os.path.join(annotation_dir, f'{level}_level', annotation_file), index_col=0)
    # Add to adata
    adata_dict[sample].obs[level] = annotations['pred_cell_type']


# --- Also load the STING signature gene lists

signature_dir = os.path.join(project_dir, 'resources', 'gene_lists', 'STING')
signature_file = 'STING_signatures.csv'

gene_df = pd.read_csv(os.path.join(signature_dir, signature_file))

# Also make a dict for singature type
signature_type_dict = {}

# Create a dict of gene lists
signature_dict = {}
for sig in gene_df['signature'].unique():
    genes = gene_df[gene_df['signature'] == sig]['gene'].tolist()
    signature_dict[sig] = genes

signature_type_dict['STING'] = gene_df['signature'].unique().tolist()

# Also add the type I interferon related signatures
signature_file = 'type_I_interferon_signatures.csv'
gene_df = pd.read_csv(os.path.join(signature_dir, signature_file))

for sig in gene_df['signature'].unique():
    genes = gene_df[gene_df['signature'] == sig]['gene'].tolist()
    signature_dict[sig] = genes

signature_type_dict['Type I IFN'] = gene_df['signature'].unique().tolist()

all_score_names = []
for sig_type in signature_type_dict.keys():
    all_score_names.extend(signature_type_dict[sig_type])


# --- Load the scoring results

score_dir = os.path.join(project_dir, 'results', '05_xenium_mechanisms', '04_score_signatures', '01_score_STING_signatures', 'csv_files')
most_recent_score_dir = find_most_recent(score_dir)

for sample in sample_list:
    score_file = f'{sample}_STING_signature_scores.csv'
    score_df = pd.read_csv(os.path.join(most_recent_score_dir, score_file), index_col=0)
    adata_dict[sample].obs[all_score_names] = score_df[all_score_names]
    # Also add niche_id info
    adata_dict[sample].obs['niche_TLS_proximity'] = score_df['niche_TLS_proximity']
    # And cell type annotation
    adata_dict[sample].obs['cell_annotation'] = score_df['subset']


# --- Subset to endothelial ACKR1+ cells only

# Subset adatas to HEV / venous niches only
for sample in sample_list:
    adata = adata_dict[sample]
    adata_HEV = adata[adata.obs['cell_annotation'] == 'Endothelial ACKR1'].copy()
    adata_dict[sample] = adata_HEV

# Extract a single adata object with all samples
adata_list = list(adata_dict.values())
all_adata = ad.concat(adata_list, merge='same')

all_adata.obs_names_make_unique()

# Filter genes to those present in the adata
for sig in signature_dict.keys():
    genes_in_adata = [gene for gene in signature_dict[sig] if gene in all_adata.var_names]
    print(f'Signature {sig}: {len(genes_in_adata)}/{len(signature_dict[sig])} genes present in the data.')
    print(f'Genes missing: {[gene for gene in signature_dict[sig] if gene not in all_adata.var_names]}')
    signature_dict[sig] = genes_in_adata


# --- Load the endothelial colors

color_dir = os.path.join(project_dir, 'resources', 'plotting_colors')
filename = 'endothelial_colors.csv'

endothelial_colors_df = pd.read_csv(os.path.join(color_dir, filename), index_col=0)
color_dict = endothelial_colors_df['color'].to_dict()


# --- Also add the HEV colors

filename = 'TLS_HEV_colors.csv'

HEV_colors_df = pd.read_csv(os.path.join(color_dir, filename), index_col=0)

# Add to color_dict
for cell_type, row in HEV_colors_df.iterrows():
    color_dict[cell_type] = row['color']

color_dict['Venous'] = color_dict['Endothelial ACKR1']
color_dict['HEV far from TLS'] = color_dict['non_TLS_HEV']
color_dict['HEV close to TLS'] = color_dict['near_TLS_HEV']

# Rename for consistency
rename_dict = {
    'HEV far from TLS': 'HEV-like (stroma)',
    'HEV close to TLS': 'HEV (TLS)',
}

# --- Prepare data for plotting

# Make ordered
plot_order_TLS = ['Venous', 'HEV far from TLS', 'HEV close to TLS']
all_adata.obs['niche_TLS_proximity'] = pd.Categorical(all_adata.obs['niche_TLS_proximity'], categories=plot_order_TLS, ordered=True)

annotation_key = 'niche_TLS_proximity'

# Prepare data for dotplot
dotplot_data = []
for sig in signature_dict.keys():
    for cell_type in all_adata.obs[annotation_key].astype('category').cat.categories:
        subset_mask = all_adata.obs[annotation_key] == cell_type
        scores = all_adata.obs.loc[subset_mask, sig]
        mean_score = scores.mean()
        median_score = scores.median()
        pct_above_0 = (scores > 0).sum() / len(scores) * 100
        dotplot_data.append({'signature': sig,
                             'cell_type': cell_type,
                             'mean_score': mean_score,
                             'median_score': median_score,
                             'pct_above_0': pct_above_0})

dotplot_df = pd.DataFrame(dotplot_data)

# Also get counts of HEV ACKR1+ per TLS proximity category
hev_counts = all_adata.obs[annotation_key].value_counts().reindex(plot_order_TLS).fillna(0).astype(int)

# Get number of genes per signature
signature_gene_counts = {sig: len(signature_dict[sig]) for sig in signature_dict.keys()}
signature_counts_df = pd.DataFrame.from_dict(signature_gene_counts, orient='index', columns=['n_genes'])

# Sort signatures by number of genes
signature_counts_df = signature_counts_df.sort_values(by='n_genes', ascending=True)
dotplot_df['signature'] = pd.Categorical(dotplot_df['signature'], categories=signature_counts_df.index.tolist(), ordered=True)

# Get vmin & vmax for mean_score
vmin = dotplot_df['mean_score'].min()
vmax = dotplot_df['mean_score'].max()

# Also get size_min & size_max
size_min = dotplot_df['pct_above_0'].min()
size_max = dotplot_df['pct_above_0'].max()

# Create a truncated RdBu_r colormap to avoid too-dark colors at extremes
cmap_full = plt.cm.RdBu_r
cmap_RdBu_truncated = LinearSegmentedColormap.from_list(
    'RdBu_r_truncated',
    cmap_full(np.linspace(0.15, 0.85, 256))
)

# Set fontsize
sc.set_figure_params(fontsize=10)

# Turn the grid off
plt.rcParams['axes.grid'] = False

# Set sizemin
min_size = size_min*0.98


# --- Plot 1: Dotplot of STING signature scores in HEV vs non-HEV niches, with count annotations on top

# Also, split the IFN scores into production & response to better see trends
signature_type_dict_zscore = {
    'STING': signature_type_dict['STING'],
    'Type I IFN\nProduction': [],
    'Type I IFN\nResponse': [],
}

for sig in signature_type_dict['Type I IFN']:
    if 'production' in sig.lower():
        signature_type_dict_zscore['Type I IFN\nProduction'].append(sig)
    else:
        signature_type_dict_zscore['Type I IFN\nResponse'].append(sig)

# Create a color dict for signature type
signature_types = signature_type_dict_zscore.keys()
color_palette = ['#f16912', '#440155', '#817dba']

signature_type_color_dict = {sig_type: color_palette[n] for n, sig_type in enumerate(signature_types)}

# Create a mapping of signature to signature type color
signature_to_type_color = {}
for sig_type, sig_list in signature_type_dict_zscore.items():
    for sig in sig_list:
        signature_to_type_color[sig] = signature_type_color_dict[sig_type]

# Make a df for z-score dotplot
z_score_dotplot_df = dotplot_df.copy()

# Z-score mean_score per signature
hue_key = 'mean_score'

for sig in signature_dict.keys():
    sig_mask = z_score_dotplot_df['signature'] == sig
    mean_scores = z_score_dotplot_df.loc[sig_mask, hue_key]
    z_scores = (mean_scores - mean_scores.mean()) / (mean_scores.std() + 1e-9)
    z_score_dotplot_df.loc[sig_mask, hue_key] = z_scores

# Pick a cmap - use truncated version to avoid too-dark colors
cmap = cmap_RdBu_truncated

vmin = z_score_dotplot_df[hue_key].min()
vmax = z_score_dotplot_df[hue_key].max()

# Set vmax & vmin equal distance from 0
abs_max = max(abs(vmin), abs(vmax))
vmax = abs_max
vmin = -abs_max

# Make the dotplot
fig, top_ax = plt.subplots(figsize=(2, 2))

axs = []

for n, (signature_type, signature_list) in enumerate(signature_type_dict_zscore.items()):
    print(f'{signature_type}: {signature_list}')
    if n == 0:
        ax = top_ax
    else:
        ax_width = ax_prev.get_position().width
        ax_height = ax_prev.get_position().height / len(ax_prev.get_yticks()) * len(signature_list)
        ax = plt.gcf().add_axes([ax_prev.get_position().x0, ax_prev.get_position().y0 - ax_height - 0.02, ax_width, ax_height])
    plot_df = z_score_dotplot_df[z_score_dotplot_df['signature'].isin(signature_list)].copy()
    # Make sure signature order matches signature_counts_df
    plot_df['signature'] = pd.Categorical(plot_df['signature'], categories=signature_counts_df.index.tolist(), ordered=True)
    # Remove unused categories
    plot_df['signature'] = plot_df['signature'].cat.remove_unused_categories()
    dotplot_ax, size_max = plotting_functions.seaborn_dotplot(plot_df, ax, x='cell_type', y='signature', hue=hue_key, size='pct_above_0',
                                                        max_size_px=250, cmap=cmap, vmin=vmin, vmax=vmax, size_max=size_max, size_min=min_size)
    if n == 0:
        # Resize ax based on xticks & yticks
        ax_pos = ax.get_position()
        ax.set_position([ax_pos.x0, ax_pos.y0, 0.15 * len(ax.get_xticks()), 0.12 * len(ax.get_yticks())])
    # Add y label
    ax.set_ylabel(signature_type, fontweight='bold')
    # Rotate yticklabels
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0, ha='right', fontsize=8)
    ax_prev = ax
    axs.append(ax)
    # Adjust xlim & ylim
    ax.set_xlim(-0.5, len(ax.get_xticks()) - 0.5)
    ax.set_ylim(-0.75, len(ax.get_yticks()) - 0.25)
    # Rotate xticklabels and rename using rename_dict
    xticklabels = [tick.get_text() for tick in ax.get_xticklabels()]
    xticklabels = [rename_dict.get(label, label) for label in xticklabels]
    ax.set_xticklabels(xticklabels, rotation=90, ha='center')

# Remove xticks from all but the bottom ax
for ax in axs[:-1]:
    ax.set_xticks([])
    # Also remove legend
    ax.get_legend().remove()

# Add xlabel
axs[-1].set_xlabel('Vessel type', fontweight='bold')

# Add countax at the top
count_ax = plt.gcf().add_axes([top_ax.get_position().x0, top_ax.get_position().y1 + 0.02, top_ax.get_position().width, top_ax.get_position().height * 0.3])
sns.barplot(x=hev_counts.index, y=hev_counts.values, color='lightgrey', edgecolor='black', ax=count_ax, hue=hev_counts.index, palette=color_dict,
               saturation=1, width=0.65)

# Beautify
plotting_functions.beautify_plot(count_ax)

# Annotate counts
for n, p in enumerate(count_ax.patches):
    # Annotate with k if over 1000
    buffer = count_ax.get_ylim()[1] * 0.05
    if p.get_height() >= 1000:
        annotate_value = f'{round(p.get_height()/1000, 1)}k'
    else:
        annotate_value = f'{int(p.get_height())}'
    count_ax.annotate(annotate_value, (p.get_x() + p.get_width() / 2., p.get_height() + buffer), ha='center', va='bottom', fontsize=8, rotation=90)

# Set top and right spine invisible
count_ax.spines['top'].set_visible(False)
count_ax.spines['right'].set_visible(False)

# Remove xticks
count_ax.set_xticks([])
count_ax.set_xlabel('')
count_ax.set_ylabel('# cells')

# Add title
count_ax.set_title('Gene signatures\n(ACKR1+ Endothelial)', fontweight='bold', fontsize=10, pad=30)

count_axs = []

# Also add a countax on the side for the number of genes per signature
for n, (signature_type, signature_list) in enumerate(signature_type_dict_zscore.items()):
    dotplot_ax = axs[n]
    gene_count_ax = plt.gcf().add_axes([dotplot_ax.get_position().x1 + 0.02, dotplot_ax.get_position().y0, 0.12*1.5, dotplot_ax.get_position().height])
    signature_counts_subset = signature_counts_df[signature_counts_df.index.isin(signature_list)]
    # Sort by the order in the dotplot
    ytick_texts = [tick.get_text() for tick in dotplot_ax.get_yticklabels()]
    signature_counts_subset = signature_counts_subset.reindex(ytick_texts, axis=0)
    # Plot horizontal barplot
    sns.barplot(y=signature_counts_subset.index, x=signature_counts_subset['n_genes'], hue=signature_counts_subset.index, palette=signature_to_type_color, ax=gene_count_ax, saturation=1)
    # Beautify
    plotting_functions.beautify_plot(gene_count_ax)
    # Remove yticks
    gene_count_ax.set_yticks([])
    gene_count_ax.set_ylabel('')
    # Set xlabel on the last one only
    if n == len(signature_type_dict_zscore) - 1:
        gene_count_ax.set_xlabel('# genes')
    else:
        gene_count_ax.set_xlabel('')
        # Also remove xticklabels
        gene_count_ax.set_xticklabels([])
    # Adjust xlim
    gene_count_ax.set_xlim(right=signature_counts_df['n_genes'].max() * 1.4)
    # Annotate bars
    annotation_offset = gene_count_ax.get_xlim()[1] * 0.05
    for n, p in enumerate(gene_count_ax.patches):
        gene_count_ax.annotate(f'{int(p.get_width())}', (p.get_width() + annotation_offset, p.get_y() + p.get_height() / 2.), ha='left', va='center', fontsize=8)
    # Adjust ylim
    gene_count_ax.set_ylim(dotplot_ax.get_ylim())
    # Set top and right spine invisible
    gene_count_ax.spines['top'].set_visible(False)
    gene_count_ax.spines['right'].set_visible(False)
    count_axs.append(gene_count_ax)

# For all but the last count_ax, keep only the second xtick
for gene_count_ax in count_axs[:-1]:
    gene_count_ax.set_xticks([gene_count_ax.get_xticks()[1]])

gene_count_ax = count_axs[-1]

# Add an ax for colorbar
cbar_ax = plt.gcf().add_axes([gene_count_ax.get_position().x1 + 0.02, gene_count_ax.get_position().y0, 0.12*0.5, 0.12*3])

# Add colorbar
norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
cbar = plt.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), cax=cbar_ax, orientation='vertical')

# Add title
hue_title = hue_key.replace('_', ' ').capitalize()
cbar_ax.set_title(f'{hue_title}\n(z-score)', fontsize=10, ha='left', position=(0, 1.05))

# Set xlim to actual vmin & vmax
vmax = z_score_dotplot_df[hue_key].max()
vmin = z_score_dotplot_df[hue_key].min()
cbar.set_ticks([vmin, 0, vmax])
cbar.set_ticklabels([f'{vmin:.2f}', '0', f'{vmax:.2f}'])
cbar_ax.set_ylim(vmin, vmax)

# Add size legend
size_legend = plotting_functions.add_size_legend(dotplot_ax, cbar_ax, size='pct_above_0', size_title='Pct. > 0', limit_circles=True,
                                                alignment='left', bbox_to_anchor=(-0.3, 1.5), loc='lower left', labelspacing=0.8, handletextpad=0.8)

# Adjust main title position
title_x = (dotplot_ax.get_position().x0 + dotplot_ax.get_position().x1) / 2
title_y = count_ax.get_position().y1 + 0.55

plt.suptitle(f'All samples - STING & type I IFN signature scores in endothelial ACKR1+ (annotation_level={level})', fontweight='bold', y=title_y, x=title_x)


# --- Save output

out_dir = os.path.join(res_dir, 'figures')
os.makedirs(out_dir, exist_ok=True)

filename = 'Fig_4C.pdf'
plt.savefig(os.path.join(out_dir, filename), bbox_inches='tight')


# --- Finish

print("Script successfully completed")