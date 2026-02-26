#!/usr/bin/env python3

# --- Aim
# Score STING-related gene signatures across cell types in Xenium data.
# Identify cell populations with high STING pathway activity.

# --- Define environment
# env: xen-hev

# --- Load packages

import argparse
import os
from datetime import datetime, date
import sys

import anndata as ad
import scanpy as sc
import pandas as pd
import numpy as np

import matplotlib as mpl
import matplotlib.pyplot as plt

import seaborn as sns

# Import custom functions
from sophiesfunctions.save_multi_image import save_multi_image

sc.set_figure_params(figsize=(8, 8), fontsize=20)
mpl.rcParams['pdf.fonttype'] = 42


# --- Argparse arguments

parser = argparse.ArgumentParser(description="Import Visium data & save as .h5ad file",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

# Required

parser.add_argument("--raw_data_dir", help="Path to raw data")
parser.add_argument("--h5ad_dir", help="Path to directory for .h5ad output files")
parser.add_argument("--sample_id", help="Sample ID")
parser.add_argument("--project_dir", help="Path to the root of the project directory")
parser.add_argument("--outs_path", help="Path to spaceranger 'outs'")
parser.add_argument("--use_data", help="Which data to use (log_norm or SCT)")
parser.add_argument("--res_dir", help="Path to results directory")
parser.add_argument("--seed", help="Seed for reproducibility", type=int)
parser.add_argument("--sample_sheet", help="Path to sample sheet CSV file", default=None)


args = vars(parser.parse_args())


# --- Arguments

h5ad_dir=args["h5ad_dir"]
sample_id=args["sample_id"]

res_dir = args['res_dir']
project_dir = args['project_dir']

# Import function script as a module
sys.path.append(os.path.join(project_dir, 'scripts', '05_xenium_mechanisms'))
import plotting_functions


# --- Extract some info from the sample sheet
# Use sample_sheet argument if provided, else default
if args.get("sample_sheet"):
    sample_sheet_path = args["sample_sheet"]
else:
    sample_dir = os.path.join(project_dir, 'config')
    sample_sheet_file = 'sample_sheet.csv'
    sample_sheet_path = os.path.join(sample_dir, sample_sheet_file)

sample_sheet = pd.read_csv(sample_sheet_path)
tissue_id = sample_sheet[sample_sheet['section_id'] == sample_id]['tissue_id'].values[0]

# Also load tissue metadata
metadata_dir = os.path.join(project_dir, 'resources',)
tissue_metadata_file = 'tissue_metadata.csv'
tissue_metadata = pd.read_csv(os.path.join(metadata_dir, tissue_metadata_file))

subtype = tissue_metadata[tissue_metadata['tissue_id'] == str(tissue_id)]['clinical_subtype'].values[0]


# --- Load data

adata = ad.read_h5ad(os.path.join(h5ad_dir, 'xenium', f'{sample_id}.h5ad'))


# --- Load the tangram annotations

# Decide which annotation level to use
level = 'major'

annotation_dir = os.path.join(project_dir, 'results', '02_cell_type_mapping', 'label_transfer', 'xenium', f'{level}_level')
annotation_file = f'{sample_id}_{level}_pred_cell_types.csv'
annotations = pd.read_csv(os.path.join(annotation_dir, annotation_file), index_col=0)

# Add to anndata object
annotation_key = 'tangram_annotation'
probability_key = 'tangram_probability'

adata.obs[annotation_key] = annotations['pred_cell_type']
adata.obs[probability_key] = annotations['pred_cell_type_prob']

# Also add the subset annotations
level = 'subset'

annotation_dir = os.path.join(project_dir, 'results', '02_cell_type_mapping', 'label_transfer', 'xenium', f'{level}_level')
annotation_file = f'{sample_id}_{level}_pred_cell_types.csv'
annotations = pd.read_csv(os.path.join(annotation_dir, annotation_file), index_col=0)
adata.obs[level] = annotations['pred_cell_type']


# --- Load plotting colors

level = 'major'
color_dir = os.path.join(project_dir, 'resources', 'plotting_colors')
color_file = f'miniatlas_{level}_colors.csv'

color_dict = pd.read_csv(os.path.join(color_dir, color_file), index_col=0).to_dict()['color']

all_cell_types = adata.obs[annotation_key].astype('category').cat.categories

# Make sure all cell types have a color
for cell_type in all_cell_types:
    if cell_type not in color_dict.keys():
        print(f'Adding color for {cell_type}')
        color_dict[cell_type] = 'lightgrey'

# Also add subset colors
miniatlas_dir = os.path.join(project_dir, 'resources', 'cell_types', 'miniatlas')
miniatlas_file = 'miniatlas_subset_classifications.csv'

miniatlas_df = pd.read_csv(os.path.join(miniatlas_dir, miniatlas_file),)

# Add the other levels to the color_dict
level = 'subset'

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
        color_dict[cell_type] = 'lightgray'

# Also define a dict with shorter names for the subset level cell types
miniatlas_sub_dict = miniatlas_df.set_index('miniatlas_subset')[['miniatlas_subset_short']].to_dict()


# --- Sort by lineage, then by major

lineage_order = ['Epithelial', 'Immune', 'Stromal']

major_cell_types = miniatlas_df[f'miniatlas_major'].dropna().unique().tolist()
minor_cell_types = miniatlas_df[f'miniatlas_minor'].dropna().unique().tolist()
subset_cell_types = miniatlas_df[f'miniatlas_subset'].dropna().unique().tolist()

# Sort based on lineage first (recorded in the miniatlas df)
def get_lineage(cell_type, level='major'):
    if cell_type in miniatlas_df[f'miniatlas_{level}'].values:
        lineage = miniatlas_df[miniatlas_df[f'miniatlas_{level}'] == cell_type]['lineage'].values[0]
    else:
        lineage = 'other'
    return lineage

major_cell_types_sorted = sorted(major_cell_types, key=lambda x: (lineage_order.index(get_lineage(x)) if get_lineage(x) in lineage_order else len(lineage_order), x))

# Now sort minor level by corresponding major level
def get_major(cell_type, level, upper_level='major'):
    level_clean = level.replace('_stepwise', '')
    if cell_type in miniatlas_df[f'miniatlas_{level_clean}'].values:
        major = miniatlas_df[miniatlas_df[f'miniatlas_{level_clean}'] == cell_type][f'miniatlas_{upper_level}'].values[0]
    else:
        major = 'other'
    return major

minor_cell_types_sorted = sorted(minor_cell_types, key=lambda x: (major_cell_types_sorted.index(get_major(x, level='minor')) if get_major(x, level='minor') in major_cell_types_sorted else len(major_cell_types_sorted), x))
subset_cell_types_sorted = sorted(subset_cell_types, key=lambda x: (major_cell_types_sorted.index(get_major(x, level='subset')) if get_major(x, level='subset') in major_cell_types_sorted else len(major_cell_types_sorted), x))

cell_type_orders = {
    'major': major_cell_types_sorted,
    'minor': minor_cell_types_sorted,
    'subset': subset_cell_types_sorted
}


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

# Also load HEV niche barcode and add to adata
filename = f'{sample_id}_HEV_niches_cell_to_niche.csv'
niche_barcode_df = pd.read_csv(os.path.join(most_recent_distance, 'csv_files', filename), index_col=0)
niche_barcode_df.dropna(subset=['niche_id'], inplace=True)
adata.obs['niche_id'] = niche_barcode_df['niche_id']

# Also load the distance calculations
distance_file = f'{sample_id}_HEV_niches_info.csv'
distance_df = pd.read_csv(os.path.join(most_recent_distance, 'csv_files', distance_file))

# Make a mapping dict of distance to TLS
distance_dict = distance_df.set_index('niche_id')['close_to_TLS'].to_dict()

# Map to adata
adata.obs['close_to_TLS'] = adata.obs['niche_id'].map(distance_dict)


# --- Load the endothelial colors

color_dir = os.path.join(project_dir, 'resources', 'plotting_colors')
filename = 'endothelial_colors.csv'

endothelial_colors_df = pd.read_csv(os.path.join(color_dir, filename), index_col=0)
color_dict_endothelial = endothelial_colors_df['color'].to_dict()


# --- Also add the HEV colors

filename = 'TLS_HEV_colors.csv'

HEV_colors_df = pd.read_csv(os.path.join(color_dir, filename), index_col=0)

# Add to color_dict
for cell_type, row in HEV_colors_df.iterrows():
    color_dict_endothelial[cell_type] = row['color']


# --- Specify venous vs HEV-like vs TLS-HEV

# Convert to either HEV or Endothelial ACKR1
adata.obs['niche_type'] = np.where(adata.obs['niche_id'].notna(), 'HEV', 'Venous')
color_dict_endothelial['Venous'] = color_dict_endothelial['Endothelial ACKR1']

# Also have a third category for HEV close to TLS vs far from TLS
adata.obs['niche_TLS_proximity'] = np.where(adata.obs[level] == 'Endothelial ACKR1', 'Venous', None)
hev_close_mask = (adata.obs['niche_id'].notna()) & (adata.obs['close_to_TLS'].values == True)
hev_far_mask = (adata.obs['niche_id'].notna()) & (adata.obs['close_to_TLS'].values == False)
adata.obs.loc[hev_close_mask, 'niche_TLS_proximity'] = 'HEV close to TLS'
adata.obs.loc[hev_far_mask, 'niche_TLS_proximity'] = 'HEV far from TLS'

# Make ordered
plot_order_TLS = ['Venous', 'HEV far from TLS', 'HEV close to TLS']
adata.obs['niche_TLS_proximity'] = pd.Categorical(adata.obs['niche_TLS_proximity'], categories=plot_order_TLS, ordered=True)

color_dict_endothelial['HEV far from TLS'] = color_dict_endothelial['non_TLS_HEV']
color_dict_endothelial['HEV close to TLS'] = color_dict_endothelial['near_TLS_HEV']


# --- Load the STING signatures

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

# Filter genes to those present in the adata
for sig in signature_dict.keys():
    genes_in_adata = [gene for gene in signature_dict[sig] if gene in adata.var_names]
    print(f'Signature {sig}: {len(genes_in_adata)}/{len(signature_dict[sig])} genes present in the data.')
    print(f'Genes missing: {[gene for gene in signature_dict[sig] if gene not in adata.var_names]}')
    signature_dict[sig] = genes_in_adata


# --- Score the STING signatures

for sig in signature_dict.keys():
    score_name = sig
    sc.tl.score_genes(adata, gene_list=signature_dict[sig], score_name=score_name)



# --- Plot 1: Violin plots of STING signature scores per cell type

# Set cell type order
adata.obs[level] = pd.Categorical(adata.obs[level], categories=cell_type_orders[level], ordered=True)

# Set the grid behind
plt.rcParams['axes.axisbelow'] = True

ncols = 1
nrows = int(np.ceil(len(signature_dict.keys()) / ncols))

plt.figure(figsize=(5, 5))

for n, sig in enumerate(signature_dict.keys()):
    ax = plt.subplot(nrows, ncols, n + 1)
    sns.violinplot(data=adata.obs, x=level, y=sig, hue=level,
                 palette=color_dict, density_norm='width',
                 ax=ax, saturation=1, cut=0)
    ax.set_title(f'{sig}\n(n={len(signature_dict[sig])})', fontweight='bold')
    ax.set_ylabel('Score')
    ax.set_xlabel('')
    # Rotate the xticklabels
    ax.set_xticklabels(ax.get_xticklabels(), rotation=90,)
    # Replace xticklabels with short names if subset level
    if level == 'subset':
        new_labels = []
        for label in ax.get_xticklabels():
            label_text = label.get_text()
            if label_text in miniatlas_sub_dict['miniatlas_subset_short']:
                new_label = miniatlas_sub_dict['miniatlas_subset_short'][label_text]
            else:
                new_label = label_text
            new_labels.append(new_label)
        ax.set_xticklabels(new_labels, rotation=90,)

axs = plt.gcf().get_axes()

# Resize plots
for ax in axs:
    ax_pos = ax.get_position()
    ax.set_position([ax_pos.x0, ax_pos.y0, len(ax.get_xticks())*0.1, 0.8])
    for collection in ax.collections: # Adjust the colors of the violin
        # Get the RGBA values of the fill color
        fill_color = collection.get_facecolor()
        darkened_color = fill_color.copy()
        darkened_color[:,:-1] *= 0.8   # multiply RGB values by 0.8 to darken
        collection.set_edgecolor(darkened_color) 
        fill_color[0][-1] = 0.8 # make the fill color a bit transparent
        collection.set_facecolor(fill_color)

# Pad the axs list with None to make it evenly divisible by ncols
padded_axs = axs + [None] * (nrows * ncols - len(axs))

# Reshape the padded list
axs_grid = np.array(padded_axs).reshape(nrows, ncols)
axs_grid_cleaned = np.array([[ax for ax in row if ax is not None] for row in axs_grid if any(ax is not None for ax in row)])

# Move plots
plotting_functions.move_plots_x_y(axs_grid, x_offset=0.15, y_offset=0.8, remove_ylabels=False,
                                  remove_xlabels=False, remove_title=False)

# Position the title
title_x = (axs_grid_cleaned[0,0].get_position().x0 + axs_grid_cleaned[0,-1].get_position().x1) / 2
title_y = axs_grid_cleaned[0,0].get_position().y1 + 0.3

plt.suptitle(f'{sample_id} ({subtype}) - STING & type I IFN signature scores per cell type (annotation_level={level})', y=title_y, x=title_x, fontweight='bold')


# --- Plot 2: Violin plots for endothelial ACKR1+ only

ncols = 3
nrows = int(np.ceil(len(signature_dict.keys()) / ncols))

plt.figure(figsize=(5, 5))

for n, sig in enumerate(signature_dict.keys()):
    ax = plt.subplot(nrows, ncols, n + 1)
    sns.violinplot(data=adata.obs, x='niche_TLS_proximity', y=sig, hue='niche_TLS_proximity',
                 palette=color_dict_endothelial, density_norm='width',
                 ax=ax, saturation=1, cut=0, legend=False)
    # Wrap the signature name if too long
    # Try to wrap on an underscore near the middle
    wrapped_sig = sig
    if len(sig) > 20:
        mid_point = len(sig) // 2
        underscore_pos = sig.rfind('_', 0, mid_point)
        if underscore_pos == -1:
            underscore_pos = sig.find('_', mid_point)
        if underscore_pos != -1:
            wrapped_sig = sig[:underscore_pos] + '\n' + sig[underscore_pos+1:]
    ax.set_title(f'{wrapped_sig}\n(n={len(signature_dict[sig])})', fontweight='bold')
    ax.set_ylabel('Score')
    ax.set_xlabel('')
    # Rotate the xticklabels
    ax.set_xticklabels(ax.get_xticklabels(), rotation=90,)

axs = plt.gcf().get_axes()

# Resize plots
for ax in axs:
    ax_pos = ax.get_position()
    ax.set_position([ax_pos.x0, ax_pos.y0, len(ax.get_xticks())*0.2, 0.8])
    for collection in ax.collections: # Adjust the colors of the violin
        # Get the RGBA values of the fill color
        fill_color = collection.get_facecolor()
        darkened_color = fill_color.copy()
        darkened_color[:,:-1] *= 0.8   # multiply RGB values by 0.8 to darken
        collection.set_edgecolor(darkened_color) 
        fill_color[0][-1] = 0.5 # make the fill color a bit transparent
        collection.set_facecolor(fill_color)

# Pad the axs list with None to make it evenly divisible by ncols
padded_axs = axs + [None] * (nrows * ncols - len(axs))

# Reshape the padded list
axs_grid = np.array(padded_axs).reshape(nrows, ncols)
axs_grid_cleaned = np.array([[ax for ax in row if ax is not None] for row in axs_grid if any(ax is not None for ax in row)])

# Move plots
plotting_functions.move_plots_x_y(axs_grid, x_offset=0.6, y_offset=0.8, remove_ylabels=False,
                                  remove_xlabels=False, remove_title=False)

# Position the title
title_x = (axs_grid_cleaned[0,0].get_position().x0 + axs_grid_cleaned[0,-1].get_position().x1) / 2
title_y = axs_grid_cleaned[0,0].get_position().y1 + 0.3

plt.suptitle(f'{sample_id} ({subtype}) - STING & type I IFN signature scores in endothelial ACKR1+ (annotation_level={level})', y=title_y, x=title_x, fontweight='bold')


# --- Plot 3: Dotplot with pct above 0 and mean score

annotation_key = 'niche_TLS_proximity'

# Prepare data for dotplot
dotplot_data = []
for sig in signature_dict.keys():
    for cell_type in adata.obs[annotation_key].astype('category').cat.categories:
        subset_mask = adata.obs[annotation_key] == cell_type
        scores = adata.obs.loc[subset_mask, sig]
        mean_score = scores.mean()
        pct_above_0 = (scores > 0).sum() / len(scores) * 100
        dotplot_data.append({'signature': sig,
                             'cell_type': cell_type,
                             'mean_score': mean_score,
                             'pct_above_0': pct_above_0})

dotplot_df = pd.DataFrame(dotplot_data)

# Also get counts of HEV ACKR1+ per TLS proximity category
hev_counts = adata.obs[annotation_key].value_counts().reindex(plot_order_TLS).fillna(0).astype(int)

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

# Set fontsize
sc.set_figure_params(fontsize=10)

# Turn the grid off
plt.rcParams['axes.grid'] = False

# Make the dotplot
fig, top_ax = plt.subplots(figsize=(2, 2))

axs = []

for n, (signature_type, signature_list) in enumerate(signature_type_dict.items()):
    print(f'{signature_type}: {signature_list}')
    if n == 0:
        ax = top_ax
    else:
        ax_width = ax_prev.get_position().width
        ax_height = ax_prev.get_position().height / len(ax_prev.get_yticks()) * len(signature_list)
        ax = plt.gcf().add_axes([ax_prev.get_position().x0, ax_prev.get_position().y0 - ax_height - 0.02, ax_width, ax_height])
    plot_df = dotplot_df[dotplot_df['signature'].isin(signature_list)].copy()
    # Remove unused categories
    plot_df['signature'] = plot_df['signature'].cat.remove_unused_categories()
    dotplot_ax, size_max = plotting_functions.seaborn_dotplot(plot_df, ax, x='cell_type', y='signature', hue='mean_score', size='pct_above_0',
                                                        max_size_px=250, cmap='RdBu_r', vmin=vmin, vmax=vmax, size_max=size_max)
    if n == 0:
        # Resize ax based on xticks & yticks
        ax_pos = ax.get_position()
        ax.set_position([ax_pos.x0, ax_pos.y0, 0.15 * len(ax.get_xticks()), 0.12 * len(ax.get_yticks())])
    # Add y label
    ax.set_ylabel(signature_type, fontweight='bold')
    # Rotate yticklabels
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0, ha='right')
    ax_prev = ax
    axs.append(ax)
    # Adjust xlim & ylim
    ax.set_xlim(-0.5, len(ax.get_xticks()) - 0.5)
    ax.set_ylim(-0.75, len(ax.get_yticks()) - 0.25)
    # Rotate xticklabels
    ax.set_xticklabels(ax.get_xticklabels(), rotation=90, ha='center')

# Remove xticks from all but the bottom ax
for ax in axs[:-1]:
    ax.set_xticks([])
    # Also remove legend
    ax.get_legend().remove()

# Add xlabel
axs[-1].set_xlabel('Vessel type', fontweight='bold')

# Add countax at the top
count_ax = plt.gcf().add_axes([top_ax.get_position().x0, top_ax.get_position().y1 + 0.02, top_ax.get_position().width, top_ax.get_position().height * 0.3])
sns.barplot(x=hev_counts.index, y=hev_counts.values, color='lightgrey', edgecolor='black', ax=count_ax, hue=hev_counts.index, palette=color_dict_endothelial,
               saturation=1)

# Beautify
plotting_functions.beautify_plot(count_ax)

# Annotate counts
for n, p in enumerate(count_ax.patches):
    # Annotate with k if over 1000
    if p.get_height() >= 1000:
        annotate_value = f'{round(p.get_height()/1000, 1)}k'
    else:
        annotate_value = f'{int(p.get_height())}'
    count_ax.annotate(annotate_value, (p.get_x() + p.get_width() / 2., p.get_height()), ha='center', va='bottom', fontsize=8)

# Adjust ylim
count_ax.set_ylim(top=hev_counts.values.max() * 1.4)

# Set top and right spine invisible
count_ax.spines['top'].set_visible(False)
count_ax.spines['right'].set_visible(False)

# Remove xticks
count_ax.set_xticks([])
count_ax.set_xlabel('')
count_ax.set_ylabel('# cells')

# Add title
count_ax.set_title('STING signatures\n(Endothelial ACKR1+)', fontweight='bold', fontsize=10)

count_axs = []

# Also add a countax on the side for the number of genes per signature
for n, (signature_type, signature_list) in enumerate(signature_type_dict.items()):
    dotplot_ax = axs[n]
    gene_count_ax = plt.gcf().add_axes([dotplot_ax.get_position().x1 + 0.02, dotplot_ax.get_position().y0, 0.12*1.5, dotplot_ax.get_position().height])
    signature_counts_subset = signature_counts_df[signature_counts_df.index.isin(signature_list)]
    # Sort by the order in the dotplot
    ytick_texts = [tick.get_text() for tick in dotplot_ax.get_yticklabels()]
    signature_counts_subset = signature_counts_subset.reindex(ytick_texts, axis=0)
    # Plot horizontal barplot
    sns.barplot(y=signature_counts_subset.index, x=signature_counts_subset['n_genes'], color='cornflowerblue', ax=gene_count_ax, saturation=1)
    # Beautify
    plotting_functions.beautify_plot(gene_count_ax)
    # Remove yticks
    gene_count_ax.set_yticks([])
    gene_count_ax.set_ylabel('')
    # Set xlabel on the last one only
    if n == len(signature_type_dict) - 1:
        gene_count_ax.set_xlabel('# genes')
    else:
        gene_count_ax.set_xlabel('')
        # Also remove xticklabels
        gene_count_ax.set_xticklabels([])
    # Annotate bars
    annotation_offset = gene_count_ax.get_xlim()[1] * 0.05
    for n, p in enumerate(gene_count_ax.patches):
        gene_count_ax.annotate(f'{int(p.get_width())}', (p.get_width() + annotation_offset, p.get_y() + p.get_height() / 2.), ha='left', va='center', fontsize=8)
    # Adjust xlim
    gene_count_ax.set_xlim(right=signature_counts_df['n_genes'].max() * 1.4)
    # Adjust ylim
    gene_count_ax.set_ylim(dotplot_ax.get_ylim())
    # Set top and right spine invisible
    gene_count_ax.spines['top'].set_visible(False)
    gene_count_ax.spines['right'].set_visible(False)
    count_axs.append(gene_count_ax)

gene_count_ax = count_axs[-1]

# Add an ax for colorbar
cbar_ax = plt.gcf().add_axes([gene_count_ax.get_position().x1 + 0.02, gene_count_ax.get_position().y0, 0.12*0.5, 0.12*3])

# Add colorbar
norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
cbar = plt.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap='RdBu_r'), cax=cbar_ax, orientation='vertical')

# Add title
cbar_ax.set_title('Mean score', fontsize=10, ha='left', position=(0, 1.05))

# Set xlim to actual vmin & vmax
cbar.set_ticks([vmin, vmax])
cbar.set_ticklabels([f'{vmin:.2f}', f'{vmax:.2f}'])
cbar_ax.set_ylim(vmin, vmax)

# Add size legend
size_legend = plotting_functions.add_size_legend(dotplot_ax, cbar_ax, size='pct_above_0', size_title='Pct. > 0', limit_circles=True,
                                                alignment='left', bbox_to_anchor=(-0.3, 1.3), loc='lower left', labelspacing=0.8, handletextpad=0.8)

# Adjust main title position
title_x = (dotplot_ax.get_position().x0 + dotplot_ax.get_position().x1) / 2
title_y = count_ax.get_position().y1 + 0.4

plt.suptitle(f'{sample_id} ({subtype}) - STING & type I IFN signature scores in endothelial ACKR1+ (annotation_level={level})', fontweight='bold', y=title_y, x=title_x)


# --- Plot 4: Same as above but normalizing scores 0-1

normalized_dotplot_df = dotplot_df.copy()
# Normalize mean_score per signature to 0-1
for sig in signature_dict.keys():
    sig_mask = normalized_dotplot_df['signature'] == sig
    min_score = normalized_dotplot_df.loc[sig_mask, 'mean_score'].min()
    max_score = normalized_dotplot_df.loc[sig_mask, 'mean_score'].max()
    normalized_dotplot_df.loc[sig_mask, 'mean_score'] = (normalized_dotplot_df.loc[sig_mask, 'mean_score'] - min_score) / (max_score - min_score + 1e-9)

vmin = 0
vmax = 1

# Set sizemin
min_size = size_min*0.98

# Make the dotplot
fig, top_ax = plt.subplots(figsize=(2, 2))

axs = []

for n, (signature_type, signature_list) in enumerate(signature_type_dict.items()):
    print(f'{signature_type}: {signature_list}')
    if n == 0:
        ax = top_ax
    else:
        ax_width = ax_prev.get_position().width
        ax_height = ax_prev.get_position().height / len(ax_prev.get_yticks()) * len(signature_list)
        ax = plt.gcf().add_axes([ax_prev.get_position().x0, ax_prev.get_position().y0 - ax_height - 0.02, ax_width, ax_height])
    plot_df = normalized_dotplot_df[normalized_dotplot_df['signature'].isin(signature_list)].copy()
    # Remove unused categories
    plot_df['signature'] = plot_df['signature'].cat.remove_unused_categories()
    dotplot_ax, size_max = plotting_functions.seaborn_dotplot(plot_df, ax, x='cell_type', y='signature', hue='mean_score', size='pct_above_0',
                                                        max_size_px=250, cmap='Reds', vmin=vmin, vmax=vmax, size_max=size_max, size_min=min_size)
    if n == 0:
        # Resize ax based on xticks & yticks
        ax_pos = ax.get_position()
        ax.set_position([ax_pos.x0, ax_pos.y0, 0.15 * len(ax.get_xticks()), 0.12 * len(ax.get_yticks())])
    # Add y label
    ax.set_ylabel(signature_type, fontweight='bold')
    # Rotate yticklabels
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0, ha='right')
    ax_prev = ax
    axs.append(ax)
    # Adjust xlim & ylim
    ax.set_xlim(-0.5, len(ax.get_xticks()) - 0.5)
    ax.set_ylim(-0.75, len(ax.get_yticks()) - 0.25)
    # Rotate xticklabels
    ax.set_xticklabels(ax.get_xticklabels(), rotation=90, ha='center')

# Remove xticks from all but the bottom ax
for ax in axs[:-1]:
    ax.set_xticks([])
    # Also remove legend
    ax.get_legend().remove()

# Add xlabel
axs[-1].set_xlabel('Vessel type', fontweight='bold')

# Add countax at the top
count_ax = plt.gcf().add_axes([top_ax.get_position().x0, top_ax.get_position().y1 + 0.02, top_ax.get_position().width, top_ax.get_position().height * 0.3])
sns.barplot(x=hev_counts.index, y=hev_counts.values, color='lightgrey', edgecolor='black', ax=count_ax, hue=hev_counts.index, palette=color_dict_endothelial,
               saturation=1)

# Beautify
plotting_functions.beautify_plot(count_ax)

# Annotate counts
for n, p in enumerate(count_ax.patches):
    # Annotate with k if over 1000
    if p.get_height() >= 1000:
        annotate_value = f'{round(p.get_height()/1000, 1)}k'
    else:
        annotate_value = f'{int(p.get_height())}'
    count_ax.annotate(annotate_value, (p.get_x() + p.get_width() / 2., p.get_height()), ha='center', va='bottom', fontsize=8)

# Adjust ylim
count_ax.set_ylim(top=hev_counts.values.max() * 1.4)

# Set top and right spine invisible
count_ax.spines['top'].set_visible(False)
count_ax.spines['right'].set_visible(False)

# Remove xticks
count_ax.set_xticks([])
count_ax.set_xlabel('')
count_ax.set_ylabel('# cells')

# Add title
count_ax.set_title('STING signatures\n(Endothelial ACKR1+)', fontweight='bold', fontsize=10)

count_axs = []

# Also add a countax on the side for the number of genes per signature
for n, (signature_type, signature_list) in enumerate(signature_type_dict.items()):
    dotplot_ax = axs[n]
    gene_count_ax = plt.gcf().add_axes([dotplot_ax.get_position().x1 + 0.02, dotplot_ax.get_position().y0, 0.12*1.5, dotplot_ax.get_position().height])
    signature_counts_subset = signature_counts_df[signature_counts_df.index.isin(signature_list)]
    # Sort by the order in the dotplot
    ytick_texts = [tick.get_text() for tick in dotplot_ax.get_yticklabels()]
    signature_counts_subset = signature_counts_subset.reindex(ytick_texts, axis=0)
    # Plot horizontal barplot
    sns.barplot(y=signature_counts_subset.index, x=signature_counts_subset['n_genes'], color='cornflowerblue', ax=gene_count_ax, saturation=1)
    # Beautify
    plotting_functions.beautify_plot(gene_count_ax)
    # Remove yticks
    gene_count_ax.set_yticks([])
    gene_count_ax.set_ylabel('')
    # Set xlabel on the last one only
    if n == len(signature_type_dict) - 1:
        gene_count_ax.set_xlabel('# genes')
    else:
        gene_count_ax.set_xlabel('')
        # Also remove xticklabels
        gene_count_ax.set_xticklabels([])
    # Annotate bars
    annotation_offset = gene_count_ax.get_xlim()[1] * 0.05
    for n, p in enumerate(gene_count_ax.patches):
        gene_count_ax.annotate(f'{int(p.get_width())}', (p.get_width() + annotation_offset, p.get_y() + p.get_height() / 2.), ha='left', va='center', fontsize=8)
    # Adjust xlim
    gene_count_ax.set_xlim(right=signature_counts_df['n_genes'].max() * 1.4)
    # Adjust ylim
    gene_count_ax.set_ylim(dotplot_ax.get_ylim())
    # Set top and right spine invisible
    gene_count_ax.spines['top'].set_visible(False)
    gene_count_ax.spines['right'].set_visible(False)
    count_axs.append(gene_count_ax)

gene_count_ax = count_axs[-1]

# Add an ax for colorbar
cbar_ax = plt.gcf().add_axes([gene_count_ax.get_position().x1 + 0.02, gene_count_ax.get_position().y0, 0.12*0.5, 0.12*3])

# Add colorbar
norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
cbar = plt.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap='Reds'), cax=cbar_ax, orientation='vertical')

# Add title
cbar_ax.set_title('Mean score\n(scaled)', fontsize=10, ha='left', position=(0, 1.05))

# Set xlim to actual vmin & vmax
cbar.set_ticks([vmin, vmax])
cbar.set_ticklabels([f'{vmin:.2f}', f'{vmax:.2f}'])
cbar_ax.set_ylim(vmin, vmax)

# Add size legend
size_legend = plotting_functions.add_size_legend(dotplot_ax, cbar_ax, size='pct_above_0', size_title='Pct. > 0', limit_circles=True,
                                                alignment='left', bbox_to_anchor=(-0.3, 1.5), loc='lower left', labelspacing=0.8, handletextpad=0.8)

# Adjust main title position
title_x = (dotplot_ax.get_position().x0 + dotplot_ax.get_position().x1) / 2
title_y = count_ax.get_position().y1 + 0.4

plt.suptitle(f'{sample_id} ({subtype}) - STING & type I IFN signature scores in endothelial ACKR1+ (annotation_level={level})', fontweight='bold', y=title_y, x=title_x)


# --- Save the output

out_dir = os.path.join(res_dir, '05_xenium_mechanisms', '04_score_signatures', '01_score_STING_signatures')

today = date.today()
today = today.strftime("%Y%m%d")

out_dir_plots = os.path.join(out_dir, 'plots', today)
os.makedirs(out_dir_plots, exist_ok=True)

out_dir_csv = os.path.join(out_dir, 'csv_files', today)
os.makedirs(out_dir_csv, exist_ok=True)

# Save the figures
filename = f'{sample_id}_STING_signatures.pdf'
save_multi_image(os.path.join(out_dir_plots, filename))

# Also save the scores
score_filename = f'{sample_id}_STING_signature_scores.csv'
to_save = adata.obs[[level, annotation_key] + [sig for sig in signature_dict.keys()]]
to_save.to_csv(os.path.join(out_dir_csv, score_filename),)

# Also save genes used per signature
genes_filename = f'{sample_id}_STING_signatures_genes_used.csv'

# Make a dataframe
genes_data = []
for sig in signature_dict.keys():
    for gene in signature_dict[sig]:
        genes_data.append({'signature': sig,
                           'gene': gene})

genes_df = pd.DataFrame(genes_data)
genes_df.to_csv(os.path.join(out_dir_csv, genes_filename), index=False)


# --- Finish

print("Script successfully completed")