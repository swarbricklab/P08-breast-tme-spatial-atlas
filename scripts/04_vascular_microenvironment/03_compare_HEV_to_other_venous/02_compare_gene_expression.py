#!/usr/bin/env python3

# --- Aim
# Compare gene expression profiles of vessels overlapping with HEVs versus non-HEV vessels.
# Analyze differences by vessel type across all samples.

# --- Define environment
# env: vis-hd

# --- Load packages
import argparse
import math
import os
from datetime import date, datetime

import anndata as ad
import matplotlib as mpl
import matplotlib.lines as mlines
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scanpy.external as sce
import scanpy.plotting.palettes as sc_palettes
from matplotlib.backends.backend_pdf import PdfPages

# Custom functions
from sophiesfunctions.auto_crop import auto_crop_scanpy
from sophiesfunctions.color_functions import generate_color_variations
from sophiesfunctions.save_multi_image import save_multi_image

sc.set_figure_params(facecolor="white", figsize=(8, 8), fontsize=20)
sc.settings.verbosity = 3

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

args = vars(parser.parse_args())

# Arguments
project_dir = args["project_dir"]
h5ad_dir = args["h5ad_dir"]
res_dir = args["res_dir"]
sample_list = args["samples"]

print(f'Samples in python: {sample_list}')


# --- Extract some info from the sample sheet

sample_dir = os.path.join(project_dir, 'config')

sample_sheet_file = 'sample_sheet.csv'

sample_sheet = pd.read_csv(os.path.join(sample_dir, sample_sheet_file))
sample_lut = sample_sheet.set_index('section_id')


# --- Also get some tissue metadata

resource_dir = os.path.join(project_dir, 'resources')
metadata_file = 'tissue_metadata.csv'

tissue_metadata = pd.read_csv(os.path.join(resource_dir, metadata_file))

# Merge the two DataFrames on the common column 'tissue_id'
merged_df = pd.merge(sample_sheet, tissue_metadata, on='tissue_id')
# Create a sample - subtype dictionary
subtype_dict = dict(zip(merged_df['section_id'], merged_df['clinical_subtype']))


# --- Load miniatlas classifications

decon_level = 'minor'

# To make it easier to visualize trends by lineage
class_dir = os.path.join(project_dir, 'resources', 'cell_types', 'miniatlas')
class_file = 'miniatlas_subset_classifications.csv'

class_df = pd.read_csv(os.path.join(class_dir, class_file))

# Get endothelial classes
endothelial = list(class_df[class_df['miniatlas_major'] == 'Endothelial'][f'miniatlas_{decon_level}'].unique())
endothelial = sorted(endothelial)


# -- Load data

use_data = 'log_norm'

resolution = 'square_016um'

# Multiple samples as a dictionary seems like an okay idea
adata_dict = {}

for sample in sample_list:
    external_id = sample_lut.loc[sample, 'external_id']
    joined_id = f'{external_id}-{sample}'
    adata_dict[sample] = ad.read_h5ad(os.path.join(h5ad_dir, resolution, use_data, f'{sample}_{use_data}.h5ad'))

# Add some info for plotting
for sample in sample_list:
    adata_dict[sample].obs['sample_id'] = sample
    adata_dict[sample].obs['subtype'] = subtype_dict[sample]
    adata_dict[sample].obs['donor_id'] = sample_lut.loc[sample, 'donor_id']


# --- Calculate crop coords for optimal presentation

coords_dict = {}

for sample in sample_list:
    coords_dict[sample] = auto_crop_scanpy(adata_dict[sample], border_size=0.1)


# --- Load the spot classifications

# Define which endothelial type want to look at
endo_of_interest = 'Endothelial ACKR1'
out_name = endo_of_interest.split(' ')[-1]

# Define threshold
threshold = 0.5
threshold_adj = str(threshold).split('.')[-1]

# Function to find the most recent results
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

TLS_group_dir = os.path.join(project_dir, 'results', 'TLS_model', 'other_signatures', 'HEV_analysis', 'per_sample', resolution, 
                             f'threshold_{threshold_adj}', 'other_vessels', 'HEV_vs_non_HEV', out_name, 'gene_expression',
                             'csv_files')


most_recent_dir_path = find_most_recent(TLS_group_dir)

# Read in the results & add to anndata
for sample in sample_list:
    TLS_file_name = f'{sample}_leiden_clusters.csv'
    TLS_results = pd.read_csv(os.path.join(most_recent_dir_path, TLS_file_name), index_col=0)
    adata_dict[sample].obs[TLS_results.columns] = TLS_results[TLS_results.columns]


# --- Get some colors for the pathology annotations

annotation_labels = {
    'epithelial': ['invasive', 'in_situ', 'benign', 'normal'],
    'immune': ['lymphocytes'],
    'stroma': ['stroma_fibrous', 'stroma_adipose'],
    'additional_features': ['necrosis', 'blood_vessel', 'LVI']
}

all_colors = sc_palettes.default_20

# Remove red to be more colorblind-friendly
del all_colors[3] 

# Extract all pathology categories
all_cats = set(val for adata in adata_dict.values() for val in adata.obs['dominant_pathology'].dropna())
all_cats = list(all_cats)

colors = sc_palettes.default_20[:len(annotation_labels)]

top_cat_colors = {}
for n, cat in enumerate(annotation_labels):
    top_cat_colors[cat] = colors[n]

# Get some color variations of the original colors (so variations of the same annotation are a similar color)
color_dict_patho = {}

# And also define a plot order
plot_order = []

for top_cat, labels in annotation_labels.items():
    base_color = top_cat_colors[top_cat]
    # Check how many of the labels are found in the object
    num_cols = sum(item in set(all_cats) for item in labels)
    # Get the matching items
    matching_items = [item for item in labels if item in set(all_cats)]
    plot_order += matching_items
    if num_cols > 4:
        step_size=0.1
    else:
        step_size=0.2
    variations = generate_color_variations(base_color, num_cols, step_size=step_size)
    for n, item in enumerate(matching_items):
        color_dict_patho[item] = variations[n]


# --- Also set colors for endothelial types

# Define some colors per endothelial type
colors_file = os.path.join(project_dir, 'resources', 'plotting_colors', 'endothelial_colors.csv')
colors_df = pd.read_csv(colors_file)

# Make a color_dict
color_dict = colors_df.set_index('endothelial_type')['color'].to_dict()

# Add a color for overlap
color_dict['overlap'] = 'indigo'
color_dict['overlap HEV'] = color_dict['overlap']


# --- Sort by subtype

# Extract all subtypes
all_subtypes = {subtype_dict[sample] for sample in sample_list if sample in subtype_dict}

# Sort according to preferred plotting order
subtype_order = ['Luminal A', 'Luminal B (ER+)', 'Luminal B (HER2+)', 'HER2+', 'TNBC']
sorted_subtypes = sorted(all_subtypes, key=lambda x: subtype_order.index(x) if x in subtype_order else float('inf'))

# Plotting order per subtype
list_all_samples = sorted(sample_list, key=lambda s: (sorted_subtypes.index(subtype_dict[s]), s))


# --- Prepare for saving

# Put in a timestamped folder to avoid overwriting older plots
qc_dir = os.path.join(res_dir, 'TLS_model')
os.makedirs(qc_dir, exist_ok=True)

today = date.today()
today = today.strftime("%Y%m%d")

out_dir = os.path.join(qc_dir, 'other_signatures', 'HEV_analysis', 'all_samples', 'square_016um', f'threshold_{threshold_adj}', 'other_vessels',
                       'HEV_vs_non_HEV', out_name, 'gene_expression')
os.makedirs(out_dir, exist_ok=True)

# Save plots
out_dir_plots = os.path.join(out_dir, 'plots', today)
os.makedirs(out_dir_plots, exist_ok=True)


# --- Plot 1a: Spatial overview of endothelial + HEV overlap

if len(adata_dict.keys()) < 5:
        ncols = len(adata_dict.keys())
else:
    ncols = 5

nrows = math.ceil((len(adata_dict.keys()))/ncols)

if nrows == 1:
    plt.figure(figsize=(8*ncols, 6*nrows))
    plt_y = 1.1
else:
    plt.figure(figsize=(8*ncols, 7*nrows))
    plt_y = 1

classification_key = 'majority_overlap'

for n, sample in enumerate(list_all_samples):
    if n == (ncols-1):
        legend_loc = 'right margin'
    else:
        legend_loc = None
    # add a new subplot iteratively
    ax = plt.subplot(nrows, ncols, n + 1)
    # filter df and plot ticker on the new subplot axis
    subtype = subtype_dict[sample]
    adata_sub = adata_dict[sample][adata_dict[sample].obs[classification_key].notna()]
    sc.pl.spatial(adata_sub, color=classification_key, bw=False, alpha_img=0.7, palette=color_dict,
                  crop_coord=coords_dict[sample], title=f'{sample} ({subtype})\n{endo_of_interest} HEV overlap', ax=ax, legend_loc=legend_loc)
    ax.set_title(ax.get_title(), fontweight='bold')

plt.suptitle(f'All samples - {endo_of_interest} HEV overlap (threshold={threshold})', y=plt_y, fontweight='bold')


# --- Plot 1b: Spatial overview pathology classifications

# Create a list of handles and labels for the legend
handles = []
labels = plot_order
for cat in labels:
    color = color_dict_patho[cat]
    # Make sure not to include non-epithelial
    handles.append(mlines.Line2D([], [], color=color, marker='o', alpha=1, lw=0, mew=0, markersize=10))

classification_key = 'dominant_pathology'

if nrows == 1:
    plt.figure(figsize=(8*ncols, 6*nrows))
    plt_y = 1.1
else:
    plt.figure(figsize=(8*ncols, 7*nrows))
    plt_y = 1

for n, sample in enumerate(list_all_samples):
    legend_loc = None
    # add a new subplot iteratively
    ax = plt.subplot(nrows, ncols, n + 1)
    # filter df and plot ticker on the new subplot axis
    subtype = subtype_dict[sample]
    adata_sub = adata_dict[sample][adata_dict[sample].obs[classification_key].notna()]
    sc.pl.spatial(adata_sub, color=classification_key, bw=False, alpha_img=0.5, palette=color_dict_patho,
                  crop_coord=coords_dict[sample], title=f'{sample} ({subtype})\n{endo_of_interest}\nDominant pathology', ax=ax, legend_loc=legend_loc)
    ax.set_title(ax.get_title(), fontweight='bold')

# Add legend to last plot on first row
plt.legend(handles=handles, labels=labels, title='Pathology classification', title_fontproperties={'weight':'bold'}, bbox_to_anchor=(1, 0.5), 
            loc='center left', ncol=1, frameon=False, alignment='left')

plt.suptitle(f'All samples - {endo_of_interest} dominant pathology (threshold={threshold})', y=plt_y, fontweight='bold')


# --- Save the spatial plots

filename = f'All_samples_{out_name}_HEV_overlap_spatial.pdf'
# save_multi_image(os.path.join(out_dir_plots, filename))

# Save with high resolution
filename_dir=os.path.join(out_dir_plots, filename)

pp = PdfPages(filename_dir)
fig_nums = plt.get_fignums()

figs = [plt.figure(n) for n in fig_nums]

for fig in figs:
    fig.savefig(pp, format='pdf', bbox_inches='tight', dpi=1000)

pp.close()

plt.close('all')


# --- Extract a single adata and subset to just the vessel spots

classification_key = 'majority_overlap'

adata_sub_dict = {}

for sample in sample_list:
    adata_sub_dict[sample] = adata_dict[sample][adata_dict[sample].obs[classification_key].notna()]

# Convert to adata list & concatenate
adata_list = list(adata_sub_dict.values())
all_adata = ad.concat(adata_list, merge='same')


# --- Get colors per subtype

color_dir = os.path.join(project_dir, 'resources', 'plotting_colors')
color_file = 'SCSubtype_colors.csv'

color_lut = pd.read_csv(os.path.join(color_dir, color_file))
color_lut = color_lut.set_index('subtype')

subtype_colors = list(color_lut.reindex(sorted_subtypes)['color'])


# --- PCA + UMAP

# seed for reproducibility
seed = 42

# Get highly variable genes
sc.pp.highly_variable_genes(all_adata, flavor="seurat", n_top_genes=2000, layer=use_data)

n_comps = 50

# create a new object
if 'scaled' not in all_adata.layers:
    # Scale on the fly
    all_adata.layers['scaled'] = sc.pp.scale(all_adata.layers['log_norm'], max_value=10, copy=True)

pca_obj = ad.AnnData(all_adata.layers['scaled'])
pca_obj.var = all_adata.var
sc.tl.pca(pca_obj, n_comps = n_comps, zero_center = True, svd_solver='arpack', random_state = seed, return_info = True, use_highly_variable = True)

all_adata.uns['pca'] = pca_obj.uns['pca']
all_adata.obsm['X_pca'] = pca_obj.obsm['X_pca']
all_adata.varm['PCs'] = pca_obj.varm['PCs']

# Adjust the PCs with harmony
sce.pp.harmony_integrate(all_adata, 'sample_id')

# Assign the new PCs to the standard slot
all_adata.obsm['X_pca'] = all_adata.obsm['X_pca_harmony']

# neighbourhood graph
n_pcs=30

sc.pp.neighbors(all_adata, n_pcs=n_pcs, random_state=seed)

# calculate umap
sc.tl.umap(all_adata, random_state=seed)

# Leiden clustering
resolution = 0.5

sc.tl.leiden(all_adata, key_added='leiden', resolution=resolution)


# --- Plot 1: Variance ratio

sc.set_figure_params(scanpy=True, fontsize=8)

sc.pl.pca_variance_ratio(all_adata, n_pcs=n_comps, log=True)
plt.suptitle(f'All samples - PCA variance (n_comps={n_comps}, n_pcs={n_pcs}) - {endo_of_interest}', y=1.1, fontweight='bold')


# --- Plot 2: Plot the genes driving each PC

sc.set_figure_params(scanpy=True, fontsize=10)

# Note: this doesn't work with the SCT data because you can't save the PCs
if use_data == 'log_norm':
    dimensions_to_check=10
    plt.rcParams['figure.figsize']=(5, 5)
    sc.pl.pca_loadings(all_adata, components = range(1, dimensions_to_check+1, 1), n_points = 10)
    plt.suptitle(f'All samples - PCA loadings (top 10) - Use data: {use_data} {endo_of_interest}', y=1, fontweight='bold')


# --- Plot 3: Plot qc metrics in UMAP

plt.rcParams["figure.figsize"] = (3, 3)

to_plot = ['sample_id', 'donor_id', 'batch', 'subtype', classification_key, 'leiden', 'vessel_size',
           'dominant_pathology']

# Make donor id type string for correct plotting
all_adata.obs['donor_id'] = all_adata.obs['donor_id'].astype('category')

# Order the subtypes
all_adata.obs['subtype'] = all_adata.obs['subtype'].astype(pd.CategoricalDtype(categories=sorted_subtypes, ordered=True))

# Also order pathology annotations
all_adata.obs['dominant_pathology'] = all_adata.obs['dominant_pathology'].astype(pd.CategoricalDtype(categories=plot_order, ordered=True))

# Also include the batch
# Create an empty dictionary to store the indices
batch_dict = {}

# Iterate over unique categories
for category in sample_lut['batch'].unique():
    key = category.replace('SciLife_FFPE_', '')
    # Get the indices for each category and convert them to a list
    indices = sample_lut.index[sample_lut['batch'] == category].tolist()
    for index in indices:
        if 'Xe1' in index:
            key = 'xenium_visium'
            if key in batch_dict:
                batch_dict[key].append(index)
            else:
                batch_dict[key] = [index]
        else:
            if key in batch_dict:
                batch_dict[key].append(index)
            else:
                batch_dict[key] = [index]

def assign_category(row):
    for category, indices in batch_dict.items():
        if row['sample_id'] in indices:
            return category
    return ''

all_adata.obs['batch'] = all_adata.obs.apply(assign_category, axis=1)

# Make plots
ncols = 3
nrows = math.ceil((len(to_plot))/ncols)

plt.figure(figsize=(4*ncols, 3*nrows))

for n, plot in enumerate(to_plot):
    # add a new subplot iteratively
    ax = plt.subplot(nrows, ncols, n + 1)
    if plot == 'subtype':
        palette = subtype_colors
    elif plot == classification_key:
        palette = color_dict
    else:
        palette = None
    if plot == 'dominant_pathology':
        palette = color_dict_patho
    if plot == 'min_dist_um':
        vmax = 200
    elif plot == 'vessel_size':
        vmax = 20
    else:
        vmax = None
    sc.pl.umap(all_adata, color=plot, ax=ax, palette=palette, vmax=vmax, size=12)
    ax.set_title(ax.get_title(), fontweight='bold')

plt.suptitle(f'All samples - UMAP qc metrics - (use_data={use_data}, leiden_resolution={resolution})', y=1.1, fontweight='bold')

plt.tight_layout()


# --- Plot 3b:  composition per cluster

leiden_groups_df = all_adata.obs[[classification_key, 'leiden']]

plot_df = leiden_groups_df.value_counts().reset_index()

# Pivot the DataFrame
pivot_df = plot_df.pivot(index=classification_key, columns='leiden', values='count').fillna(0)

# Calculate rowsums
rowsums = pd.DataFrame(pivot_df.sum(axis=1))

# Also make a proportion df
proportion_df = pivot_df.div(pivot_df.sum(axis=1), axis=0)

nrows = 1
ncols = 2

# Make a function for ease
def plot_composition(pivot_df, ax, ylabel=None, title=None, annotate=True):
    pivot_df.plot(kind='bar', stacked=True, ax=ax, grid=False, width=0.8, linewidth=0.2, edgecolor='black')
    # Annotate the bars
    if annotate:
        for index, value in enumerate(rowsums[0]):
            ax.annotate(str(int(value)), xy=(index, value), ha='center', va='bottom', fontsize=8)
    # Rotate xlabels
    xlabels = ax.get_xticklabels()
    ax.set_xticklabels(xlabels, rotation=45, ha='right')
    # Set axis labels
    ax.set_xlabel('')
    ax.set_ylabel(ylabel)
    # Set title
    ax.set_title(title, fontweight='bold')
    # Move legend
    ax.legend(loc='upper left', bbox_to_anchor=(1.0, 1.0), fancybox=False, frameon=False, handlelength=1.5)

sc.set_figure_params(fontsize=10)

fig, axs = plt.subplots(nrows, ncols, figsize=(2*ncols, 2*nrows), gridspec_kw={'hspace':0.5})

plot_composition(pivot_df, ax=axs[0], ylabel='Number of spots', title='Composition by cluster')
axs[0].set_ylim(axs[0].get_ylim()[0], axs[0].get_ylim()[1]*1.05)
plot_composition(proportion_df, ax=axs[1], ylabel='Proportion of spots', title='Composition by cluster', annotate=False)

plt.subplots_adjust(wspace=2)

plt.suptitle(f'All samples - Labels per cluster - (use_data={use_data}, leiden_resolution={resolution})', y=1.2, fontweight='bold')


# --- Extract marker genes

rank_genes_key = 'rank_genes_leiden'

all_keys = []
all_keys.append(rank_genes_key)

# Remove clusters with a single spot
single_spot_clusters = all_adata.obs['leiden'].value_counts()[all_adata.obs['leiden'].value_counts() == 1].index

single_spot_clusters = list(single_spot_clusters.unique())

all_adata.obs['leiden_2'] = all_adata.obs['leiden']
all_adata.obs['leiden_2'] = np.where(all_adata.obs['leiden_2'].isin(single_spot_clusters), None, all_adata.obs['leiden_2'])

sc.tl.rank_genes_groups(all_adata, 'leiden_2', use_raw = False, layer=use_data, method='wilcoxon', key_added=rank_genes_key,
                        pts=True)
sc.tl.dendrogram(all_adata, 'leiden_2')


# --- Plot 7: Dot plots for top 20 genes per cluster

sc.set_figure_params(scanpy=True, fontsize=12)

sc.pl.rank_genes_groups_dotplot(all_adata, n_genes=20, groupby='leiden_2', key=rank_genes_key, layer=use_data)

plt.suptitle(f'All samples - leiden clusters {endo_of_interest}, top 20 genes per cluster (use_data={use_data})', y=1.1, fontweight='bold')


# --- Plot 8: 1 vs rest plots

sc.pl.rank_genes_groups(all_adata, n_genes=25, groupby='leiden_2', key=rank_genes_key, layer=use_data)

plt.suptitle(f'All samples - leiden clusters {endo_of_interest} 1 vs rest plots, 25 genes (use_data={use_data})', y=1.1, fontweight='bold')


# --- Save plots so far

filename = f'All_samples_{out_name}_HEV_overlap_umap.pdf'
save_multi_image(os.path.join(out_dir_plots, filename))

plt.close('all')


# --- Load the HEV genes too to check they're similar between the groups

gene_dir = os.path.join(project_dir, 'resources', 'gene_lists', 'benchmarking_cohort_signatures')
filename =  'other_scores/HEV_vs_endothelial_positive_LFC_filtered.csv'

gene_table = pd.read_csv(os.path.join(gene_dir, filename))

# Extract the genes
if 'gene' in gene_table.columns:
    gene_list = list(gene_table.gene)
elif 'genes' in gene_table.columns:
    gene_list = list(gene_table.genes)
else:
    gene_list = list(gene_table.iloc[:, 0]) # Assuming the genes are in the first column of this table

# Remove any potential nan values
gene_list = [x for x in gene_list if isinstance(x, str) or not np.isnan(x)]

# Plot the top x
top = 50
to_plot = gene_list[:top]


# --- Now find marker genes for HEV classes

rank_genes_key = 'rank_genes_HEV_overlap'
all_keys.append(rank_genes_key)

sc.tl.rank_genes_groups(all_adata, classification_key, use_raw = False, layer=use_data, method='wilcoxon', key_added=rank_genes_key,
                        pts=True)


# --- Plot 9: Plot top genes

sc.pl.rank_genes_groups_dotplot(all_adata, n_genes=50, groupby=classification_key, key=rank_genes_key, layer=use_data)

plt.suptitle(f'All samples - {endo_of_interest} HEV overlap DEGs, top 50 genes per cluster (use_data={use_data})', y=1.1, fontweight='bold')


# --- Plot 9a: Plot top genes HEV score

sc.pl.dotplot(all_adata, to_plot, groupby=classification_key, layer=use_data,)

plt.suptitle(f'All samples - HEV gene signature top {top} genes - per class (use_data={use_data})', y=1.1, fontweight='bold')


# --- Plot 9b: Same as a bove but fewer genes

top = 30
to_plot = gene_list[:top]

sc.pl.dotplot(all_adata, to_plot, groupby=classification_key, layer=use_data,)

plt.suptitle(f'All samples - HEV gene signature top {top} genes - per class (use_data={use_data})', y=1.1, fontweight='bold')


# --- Plot 10: 1 vs rest plots

sc.pl.rank_genes_groups(all_adata, n_genes=25, groupby=classification_key, key=rank_genes_key, layer=use_data)

plt.suptitle(f'All samples - {endo_of_interest} HEV overlap DEGs 1 vs rest plots, 25 genes (use_data={use_data})', y=1.1, fontweight='bold')


# --- Save plots so far

filename = f'All_samples_{out_name}_HEV_overlap_DEGs.pdf'
save_multi_image(os.path.join(out_dir_plots, filename))

plt.close('all')


# --- Save output

# Save the clustering
out_dir_csv = os.path.join(out_dir, 'csv_files', today)
os.makedirs(out_dir_csv, exist_ok=True)
clusters_df = all_adata.obs[['sample_id', 'donor_id', 'subtype', 'leiden', classification_key, 'vessel_size', 'dominant_pathology']]

clusters_df.to_csv(os.path.join(out_dir_csv, f'All_samples_{out_name}_leiden_clusters.csv'))

# Also save the differential genes

for rank_genes_key in all_keys:

    dfs = []

    result = all_adata.uns[rank_genes_key]
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
        # Append the DataFrame to the list
        dfs.append(df)

    # Concatenate the DataFrames from all groups into one DataFrame
    final_df = pd.concat(dfs, ignore_index=True)

    filename = f'All_samples_{rank_genes_key}.csv'
    final_df.to_csv(os.path.join(out_dir_csv, filename), index=False)


# Save the whole object as well to have all the umap stuff
ad.AnnData.write(all_adata, filename=os.path.join(out_dir_csv, f"All_samples_HEV_umap.h5ad"))


# --- Finish

print("Script successfully completed")