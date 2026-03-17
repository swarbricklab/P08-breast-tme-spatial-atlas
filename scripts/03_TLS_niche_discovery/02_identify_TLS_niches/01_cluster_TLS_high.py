#!/usr/bin/env python3

# --- Aim
# Integrate spatial transcriptomics samples to identify TLS clusters:
# - Load and merge sample/tissue metadata
# - Perform PCA, Harmony integration, UMAP, Leiden clustering across samples
# - Visualize QC metrics, subtype/sample/donor distributions
# - Extract and plot marker genes for TLS clusters
# - Save clustering results, marker gene lists, and integrated AnnData for downstream analysis

# --- Define environment
# env: spatial-pilot

# --- Load packages
import argparse
import scanpy as sc
import scanpy.external as sce
import anndata as ad
import os
import sys
import pandas as pd
import matplotlib.pyplot as plt
from datetime import date
from datetime import datetime
import numpy as np
import matplotlib.patches as patches

import math
from matplotlib.colors import LinearSegmentedColormap
from mycolorpy import colorlist as mcp
import scanpy.plotting.palettes as sc_palettes

sc.set_figure_params(facecolor="white", figsize=(8, 8), fontsize=10)
sc.settings.verbosity = 3


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

print(f'Samples in python: {sample_list}')

sys.path.append(os.path.join(project_dir, 'scripts'))
from helper_functions import auto_crop_scanpy, save_multi_image


# --- Extract some info from the sample sheet
# Use sample_sheet argument if provided, else default
if args.get("sample_sheet"):
    sample_sheet_path = args["sample_sheet"]
else:
    sample_sheet_path = os.path.join(project_dir, 'config', 'sample_sheet_visium.csv')

sample_sheet = pd.read_csv(sample_sheet_path)
sample_lut = sample_sheet.set_index('section_id')


# --- Also get some tissue metadata

resource_dir = os.path.join(project_dir, 'resources')
metadata_file = 'tissue_metadata.csv'

tissue_metadata = pd.read_csv(os.path.join(resource_dir, metadata_file))

# Merge the two DataFrames on the common column 'tissue_id'
merged_df = pd.merge(sample_sheet, tissue_metadata, on='tissue_id')

# Create a sample - subtype dictionary
subtype_dict = dict(zip(merged_df['section_id'], merged_df['clinical_subtype']))

# And make a sample patient LUT
sample_donor_lut = {}

for sample in sample_list:
    sample_donor_lut[sample] = sample_lut.loc[sample, 'donor_id']


# -- Load data
# Multiple samples as a dictionary seems like an okay idea
adata_dict = {}

use_data = 'log_norm'
use_layer = 'log_norm'

for sample in sample_list:
    adata_dict[sample] = ad.read_h5ad(os.path.join(h5ad_dir, 'visium', f'{sample}.h5ad'))

# Add some info for plotting

for sample in sample_list:
    adata_dict[sample].obs['sample_id'] = sample
    adata_dict[sample].obs['subtype'] = subtype_dict[sample]
    adata_dict[sample].obs['donor_id'] = sample_lut.loc[sample, 'donor_id']


# --- Calculate crop coords for optimal presentation

coords = {}
for sample in adata_dict.keys():
    coords[sample] = auto_crop_scanpy(adata_dict[sample], border_size=0.1)


# --- Order based on sample ID

ordered_keys = sorted(adata_dict.keys())


# --- Set spot size for large & small capture area
# Plotting the spots larger to make it easier to see

spot_size={}

for sample in adata_dict.keys():
    if sample_lut.loc[sample, 'slide_size'] == '11 mm':
        spot_size[sample]=1.8
    else:
        spot_size[sample]=1.6


# --- Load the TLS model results for each sample

TLS_res_dir = os.path.join(res_dir, '03_TLS_niche_discovery', '01_identify_TLS', '01_run_TLS_model')

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
    else:
        print("No valid directories found.")
    return most_recent_dir_path

most_recent_dir_path = find_most_recent(TLS_res_dir)

# Read in the results & add to anndata
for sample in adata_dict.keys():
    TLS_file_name = f'{sample}_TLS_model.csv'
    TLS_results = pd.read_csv(os.path.join(most_recent_dir_path, TLS_file_name), index_col=0)
    adata_dict[sample].obs['TLS_model'] = TLS_results


# --- Make a custom color map for better visualization

# Set everything under 0 as the lowest color of a color map

cmap = "YlOrRd"

color_list = mcp.gen_color(cmap=cmap, n=11)

color_list_extracted = ['lightgray'] + color_list[2:]

n_bins = 100
cmap_name = 'my_list'
cmap_adj = LinearSegmentedColormap.from_list(cmap_name, color_list_extracted, N=n_bins)

cmap_adj.set_under(color='lightgray', alpha=0.4)


# --- Also load the TLS groupings

threshold = 0.01
threshold_adj = str(threshold).split('.')[-1]

TLS_group_dir = os.path.join(res_dir, '03_TLS_niche_discovery', '01_identify_TLS', '02_group_TLS_high', 
                                f'threshold_{threshold_adj}', 'csv_files')

most_recent_dir_path = find_most_recent(TLS_group_dir)

# Read in the results & add to anndata
group_key = 'group_singleton_removed'

for sample in sample_list:
    TLS_file_name = f'{sample}_TLS_groups.csv'
    TLS_results = pd.read_csv(os.path.join(most_recent_dir_path, TLS_file_name), index_col=0)
    adata_dict[sample].obs[group_key] = TLS_results[group_key]
    # Set singletons to NA to exclude them
    # adata_dict[sample].obs[group_key] = np.where(adata_dict[sample].obs[group_key] == 'singleton', None, adata_dict[sample].obs[group_key])

# Subset adata to include only the TLS  groups
for sample, adata in adata_dict.items():
    adata_dict[sample] = adata[~adata.obs[group_key].isna()]


# --- Prepare for plotting the result per subtype: define a plotting order based on subtype

# Convert to adata list & concatenate
adata_list = list(adata_dict.values())
all_adata = ad.concat(adata_list, merge='same')

# all_adata.obs_names_make_unique()

# Extract all subtypes
all_subtypes = all_adata.obs.subtype.unique()

# Sort according to preferred plotting order
subtype_order = ['Luminal A', 'Luminal B', 'Luminal B (HER2+)', 'HER2+', 'TNBC']
sorted_subtypes = sorted(all_subtypes, key=lambda x: subtype_order.index(x) if x in subtype_order else float('inf'))


# -- Define a plotting order

# Extract all sample ids per subtype
adata_per_subtype_list = {}

for subtype in sorted_subtypes:
    adata_per_subtype_list[subtype] = []
    for sample in adata_dict.keys():
        if subtype_dict[sample] == subtype:
            adata_per_subtype_list[subtype] += [sample]

# Sort each list of sample IDs based
for list_samples in adata_per_subtype_list.values():
    list_samples.sort()


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
sc.pp.highly_variable_genes(all_adata, flavor="seurat", n_top_genes=2000, layer=use_layer)

n_comps = 50

# create a new object
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
resolution = 1.5

sc.tl.leiden(all_adata, key_added='leiden', resolution=resolution)


# --- Plot 1: Variance ratio

sc.set_figure_params(scanpy=True, fontsize=8)

sc.pl.pca_variance_ratio(all_adata, n_pcs=n_comps, log=True)
plt.suptitle(f'All samples - PCA variance (n_comps={n_comps}, n_pcs={n_pcs}) - TLS groups threshold={threshold}', y=1.1, fontweight='bold')


# --- Plot 2: Genes driving each PC

sc.set_figure_params(scanpy=True, fontsize=10)

# Note: this doesn't work with the SCT data because you can't save the PCs
if use_data == 'log_norm':
    dimensions_to_check=10
    plt.rcParams['figure.figsize']=(5, 5)
    sc.pl.pca_loadings(all_adata, components = range(1, dimensions_to_check+1, 1), n_points = 10)
    plt.suptitle(f'All samples - PCA loadings (top 10) - Use data: {use_data} TLS groups threshold={threshold}', y=1, fontweight='bold')


# --- Plot 3: QC metrics in UMAP

plt.rcParams["figure.figsize"] = (4, 4)

to_plot = ['sample_id', 'donor_id', 'subtype', 'leiden']

# Make donor id type string for correct plotting
all_adata.obs['donor_id'] = all_adata.obs['donor_id'].astype('category')

# Order the subtypes
all_adata.obs['subtype'] = all_adata.obs['subtype'].astype(pd.CategoricalDtype(categories=sorted_subtypes, ordered=True))

# Make plots
ncols = 3
nrows = math.ceil((len(to_plot))/ncols)

plt.figure(figsize=(7*ncols, 4*nrows))

for n, plot in enumerate(to_plot):
    # add a new subplot iteratively
    ax = plt.subplot(nrows, ncols, n + 1)
    if plot == 'subtype':
        palette = subtype_colors
    else:
        palette = None
    sc.pl.umap(all_adata, color=plot, ax=ax, palette=palette)
    ax.set_title(ax.get_title(), fontweight='bold')

plt.suptitle(f'All samples - UMAP qc metrics - (use_data={use_data}, resolution={resolution})', y=1.1, fontweight='bold')

plt.tight_layout()


# --- Plot 4a: UMAP per subtype

# Function to avoid repetition
def plot_separate_umaps(adata, group_key, group_list, color_list, ncols=3, add_subtype=False, subtype_lut=None):
    nrows = math.ceil((len(group_list) + 2)/ncols)
    plt.figure(figsize=(3*ncols, 3*nrows))
    # Then plot the groups
    for n, group in enumerate(group_list):
        # add a new subplot iteratively
        ax = plt.subplot(nrows, ncols, n + 1)
        if add_subtype:
            subtype = adata[adata.obs[group_key] == group].obs['subtype'].unique()[0]
            sc.pl.umap(adata, color=group_key, groups=[group], ax=ax, palette=[color_list[n]], na_color='lightgray', title=f'{group} ({subtype})')
        else:
            # filter df and plot ticker on the new subplot axis
            sc.pl.umap(adata, color=group_key, groups=[group], ax=ax, palette=[color_list[n]], na_color='lightgray', title=group)
        ax.set_title(ax.get_title(), fontweight='bold')
        ax.get_legend().remove()


# Plot separate umaps for subtype
group_key = 'subtype'

plot_separate_umaps(all_adata, group_key, sorted_subtypes, color_list=subtype_colors)

plt.suptitle(f'All samples - UMAP per subtype TLS groups (threshold={threshold})', y=1, fontweight='bold')

plt.tight_layout()


# --- Plot 4b: UMAP per sample

group_key = 'sample_id'

sorted_samples = adata_per_subtype_list.values()
sorted_samples = [item for sublist in sorted_samples for item in sublist]

# Check that there aren't any samples that didn't make it because they don't have aggregates
samples_w_aggregates = list(all_adata.obs['sample_id'].cat.categories)

cleaned_list = [item for item in sorted_samples if item in samples_w_aggregates]

group_list = cleaned_list

# Extract colors
if len(group_list) < 21:
    color_list = sc_palettes.default_20[:len(group_list)]
else:
    color_list = sc_palettes.godsnot_102[:len(group_list)]

plot_separate_umaps(all_adata, group_key, group_list, color_list=color_list, add_subtype=True, ncols=5)

plt.suptitle(f'All samples - UMAP per sample TLS groups (threshold={threshold})', y=1, fontweight='bold')

plt.tight_layout()


# --- Plot 4c: UMAP per donor

group_key = 'donor_id'

# Create a list of donors corresponding to the samples
donor_list = [sample_donor_lut[sample] for sample in cleaned_list]

# Remove duplicates while preserving order
unique_list = []
seen = set()
for item in donor_list:
    if item not in seen:
        unique_list.append(item)
        seen.add(item)

group_list = donor_list

# Extract colors
if len(group_list) < 21:
    color_list = sc_palettes.default_20[:len(group_list)]
else:
    color_list = sc_palettes.godsnot_102[:len(group_list)]

plot_separate_umaps(all_adata, group_key, group_list, color_list=color_list, add_subtype=True, ncols=5)

plt.suptitle(f'All samples - UMAP per patient TLS groups (threshold={threshold})', y=1, fontweight='bold')

plt.tight_layout()


# --- Extract marker genes

rank_genes_key = 'rank_genes_leiden'

# Remove clusters with a single spot
single_spot_clusters = all_adata.obs['leiden'].value_counts()[all_adata.obs['leiden'].value_counts() == 1].index

single_spot_clusters = list(single_spot_clusters.unique())

all_adata.obs['leiden_2'] = all_adata.obs['leiden']
all_adata.obs['leiden_2'] = np.where(all_adata.obs['leiden_2'].isin(single_spot_clusters), None, all_adata.obs['leiden_2'])

sc.tl.rank_genes_groups(all_adata, 'leiden_2', use_raw = False, layer=use_layer, method='wilcoxon', key_added=rank_genes_key,
                        pts=True)
sc.tl.dendrogram(all_adata, 'leiden_2')


# --- Plot 5: Dot plots for top 10 genes per cluster

sc.pl.rank_genes_groups_dotplot(all_adata, n_genes=10, groupby='leiden_2', key=rank_genes_key, layer=use_layer)

plt.suptitle(f'All samples - leiden clusters TLS groups (threshold={threshold}), top 10 genes per cluster (use_data={use_data})', y=1.1, fontweight='bold')


# --- Plot 6: 1 vs rest plots

sc.pl.rank_genes_groups(all_adata, n_genes=25, groupby='leiden_2', key=rank_genes_key, layer=use_layer)

plt.suptitle(f'All samples - leiden clusters TLS groups (threshold={threshold}) 1 vs rest plots, 25 genes (use_data={use_data})', y=1.1, fontweight='bold')


# --- Save output

# Put in a timestamped folder to avoid overwriting older plots
qc_dir = os.path.join(res_dir, '03_TLS_niche_discovery', '02_identify_TLS_niches', '01_cluster_TLS_high')
os.makedirs(qc_dir, exist_ok=True)

today = date.today()
today = today.strftime("%Y%m%d")

out_dir = os.path.join(qc_dir, today)
os.makedirs(out_dir, exist_ok=True)

# Save plots
out_dir_plots = os.path.join(out_dir, 'plots')
os.makedirs(out_dir_plots, exist_ok=True)

filename = f'All_samples_TLS_groups_umap.pdf'

save_multi_image(os.path.join(out_dir_plots, filename))

# Save the clustering
out_dir_csv = os.path.join(out_dir, 'csv_files')
os.makedirs(out_dir_csv, exist_ok=True)
clusters_df = all_adata.obs[['sample_id', 'donor_id', 'subtype', 'leiden']]

clusters_df.to_csv(os.path.join(out_dir_csv, f'All_samples_leiden_clusters.csv'))

# Also save the differential genes
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

filename = f'All_samples_leiden_clusters_DEGS.csv'
final_df.to_csv(os.path.join(out_dir_csv, filename), index=False)


# Save the whole object as well to have all the umap stuff
ad.AnnData.write(all_adata, filename=os.path.join(out_dir_csv, f"All_samples_TLS_umap.h5ad"))


# --- Finish

print("Script successfully completed")