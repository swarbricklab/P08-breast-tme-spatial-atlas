#!/usr/bin/env python3

# --- Aim
# Identify robust marker genes for TLS clusters across samples:
# - Load integrated clustering results and DEGs
# - Filter genes by expression specificity (pct expressing in cluster vs rest)
# - Visualize top filtered marker genes per cluster
# - Save filtered gene lists and expression statistics for downstream analysis

# --- Define environment
# env: spatial-pilot

# --- Load packages
import argparse
import os
from datetime import datetime

import scanpy as sc
import anndata as ad
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


# Custom functions
from sophiesfunctions.save_multi_image import save_multi_image

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


# --- Load the sample sheet

if args.get("sample_sheet"):
    sample_sheet_path = args["sample_sheet"]
else:
    sample_sheet_path = os.path.join(project_dir, 'config', 'sample_sheet.csv')

sample_sheet = pd.read_csv(sample_sheet_path)

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

harmony_dir = os.path.join(project_dir, 'results', '03_TLS_niche_discovery', '02_identify_TLS_niches', '01_cluster_TLS_high')
most_recent = find_most_recent(harmony_dir)

# Load the adata containing results
use_layer = 'log_norm'
adata = ad.read_h5ad(os.path.join(most_recent, 'csv_files', 'All_samples_TLS_umap.h5ad'))


# --- Also load top genes

filename = 'All_samples_leiden_clusters_DEGS.csv'
top_genes_table = pd.read_csv(os.path.join(most_recent, 'csv_files', filename))

all_groups = list(top_genes_table['group'].unique())
# Decide if want to limit the number of genes considered
num_top = None

dfs = []

for group in all_groups:
    # Get the list of genes for this group
    genes = top_genes_table[top_genes_table['group'] == group]
    # Sort by LFC
    top_genes = genes[genes['pvals_adj'] < 0.05]
    # Sort of LFC to get a very specific GC score
    top_genes_alt = top_genes.sort_values(by='logfoldchanges', ascending=False)
    if num_top is None:
        num_top = top_genes_alt.shape[0]
    gene_list = list(top_genes_alt['genes'][:num_top])
    adata_sub = adata[adata.obs['leiden_2'] == str(group)]
    adata_rest = adata[~adata.obs['leiden_2'].isna()]
    adata_rest = adata_rest[adata_rest.obs['leiden_2'] != str(group)]
    # Convert to DataFrame once and calculate percentage expressing for all genes at once
    df_sub = adata_sub.to_df(layer=use_layer)[gene_list]
    pct_expressing_list = (df_sub > 0).mean(axis=0).tolist()
    # Same for rest
    df_rest = adata_rest.to_df(layer=use_layer)[gene_list]
    pct_expressing_rest_list = (df_rest > 0).mean(axis=0).tolist()
    tmp_df = top_genes_alt.iloc[:num_top,]
    tmp_df['pct_expressing'] = pct_expressing_list
    tmp_df['pct_expressing_rest'] = pct_expressing_rest_list
    # Append the DataFrame to the list
    dfs.append(tmp_df)

# Concatenate the DataFrames from all groups into one DataFrame
final_df = pd.concat(dfs, ignore_index=True)


# --- Filter the top genes

# Set some thresholds
min_expressed = 0.3
max_expressed_rest = 0.75

# Also note in final df which genes filtered
final_df['filter'] = np.where(final_df['pct_expressing'] < min_expressed, 'too_low', 'pass')
final_df['filter'] = np.where(final_df['pct_expressing_rest'] > max_expressed_rest, 'too_general', final_df['filter'])

# Filter genes expressed in less than 30% of spots
final_df_sub = final_df[final_df['pct_expressing'] > min_expressed]
final_df_sub = final_df_sub[final_df_sub['pct_expressing_rest'] < max_expressed_rest]


# --- Plot 1: Dot plots for top 30 genes per cluster - separated by cluster

# Prepare for saving
out_dir_plots = os.path.join(most_recent, 'top_genes_dotplot', 'post_filtering')
os.makedirs(out_dir_plots, exist_ok=True)

# Slash to be only the top x per cluster
# Keep only the first x rows for each group

num_tops = [20, 25, 30]

for num_top in num_tops:
    top_df = final_df_sub.groupby('group').head(num_top)
    # Plotting per cluster
    for cluster in all_groups:
        top_genes = top_df[top_df['group'] == cluster]['genes'].values
        # Check that there are top genes
        if len(top_genes) > 0:
            dp = sc.pl.dotplot(adata, var_names=top_genes, groupby=['leiden'], layer=use_layer, swap_axes=False, return_fig=True)
            dp.add_totals().style(dot_edge_color='black', dot_edge_lw=0.5, cmap='Reds').show()
            axs = plt.gcf().get_axes()
            # Get the legend ax & adjust the legend
            dots_ax = axs[3]
            dots_ax.set_title('Fraction of spots\nin group (%)', y=0.5, fontsize=8.5)
            # Get the top_ax to add a title
            top_ax = axs[0]
            top_ax.set_title(f'Top {num_top} genes cluster {cluster} (post filtering)', fontweight='bold', fontsize=9.5)
            plt.suptitle(f'All samples - leiden clusters TLS groups (top {num_top} genes (by LFC) for cluster {cluster} (use_data={use_layer}, min_expressed={min_expressed}, max_expressed_rest={max_expressed_rest})', y=1, fontweight='bold')
    # Save plots
    filename = f'All_samples_TLS_umap_marker_genes_top_{num_top}_LFC_filtered.pdf'
    save_multi_image(os.path.join(out_dir_plots, filename))
    plt.close('all')


# --- Save some csv files

# Prepare for saving
out_dir_csv = os.path.join(out_dir_plots, 'csv_files')
os.makedirs(out_dir_csv, exist_ok=True)

filename = 'All_clusters_pct_expressing.csv'
final_df.to_csv(os.path.join(out_dir_csv, filename), index=False)
final_df = pd.read_csv(os.path.join(out_dir_csv, filename))

# Also save filtered list per cluster
for cluster in all_groups:
    top_genes = final_df_sub[final_df_sub['group'] == cluster]
    filename = f'Cluster_{cluster}_top_genes_filtered.csv'
    top_genes.to_csv(os.path.join(out_dir_csv, filename), index=False)


# --- Also identify the clusters of interest

# Cluster numbers may shift a bit but should be pretty clear from gene expression which is which
top_marker_dict = {
    'T cell zone': ['CCL19', 'LAMP3'],
    'B cell zone': ['CXCL13', 'MS4A1'],
    'HEVs': ['CCL21', 'ACKR1'],
    'Germinal center': ['RGS13', 'MEF2B'],
}

# Extract cluster numbers of interest based on marker gene expression
clusters_of_interest = {}

# Ensure clusters are sorted numerically (as strings)
sorted_clusters = sorted(all_groups, key=lambda x: int(x) if str(x).isdigit() else x)
for region, markers in top_marker_dict.items():
    for cluster in sorted_clusters:
        # Get the top 20 genes for this cluster
        cluster_genes = final_df[final_df['group'] == cluster]['genes'].values[:20]
        if all(marker in cluster_genes for marker in markers):
            clusters_of_interest[region] = cluster
            break

print("Clusters of interest based on marker gene expression:")
for region, cluster in clusters_of_interest.items():
    print(f"{region}: Cluster {cluster}")

# Save the filtered gene lists for the clusters of interest in a separate directory for downstream analysis
signature_dir = os.path.join(project_dir, 'results', '03_TLS_niche_discovery', '03_extract_niche_signatures', 'niche_signatures')
os.makedirs(signature_dir, exist_ok=True)

for region, cluster in clusters_of_interest.items():
    region_adj = region.replace(' ', '_').lower()
    top_genes = final_df_sub[final_df_sub['group'] == cluster]
    filename = f'{region_adj}_top_genes_filtered.csv'
    top_genes.to_csv(os.path.join(signature_dir, filename), index=False)


# --- Finish

print("Script successfully completed")