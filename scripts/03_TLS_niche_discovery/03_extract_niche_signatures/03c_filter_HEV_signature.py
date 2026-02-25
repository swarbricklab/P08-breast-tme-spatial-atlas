#!/usr/bin/env python3

# --- Aim
# Filter HEV marker genes for specificity and expression in spatial transcriptomics data.
# Remove genes with low expression in HEV spots or high expression in general endothelial spots.
# Visualize filtered gene sets and save results for downstream marker selection and validation.

# --- Define environment
# env: spatial-pilot

# --- Load packages
import argparse
import scanpy as sc
import anndata as ad
import os
import sys
import pandas as pd
import matplotlib.pyplot as plt
from datetime import date
from datetime import datetime
import numpy as np

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
from helper_functions import save_multi_image

# --- Extract some info from the sample sheet

# Use sample_sheet argument if provided, else default
if args.get("sample_sheet"):
    sample_sheet_path = args["sample_sheet"]
else:
    sample_sheet_path = os.path.join(project_dir, 'config', 'sample_sheet.csv')

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
# Load the integrated anndata object

harmony_dir = os.path.join(res_dir, '03_TLS_niche_discovery', '03_extract_niche_signatures', '03_filter_HEV_signature', 'create_endothelial_object')

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

filename = 'All_samples_endothelial_object.h5ad'

most_recent = find_most_recent(harmony_dir)

# Load the adata containing results
adata = ad.read_h5ad(os.path.join(most_recent, 'csv_files', filename))

use_data = 'log_norm'
use_layer = 'log_norm'

classification_key = 'endo_sub'


# --- Also load top genes

filename = 'All_samples_HEV_vs_endothelial_positive_LFC.csv'
top_genes_table = pd.read_csv(os.path.join(most_recent, 'csv_files', filename))

# Limit number of genes if desired
num_top = None

if num_top is None:
    num_top = top_genes_table.shape[0]

group = 'HEV'

dfs = []

# Get the list of genes
gene_list = list(top_genes_table['genes'][:num_top])
adata_sub = adata[adata.obs[classification_key] == group]

adata_rest = adata[~adata.obs[classification_key].isna()]
adata_rest = adata_rest[adata_rest.obs[classification_key] != str(group)]

# Convert to DataFrame once and calculate percentage expressing for all genes at once
df_sub = adata_sub.to_df(layer=use_layer)[gene_list]
pct_expressing_list = (df_sub > 0).mean(axis=0).tolist()
# Same for rest
df_rest = adata_rest.to_df(layer=use_layer)[gene_list]
pct_expressing_rest_list = (df_rest > 0).mean(axis=0).tolist()

final_df = top_genes_table.iloc[:num_top,]
final_df['pct_expressing'] = pct_expressing_list
final_df['pct_expressing_rest'] = pct_expressing_rest_list


# --- Filter the top genes

# Set some thresholds
min_expressed = 0.15
max_expressed_rest = 0.75

# Also note in final df which genes filtered
final_df['filter'] = np.where(final_df['pct_expressing'] < min_expressed, 'too_low', 'pass')
final_df['filter'] = np.where(final_df['pct_expressing_rest'] > max_expressed_rest, 'too_general', final_df['filter'])

# Filter genes expressed in less than 30% of spots
final_df_sub = final_df[final_df['pct_expressing'] > min_expressed]
final_df_sub = final_df_sub[final_df_sub['pct_expressing_rest'] < max_expressed_rest]


# --- Plot 1: Dot plots for top 30 genes per cluster - separated by cluster

base_results = os.path.join(project_dir, 'results', '03_TLS_niche_discovery', '03_extract_niche_signatures', '03_filter_HEV_signature', 'filter_HEV_signature',)
today = date.today().strftime("%Y%m%d")
out_dir_plots = os.path.join(base_results, today, 'plots')
os.makedirs(out_dir_plots, exist_ok=True)

# Slash to be only the top x per cluster
# Keep only the first x rows for each group

num_tops = [50]

for num_top in num_tops:
    top_df = final_df_sub.head(num_top)
    # Plotting per cluster
    top_genes = top_df['genes'].to_list()
    # Check that there are top genes
    if len(top_genes) > 0:
        dp = sc.pl.dotplot(adata, var_names=top_genes, groupby=[classification_key], layer=use_layer, swap_axes=False, return_fig=True)
        dp.add_totals().style(dot_edge_color='black', dot_edge_lw=0.5, cmap='Reds').show()
        axs = plt.gcf().get_axes()
        # Get the legend ax & adjust the legend
        dots_ax = axs[3]
        dots_ax.set_title('Fraction of spots\nin group (%)', y=0.5, fontsize=8.5)
        # Get the top_ax to add a title
        top_ax = axs[0]
        top_ax.set_title(f'Top {num_top} genes HEVs (post filtering)', fontweight='bold', fontsize=9.5)
        plt.suptitle(f'All samples - top genes HEV after filtering based on expression (top {num_top} genes (by LFC) (use_data={use_data}, min_expressed={min_expressed}, max_expressed_rest={max_expressed_rest}))', y=1, fontweight='bold')
    # Save plots
    filename = f'HEV_signature_top_{num_top}_LFC_filtered.pdf'
    save_multi_image(os.path.join(out_dir_plots, filename))
    plt.close('all')


# --- Save some csv files

out_dir_csv = os.path.join(base_results, today, 'csv_files')
os.makedirs(out_dir_csv, exist_ok=True)

filename = 'hevs_filtering_documentation.csv'
final_df.to_csv(os.path.join(out_dir_csv, filename), index=False)
final_df = pd.read_csv(os.path.join(out_dir_csv, filename))

# Also save filtered list
top_genes = final_df_sub
filename = 'hevs_filtered_endothelial_genes.csv'
top_genes.to_csv(os.path.join(out_dir_csv, filename), index=False)


# --- Finish

print("Script successfully completed")
