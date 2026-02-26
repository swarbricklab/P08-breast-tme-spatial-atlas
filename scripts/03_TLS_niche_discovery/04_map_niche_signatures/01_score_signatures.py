#!/usr/bin/env python3

# --- Aim
# Score spatial transcriptomics spots for germinal center and other immune niche gene signatures.
# Test multiple gene cutoffs.
# Save scoring results and gene usage for downstream analysis.

# --- Define environment
# env: spatial-pilot

# --- Load packages
import argparse
import scanpy as sc
import anndata as ad
import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from datetime import date
from datetime import datetime

from matplotlib.colors import LinearSegmentedColormap
from mycolorpy import colorlist as mcp

sc.set_figure_params(facecolor="white", figsize=(8, 8))
sc.settings.verbosity = 3

# --- Argparse arguments

parser = argparse.ArgumentParser(description="Import Visium data & do basic QC",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

# Required


parser.add_argument("--res_dir", help="Path to results directory")
parser.add_argument("--h5ad_dir", help="Path to directory for .h5ad output files")
parser.add_argument("--project_dir", help="Project directory (root)")
parser.add_argument("--sample_id", help="Sample ID")
parser.add_argument("--sample_sheet", help="Path to sample sheet CSV file", default=None)


args = vars(parser.parse_args())


# --- Arguments

res_dir = args["res_dir"]
h5ad_dir=args["h5ad_dir"]
sample_id=args["sample_id"]
project_dir=args["project_dir"]

sys.path.append(os.path.join(project_dir, 'scripts'))
from helper_functions import auto_crop_scanpy, save_multi_image


# --- Get some info from the sample sheet
# Use sample_sheet argument if provided, else default
if args.get("sample_sheet"):
    sample_sheet_path = args["sample_sheet"]
else:
    sample_dir = os.path.join(project_dir, 'config')
    sample_sheet_file = 'sample_sheet.csv'
    sample_sheet_path = os.path.join(sample_dir, sample_sheet_file)

sample_sheet = pd.read_csv(sample_sheet_path)
tissue_id = sample_sheet[sample_sheet['section_id'] == sample_id]['tissue_id'].values[0]
slide_size = sample_sheet[sample_sheet['section_id'] == sample_id]['slide_size'].values[0]

# Check if tissue_type is Lymph node
tissue_type = sample_sheet.loc[sample_sheet['section_id'] == sample_id, 'tissue_type'].values[0]
run_lymph_node = False
if tissue_type == 'Lymph node':
    run_lymph_node = True


# --- Also get some tissue metadata

resource_dir = os.path.join(project_dir, 'resources')
metadata_file = 'tissue_metadata.csv'

tissue_metadata = pd.read_csv(os.path.join(resource_dir, metadata_file))

# Extract the subtype
if run_lymph_node:
    subtype = 'lymph_node'
else:
    subtype = tissue_metadata[tissue_metadata['tissue_id'] == str(tissue_id)]['clinical_subtype'].values[0]


# --- Load the adata object

adata = ad.read_h5ad(os.path.join(h5ad_dir, 'visium', f'{sample_id}.h5ad'))

# Specify which layer to use
use_layer = 'log_norm'


# --- Calculate crop coords for optimal presentation
coords = auto_crop_scanpy(adata, border_size=0.1)


# --- Set spot size for large & small capture area
# Plotting the spots larger to make it easier to see

if slide_size == '11 mm':
    spot_size = 1.8
else:
    spot_size = 1.3


# --- Load the TLS model results

TLS_res_dir = os.path.join(res_dir, '03_TLS_niche_discovery', '01_identify_TLS', '01_run_TLS_model')

if run_lymph_node:
    TLS_res_dir = os.path.join(TLS_res_dir, 'lymph_nodes')

# Find most recent results
subdirs = [d for d in os.listdir(TLS_res_dir) if os.path.isdir(os.path.join(TLS_res_dir, d))]

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
    most_recent_dir_path = os.path.join(TLS_res_dir, most_recent_dir)
    print("Most recent results:", most_recent_dir_path)
else:
    print("No valid directories found.")

# Read in the results & add to anndata

TLS_file_name = f'{sample_id}_TLS_model.csv'
TLS_results = pd.read_csv(os.path.join(most_recent_dir_path, TLS_file_name), index_col=0)
adata.obs['TLS_model'] = TLS_results


# --- Also load the group results

threshold = 0.01
print(f'Using threshold {threshold}')

threshold_adj = str(threshold).split('.')[-1]

TLS_group_dir = os.path.join(res_dir, '03_TLS_niche_discovery', '01_identify_TLS', '02_group_TLS_high', f'threshold_{threshold_adj}', 'csv_files')

if run_lymph_node:
    TLS_group_dir = os.path.join(res_dir, '03_TLS_niche_discovery', '01_identify_TLS', '02_group_TLS_high', f'threshold_{threshold_adj}', 'lymph_nodes', 'csv_files')

# Find most recent results
subdirs = [d for d in os.listdir(TLS_group_dir) if os.path.isdir(os.path.join(TLS_group_dir, d))]

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
    most_recent_dir_path = os.path.join(TLS_group_dir, most_recent_dir)
    print("Most recent results:", most_recent_dir_path)
else:
    print("No valid directories found.")

# Read in the results & add to anndata
TLS_file_name = f'{sample_id}_TLS_groups.csv'
TLS_results = pd.read_csv(os.path.join(most_recent_dir_path, TLS_file_name), index_col=0)
adata.obs[TLS_results.columns] = TLS_results


# --- Load the genes from a file

# Choose between resource signatures or newly generated results based on argument
use_resource_signatures = False

if use_resource_signatures:
    print("Using niche signatures from resources directory")
    # Load from resources
    niche_sig_file = os.path.join(resource_dir, 'gene_lists', 'TLS_niche_signatures', 'niche_signatures.csv')
    niche_signatures_df = pd.read_csv(niche_sig_file)
    
    # Define mapping between score names and niche names in CSV
    niche_mapping = {
        'GC': 'Germinal center',
        'B_cell_zone': 'B cell zone',
        'T_cell_zone': 'T cell zone',
    }
    
    # Create gene lists dict from the niche signatures file
    gene_lists_dict = {}
    for score_name, niche_name in niche_mapping.items():
        niche_genes = niche_signatures_df[niche_signatures_df['niche'] == niche_name].copy()
        # Sort by rank to get genes in order
        niche_genes = niche_genes.sort_values('rank')
        gene_lists_dict[score_name] = niche_genes['gene_symbol'].tolist()
    
else:
    print("Using niche signatures from generated results")
    # Load from results directory (original behavior)
    gene_dir = os.path.join(res_dir, '03_TLS_niche_discovery', '03_extract_niche_signatures', '04_final_niche_signature_files')
    
    clusters_of_i = {
        'GC': 'germinal_center_filtered_epithelial_genes.csv',
        'B_cell_zone': 'b_cell_zone_top_genes_filtered.csv',
        'T_cell_zone': 't_cell_zone_top_genes_filtered.csv',
    }
    
    gene_lists_dict = {}
    for score_name, filename in clusters_of_i.items():
        gene_table = pd.read_csv(os.path.join(gene_dir, filename))
        
        # Extract the genes
        if 'gene' in gene_table.columns:
            gene_list = list(gene_table.gene)
        elif 'genes' in gene_table.columns:
            gene_list = list(gene_table.genes)
        else:
            gene_list = list(gene_table.iloc[:, 0])  # Assuming the genes are in the first column of this table
        
        # Remove any potential nan values
        gene_list = [x for x in gene_list if isinstance(x, str) or not np.isnan(x)]
        gene_lists_dict[score_name] = gene_list


# --- Assess multiple scores

for score_name in gene_lists_dict.keys():

    gene_list = gene_lists_dict[score_name]
    
    print(f'\nProcessing {score_name} signature with {len(gene_list)} genes')

    # Select only the genes that are present in the anndata objects
    to_plot = np.intersect1d(gene_list , adata.var.index).tolist()

    missing = list(set(gene_list) - set(to_plot))

    print(f'Genes not found: {missing}') 

    gene_dict = {f'score_genes': gene_list}


    # --- Test different cutoffs

    # How many genes is optimal?
    cutoffs = [20, 50, 100]

    all_genes = [value for sublist in gene_dict.values() for value in sublist]

    # Set seed
    seed = 42

    # Extract the layer you want to score on

    score_obj = ad.AnnData(adata.X)
    score_obj.obs_names = adata.obs_names
    score_obj.var_names = adata.var_names

    for cutoff in cutoffs:
        # Get top x genes
        to_score = list(all_genes[0:cutoff])
        print(f'Top {cutoff} genes: {to_score}')
        # Save per cutoff
        score_key = f'{score_name}_score_{cutoff}'
        sc.tl.score_genes(score_obj, to_score, random_state=seed, score_name=score_key)
        # Put the score back in the original adata
        adata.obs[score_key] = score_obj.obs[score_key]


    # --- Prepare for plotting

    # Subset anndata to TLS score above threshold
    adata_sub = adata[adata.obs['TLS_model'] > threshold]

    adata.obs['group_singleton_removed'] = adata_sub.obs['group_singleton_removed']

    all_groups = adata_sub.obs['group_singleton_removed'].astype('category').unique()
    all_groups = all_groups.categories

    # Drop 'singleton' from the list
    all_groups = [x for x in all_groups if x != 'singleton']

    # Sort the groups by number
    sorted_groups = sorted(all_groups, key=lambda x: int(x[5:]))  # Sort by the number part of the group name
    print(sorted_groups)

    adata_sub.obs['group_singleton_removed'] = adata_sub.obs['group_singleton_removed'].astype('category').cat.reorder_categories(sorted_groups, ordered=True)


    # --- Plot 1: TLS score & groups

    # Define custom colormap for plotting TLS model results
    # Set everything under certain threshold as the lowest color of a color map

    cmap = "YlOrRd"
    color_list = mcp.gen_color(cmap=cmap, n=11)
    color_list_extracted = ['lightgray'] + color_list[2:]

    n_bins = 100
    cmap_name = 'my_list'
    cmap_adj = LinearSegmentedColormap.from_list(cmap_name, color_list_extracted, N=n_bins)

    cmap_adj.set_under(color='lightgray', alpha=0.4)

    # Plotting
    sc.set_figure_params(fontsize=16)

    nrows = 1
    ncols = 2

    fig, axs = plt.subplots(nrows, ncols, figsize=(8*ncols, 6.5*nrows))

    # TLS model results
    sc.pl.spatial(adata, color='TLS_model', size=spot_size, cmap=cmap_adj, bw=False, alpha_img=1, crop_coord=coords, title=f'TLS score', vmin=0, ax=axs[0])

    # Spots above threshold
    sc.pl.spatial(adata_sub, color='group_singleton_removed', size=spot_size, bw=False, alpha_img=1, crop_coord=coords, title=f'Spots >{threshold}', vmin=0, 
                ax=axs[1])

    plt.suptitle(f'{sample_id} ({subtype}) - TLS score + TLS-high spots grouped on spatial continuity', y=1.05, fontweight='bold', x=0.5)


    # --- Plot 2: Plot the score

    # Create a single plot with a subplot for each cutoff
    nrows = 1
    ncols = len(cutoffs)

    fig, axs = plt.subplots(nrows, ncols, figsize=(8*ncols, 6.5*nrows))
    if ncols == 1:
        axs = [axs]

    for i, cutoff in enumerate(cutoffs):
        score_key = f'{score_name}_score_{cutoff}'
        sc.pl.spatial(
            adata,
            color=score_key,
            size=spot_size,
            cmap=cmap_adj,
            bw=False,
            alpha_img=1,
            crop_coord=coords,
            title=f'{score_key}',
            vmin=0,
            ax=axs[i]
        )


    # --- Save output

    out_name = f'{score_name}_score'

    # Put in a timestamped folder to avoid overwriting older plots
    qc_dir = os.path.join(res_dir, '03_TLS_niche_discovery', '04_map_niche_signatures', '01_score_signatures')

    if run_lymph_node:
        qc_dir = os.path.join(qc_dir, 'lymph_nodes')

    os.makedirs(qc_dir, exist_ok=True)

    filename = f"{sample_id}_{out_name}.pdf"
    filename = filename.replace('/', '_')  # Replace / with underscores for file naming

    today = date.today()
    today = today.strftime("%Y%m%d")

    out_dir = os.path.join(qc_dir, out_name, today, 'plots')
    os.makedirs(out_dir, exist_ok=True)

    print(f'Saving plot as {filename} in {out_dir}')
    save_multi_image(os.path.join(out_dir, filename))
    plt.close('all')

    # Also save the scores + which genes were used
    out_dir = os.path.join(qc_dir, out_name, today, 'csv_files')
    os.makedirs(out_dir, exist_ok=True)

    score_keys = [f'{score_name}_score_' + str(item) for item in cutoffs]
    save_score = adata.obs[score_keys]
    filename = f'{sample_id}_{out_name}.csv'
    filename = filename.replace('/', '_')  # Replace / with underscores for file naming

    save_score.to_csv(os.path.join(out_dir, filename))

    # Save genes used
    # Also save the pandas dataframe as a csv file (to be able to check which genes were used in each sample)
    genes_used = pd.DataFrame.from_dict(gene_dict, orient='index').transpose()
    genes_used.to_csv(os.path.join(out_dir, f'{sample_id}_genes_used.csv'), index=False)


print("Script succesfully completed")
