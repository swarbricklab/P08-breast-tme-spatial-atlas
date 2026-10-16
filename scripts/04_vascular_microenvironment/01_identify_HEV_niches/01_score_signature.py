#!/usr/bin/env python3

# --- Aim
# Score Visium HD bins for HEV and TLS zone signatures derived from the Visium cohort.
# Assess the performance of these signatures at different spatial resolutions.

# --- Define environment
# env: vis-hd

# --- Load packages
import argparse
import anndata as ad
import os
import pandas as pd
import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
import squidpy as sq
import math
import sys


# --- Argparse arguments

parser = argparse.ArgumentParser(description="Score Visium HD bins for the HEV and TLS zone signatures",
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
from helper_functions import auto_crop, save_multi_image


# --- Get some info from the sample sheet
# Use sample_sheet argument if provided, else default
if args.get("sample_sheet"):
    sample_sheet_path = args["sample_sheet"]
else:
    sample_dir = os.path.join(project_dir, 'config')
    sample_sheet_file = 'sample_sheet_visium_hd.csv'
    sample_sheet_path = os.path.join(sample_dir, sample_sheet_file)

sample_sheet = pd.read_csv(sample_sheet_path)
tissue_id = sample_sheet[sample_sheet['section_id'] == sample_id]['tissue_id'].values[0]


# --- Also get some tissue metadata

resource_dir = os.path.join(project_dir, 'resources')
metadata_file = 'tissue_metadata.csv'

tissue_metadata = pd.read_csv(os.path.join(resource_dir, metadata_file))

subtype = tissue_metadata[tissue_metadata['tissue_id'] == tissue_id]['clinical_subtype'].values[0]


# --- Load the adata object

# Decide which resolutions to run on
resolutions = ['square_048um', 'square_016um', 'square_008um']

use_layer = 'log_norm'

sc.set_figure_params(fontsize=15)

# Choose between resource signatures or newly generated results
# Set to True to use curated signatures from resources directory
use_resource_signatures = True   # published signatures (resources/gene_lists/TLS_niche_signatures)

for resolution in resolutions:

    # --- Load anndata

    adata = ad.read_h5ad(os.path.join(h5ad_dir, 'visium_hd', resolution, f'{sample_id}.h5ad'))

    # Calculate crop coords for optimal presentation
    coords = auto_crop(adata, border_size=0.1)


    # --- Load the genes from a file

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
            'HEVs': 'HEVs',
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
            'HEVs': 'hevs_filtered_endothelial_genes.csv'
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

        # Ensure cutoffs do not exceed the number of available genes
        cutoffs = [c for c in cutoffs if c <= len(all_genes)]
        if not cutoffs or (cutoffs[-1] < len(all_genes)):  
            cutoffs.append(len(all_genes))  # Ensure the full list is included

        # Set seed
        seed = 42

        # Extract the layer you want to score on
        score_obj = ad.AnnData(adata.layers[use_layer])
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


        # --- Plot 1: Plot the score results for different cutoffs

        ncols = 3
        nrows = math.ceil((len(cutoffs))/ncols)

        plt.figure(figsize=(8*ncols, 6.5*nrows))

        # Scoring results
        for n, cutoff in enumerate(cutoffs):
            score_key = f'{score_name}_score_{cutoff}'
            ax = plt.subplot(nrows, ncols, n + 1)
            # Scoring results per cutoff
            sq.pl.spatial_scatter(adata, color=score_key, ax=ax, img_alpha=0.4, cmap='viridis', 
                crop_coord=coords, vmin=0, shape='square', )

        plt.suptitle(f'{sample_id} ({subtype}) - {score_name} score with different top gene cutoffs ({cutoffs}, resolution={resolution})', y=1.05, fontweight='bold', x=0.5)


        # --- Plot 2: Different vmax

        vmax = 0.8

        ncols = 3
        nrows = math.ceil((len(cutoffs))/ncols)

        plt.figure(figsize=(8*ncols, 6.5*nrows))

        # Scoring results
        for n, cutoff in enumerate(cutoffs):
            score_key = f'{score_name}_score_{cutoff}'
            ax = plt.subplot(nrows, ncols, n + 1)
            # Scoring results per cutoff
            sq.pl.spatial_scatter(adata, color=score_key, ax=ax, img_alpha=0.4, cmap='viridis', vmax=0.8,
                crop_coord=coords, vmin=0, shape='square', )

        plt.suptitle(f'{sample_id} ({subtype}) - {score_name} score with different top gene cutoffs ({cutoffs}) (vmax={vmax}, resolution={resolution})', y=1.05, fontweight='bold', x=0.5)


        # --- Save the results

        out_name = f'{score_name}_score'

        qc_dir = os.path.join(res_dir, '04_vascular_microenvironment', '01_identify_HEV_niches', '01_score_signature', out_name)
        os.makedirs(qc_dir, exist_ok=True)

        out_dir = os.path.join(qc_dir, resolution, 'plots')
        os.makedirs(out_dir, exist_ok=True)

        filename = f'{sample_id}_{out_name}_{resolution}.pdf'

        print(f'Saving plot as {filename} in {out_dir}')
        save_multi_image(os.path.join(out_dir, filename))
        plt.close('all')

        # Also save the scores + which genes were used
        out_dir_csv = os.path.join(qc_dir, resolution, 'csv_files')
        os.makedirs(out_dir_csv, exist_ok=True)

        score_keys = [f'{score_name}_score_' + str(item) for item in cutoffs]

        save_score = adata.obs[score_keys]
        filename = f'{sample_id}_{out_name}.csv'

        save_score.to_csv(os.path.join(out_dir_csv, filename))

        # Save genes used
        # Also save the pandas dataframe as a csv file (to be able to check which genes were used in each sample)
        genes_used = gene_list[0:cutoff]
        df = pd.DataFrame(genes_used, columns=["genes_used"])
        df.to_csv(os.path.join(out_dir_csv, f'{sample_id}_genes_used.csv'), index=False)

        print(f'Finished scoring for resolution {resolution}')


# --- Finish

print("Script successfully completed")
