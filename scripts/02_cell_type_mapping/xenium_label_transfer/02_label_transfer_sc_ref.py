#!/usr/bin/env python3

# --- Aim
# Try out label transfer from miniatlas for quick annotation of xenium data
# Use tangram for label transfer


# --- Define environment
# env: tangram

# --- Load packages
import argparse
import os

import scanpy as sc
import pandas as pd
import numpy as np

import torch
import tangram as tg
import anndata as ad


# --- Argparse arguments

parser = argparse.ArgumentParser(description="Import Visium data & save as .h5ad file",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

# Required
parser.add_argument("--raw_data_dir", help="Path to raw data")
parser.add_argument("--h5ad_dir", help="Path to directory for .h5ad output files")
parser.add_argument("--sample_id", help="Sample ID")
parser.add_argument("--project_dir", help="Path to the root of the project directory")
parser.add_argument("--min_count", help="Minimal number of counts per spot", type=int)
parser.add_argument("--min_gene", help="Minimal number of genes per spot", type=int)
parser.add_argument("--outs_path", help="Path to spaceranger 'outs'")
parser.add_argument("--use_data", help="Which data to use (log_norm or SCT)")


args = vars(parser.parse_args())


# --- Arguments
raw_data_dir = args["raw_data_dir"]
h5ad_dir=args["h5ad_dir"]
sample_id=args["sample_id"]

min_count=args["min_count"]
min_gene=args["min_gene"]

outs_path=args["outs_path"]
project_dir=args['project_dir']


# --- Load data

# Load the log_norm file
adata = ad.read_h5ad(os.path.join(h5ad_dir, 'log_norm', f'{sample_id}_log_norm.h5ad'))
to_keep = adata.obs_names

# Load the cell feature matrix as an anndata object
raw_data_dir_sample = os.path.join(raw_data_dir, sample_id)

# Load the cell feature matrix as an anndata object
cell_feature_matrix = 'cell_feature_matrix.h5'

adata_raw = sc.read_10x_h5(
    filename=os.path.join(raw_data_dir_sample, cell_feature_matrix)
)

adata.X = adata_raw.to_df().loc[list(to_keep),].values


# --- Load pre-processed reference and markers
# Note: Run 01_prepare_sc_ref_for_tangram.py first to generate these files

level = 'subset'
cell_type_key = f'celltype_{level}'

ref_dir = os.path.join(project_dir, 'results', 'annotation', 'tangram', 'reference')

# Load the processed reference
print(f"Loading pre-processed reference...")
ref = ad.read_h5ad(os.path.join(ref_dir, 'miniatlas_GEO_raw_processed.h5ad'))

# Load the markers
print(f"Loading markers...")
markers_df = pd.read_csv(os.path.join(ref_dir, 'markers', f'{level}_cell_type_markers.csv'))
markers = list(np.unique(markers_df.melt().value.values))
print(f'Total unique markers: {len(markers)}')

# Convert to ensembl ids using shared_genes.csv
shared_genes = pd.read_csv(os.path.join(ref_dir, 'shared_genes.csv'))
gene_to_ensembl = dict(zip(shared_genes['gene_symbol'], shared_genes['ensembl_id']))
markers = [gene_to_ensembl.get(gene, gene) for gene in markers]
markers_df = markers_df.applymap(lambda gene: gene_to_ensembl.get(gene, gene))

# Set adata var names to ensembl ids
adata.var['gene_symbol'] = adata.var_names
adata.var_names = adata.var['gene_ids']

# --- Normalize the data

# But do not log transform
sc.pp.normalize_total(adata, target_sum=1e4)

# For ref: set the norm layer as the main layer
ref.X = ref.layers['norm']

# Run scanpy preprocessing
tg.pp_adatas(ref, adata, genes=markers)


# --- Prepare for saving

out_dir = os.path.join(project_dir, 'results', 'annotation', 'tangram', sample_id, f'{level}_level')
os.makedirs(out_dir, exist_ok=True)


# --- Cell type mapping

# Set the random seed for reproducibility
seed = 42

# Run Tangram mapping
ad_map = tg.map_cells_to_space(
    adata_sc=ref,
    adata_sp=adata,
    mode='clusters',  # or 'clusters'
    cluster_label=cell_type_key,
    device="cuda:0" if tg.torch.cuda.is_available() else "cpu",
    num_epochs=1000,
    random_state=seed,
)

# Change the obs names to the cell types
cell_types = list(ad_map.obs[cell_type_key])
ad_map.obs_names = cell_types

# Extract a df with cluster density
cluster_density = ad_map.obs[[cell_type_key, 'cluster_density']]
# Save csv
cluster_density.to_csv(os.path.join(out_dir, f"{sample_id}_cluster_density.csv"), index=False)

# Get the cell type probabilities
prob_df = ad_map.to_df()

# Extract a df with for each cell the label with max probability + probability
# Cells are columns
# Cell type labels are index
annotation_df = prob_df.idxmax(axis=0).to_frame(name='pred_cell_type')
annotation_df['pred_cell_type_prob'] = prob_df.max(axis=0)


# --- Save the output

# Save probabilities df to gzipped csv
prob_df = prob_df.T  # Transpose to have cell types as index
prob_df.to_csv(
    os.path.join(out_dir, f"{sample_id}_pred_cell_type_probabilities.csv.gz"),
    compression='gzip'
)

# Save the annotation DataFrame to a csv
annotation_df.to_csv(
    os.path.join(out_dir, f"{sample_id}_pred_cell_types.csv"),
)

print("Script succesfully completed")