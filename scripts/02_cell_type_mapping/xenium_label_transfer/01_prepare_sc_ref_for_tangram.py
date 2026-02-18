#!/usr/bin/env python3

# --- Aim
# Prepare single-cell reference for Tangram label transfer
# - Subset to genes present in Xenium data
# - Find marker genes per cell type
# - Save processed reference and markers

# --- Define environment
# env: tangram

# --- Load packages
import argparse
import os

import scanpy as sc
import pandas as pd
import numpy as np
import anndata as ad


# --- Argparse arguments

parser = argparse.ArgumentParser(description="Prepare sc reference for Tangram",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

# Required
parser.add_argument("--project_dir", help="Path to the root of the project directory")
parser.add_argument("--annotation_level", help="Annotation level to use", default='subset')
parser.add_argument("--sample_id", help="Sample ID (optional, for output organization)", default=None)

args = vars(parser.parse_args())


# --- Arguments

project_dir = args['project_dir']


# --- Load data

# Load probe set files to get gene list
probe_set_dir = os.path.join(project_dir, 'resources', 'probe_set')

print("Loading probe set files...")
# Load the 5k panel metadata
panel_5k = pd.read_csv(os.path.join(probe_set_dir, 'XeniumPrimeHuman5Kpan_tissue_pathways_metadata.csv'))
# Load the 100 gene breast panel
panel_100g = pd.read_csv(os.path.join(probe_set_dir, '23Q3D7_hBreast_100g_gene_list.csv'))

# Extract Ensembl IDs from both files
genes_5k = panel_5k['gene_id'].values
genes_100g = panel_100g['Ensembl ID'].values

# Combine and get unique genes
genes_in_adata = np.unique(np.concatenate([genes_5k, genes_100g]))

print(f"Total unique genes from probe sets: {len(genes_in_adata)}")

# Also make an ensembl to gene symbol dict based on the csv files
ensembl_to_symbol_5k = dict(zip(panel_5k['gene_id'], panel_5k['gene_name']))
ensembl_to_symbol_100g = dict(zip(panel_100g['Ensembl ID'], panel_100g['Gene']))
ensembl_to_symbol = {**ensembl_to_symbol_5k, **ensembl_to_symbol_100g}

# Load the sc_ref
sc_ref_dir = os.path.join(project_dir, 'resources', 'sc_ref')
ref_file = 'miniatlas_GEO_raw.h5ad'

print(f"Loading reference from {os.path.join(sc_ref_dir, ref_file)}")
ref = ad.read_h5ad(os.path.join(sc_ref_dir, ref_file))


# --- Prepare reference

print("Subsetting reference to Xenium genes...")

# Set raw layer to .X BEFORE subsetting (more efficient)
ref.X = ref.raw.X.copy()
del ref.raw  # Remove raw to save memory

# Subset using boolean mask (much faster than label-based subsetting)
gene_mask = ref.var_names.isin(genes_in_adata)
ref = ref[:, gene_mask].copy()

print(f"Reference subsetted to {ref.n_vars} genes")


# --- Find marker genes

def map_ensembl(gene):
    return ensembl_to_symbol.get(gene, gene)

levels = ['major', 'minor', 'subset']

print("Finding marker genes for all annotation levels...")

# Log norm to find markers
ref.layers['norm'] = sc.pp.normalize_total(ref, target_sum=1e4, inplace=False)['X']

# log transform
ref.layers['log_norm'] = sc.pp.log1p(ref.layers['norm'], copy=True)

# Create output directory
out_dir = os.path.join(project_dir, 'results', 'annotation', 'tangram', 'reference')
os.makedirs(out_dir, exist_ok=True)

# Loop through each annotation level
for current_level in levels:
    print(f"\nProcessing {current_level} level...")
    
    # Set cell type key
    cell_type_key = f'celltype_{current_level}'
    
    # Find markers
    sc.tl.rank_genes_groups(ref, groupby=cell_type_key, method='wilcoxon', use_raw=False, layer='log_norm')
    
    # Get the top markers
    # Calculate how many markers per cell type to use
    # Aim for a total of ~1500 markers
    markers_per_cell_type = 1500 // ref.obs[cell_type_key].nunique()
    markers_df = pd.DataFrame(ref.uns["rank_genes_groups"]["names"]).iloc[0:markers_per_cell_type, :]
    markers = list(np.unique(markers_df.melt().value.values))
    
    print(f'  Total unique markers: {len(markers)}')
    print(f'  Number of markers per cell type: {markers_per_cell_type}')

    # Map the gene symbols to ensembl IDs for saving
    # More interpretable
    markers_df_mapped = markers_df.applymap(map_ensembl)
    
    # Save the markers per cell type
    markers_out_path = os.path.join(out_dir, 'markers', f'{current_level}_cell_type_markers.csv')
    os.makedirs(os.path.dirname(markers_out_path), exist_ok=True)
    print(f'  Saving markers to {markers_out_path}')
    markers_df_mapped.to_csv(markers_out_path, index=False)


# --- Save outputs

# Save the processed reference
ref_name = ref_file.replace('.h5ad', '_processed.h5ad')
ref_out_path = os.path.join(out_dir, ref_name)
print(f"\nSaving processed reference to {ref_out_path}")
ref.write_h5ad(ref_out_path)

# Also save a csv of shared genes + ensembl id
shared_genes = pd.DataFrame({
    'ensembl_id': ref.var_names,
    'gene_symbol': [map_ensembl(gene) for gene in ref.var_names]
})

shared_genes_out_path = os.path.join(out_dir, 'shared_genes.csv')
print(f"Saving shared genes to {shared_genes_out_path}")
shared_genes.to_csv(shared_genes_out_path, index=False)

print("\nReference preparation completed successfully")
