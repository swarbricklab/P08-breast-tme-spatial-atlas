#!/usr/bin/env python3

# --- Aim
# Create an object with raw counts from miniatlas GEO dataset
# Save results as .h5ad file

# --- Define environment
# env: spatial-pilot

# --- Load packages
import argparse
import os

import scanpy as sc
import pandas as pd
import anndata as ad


# --- Argparse arguments

parser = argparse.ArgumentParser(description="Import Visium data & save as .h5ad file",
								 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

# Required
parser.add_argument("--raw_data_dir", help="Path to raw data")
parser.add_argument("--h5ad_dir", help="Path to directory for .h5ad output files")
parser.add_argument("--sample_id", help="Sample ID")
parser.add_argument("--project_dir", help="Path to the root of the project directory")

args = vars(parser.parse_args())


# --- Arguments
raw_data_dir = args["raw_data_dir"]
h5ad_dir=args["h5ad_dir"]
sample_id=args["sample_id"]
project_dir=args["project_dir"]


# --- Path to miniatlas count matrix

miniatlas_path = os.path.join(project_dir, 'resources', 'cell2location', 'miniatlas_GEO')
subdir = 'Wu_etal_2021_BRCA_scRNASeq'
count_matrix_file = 'count_matrix_sparse.mtx'
genes_file = 'count_matrix_genes.tsv'
cells_file = 'count_matrix_barcodes.tsv'
metadata_file = 'metadata.csv'

# Load data
counts = sc.read_mtx(os.path.join(miniatlas_path, subdir, count_matrix_file)).T
genes = pd.read_csv(os.path.join(miniatlas_path, subdir, genes_file), header=None, sep='\t')
cells = pd.read_csv(os.path.join(miniatlas_path, subdir, cells_file), header=None, sep='\t')
metadata = pd.read_csv(os.path.join(miniatlas_path, subdir, metadata_file), header=0, index_col=0)

# Make anndata object
adata = ad.AnnData(X=counts.X)
adata.var_names = genes[0].values
adata.obs_names = cells[0].values
adata.obs = metadata.set_index(adata.obs_names)


# --- Load mapping file
# John's mapping file for mapping gene symbol to ensembl ids
subdir = 'ensembl_mapping'
mapping_file = 'map.csv'
mapping = pd.read_csv(os.path.join(miniatlas_path, subdir, mapping_file))

# Create a dictionary from gene symbol to ensembl id
symbol_to_ensembl = dict(zip(mapping['dedup_name'], mapping['ensembl_id']))

# Map gene symbols to ensembl ids
adata.var['ensembl_id'] = adata.var_names.map(symbol_to_ensembl)

# Also map the original names
ensembl_to_original = dict(zip(mapping['ensembl_id'], mapping['original_name']))
adata.var['original_name'] = adata.var['ensembl_id'].map(ensembl_to_original)

# Check every gene has an ensembl id
if adata.var['ensembl_id'].isnull().sum() > 0:
	missing_genes = adata.var_names[adata.var['ensembl_id'].isnull()]
	print("!! The following genes are missing ensembl ids:", missing_genes)

# For the ones that are missing, get them straight from genes.tsv
genes_file = pd.read_csv(os.path.join(miniatlas_path, subdir, 'genes.tsv'), header=None, sep='\t', names=['ensembl_id', 'gene_symbol'])

# Replace __ by --
genes_file['gene_symbol'] = genes_file['gene_symbol'].str.replace('__', '--')

# Adjust XXyac-YX65C7_A.2 to XXyac-YX65C7-A.2
genes_file['gene_symbol'] = genes_file['gene_symbol'].str.replace(r'_([A-Z])\.', r'-\1.', regex=True)

# For the missing genes, identify the ensembl ids from genes.tsv
for gene in adata.var_names[adata.var['ensembl_id'].isnull()]:
	ensembl_id = genes_file.loc[genes_file['gene_symbol'] == gene, 'ensembl_id']
	if len(ensembl_id) == 1:
		adata.var.loc[gene, 'ensembl_id'] = ensembl_id.values[0]
	else:
		print(f"!! Gene {gene} not found in genes.tsv")

# Check if every gene has an ensembl id now
if adata.var['ensembl_id'].isnull().sum() > 0:
	missing_genes = adata.var_names[adata.var['ensembl_id'].isnull()]
	print("!! The following genes are still missing ensembl ids:", missing_genes)
else:
	print("All genes have ensembl ids now.")

# Set ensembl id as var_names
adata.var_names = adata.var['ensembl_id']

# Also save the dedup name
adata.var['dedup_name'] = adata.var_names.map({v: k for k, v in symbol_to_ensembl.items()})

# Check if any values are missing
if adata.var['dedup_name'].isnull().sum() > 0:
	missing_ensembls = adata.var_names[adata.var['dedup_name'].isnull()]
	print("!! The following ensembl ids are missing dedup names:", missing_ensembls)


# --- Do some filtering similar to what was done in the paper

# Annotate mitochondrial genes
adata.var['mt'] = adata.var['dedup_name'].str.startswith('MT-', na=False)
print(f"Mitochondrial genes: {adata.var['mt'].sum()}")

# Calculate QC metrics
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

# Filter cells:
# - gene count > 200
# - UMI count > 250
# - mitochondrial percentage < 20%
print(f"Cells before filtering: {adata.n_obs}")
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_cells(adata, min_counts=250)
adata = adata[adata.obs['pct_counts_mt'] < 20, :].copy()
print(f"Cells after filtering: {adata.n_obs}")


# --- Store raw counts and normalize

# Store raw counts in .raw slot
adata.raw = adata

# Normalize and log-transform, store in .X
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)


# --- Save output

out_file = 'miniatlas_GEO_raw.h5ad'
adata.write(os.path.join(miniatlas_path, out_file))

print("Script succesfully completed")
