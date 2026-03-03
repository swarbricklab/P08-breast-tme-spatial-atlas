#!/usr/bin/env python3

# --- Aim
# Use cell2location for deconvolution
# Recycle a pre-built model and do the cell type mapping only

# This script: run on all samples as a single job, splitting the object into chunks
# this is necessary to allow cell2loc to run on objects with large number of spots (e.g. higher res Vis HD)
# without splitting the GPU runs out of memory

# --- Define environment
# env: cell2location

# --- Load packages
import argparse
import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
import json

import cell2location
import pyro

# Set seed for reproducibility
import scvi
seed = 2023
scvi.settings.seed = seed
np.random.seed(seed)

from matplotlib import rcParams
rcParams['pdf.fonttype'] = 42 # enables correct plotting of text for PDFs

import pandas as pd
import os

import anndata as ad


# --- Argparse arguments

parser = argparse.ArgumentParser(description="Import Visium data & save as .h5ad file",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument("--raw_data_dir", help="Path to raw data")
parser.add_argument("--h5ad_dir", help="Path to directory for .h5ad output files")
parser.add_argument("--sample_id", help="Sample ID")
parser.add_argument("--project_dir", help="Path to the root of the project directory")
parser.add_argument("--sample_sheet", help="Path to sample sheet CSV file", default=None)
parser.add_argument("--res_dir", help="Path to results output directory", default=None)

args = vars(parser.parse_args())


# --- Arguments
sample_id = args['sample_id']
project_dir = args['project_dir']
h5ad_dir = args['h5ad_dir']


# --- Extract some info from the sample sheet

# Determine sample sheet path
if args["sample_sheet"] is not None:
    sample_sheet_path = args["sample_sheet"]
else:
    sample_sheet_path = os.path.join(project_dir, 'config', 'sample_sheet.csv')

sample_sheet = pd.read_csv(sample_sheet_path)


# --- Load the model

# Decide which sc reference to use
sc_ref = 'miniatlas_GEO'

# Decide which level to run
level = 'major'

model_dir = os.path.join(project_dir, 'results', '02_cell_type_mapping', 'deconvolution', 'model', f'{level}_level')

adata_file = f"reference_signatures/sc.h5ad"
adata_ref = sc.read_h5ad(os.path.join(model_dir, adata_file))
mod = cell2location.models.RegressionModel.load(os.path.join(model_dir, 'reference_signatures'), adata_ref)

# Export estimated expression in each cluster
if 'means_per_cluster_mu_fg' in adata_ref.varm.keys():
    inf_aver = adata_ref.varm['means_per_cluster_mu_fg'][[f'means_per_cluster_mu_fg_{i}'
                                    for i in adata_ref.uns['mod']['factor_names']]].copy()
else:
    inf_aver = adata_ref.var[[f'means_per_cluster_mu_fg_{i}'
                                    for i in adata_ref.uns['mod']['factor_names']]].copy()

inf_aver.columns = adata_ref.uns['mod']['factor_names']
inf_aver.iloc[0:5, 0:5]


# --- Load the visium data

# Decide which resolution to use
resolution = 'square_016um'

input_adata_dir = os.path.join(h5ad_dir, 'visium_hd', resolution)
adata_vis = ad.read_h5ad(os.path.join(input_adata_dir, f'{sample_id}.h5ad'))

# Note: .X contains raw counts (cell2location requires raw counts)

# --- Visium preprocessing

# Switch out gene_ids for ensembl ids
adata_vis.var['SYMBOL'] = adata_vis.var_names
adata_vis.var.set_index('gene_ids', drop=True, inplace=True)

# find mitochondria-encoded (MT) genes
adata_vis.var['MT_gene'] = [gene.startswith('MT-') for gene in adata_vis.var['SYMBOL']]

# remove MT genes for spatial mapping (keeping their counts in the object)
adata_vis.obsm['MT'] = adata_vis[:, adata_vis.var['MT_gene'].values].X.toarray()
adata_vis = adata_vis[:, ~adata_vis.var['MT_gene'].values]

# find shared genes and subset both anndata and reference signatures
intersect = np.intersect1d(adata_vis.var_names, inf_aver.index)
adata_vis = adata_vis[:, intersect].copy()
inf_aver = inf_aver.loc[intersect, :].copy()


# --- Prepare for saving

out_dir = os.path.join(project_dir, 'results', '02_cell_type_mapping', 'deconvolution', 'visium_hd', f'{level}_level', resolution, sample_id)
os.makedirs(out_dir, exist_ok=True)

# Set ncells per location based on resolution
cells_dict = {
    'square_048um': 20,
    'square_016um': 3,
}

if resolution in cells_dict:
    N_cells_per_location = cells_dict[resolution]
else:
    N_cells_per_location = 10

# Also make an output directory for the plots
out_dir_plots = os.path.join(out_dir, 'plots')
os.makedirs(out_dir_plots, exist_ok=True)


# --- Split visium data into chunks

chunk_size = 30_000
chunks = [i for i in range(int(np.ceil(adata_vis.n_obs / chunk_size)))]

# Initialize training_batch
adata_vis.obs['training_batch'] = 0  

# Assign random chunks directly (since there's only one unique sample)
adata_vis.obs['training_batch'] = np.random.choice(
    chunks, size=adata_vis.n_obs, replace=True, p=None
)

# Save the training batch assignments
adata_vis.obs['training_batch'].to_csv(os.path.join(out_dir, f'{sample_id}_training_batches.csv'))

# Create a copy
adata_vis_full = adata_vis.copy()

# Initialize obsm matrices
for k in ['means', 'stds', 'q05', 'q95']:
    adata_vis_full.obsm[f"{k}_cell_abundance_w_sf"] = np.zeros((adata_vis_full.n_obs, inf_aver.shape[1]))

# Check batch assignments
print(f'Training batch sizes: {adata_vis.obs["training_batch"].value_counts()}')    


# --- Record parameters

detection_alpha=20
batch_size=None
max_epochs=30000
train_size=1
use_gpu=True

parameters_log = {'N_cells_per_location': N_cells_per_location, 'detection_alpha': detection_alpha,'path_to_model': model_dir, 'max_epochs': max_epochs, 
                  'batch_size': str(batch_size), 'train_size':train_size, 'use_gpu': use_gpu, 'subset_level': level, 'sc_ref': sc_ref, 'chunk_size': chunk_size,}

chunk_sizes = adata_vis.obs["training_batch"].value_counts().to_dict()

for key, value in chunk_sizes.items():
    parameters_log[f'chunk_{key}'] = value

# Record parameters
with open(os.path.join(out_dir, f'cell2location_parameters.json'), 'w') as json_file:
    json.dump(parameters_log, json_file, indent=4)

# Run cell2location
print(f'Running sample {sample_id}')


# --- Spatial mapping

run_name_global = sample_id
os.makedirs(os.path.join(out_dir, 'model_output'), exist_ok=True)

# Save cell abundance
abundance_df_dict = {}

# Submit chunk as separate jobs
for batch in sorted(adata_vis.obs['training_batch'].unique()):

    # create and train the model
    scvi_run_name = f'{run_name_global}_batch{batch}_seed{seed}'
    print(scvi_run_name)

    training_batch_index = adata_vis_full.obs['training_batch'].isin([batch])
    adata_vis = adata_vis_full[training_batch_index, :].copy()
    
    # prepare anndata for scVI model
    cell2location.models.Cell2location.setup_anndata(
        adata=adata_vis, batch_key=None
    )

    # train as normal
    # create and train the model
    mod = cell2location.models.Cell2location(
        adata_vis, cell_state_df=inf_aver,
        # the expected average cell abundance: tissue-dependent
        # hyper-prior which can be estimated from paired histology:
        N_cells_per_location=N_cells_per_location,
        # hyperparameter controlling normalisation of
        # within-experiment variation in RNA detection:
        detection_alpha=detection_alpha
    )

    mod.view_anndata_setup()

    # --- Train cell2location
    mod.train(max_epochs=max_epochs,
            # train using full data (batch_size=None)
            batch_size=batch_size,
            # use all data points in training because
            # we need to estimate cell abundance at all locations
            train_size=train_size,
            use_gpu=use_gpu,
            )

    # plot ELBO loss history during training, removing first 100 epochs from the plot
    mod.plot_history(1000)
    plt.legend(labels=['full data training'])
    plt.savefig(os.path.join(out_dir_plots, f'training_elbo_loss_{scvi_run_name}.pdf'))
    plt.close()

    # export posterior
    # In this section, we export the estimated cell abundance (summary of the posterior distribution).
    adata_vis = mod.export_posterior(
        adata_vis, sample_kwargs={
            'batch_size': int(np.ceil(adata_vis.n_obs / 4)), 'use_gpu': True,
        },
    )

    # --- Assess mapping quality
    mod.plot_QC()
    plt.savefig(os.path.join(out_dir_plots, f'mapping_quality_{scvi_run_name}.pdf'))
    plt.close()

    # -- Save cell type mapping
    # 5% quantile of the posterior distribution, representing the value of cell abundance that the model has high confidence in (aka ‘at least this amount is present’)
    adata_vis.obs[adata_vis.uns['mod']['factor_names']] = adata_vis.obsm['q05_cell_abundance_w_sf']
    abundance_df = adata_vis.obs[adata_vis.uns['mod']['factor_names']]
    # abundance_df.to_csv(os.path.join(out_dir, 'model_output', f'q05_cell_abundance_w_sf_{scvi_run_name}.csv'))
    abundance_df_dict[scvi_run_name] = abundance_df

    # -- Copy cell2location results to the main object
    for k in adata_vis_full.obsm.keys():
        adata_vis_full.obsm[k][training_batch_index, :] = adata_vis.obsm[k].copy()
    
    adata_vis_full.uns[f'mod_{batch}'] = adata_vis.uns['mod'].copy()


# --- Complete the full object with cell abundance estimates from batched analysis

# Save anndata object with results
adata_file = f"{sample_id}_cell2location.h5ad"
# Save with lzf to stop dvc from truncating the file
adata_vis_full.write(os.path.join(out_dir, 'model_output', adata_file), compression='lzf')

# adata_vis_full = sc.read_h5ad(os.path.join(out_dir, 'model_output', adata_file))


# --- Export estimated cell abundance

# Also save the cell type mapping as a csv file for ease of access
# 5% quantile of the posterior distribution, representing the value of cell abundance that the model has high confidence in (aka ‘at least this amount is present’)

# Get factor names from the first batch (all batches should have the same factor names)
first_batch = sorted(adata_vis_full.obs['training_batch'].unique())[0]
factor_names = adata_vis_full.uns[f'mod_{first_batch}']['factor_names']

adata_vis_full.obs[factor_names] = adata_vis_full.obsm['q05_cell_abundance_w_sf']
abundance_df = adata_vis_full.obs[factor_names]

# Save one dir up
out_dir_cell_abundance = os.path.dirname(out_dir)
abundance_df.to_csv(os.path.join(out_dir_cell_abundance, f'{sample_id}_{level}_q05_cell_abundance_w_sf.csv'))


# --- Finish

print("Script succesfully completed")