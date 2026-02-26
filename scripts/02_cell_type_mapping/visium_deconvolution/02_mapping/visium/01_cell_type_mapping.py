#!/usr/bin/env python3

# --- Aim
# Use cell2location for deconvolution
# Recycle a pre-built model and do the cell type mapping only

# This script: run on all samples as a single job

# --- Define environment
# env: cell2location

# --- Load packages
import argparse
import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
import json

from datetime import date
from datetime import datetime

import cell2location
import json

# Set seed for reproducibility
import scvi
scvi.settings.seed = 2023

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
level = 'subset'

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

# Apply the function
if sc_ref == 'miniatlas_GEO':
    model_dir = os.path.join(project_dir, 'results', 'cell2location', sc_ref, f'{level}_level', 'model')
else:
    model_dir = os.path.join(project_dir, 'resources', 'cell2location', 'pretrained_models', sc_ref,)
    if 'miniatlas' in sc_ref:
        model_dir = os.path.join(model_dir, f'{level}_level',)

most_recent_dir_path = find_most_recent(model_dir)

adata_file = f"reference_signatures/sc.h5ad"
adata_ref = sc.read_h5ad(os.path.join(most_recent_dir_path, adata_file))
mod = cell2location.models.RegressionModel.load(os.path.join(most_recent_dir_path, 'reference_signatures'), adata_ref)

# Record parameters
N_cells_per_location=20
detection_alpha=20
batch_size=None
max_epochs=30000
train_size=1
use_gpu=True

parameters_log = {'N_cells_per_location': N_cells_per_location, 'detection_alpha': detection_alpha,'path_to_model': most_recent_dir_path, 'max_epochs': max_epochs, 
                  'batch_size': str(batch_size), 'train_size':train_size, 'use_gpu': use_gpu, 'subset_level': level, 'sc_ref': sc_ref}


# Run cell2location
print(f'Running sample {sample_id}')

# --- Load the visium data
input_adata_dir = os.path.join(h5ad_dir, 'raw')
adata_vis = ad.read(os.path.join(input_adata_dir, f'{sample_id}_raw.h5ad'))

# --- Visium preprocessing

# Switch out gene_ids for ensembl ids
adata_vis.var['SYMBOL'] = adata_vis.var_names

if sc_ref == 'miniatlas_GEO':
    adata_vis.var.set_index('gene_ids', drop=True, inplace=True)

# find mitochondria-encoded (MT) genes
adata_vis.var['MT_gene'] = [gene.startswith('MT-') for gene in adata_vis.var['SYMBOL']]

# remove MT genes for spatial mapping (keeping their counts in the object)
adata_vis.obsm['MT'] = adata_vis[:, adata_vis.var['MT_gene'].values].X.toarray()
adata_vis = adata_vis[:, ~adata_vis.var['MT_gene'].values]

# --- Prepare for saving

res_dir = os.path.join(project_dir, 'results')
qc_dir = os.path.join(res_dir, 'cell2location', sc_ref,)
if 'miniatlas' in sc_ref:
    qc_dir = os.path.join(qc_dir, f'{level}_level',)

qc_dir = os.path.join(qc_dir, sample_id)
os.makedirs(qc_dir, exist_ok=True)

out_dir = os.path.join(qc_dir, f'n_cells_{N_cells_per_location}')
os.makedirs(out_dir, exist_ok=True)

# Also make an output directory for the plots
out_dir_plots = os.path.join(out_dir, 'plots')
os.makedirs(out_dir_plots, exist_ok=True)

# Record parameters
with open(os.path.join(out_dir, f'cell2location_parameters.json'), 'w') as json_file:
    json.dump(parameters_log, json_file, indent=4)


# --- Export cell type signatures

# export estimated expression in each cluster
if 'means_per_cluster_mu_fg' in adata_ref.varm.keys():
    inf_aver = adata_ref.varm['means_per_cluster_mu_fg'][[f'means_per_cluster_mu_fg_{i}'
                                    for i in adata_ref.uns['mod']['factor_names']]].copy()
else:
    inf_aver = adata_ref.var[[f'means_per_cluster_mu_fg_{i}'
                                    for i in adata_ref.uns['mod']['factor_names']]].copy()

inf_aver.columns = adata_ref.uns['mod']['factor_names']
inf_aver.iloc[0:5, 0:5]

# --- Spatial mapping

# find shared genes and subset both anndata and reference signatures
intersect = np.intersect1d(adata_vis.var_names, inf_aver.index)
adata_vis = adata_vis[:, intersect].copy()
inf_aver = inf_aver.loc[intersect, :].copy()

# prepare anndata for cell2location model
cell2location.models.Cell2location.setup_anndata(adata=adata_vis, batch_key=None)

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
plt.legend(labels=['full data training']);
plt.savefig(os.path.join(out_dir_plots, 'training_elbo_loss.pdf'))
plt.close()

# --- Export estimated cell abundance

# In this section, we export the estimated cell abundance (summary of the posterior distribution).
adata_vis = mod.export_posterior(
    adata_vis, sample_kwargs={'num_samples': 1000, 'batch_size': mod.adata.n_obs, 'use_gpu': True}
)
# Save model
mod.save(os.path.join(out_dir, 'model_output'), overwrite=True)

# Save anndata object with results
adata_file = f"{sample_id}_cell2location.h5ad"
adata_vis.write(os.path.join(out_dir, 'model_output', adata_file))

# Can load again like this
# adata_vis = sc.read_h5ad(os.path.join(out_dir, 'model_output', adata_file))
# mod = cell2location.models.Cell2location.load(os.path.join(out_dir, 'model_output'), adata_vis)

# Also save the cell type mapping as a csv file for ease of access
# 5% quantile of the posterior distribution, representing the value of cell abundance that the model has high confidence in (aka ‘at least this amount is present’)
adata_vis.obs[adata_vis.uns['mod']['factor_names']] = adata_vis.obsm['q05_cell_abundance_w_sf']
abundance_df = adata_vis.obs[adata_vis.uns['mod']['factor_names']]
abundance_df.to_csv(os.path.join(out_dir, 'model_output', 'q05_cell_abundance_w_sf.csv'))

# --- Assess mapping quality
mod.plot_QC()
plt.savefig(os.path.join(out_dir_plots, 'mapping_quality.pdf'))
plt.close()


# --- Finish

print("Script succesfully completed")