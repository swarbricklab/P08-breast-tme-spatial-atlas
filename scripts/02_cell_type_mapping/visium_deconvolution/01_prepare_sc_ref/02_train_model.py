#!/usr/bin/env python3

# --- Aim
# Train cell2location model to be used on many samples


# --- Define environment
# env: cell2location

# --- Load packages
import argparse
import os

import scanpy as sc
import matplotlib.pyplot as plt

# this line should go before importing cell2location
os.environ["THEANO_FLAGS"] = 'device=cuda,floatX=float32,force_device=True'
import cell2location

# Set seed for reproducibility
import scvi
scvi.settings.seed = 2023

from matplotlib import rcParams
rcParams['pdf.fonttype'] = 42 # enables correct plotting of text for PDFs

# Some imports for qc plots
import numpy as np
import matplotlib
from scvi import REGISTRY_KEYS
from scipy.sparse import issparse


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
sample_id=args["sample_id"]
project_dir=args["project_dir"]
h5ad_dir=args["h5ad_dir"]


# --- Prepare for saving

# Decide which level to run
level = 'major'

out_dir = os.path.join(project_dir, 'results', '02_cell_type_mapping', 'deconvolution', 'model', f'{level}_level')
os.makedirs(out_dir, exist_ok=True)

# Also make an output directory for the plots
out_dir_plots = os.path.join(out_dir, 'plots')
os.makedirs(out_dir_plots, exist_ok=True)


# --- Load the single-cell ref

# Use the miniatlas downloaded from cellxgene
resource_dir = os.path.join(h5ad_dir, 'sc_ref')
reference = 'miniatlas_GEO'

ref_file = f'{reference}.h5ad'

adata_ref = sc.read(os.path.join(resource_dir, ref_file))

# Extract raw counts from the counts layer (cell2location requires raw counts)
if 'counts' in adata_ref.layers.keys():
    adata_ref.X = adata_ref.layers['counts'].copy()
elif adata_ref.raw is not None:
    # Fallback to .raw if counts layer doesn't exist
    adata_ref = adata_ref.raw.to_adata()
else:
    raise ValueError("No raw counts found. Expected 'counts' layer or .raw attribute.")


# --- Filter genes

from cell2location.utils.filtering import filter_genes
selected = filter_genes(adata_ref, cell_count_cutoff=5, cell_percentage_cutoff2=0.03, nonz_mean_cutoff=1.12)

# filter the object
adata_ref = adata_ref[:, selected].copy()

plt.savefig(os.path.join(out_dir_plots, 'genes_filtered.pdf'))

plt.close()

# --- Prepare for regression model

ref_keys_dict = {
    'miniatlas_GEO': {
        'labels_key': f'celltype_{level}',
        'batch_key': 'orig.ident',
        'covariate_keys': [],
    }
}

# Print cell type numbers
labels_key = ref_keys_dict[reference]['labels_key']
print(f'Number of cells: {adata_ref.obs[labels_key].value_counts()}')

batch_key = ref_keys_dict[reference]['batch_key']
covariate_keys = ref_keys_dict[reference]['covariate_keys']

# prepare anndata for the regression model
cell2location.models.RegressionModel.setup_anndata(adata=adata_ref,
                        # 10X reaction / sample / batch
                        batch_key=batch_key,
                        # cell type, covariate used for constructing signatures
                        labels_key=labels_key,
                        # multiplicative technical effects (platform, 3' vs 5', donor effect)
                        categorical_covariate_keys=covariate_keys,
                       )


# create the regression model
# Note: ref needs to be raw (i.e. untransformed, unnormalized) for this to work
# If it's not raw this step will throw an error
from cell2location.models import RegressionModel
mod = RegressionModel(adata_ref)

# view anndata_setup as a sanity check
mod.view_anndata_setup()


# --- Training the model

mod.train(max_epochs=250, use_gpu=True)

# Determine if model needs more training

mod.plot_history(20)

plt.savefig(os.path.join(out_dir_plots, 'model_training.pdf'))

plt.close()

# In this section, we export the estimated cell abundance (summary of the posterior distribution).
adata_ref = mod.export_posterior(
    adata_ref, sample_kwargs={'num_samples': 1000, 'batch_size': 2500, 'use_gpu': True}
)

# Save model
mod.save(os.path.join(out_dir, f"reference_signatures"), overwrite=True)

# Save anndata object with results
adata_file = f"reference_signatures/sc.h5ad"
adata_ref.write(os.path.join(out_dir, adata_file))


# --- Examine QC plots

# cell2location's built-in mod.plot_QC() overlays its two QC panels on a single
# axis; the panels are re-plotted separately below so each is legible.

# Some default params
summary_name = 'means'
use_n_obs = 1000
scale_average_detection = True

# Plot 1: Reconstruction accuracy
# Should be roughly diagonal

if getattr(mod, "samples", False) is False:
    raise RuntimeError("self.samples is missing, please run self.export_posterior() first")

if use_n_obs is not None:
    # --- Stratified Sampling to ensure all categories are present ---
    print(f"Performing stratified sampling for {use_n_obs} observations...")
    # 1. Identify all categorical keys used by the model
    all_categorical_keys = [batch_key] + covariate_keys
    # 2. Get the integer positions of one representative for each category
    guaranteed_indices_pos = set()
    for key in all_categorical_keys:
        if key in adata_ref.obs.columns:
            # .drop_duplicates() is a fast way to get one index per category
            unique_indices_str = adata_ref.obs.drop_duplicates(subset=[key]).index
            # Convert string indices (e.g., 'barcode_AAAC...') to integer positions
            indices_pos = adata_ref.obs.index.get_indexer(unique_indices_str)
            guaranteed_indices_pos.update(indices_pos)
    n_guaranteed = len(guaranteed_indices_pos)
    print(f"  Found {n_guaranteed} guaranteed observations to cover all categories.")
    # 3. Calculate how many more random samples we need
    total_sample_size = np.min((use_n_obs, mod.adata.n_obs))
    n_to_sample_randomly = total_sample_size - n_guaranteed
    if n_to_sample_randomly > 0:
        # 4. Get the remaining indices (positions) that were not already selected
        all_indices_pos = set(range(mod.adata.n_obs))
        remaining_indices_pos = list(all_indices_pos - guaranteed_indices_pos)
        # 5. Randomly sample from the remaining indices
        random_indices_pos = np.random.choice(
            remaining_indices_pos, 
            n_to_sample_randomly, 
            replace=False
        )
        # 6. Combine guaranteed and random indices
        ind_x = np.array(list(guaranteed_indices_pos) + list(random_indices_pos))
    else:
        # The guaranteed set is already larger than use_n_obs
        # Just take a random sample from the guaranteed set
        print(f"  Guaranteed set ({n_guaranteed}) is larger than/equal to 'use_n_obs'. Subsampling from this set.")
        ind_x = np.random.choice(
            list(guaranteed_indices_pos),
            total_sample_size,
            replace=False
        )
    print(f"  Final sample size for QC: {len(ind_x)}")
else:
    ind_x = None

mod.expected_nb_param = mod.module.model.compute_expected(
    mod.samples[f"post_sample_{summary_name}"], mod.adata_manager, ind_x=ind_x
)
x_data = mod.adata_manager.get_from_registry(REGISTRY_KEYS.X_KEY)[ind_x, :]

if issparse(x_data):
    x_data = np.asarray(x_data.toarray())

mod.plot_posterior_mu_vs_data(mod.expected_nb_param["mu"], x_data)

plt.savefig(os.path.join(out_dir_plots, 'model_qc-1.pdf'))

plt.close()


# Plot 2: Estimated reference expression signatures comp. to average expression in cluster
# Should look like a perfect diagonal

inf_aver = mod.samples[f"post_sample_{summary_name}"]["per_cluster_mu_fg"].T
if scale_average_detection and ("detection_y_c" in list(mod.samples[f"post_sample_{summary_name}"].keys())):
    inf_aver = inf_aver * mod.samples[f"post_sample_{summary_name}"]["detection_y_c"].mean()

aver = mod._compute_cluster_averages(key=REGISTRY_KEYS.LABELS_KEY)
aver = aver[mod.factor_names_]

plt.hist2d(
    np.log10(aver.values.flatten() + 1),
    np.log10(inf_aver.flatten() + 1),
    bins=50,
    norm=matplotlib.colors.LogNorm(),
)
plt.xlabel("Mean expression for every gene in every cluster")
plt.ylabel("Estimated expression for every gene in every cluster")
plt.show()

plt.savefig(os.path.join(out_dir_plots, 'model_qc-2.pdf'))

plt.close()

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

# Save this as a table
inf_aver.to_csv(os.path.join(out_dir, 'reference_signatures', 'cell_type_signatures.csv'))

# --- Finish

print("Script succesfully completed")