#!/usr/bin/env python3

# --- Aim
# Investigate gene expression correlations with SPP1 or TMEM173 (STING) expression in HEV niches.
# Identify gene-gene relationships.

# --- Define environment
# env: xen-hev

# --- Load packages
import argparse
import os
import sys
from datetime import datetime

import anndata as ad

import matplotlib as mpl
import matplotlib.pyplot as plt

from scipy.stats import pearsonr

import scanpy as sc

import pandas as pd
import numpy as np

import seaborn as sns

mpl.rcParams['pdf.fonttype'] = 42


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

# Import custom functions
sys.path.append(os.path.join(project_dir, 'scripts'))
from helper_functions import save_multi_image


# --- Extract some info from the sample sheet

# Use sample_sheet argument if provided, else default
if args.get("sample_sheet"):
    sample_sheet_path = args["sample_sheet"]
else:
    sample_dir = os.path.join(project_dir, 'config')
    sample_sheet_file = 'sample_sheet_xenium.csv'
    sample_sheet_path = os.path.join(sample_dir, sample_sheet_file)

sample_sheet = pd.read_csv(sample_sheet_path)

# Also load tissue metadata
metadata_dir = os.path.join(project_dir, 'resources',)
tissue_metadata_file = 'tissue_metadata.csv'
tissue_metadata = pd.read_csv(os.path.join(metadata_dir, tissue_metadata_file))

# Make a section_id to subtype mapping
section_to_subtype = tissue_metadata.set_index('tissue_id')['clinical_subtype'].to_dict()
# Map section_id to subtype
sample_sheet['subtype'] = sample_sheet['tissue_id'].map(section_to_subtype)

subtype_dict = sample_sheet.set_index('section_id')['subtype'].to_dict()


# --- Load the anndata objects

adata_dict = {}

for sample in sample_list:
    adata = ad.read_h5ad(os.path.join(h5ad_dir, 'xenium', f'{sample}.h5ad'))
    adata_dict[sample] = adata


# --- Load the HEV annotations

# Find most recent results
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
        return most_recent_dir_path
    else:
        print("No valid directories found.")

distance_dir = os.path.join(res_dir, '05_xenium_mechanisms', '02_identify_HEV_regions',)

most_recent_distance = find_most_recent(distance_dir)

for sample in sample_list:
    # Also load HEV niche barcode and add to adata
    filename = f'{sample}_HEV_niches_cell_to_niche.csv'
    niche_barcode_df = pd.read_csv(os.path.join(most_recent_distance, 'csv_files', filename), index_col=0)
    niche_barcode_df.dropna(subset=['niche_id'], inplace=True)
    # Prepend sample name to niche_id
    niche_barcode_df['niche_id'] = sample + '_' + niche_barcode_df['niche_id'].astype(str)
    adata_dict[sample].obs['niche_id'] = niche_barcode_df['niche_id']
    # Also load the distance calculations
    distance_file = f'{sample}_HEV_niches_info.csv'
    distance_df = pd.read_csv(os.path.join(most_recent_distance, 'csv_files', distance_file))
    # Prepend sample name to niche_id
    distance_df['niche_id'] = sample + '_' + distance_df['niche_id'].astype(str)
    # Make a mapping dict of distance to TLS
    distance_dict = distance_df.set_index('niche_id')['close_to_TLS'].to_dict()
    # Map to adata
    adata_dict[sample].obs['close_to_TLS'] = adata_dict[sample].obs['niche_id'].map(distance_dict)


# --- Look at gene expression in HEV niches

adata_dict_HEV = {}

# Subset adatas to HEV only
for sample in sample_list:
    adata = adata_dict[sample]
    adata_HEV = adata[adata.obs['niche_id'].notna()].copy()
    adata_dict_HEV[sample] = adata_HEV

# Extract a single adata object with all samples
adata_list = list(adata_dict_HEV.values())
all_adata = ad.concat(adata_list, merge='same')


# --- Create an adata with mean expression per niche

mean_adata = all_adata.to_df().groupby(all_adata.obs['niche_id']).mean()

mean_adata = ad.AnnData(mean_adata)

mean_adata.obs['niche_id'] = mean_adata.obs.index


# --- Calculate correlations between STING expression and all other genes

# Pick a gene of interest
gene_of_i = 'TMEM173'  # STING gene symbol

correlation_df = pd.DataFrame(columns=['gene', 'correlation', 'p_value'])

for gene in mean_adata.var_names:
    if gene == gene_of_i:
        continue
    # Calculate Pearson correlation
    corr, p_val = pearsonr(mean_adata[:, gene_of_i].X.flatten(), mean_adata[:, gene].X.flatten())
    correlation_df = correlation_df._append({'gene': gene, 'correlation': corr, 'p_value': p_val}, ignore_index=True)

# Sort by correlation
correlation_df = correlation_df.sort_values(by='correlation', ascending=False).reset_index(drop=True)

# Drop any NaN values
correlation_df = correlation_df.dropna(subset=['correlation', 'p_value'])
correlation_df = correlation_df.set_index('gene')

# Add rank
correlation_df['rank'] = range(1, len(correlation_df) + 1)

# Filter significant correlations
pval_threshold = 0.05
significant_corr_df = correlation_df[correlation_df['p_value'] < pval_threshold].copy()
correlation_df['significant'] = correlation_df['p_value'] < pval_threshold


# --- Plot 1: Plot the correlations

# Add grid behind
plt.rcParams['axes.axisbelow'] = True

fig, ax = plt.subplots(figsize=(8, 4))

sns.scatterplot(data=correlation_df, x='rank', y='correlation', color='lightgrey', edgecolor=None, s=2, ax=ax)
sns.scatterplot(data=significant_corr_df, x='rank', y='correlation', color='red', edgecolor=None, s=2, ax=ax)

# Annotate some genes
annotate_genes = ['SELE', 'SELP', 'CXCL10', 'CXCL9', 'CCL5', 'CCL2', 'CCL19', 'CCL21', 'CXCL13',
                    'MCAM', 'ICAM1', 'ICAM2', 'VCAM1', 'BCAM', 'PECAM1']

annotated_positions = []
min_dist = 15  # minimum pixel distance between annotations

for gene in annotate_genes:
    if gene not in correlation_df.index:
        continue
    x = correlation_df.loc[gene, 'rank']
    y = correlation_df.loc[gene, 'correlation']
    xytext = np.array([5, 5])
    # Check for overlap with previous annotations
    for prev_x, prev_y, prev_xytext in annotated_positions:
        # Calculate screen position of previous annotation
        prev_screen = ax.transData.transform((prev_x, prev_y)) + prev_xytext
        curr_screen = ax.transData.transform((x, y)) + xytext
        dist = np.linalg.norm(curr_screen - prev_screen)
        if dist < min_dist:
            # If too close, shift the annotation further
            xytext += np.array([min_dist, min_dist])
    ax.annotate(
        gene,
        xy=(x, y),
        xytext=tuple(xytext),
        textcoords='offset points',
        fontsize=8,
        arrowprops=dict(arrowstyle='->', color='black', lw=0.7)
    )
    annotated_positions.append((x, y, xytext))

# Add grid behind plot
ax.grid(which='major', linestyle='-', color='lightgrey')

# Also add a line at the middle rank
mid_rank = len(correlation_df) / 2
ax.axvline(mid_rank, color='lightgray', lw=3, ls='-', zorder=-1)

# Adjust labels
ax.set_xlabel(f'Gene rank (by correlation with {gene_of_i} expression)',)
ax.set_ylabel('Pearson correlation coefficient', )

# Add a legend for pvalue
ax.legend(title=f'p<{pval_threshold}', labels=['False', 'True'], loc='center left',  frameon=False, markerscale=3, bbox_to_anchor=(1, 0.5))

# Add a title
ax.set_title(f'Correlation mean {gene_of_i} expression & other genes per HEV', fontweight='bold')

plt.suptitle(f'All samples - Correlation between mean {gene_of_i} expression and other genes per HEV niche', fontweight='bold', y=1.02)


# --- Save the output

out_dir = os.path.join(res_dir, '05_xenium_mechanisms', '05_gene_gene_correlations', '01_investigate_gene_correlations', f'{gene_of_i}_correlations')

today = datetime.now().strftime('%Y%m%d')
out_dir_plots = os.path.join(out_dir, today, 'plots')
os.makedirs(out_dir_plots, exist_ok=True)

# Save the figures
filename = f'All_samples_HEV_{gene_of_i}_gene_correlations.pdf'
save_multi_image(os.path.join(out_dir_plots, filename))

# Also save the correlation data
out_dir_csv = os.path.join(out_dir, today, 'csv_files')
os.makedirs(out_dir_csv, exist_ok=True)

correlation_df.to_csv(os.path.join(out_dir_csv, f'All_samples_HEV_{gene_of_i}_gene_correlations.csv'))


# --- Finish

print("Script successfully completed")