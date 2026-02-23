#!/usr/bin/env python3

# --- Aim
# Perform gene set enrichment analysis (GSEA) on genes correlated with SPP1 or TMEM173 (STING) expression in HEV niches.
# Identify pathways and biological processes enriched among correlated genes.

# --- Define environment
# env: gseapy

# --- Load packages

import pandas as pd
import os
import gseapy as gp

from datetime import date
from datetime import datetime

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

import argparse

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


# --- Load the correlation results

gene_of_i = 'SPP1'  # STING gene symbol

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

correlation_dir = os.path.join(res_dir, '05_xenium_mechanisms', '05_gene_gene_correlations', '01_investigate_gene_correlations', f'{gene_of_i}_correlations')
most_recent = find_most_recent(correlation_dir)

correlation_file = f'All_samples_HEV_{gene_of_i}_gene_correlations.csv'
correlation_df = pd.read_csv(os.path.join(most_recent, 'csv_files', correlation_file), index_col=0)


# --- Prepare for saving

# Decide which gene set to use
gene_set = 'GO_Biological_Process_2025'
out_dir = os.path.join(res_dir, '05_xenium_mechanisms', '05_gene_gene_correlations', '02_GSEA_analysis', f'{gene_of_i}_correlations')

out_dir = os.path.join(out_dir, gene_set)

today = date.today()
today = today.strftime("%Y%m%d")

out_dir_plots = os.path.join(out_dir, 'plots', today)
os.makedirs(out_dir_plots, exist_ok=True)

out_dir_csv = os.path.join(out_dir, 'csv_files', today)
os.makedirs(out_dir_csv, exist_ok=True)


# --- Filter based on pval & subset up or down-regulated genes

# Extract significant genes
pval_threshold = 0.05
significant_corr_df = correlation_df[correlation_df['p_value'] < pval_threshold].copy()


# --- Go term analysis

# Load from file (server connection unstable)
# Downloaded from https://maayanlab.cloud/Enrichr/#libraries
gene_set_dir = os.path.join(project_dir, 'resources', 'gene_lists', 'enrichment_analysis')

# Find the corresponding file
file_list = os.listdir(gene_set_dir)
filename = [f for f in file_list if gene_set.lower() in f.lower() and f.endswith('.gmt')][0]

go_bio = gp.read_gmt(os.path.join(gene_set_dir, filename))

# Filter for genes that are actually in the panel
all_genes = correlation_df.index.tolist()

go_bio_filtered = {}

for term in go_bio:
    genes_in_term = [gene for gene in go_bio[term] if gene in all_genes]
    if len(genes_in_term) >= 5:  # Keep only terms with at least 2 genes present
        go_bio_filtered[term] = genes_in_term


# --- Run GSEA preranked analysis

ranked_gene_list = correlation_df['correlation']

pre_res = gp.prerank(
    rnk=ranked_gene_list,
    gene_sets=go_bio_filtered,  # Use your filtered dictionary
    min_size=5,                 # gseapy will apply this *after* loading
    max_size=500,
)


# --- Examine preranked results

results_df = pre_res.res2d

# Filter for top terms at the 'beginning' of your list
# These have a positive NES and a significant FDR (e.g., < 0.25 or < 0.05)
top_positive_terms = results_df[
    (results_df['NES'] > 0) & 
    (results_df['FDR q-val'] < 0.25)
]

print("Top positively enriched terms:")
print(top_positive_terms[['Term', 'NES', 'FDR q-val']])

# Sort by NES to see the most enriched first
top_positive_terms = top_positive_terms.sort_values(by='NES', ascending=False)

# Check if you have any significant terms
if not top_positive_terms.empty:
    # Create the bar plot for the top 10 terms
    gp.barplot(
        top_positive_terms.head(10),
        title='Top 10 Enriched GO Terms (Positive Correlation)',
        column='FDR q-val',
        cutoff=0.25,
        ofname=os.path.join(out_dir_plots, 'Top_10_positive_enriched_terms_barplot.pdf')
    )
else:
    print("No significant positively enriched terms found.")


# --- Also save results as csv

# Sort results by NES
results_df = results_df.sort_values(by='NES', ascending=False)
results_df.to_csv(os.path.join(out_dir_csv, f'All_samples_HEV_{gene_of_i}_preranked_GSEA_results.csv'), index=False)


# --- Finish

print("Script successfully completed")
