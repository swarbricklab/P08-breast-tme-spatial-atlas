#!/usr/bin/env python3

# --- Aim
# Plot Figure 4D (part 2): GSEA pathway enrichment for STING-correlated genes
# Shows biological pathways enriched in genes that correlate with STING in HEV niches

# --- Define environment
# env: gseapy

# --- Load packages

import os
import argparse

import pandas as pd
import gseapy as gp

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

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


# --- Load the gene set file

# Load from file (server connection unstable)
# Downloaded from https://maayanlab.cloud/Enrichr/#libraries
gene_set_dir = os.path.join(project_dir, 'resources', 'gene_lists', 'enrichment_analysis')

# Find the corresponding file
file_list = os.listdir(gene_set_dir)
gene_set = 'GO_Molecular_Function_2025'
filename = [f for f in file_list if gene_set.lower() in f.lower() and f.endswith('.gmt')][0]

go_bio = gp.read_gmt(os.path.join(gene_set_dir, filename))


# --- Load the correlation results

gene_of_i = 'TMEM173'  # STING gene symbol

correlation_dir = os.path.join(res_dir, '05_xenium_mechanisms', '05_gene_gene_correlations', '02_GSEA_analysis', f'{gene_of_i}_correlations', gene_set, 'csv_files')
most_recent = correlation_dir

correlation_file = f'All_samples_HEV_{gene_of_i}_preranked_GSEA_results.csv'
results_df = pd.read_csv(os.path.join(most_recent, correlation_file), index_col=0)

# Also load the ranked list
correlation_dir = os.path.join(res_dir, '05_xenium_mechanisms', '05_gene_gene_correlations', '01_investigate_gene_correlations', f'{gene_of_i}_correlations')
most_recent = correlation_dir

correlation_file = f'All_samples_HEV_{gene_of_i}_gene_correlations.csv'
correlation_df = pd.read_csv(os.path.join(most_recent, 'csv_files', correlation_file), index_col=0)

ranked_list = correlation_df

# Filter for genes that are actually in the panel
all_genes = correlation_df.index.tolist()

go_bio_filtered = {}

for term in go_bio:
    genes_in_term = [gene for gene in go_bio[term] if gene in all_genes]
    if len(genes_in_term) >= 5:  # Keep only terms with at least 5 genes present
        go_bio_filtered[term] = genes_in_term


# --- Examine preranked results

terms_of_interest = ['cell adhesion mediator', 'chemokine activity']

matching_terms = [
    term for term in results_df['Term']
    if any(term.lower().startswith(substring.lower()) for substring in terms_of_interest)
]


# --- Plot the running score for the top terms

# Prepare the ranked list for gseapy (needs to be a pd.Series with gene names as index and correlation as values)
rnk = ranked_list.sort_values(ascending=False, by='correlation')['correlation']

# Prepare for saving in a single file
out_dir = os.path.join(res_dir, 'figures')
os.makedirs(out_dir, exist_ok=True)

filename = os.path.join(out_dir, 'Fig_4D-2.pdf')

with PdfPages(filename) as pdf:

    for term_name in matching_terms:

        print(f"Plotting enrichment for: {term_name}")

        # Get the list of genes in the ranked list (in the correct order)
        term_genes = go_bio_filtered[term_name]

        # Run a "Dry Run" GSEA to get the curve data (RES)
        # We pass the single term as a dictionary so we don't need a GMT file
        # permutation_num=0 makes this run in milliseconds
        pre_res = gp.prerank(
            rnk=rnk,
            gene_sets={term_name: term_genes}, 
            permutation_num=0, 
            seed=42,
            min_size=5,                 # gseapy will apply this *after* loading
            max_size=500,
        )

        # Manually overwrite the stats in this dictionary
        # Inject the saved values 'results_df' into the dictionary
        plot_data = pre_res.results[term_name]
        plot_data['nes'] = results_df.loc[results_df['Term'] == term_name, 'NES'].values[0]
        plot_data['pval'] = results_df.loc[results_df['Term'] == term_name, 'NOM p-val'].values[0]
        plot_data['fdr'] = results_df.loc[results_df['Term'] == term_name, 'FDR q-val'].values[0]

        # Now we can plot the curve for this term using the gseaplot function, which will use the injected stats for the title and annotations
        axes = gp.gseaplot(
            rank_metric=pre_res.ranking, # The sorted ranking
            term=term_name,
            ofname=None,
            **plot_data # Unpacks 'hits', 'RES', 'nes', etc. automatically
        )

        # Capture the current figure and save to PDF
        # gseaplot creates a matplotlib figure; we grab it with gcf()
        fig = axes[0].figure
        pdf.savefig(fig, bbox_inches='tight')
        
        # Close the figure to free memory
        plt.close(fig)

print("PDF saved successfully.")


# --- Finish

print("Script successfully completed")
