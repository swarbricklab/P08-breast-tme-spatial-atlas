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

args = vars(parser.parse_args())

# Arguments
project_dir = args["project_dir"]
h5ad_dir = args["h5ad_dir"]
res_dir = args["res_dir"]
sample_list = args["samples"]


# --- Load the correlation results

gene_of_i = 'TMEM173'  # STING gene symbol

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

# Decide whether to filter for endothelial ACKR1 only
filter_endothelial_ACKR1 = False

correlation_dir = os.path.join(res_dir, 'TLS_analysis', 'HEV_regions', 'all_samples', 'gene_gene_correlations', f'HEV_{gene_of_i}_correlations')

if filter_endothelial_ACKR1:
    correlation_dir = os.path.join(correlation_dir, 'filter_endothelial_ACKR1')

most_recent = find_most_recent(correlation_dir)

correlation_file = f'All_samples_HEV_{gene_of_i}_gene_correlations.csv'
correlation_df = pd.read_csv(os.path.join(most_recent, 'csv_files', correlation_file), index_col=0)


# --- Prepare for saving

# Decide which gene set to use
gene_set = 'GO_Molecular_Function_2025'

out_dir = os.path.join(res_dir, 'TLS_analysis', 'HEV_regions', 'all_samples', 'gene_gene_correlations', f'HEV_{gene_of_i}_correlations',)

if filter_endothelial_ACKR1:
    out_dir = os.path.join(out_dir, 'filter_endothelial_ACKR1')

out_dir = os.path.join(out_dir, 'enrichment_analysis', gene_set)

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

# Do a few cutoffs
cutoffs = [50, 100, 200, 500, 1000]

# Upregulated genes
enr_up_dict = {}

for cutoff in cutoffs:
    enr_up_dict[cutoff] = {}
    # Get top genes
    genes_up = significant_corr_df[significant_corr_df['correlation'] > 0].head(cutoff).reset_index()
    if not genes_up.empty:
        # Enricr API
        enr_up = gp.enrichr(genes_up['gene'],
                            gene_sets=go_bio_filtered,
                            outdir=None)
        enr_up.res2d.Term = enr_up.res2d.Term.str.split(" \(GO").str[0]
        enr_up_dict[cutoff] = enr_up


# --- Plot 1 - gsea dotplots & gsea barplots

# Define cutoff for considering significance
cutoff = 0.05

for layer in enr_up_dict:
    enr_up = enr_up_dict[layer]
    # Check if there are any terms with adjust p value < 0.05
    if sum(enr_up.res2d['Adjusted P-value'] < cutoff) > 0:
        # Make dotplot
            filename = f'All_samples_{layer}_gene_enrichment_dotplot.pdf'
            # Make sure there are no slashes or spaces in the filename
            filename = filename.replace('/', '_')
            filename = filename.replace(' ', '_')
            gp.dotplot(enr_up.res2d, figsize=(3,5), title=f'Top {layer} Gene Enrichment', show_ring=True, cmap=plt.cm.viridis_r, ofname=os.path.join(out_dir_plots, filename))
            # Make barplot
            filename = f'All_samples_top_{layer}_gene_enrichment_barplot.pdf'
            # Make sure there are no slashes or spaces in the filename
            filename = filename.replace('/', '_')
            filename = filename.replace(' ', '_')
            gp.barplot(enr_up.res2d, title=f'Top {layer} Gene Enrichment', color = 'cornflowerblue', ofname=os.path.join(out_dir_plots, filename))
    else:
        print(f"No significant enriched terms found for top {layer} genes.")


# --- Save the dataframes - can use these later to make more plots

for layer in enr_up_dict:
    df = enr_up_dict[layer].res2d
    # Sort by Adjusted P-value
    df = df.sort_values(by='Adjusted P-value', ascending=True)
    filename = f'All_samples_top_{layer}_gene_enrichment.csv'
    # Make sure there are no slashes or spaces in the filename
    filename = filename.replace('/', '_')
    filename = filename.replace(' ', '_')
    df.to_csv(os.path.join(out_dir_csv, filename), index=False)


# --- Also run GSEA preranked analysis

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

# Prepare output directory for preranked plots
out_dir_ranked_plots = os.path.join(out_dir_plots, 'preranked_plots')
os.makedirs(out_dir_ranked_plots, exist_ok=True)

# Check if you have any significant terms
if not top_positive_terms.empty:
    # Create the bar plot for the top 10 terms
    gp.barplot(
        top_positive_terms.head(10),
        title='Top 10 Enriched GO Terms (Positive Correlation)',
        column='FDR q-val',
        cutoff=0.25,
        ofname=os.path.join(out_dir_ranked_plots, 'Top_10_positive_enriched_terms_barplot.pdf')
    )
else:
    print("No significant positively enriched terms found.")



# --- Plot 2: Plot all top positive terms

plt.close('all')  # Close previous plots

top_terms = top_positive_terms['Term'].tolist()

if len(top_terms) > 0:
    if len(top_terms) > 20:
        top_terms = top_terms[:20]  # Limit to top 20 terms
    for n, top_term_name in enumerate(top_terms):
        print(f"Plotting enrichment for: {top_term_name}")
        # Plot the running score for that specific term
        # 'term' needs the exact name from the results index
        gp.gseaplot(
            rank_metric=pre_res.ranking,      # The ranked list you used
            term=top_term_name,          # The term you want to plot
            **pre_res.results[top_term_name]  # This passes the term's specific stats
        )
        filename = os.path.join(out_dir_ranked_plots, 'All_top_positive_enriched_terms_gseaplots.pdf')
        # Save all in a single figure
        pp = PdfPages(filename)
        fig_nums = plt.get_fignums()
        figs = [plt.figure(n) for n in fig_nums]
        for fig in figs:
            fig.savefig(pp, format='pdf', bbox_inches='tight')
        pp.close()

plt.close('all')  # Close previous plots


# --- If molecular function, also specifically plot terms related to adhesion

if gene_set == 'GO_Molecular_Function_2025':
    adhesion_terms = [term for term in results_df['Term'] if 'adhesion' in term.lower()]
    for n, adhesion_term_name in enumerate(adhesion_terms):
        print(f"Plotting enrichment for adhesion term: {adhesion_term_name}")
        # Plot the running score for that specific term
        gp.gseaplot(
            rank_metric=pre_res.ranking,      # The ranked list you used
            term=adhesion_term_name,          # The term you want to plot
            **pre_res.results[adhesion_term_name]  # This passes the term's specific stats
        )
    filename = os.path.join(out_dir_ranked_plots, 'Adhesion_related_enriched_terms_gseaplots.pdf')
    # Save all in a single figure
    pp = PdfPages(filename)
    fig_nums = plt.get_fignums()
    figs = [plt.figure(n) for n in fig_nums]
    for fig in figs:
        fig.savefig(pp, format='pdf', bbox_inches='tight')
    pp.close()

plt.close('all')  # Close previous plots


# --- Also plot specific terms of interest

keywords = ['cytokine']

terms_of_interest = [term for term in results_df['Term'] if any(keyword in term.lower() for keyword in keywords)]

# Make sure not to get anything with cytokinesis
terms_of_interest = [term for term in terms_of_interest if 'cytokinesis' not in term.lower() and 'cytokinetic' not in term.lower()]

for n, term_name in enumerate(terms_of_interest):
    print(f"Plotting enrichment for term of interest: {term_name}")
    # Plot the running score for that specific term
    gp.gseaplot(
        rank_metric=pre_res.ranking,      # The ranked list you used
        term=term_name,          # The term you want to plot
        **pre_res.results[term_name]  # This passes the term's specific stats
    )

filename = os.path.join(out_dir_ranked_plots, 'Cytokine_related_enriched_terms_gseaplots.pdf')

# Save all in a single figure
pp = PdfPages(filename)
fig_nums = plt.get_fignums()

figs = [plt.figure(n) for n in fig_nums]

for fig in figs:
    fig.savefig(pp, format='pdf', bbox_inches='tight')

pp.close()

plt.close('all')  # Close previous plots


# --- Same for chemokine related terms

keywords = ['chemokine']

terms_of_interest = [term for term in results_df['Term'] if any(keyword in term.lower() for keyword in keywords)]

for n, term_name in enumerate(terms_of_interest):
    print(f"Plotting enrichment for term of interest: {term_name}")
    # Plot the running score for that specific term
    gp.gseaplot(
        rank_metric=pre_res.ranking,      # The ranked list you used
        term=term_name,          # The term you want to plot
        **pre_res.results[term_name]  # This passes the term's specific stats
    )

filename = os.path.join(out_dir_ranked_plots, 'Chemokine_related_enriched_terms_gseaplots.pdf')

# Save all in a single figure
pp = PdfPages(filename)
fig_nums = plt.get_fignums()

figs = [plt.figure(n) for n in fig_nums]

for fig in figs:
    fig.savefig(pp, format='pdf', bbox_inches='tight')

pp.close()

plt.close('all')  # Close previous plots


# --- Also save results as csv

# Sort results by NES
results_df = results_df.sort_values(by='NES', ascending=False)
results_df.to_csv(os.path.join(out_dir_csv, f'All_samples_HEV_{gene_of_i}_preranked_GSEA_results.csv'), index=False)


# --- Finish

print("Script successfully completed")
