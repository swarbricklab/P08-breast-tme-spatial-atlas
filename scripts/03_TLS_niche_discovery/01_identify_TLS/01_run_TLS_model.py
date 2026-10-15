#!/usr/bin/env python3

# --- Aim
# Apply the TLS scoring model from Andersson 2021 to spatial transcriptomics data (Visium).
# Output TLS scores per spot.
# Modified their original code (found here: https://github.com/almaan/her2st/blob/master/scripts/predict.py)


# --- Define environment
# env: spatial-pilot

# --- Load packages
import argparse
import anndata as ad
import os
import pandas as pd
import numpy as np


# --- Argparse arguments

parser = argparse.ArgumentParser(description="Import Visium data & do basic QC",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

# Required
parser.add_argument("--res_dir", help="Path to results directory")
parser.add_argument("--h5ad_dir", help="Path to directory for .h5ad output files")
parser.add_argument("--project_dir", help="Project directory (root)")
parser.add_argument("--sample_id", help="Sample ID")
parser.add_argument("--sample_sheet", help="Path to sample sheet CSV file")

args = vars(parser.parse_args())


# --- Arguments

res_dir = args["res_dir"]
h5ad_dir=args["h5ad_dir"]
sample_id=args["sample_id"]
project_dir=args["project_dir"]
res_dir = res_dir or os.path.join(project_dir, "results")


# --- Extract some info from the sample sheet

# Use sample_sheet argument if provided, else default
if args.get("sample_sheet"):
    sample_sheet_path = args["sample_sheet"]
else:
    sample_dir = os.path.join(project_dir, 'config')
    sample_sheet_path = os.path.join(sample_dir, 'sample_sheet_visium.csv')

sample_sheet = pd.read_csv(sample_sheet_path)

# Check if tissue_type is Lymph node
tissue_type = sample_sheet.loc[sample_sheet['section_id'] == sample_id, 'tissue_type'].values[0]
run_lymph_node = False
if tissue_type == 'Lymph node':
    run_lymph_node = True


# --- Load the adata object

adata = ad.read_h5ad(os.path.join(h5ad_dir, 'visium', f'{sample_id}.h5ad'))


# --- Load the coefficient file from andersson 2021
# Note: gene names in this coefficient file have been mapped to the Visium probe set names

gene_list_dir = os.path.join(project_dir, 'resources/gene_lists', 'andersson_2021')

coef_name = "andersson_2021_coef_full_visium.csv"

coef_path = os.path.join(gene_list_dir, coef_name)

coefs =  pd.read_csv(coef_path, index_col = 0, header = None, names = [0]) # names = [0] to make the name of the first column (coefficients) 0, necessary for downstream script

# Extract the relevant layer as count object 

cnt = ad.AnnData(adata.layers['log_norm'])
cnt.var_names = adata.var_names
cnt.obs_names = adata.obs_names

# Convert adata to count matrix

cnt = cnt.to_df()

# Other arguments

maxV = -np.inf
minV = np.inf
valS = []
crdS = []
threshold = None # In case setting a threshold is desired

to_score = {} # Save which genes are used in which samples

to_score['full gene list'] = list(coefs.index)

# --- Run their code to apply the model

cnt = pd.DataFrame(cnt.values / cnt.values.sum(axis = 1, keepdims = True), index = cnt.index, columns = cnt.columns,)
std = cnt.values.std(axis=0,keepdims = True)
cnt = pd.DataFrame(np.divide(cnt.values, std, where = std > 0), index = cnt.index, columns = cnt.columns,)

# Check if intercept is included in the coef file
if 'intercept' in coefs.index:
    print("intercept included in prediction")
    cnt['intercept'] = np.ones(cnt.shape[0])

cnt.columns = pd.Index([x.upper().replace(" ",'') for x in cnt.columns.values])
coefs.index = pd.Index([x.upper().replace(" ","") for x in coefs.index.values])
uniq_genes,gene_occ = np.unique(cnt.columns.values,return_counts =True)
cnt = cnt.loc[:,uniq_genes[gene_occ == 1]]

inter = pd.Index(set(coefs.index.intersection(cnt.columns)))
to_score[sample_id] = inter # Save which genes were used for future reference
print("{} / {} genes present in count data".format(inter.shape[0], coefs.shape[0],))

coefs = coefs.loc[inter,:]
cnt = cnt.loc[:,inter]
uu,cc = np.unique(cnt.columns.values,return_counts = True)
diff = coefs.index.difference(cnt.columns)
val = np.dot(cnt.values,coefs.values)
val = val.flatten()
print("Dynamic range of values : MIN({:0.4f}) , MAX({:0.4f})".format(val.min(),val.max()))

# Apply threshold (if there is one)
if threshold is not None:
    print("using cutoff : {}".format(threshold))
    val[val < threshold] = 0.0

# Copying because I got some kind of warning
val_new = val.copy()
adata.obs['TLS_model'] = val_new
print("Successfully processed : {}".format(sample_id))

print("Data read : COMPLETE")

# Save which genes were used per sample as a dataframe

to_score_df = pd.DataFrame.from_dict(to_score, orient='index').T


# --- Save the plots

# Output to results/03_TLS_niche_discovery/01_identify_TLS/01_run_TLS_model
base_results = os.path.join(res_dir, '03_TLS_niche_discovery', '01_identify_TLS', '01_run_TLS_model')

if run_lymph_node:
    base_results = os.path.join(base_results, 'lymph_nodes')

out_dir = base_results
os.makedirs(out_dir, exist_ok=True)

# Save the TLS model output as a dataframe
TLS_model_result = adata.obs['TLS_model']
filename = f'{sample_id}_TLS_model.csv'
TLS_model_result.to_csv(os.path.join(out_dir, filename), index=True)


print("Script succesfully completed")
