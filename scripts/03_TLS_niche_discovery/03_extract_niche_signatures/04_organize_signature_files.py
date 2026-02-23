#!/usr/bin/env python3

# --- Aim
# For convenience, copy the filtered gene lists to a new directory for easier access in downstream analyses.
# Note that the files are filtered to the top 300 genes (if there are that many) to make them more manageable for downstream analyses. The original, unfiltered files are still available in the respective subdirectories if needed.

# --- Define environment
# env: spatial-pilot

# --- Load packages
import argparse
import os
from datetime import datetime

import pandas as pd


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

print(f'Samples in python: {sample_list}')

# --- Prepare for copying

# Define source and destination directories
source_dir = os.path.join(res_dir, '03_TLS_niche_discovery', '03_extract_niche_signatures',)
destination_dir = os.path.join(res_dir, '03_TLS_niche_discovery', '03_extract_niche_signatures', '04_final_niche_signature_files')

# Create destination directory if it doesn't exist
os.makedirs(destination_dir, exist_ok=True)


# --- Identify files

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

def copy_top_300_csv(source_file, destination_file):
    if os.path.exists(source_file):
        try:
            df = pd.read_csv(source_file)
            df_top = df.head(300)
            df_top.to_csv(destination_file, index=False)
            print(f'Copied top 300 rows of {os.path.basename(source_file)} to {destination_file}')
        except Exception as e:
            print(f'Error processing {source_file}: {e}')
    else:
        print(f'File {os.path.basename(source_file)} not found in {os.path.dirname(source_file)}')

# First copy T cell zone & B cell zone signatures - these didn't get additional / custom filtering
sig_dir = os.path.join(source_dir, '01_extract_niche_signatures')
files = ['b_cell_zone_top_genes_filtered.csv', 't_cell_zone_top_genes_filtered.csv']
for file in files:
    source_file = os.path.join(sig_dir, file)
    destination_file = os.path.join(destination_dir, file)
    copy_top_300_csv(source_file, destination_file)

# Now copy gc post invasive epithelial gene filtering
sig_dir = os.path.join(source_dir, '02_filter_GC_signature')
filename = 'germinal_center_filtered_epithelial_genes.csv'
source_file = os.path.join(sig_dir, filename)
destination_file = os.path.join(destination_dir, filename)
copy_top_300_csv(source_file, destination_file)

# Now copy hev signature post general endothelial filtering
sig_dir = os.path.join(source_dir, '03_filter_HEV_signature', 'filter_HEV_signature')
most_recent_dir = find_most_recent(sig_dir)
filename = 'hevs_filtered_endothelial_genes.csv'
if most_recent_dir:
    source_file = os.path.join(most_recent_dir, 'csv_files', filename)
    destination_file = os.path.join(destination_dir, filename)
    copy_top_300_csv(source_file, destination_file)


# --- Finish

print("Script successfully completed")
