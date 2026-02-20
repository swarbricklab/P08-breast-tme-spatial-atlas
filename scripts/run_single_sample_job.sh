#!/bin/bash

# --- Aim: Run a job for a single sample (per-sample processing)

# print the date & time for the log file
now=$(date)
echo
echo "Start time: $now"

# arguments from submission file
function_script_file=$1
h5ad_dir=$2
project_dir=$3
sample_id=$4
environment=$5
min_gene=$6
min_count=$7
min_spots=$8
res_dir=$9
visium_sample_sheet=${10}

# Initialize Conda environment
if [ -f /g/data/a56/conda/miniconda3/etc/profile.d/conda.sh ]; then
  source /g/data/a56/conda/miniconda3/etc/profile.d/conda.sh
fi

# load environment
conda activate $environment
echo "Run script $function_script_file"

# script file 
$function_script_file \
  --h5ad_dir $h5ad_dir \
  --sample_id $sample_id \
  --project_dir $project_dir \
  --min_gene $min_gene \
  --min_count $min_count \
  --min_spots $min_spots \
  --res_dir $res_dir \
  --sample_sheet $visium_sample_sheet \
