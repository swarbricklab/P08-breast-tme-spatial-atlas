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
res_dir=$6
sample_sheet=$7

# Initialize Conda environment
if [ -f /g/data/a56/conda/miniconda/etc/profile.d/conda.sh ]; then
  source /g/data/a56/conda/miniconda/etc/profile.d/conda.sh
fi

# load environment
conda activate $environment
echo "Run script $function_script_file"

# script file 
$function_script_file \
  --h5ad_dir $h5ad_dir \
  --sample_id $sample_id \
  --project_dir $project_dir \
  --res_dir $res_dir \
  --sample_sheet $sample_sheet \
