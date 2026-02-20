#!/bin/bash

# --- Aim: Run a single job for multiple samples (batch processing)

# print the date & time for the log file
now=$(date)
echo
echo "Start time: $now"

# arguments from submission file
function_script_file=$1
h5ad_dir=$2
project_dir=$3
environment=$4
res_dir=$5
visium_sample_sheet=$6
shift 6 # need to shift so that the first 6 arguments are gone and I get only the array of samples now
samples=$@

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
    --project_dir $project_dir \
    --res_dir $res_dir \
    --sample_sheet $visium_sample_sheet \
    --samples $samples \
