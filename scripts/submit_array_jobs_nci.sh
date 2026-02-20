#!/bin/bash

set -euxo pipefail

# --- Aim: Submit one job per sample for a given script (array jobs) on NCI HPC

while getopts s:p: flag
do
    case "${flag}" in
        s) file_name=${OPTARG};; # name of the script you want to run (inlcuding the path to it from the current directory)
    esac
done

echo "Run script $file_name"

# --- Load parameters

source "../config/params.sh"

# --- Paths to donwstream scripts

run_script_file="$cwd/run_single_sample_job.sh" # downstream bash script that will run the python script
function_script_file="$cwd/$file_name" # actual python script to be run


# --- Identify samples
# Identify all from the sample sheet
# Extract the section_id column
samples=$(awk -F',' 'NR > 1 {print $3}' $visium_sample_sheet) 

echo "Identified samples: $samples"


# --- Submit jobs

n_cores=4

# Use script name as name for job

s=$file_name
tmp=${s##*/}   # trim the directory path (split on the last '/' and take the last element)
tmp=${tmp##*_}   # also trim the digit at the start
name=${tmp%%.*}  # get rid of .py (split on '.' and take first element)

for s in ${samples[@]}
do 
    cd $project_dir
    echo $s 

    qsub -N "${name}_${s}" \
        -S /bin/bash \
        -P a56 \
        -q normal \
        -o $log_dir/$s.out \
        -e $log_dir/$s.err \
        -m a \
        -l walltime=01:00:00 \
        -l ncpus=$n_cores \
        -l mem=100g \
        -l storage=gdata/a56+gdata/px14+scratch/px14 \
        -V \
        -- $run_script_file $function_script_file $h5ad_dir $project_dir $s $environment $min_gene $min_count $min_spots $res_dir $visium_sample_sheet
        
done
