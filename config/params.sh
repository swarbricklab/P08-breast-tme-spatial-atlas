# Parameters for submission scripts (submit_wolfpack.sh)

# --- Directories

cwd=$(pwd)

# Paths to downstream scripts

# Get the project directory (root)
# Use dirname to get one directory up (assuming the scripts folder is one down from the root)
project_dir=$(dirname "$cwd")

res_dir="$project_dir/results"
h5ad_dir="$project_dir/data/anndata_objects" # Path to the h5ad files (anndata objects)

log_dir="$project_dir/logs"
mkdir -p $log_dir

# Environment to use
environment="spatial-pilot"

# Visium sample sheet: contains the names of the visium samples
config_dir="$project_dir/config"
visium_sample_sheet="$config_dir/sample_sheet_visium.csv"