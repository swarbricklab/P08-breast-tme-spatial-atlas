# Analysis Scripts

This directory contains all computational analysis scripts for the breast cancer TLS / HEV project, organized into sequential analysis workflows and figure generation scripts.

## Directory Structure

### Analysis Workflows

The main analysis is organized into five sequential steps, each containing multiple Python scripts:

- **[01_preprocessing/](01_preprocessing/)** - Data import, quality control, and preprocessing for Visium, Visium HD, and Xenium datasets. Generates AnnData objects for downstream analysis.

- **[02_cell_type_mapping/](02_cell_type_mapping/)** - Cell type deconvolution and annotation workflows. Maps single-cell reference signatures onto spatial data to infer cell type composition at each spatial location.

- **[03_TLS_niche_discovery/](03_TLS_niche_discovery/)** - Identification and characterization of tertiary lymphoid structure (TLS) niches in Visium data. Includes TLS scoring, spatial segmentation, internal architecture mapping (T cell zone, B cell zone, germinal centers, HEV), and extraction of niche-specific gene signatures.

- **[04_vascular_microenvironment/](04_vascular_microenvironment/)** - High-resolution analysis of vascular niches using Visium HD (16 µm resolution). Includes HEV identification and segmentation, vessel classification (arterial, capillary, venous, lymphatic), comparative analysis of HEV vs. canonical venous structures, and ligand-receptor interaction analysis.

- **[05_xenium_mechanisms/](05_xenium_mechanisms/)** - Single-cell resolution analysis of immune and vascular mechanisms using Xenium data. Includes TLS and HEV identification, distance-based cellular enrichment analysis, pathway enrichment, and investigation of specific signaling axes.

Each subdirectory contains its own detailed README documenting the specific workflow steps and individual script functions.

### Figure Generation

- **[figures/](figures/)** - Scripts for generating manuscript figure panels, organized by figure number (figure1/ through figure5/). Each script is named according to its corresponding panel (e.g., `3B.py` generates Figure 3B). These scripts read preprocessed data from the analysis workflows and produce publication-ready plots.

### Utility Scripts

- **`helper_functions.py`** - Shared utility functions used across multiple analysis scripts.

## Job Submission Scripts

The repository includes bash scripts for running analyses on NCI Gadi HPC using the PBS job scheduler:

### Array Job Submission (One Job Per Sample)

Use these scripts when each sample needs to be processed independently in parallel:

- **`submit_array_jobs_nci.sh`** - Main submission script that creates one PBS job per sample
  ```bash
  ./submit_array_jobs_nci.sh -s <path/to/script.py>
  ```
  - Reads sample IDs from the configured sample sheet
  - Submits parallel jobs for each sample with appropriate resource allocation
  - Each job processes a single sample via `run_single_sample_job.sh`

- **`submit_array_jobs_nci_gpu.sh`** - GPU-enabled array job submission for computationally intensive tasks
  ```bash
  ./submit_array_jobs_nci_gpu.sh -s <path/to/script.py>
  ```
  - Same as `submit_array_jobs_nci.sh` but requests GPU resources (gpuvolta queue, 1 GPU per job)
  - Use for scripts requiring GPU acceleration (e.g., deconvolution with cell2location, label transfer with tangram)
  - Configured with 12 CPUs, 64GB RAM, 4-hour walltime per default

- **`run_single_sample_job.sh`** - Wrapper script executed by each array job
  - Activates the specified conda environment
  - Runs the Python script for a single sample
  - Passes required arguments (project directory, sample ID, paths, etc.)

### Batch Job Submission (Multiple Samples in One Job)

Use these scripts when all samples should be processed together (e.g., for cross-sample integration):

- **`submit_multi_sample_job_nci.sh`** - Submits a single PBS job for all samples
  ```bash
  ./submit_multi_sample_job_nci.sh -s <path/to/script.py>
  ```
  - Collects all sample IDs from the sample sheet
  - Submits one job that processes all samples sequentially or in batch
  - Calls `run_multi_sample_job.sh` with the full sample list

- **`run_multi_sample_job.sh`** - Wrapper script for batch processing
  - Activates the conda environment
  - Executes the Python script with all samples passed as arguments
  - Useful for integration, cross-sample normalization, or comparative analyses

**Note:** You can determine whether a script is designed for single-sample or multi-sample processing by checking its argparse arguments:
- **Single-sample scripts** accept `--sample_id` (use with `submit_array_jobs_nci.sh`)
- **Multi-sample scripts** accept `--samples` (use with `submit_multi_sample_job_nci.sh`)

### Usage Notes

- Both submission scripts source parameters from `../config/params.sh`, which must be configured before running:
  - `sample_sheet` - Determines which samples will be processed (scripts read sample IDs from this file)
  - `environment` - Specifies the conda environment to activate for the analysis
- Job names are automatically derived from the script filename
- Logs are written to `$project_dir/logs/`
- Resource requirements (walltime, memory, CPUs) are configured within each submission script and should be adjusted based on the specific analysis
- The scripts are configured for the NCI Gadi HPC environment and may need modification for other HPC systems

## Workflow Dependencies

Most analysis scripts expect the following structure:

**Input:**
- Preprocessed AnnData objects in `$h5ad_dir` (generated by 01_preprocessing)
- Configuration files in `../config/` (sample sheets, parameters)
- Reference data in `../resources/` (gene signatures, sample metadata)
- See [../docs/tree.txt](../docs/tree.txt) for a detailed overview of the expected data organization

**Output:**
- Results files (e.g., CSVs) in `$res_dir`
- Plots and figures in designated output directories

## Reproducibility

To reproduce the analysis:

1. **Run preprocessing** (01_preprocessing/) for your data platform
2. **Perform cell type mapping** (02_cell_type_mapping/) 
3. **Execute sequential analysis workflows** (03-05) as needed
4. **Generate figures** using scripts in figures/ directory

Each script includes:
- Argparse-based command-line interface for flexible execution
- Environment specifications in script headers (e.g., `# env: xen-hev`)
- Conda environment YAML files can be found in [../environment/](../environment/)

For specific analysis details and script execution order within each workflow, refer to the README files in each subdirectory.