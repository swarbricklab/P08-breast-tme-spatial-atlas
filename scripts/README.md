# Analysis Scripts

The analysis runs in five sequential steps, followed by the figure scripts. Each step folder has a README with the scripts in run order, their conda environment and what they produce.

| Step | Data | What it does |
|---|---|---|
| [01_preprocessing/](01_preprocessing/) | Visium, Visium HD, Xenium | Import, basic QC and normalization; builds the AnnData and SpatialData objects |
| [02_cell_type_mapping/](02_cell_type_mapping/) | all | Cell type deconvolution (cell2location; Visium, Visium HD) and label transfer (Tangram; Xenium) from a scRNA-seq reference |
| [03_TLS_niche_discovery/](03_TLS_niche_discovery/) | Visium | TLS scoring, TLS niche clustering and the niche gene signatures |
| [04_vascular_microenvironment/](04_vascular_microenvironment/) | Visium HD | HEV and vessel segmentation, HEV vs other venous vessels, ligand-receptor analysis |
| [05_xenium_mechanisms/](05_xenium_mechanisms/) | Xenium | TLS and HEV regions at single-cell resolution, TLS zonation, chemokine sources, SPP1 |
| [figures/](figures/) | | One script per figure panel, named after the panel (e.g. `figure3/3B.py` makes Fig. 3b) |

`helper_functions.py` holds functions shared across scripts.

## Running the scripts

Steps 03–05 and the figure scripts are run by the Snakemake workflow, starting from the published inputs; see [workflow/README.md](../workflow/README.md). The outputs of steps 01 and 02 are published, so these two steps only need to be run to regenerate them; the workflow also runs step 01 for any object that has not been pulled from Zenodo.

Every script has a command-line interface (`python <script> --help`) and names its conda environment in its header (`# env: xen-hev`; see [environment/](../environment/)). Scripts take `--project_dir`, `--h5ad_dir`, `--res_dir` and `--sample_sheet`, plus either `--sample_id` (one sample per run) or `--samples` (all samples in one run).

### Job submission on NCI Gadi

The bash scripts in this folder submit a script to the PBS scheduler on NCI Gadi; on other systems, adapt the `qsub` calls or run the scripts directly. Run them from this folder after setting `environment` and `sample_sheet` in [config/params.sh](../config/params.sh):

```bash
./submit_array_jobs_nci.sh -s <path/to/script.py>        # one job per sample (scripts with --sample_id)
./submit_array_jobs_nci_gpu.sh -s <path/to/script.py>    # the same on a GPU node (cell2location, Tangram)
./submit_multi_sample_job_nci.sh -s <path/to/script.py>  # one job for all samples (scripts with --samples)
```

Each job runs the script through `run_single_sample_job.sh` or `run_multi_sample_job.sh`, which activate the conda environment and pass the arguments above. Walltime, memory and CPUs are set in the submission scripts; logs go to `logs/`.
