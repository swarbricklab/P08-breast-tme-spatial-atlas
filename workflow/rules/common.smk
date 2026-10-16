# Shared settings for all rule files: sample lists, absolute project path, default script args.
import os
import pandas as pd

PROJECT_DIR = os.getcwd()            # snakemake is run from the repository root
RES = config["res_dir"]
# Rule logs per results folder, so runs with different res_dir never overwrite each other's logs
LOGW = "logs/workflow" if RES.rstrip("/") == "results" else "logs/workflow_" + RES.strip("/").replace("/", "_")
H5AD = config["h5ad_dir"]
SHEETS = config["sample_sheets"]

VIS_SAMPLES = pd.read_csv(SHEETS["visium"])["section_id"].tolist()
HD_SAMPLES = pd.read_csv(SHEETS["visium_hd"])["section_id"].tolist()
LN_SAMPLES = pd.read_csv(SHEETS["visium_lymph_nodes"])["section_id"].tolist()
XEN_SAMPLES = pd.read_csv(SHEETS["xenium"])["section_id"].tolist()

# NB: no f-strings in .smk files: Snakemake 7.32 on Python 3.12 mangles them (pads the fields with
# spaces, breaks {{...}}). Use concatenation, and expand(..., allow_missing=True) for wildcards.

# Arguments every analysis script accepts (see scripts/run_single_sample_job.sh)
def base_args(sheet):
    """--h5ad_dir/--project_dir/--res_dir/--sample_sheet for a sample sheet key of the config."""
    return ("--h5ad_dir " + PROJECT_DIR + "/" + H5AD + " --project_dir " + PROJECT_DIR
            + " --res_dir " + PROJECT_DIR + "/" + RES + " --sample_sheet " + PROJECT_DIR + "/" + SHEETS[sheet])

# Target list filled by the rule files; rule all collects it.
TARGETS = []
