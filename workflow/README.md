# Analysis workflow

A Snakemake workflow that runs the analysis scripts of steps 03–05 and the figure scripts in dependency order, starting from the published inputs; it also builds the scRNA-seq reference object from its GEO archive. Preprocessing (step 01) and the cell type mapping models (step 02) are not part of it; their outputs are pulled from Zenodo. To rerun those steps, see [scripts/README.md](../scripts/README.md).

## Requirements

- Snakemake 7.32.4: `conda env create -f environment/EDF/snakemake.yaml`, then `conda activate snakemake` and run from the repository root. Activate the environment rather than only putting `snakemake` on the `PATH`: the job submission script imports `snakemake` and `yaml` from it.
- The six conda environments the rules activate by name (`spatial-pilot`, `vis-hd`, `xen-hev`, `spatdat`, `shapely`, `gseapy`), so they must appear in `conda env list`: `conda env create -f environment/EDF/<env>.yaml`. Exact versions are in `environment/ESF/`.
- The inputs, pulled with DVC (see [Getting the data](../README.md#getting-the-data)):
  - from Zenodo (~50 GB): `data/anndata_objects`, `data/spatialdata_objects`, `data/pathology_annotations`, `results/02_cell_type_mapping`;
  - from GEO (533 MB): `resources/sc_ref/GSE176078_Wu_etal_2021_BRCA_scRNASeq.tar.gz`, from which the workflow builds the scRNA-seq reference;
  - the aligned H&E images, read only by the panels and zoom plots drawn on the H&E (leave them out with the `all_no_images` target). The workflow reads 15 of them (~133 GB). Pull the three folders in `data/images/he/aligned_images/` with DVC (~400 GB, all sections), or download just these files (each `<name>.ome.tif`) from the BioImage Archive ([S-BIAD3080](https://doi.org/10.6019/S-BIAD3080), `images/he/aligned_images/`) into the same paths under `data/`:
    - `visium/`: `P12_LN_Vis`, `P18_Vis`, `P22_Vis`, `P24_Vis`, `P29_Vis_Rep1`
    - `visium_hd/`: `P26_HD`, `P27_HD`, `P29_HD`, `P32_HD`
    - `xenium/`: `P17_Xe5k`, `P22_Xe5k`, `P24_Xe5k`, `P26_Xe5k`, `P29_Xe5k`, `P32_Xe5k`

## Running

```bash
snakemake --profile config/profiles/nci_a56 -n              # dry run
snakemake --profile config/profiles/nci_a56                 # everything except the ligand-receptor analysis
snakemake --profile config/profiles/nci_a56 all_long        # everything
snakemake --profile config/profiles/nci_a56 results/figures/Fig_2I.pdf --forcerun fig_2I   # one panel, remade
```

| Target | What it runs |
|---|---|
| `all` | all analyses and figures, except the Visium HD ligand-receptor analysis |
| `all_long` | `all` plus the ligand-receptor analysis (stLearn, up to ~4 h per sample on 16 CPUs) and the panels built on it |
| `all_no_images` | `all` without the steps that need the aligned H&E images |
| `all_xenium` | Xenium analyses and figure panels only |
| `all_visium` / `all_visium_long` | Visium and Visium HD analyses and figures, without / with the ligand-receptor analysis |

Any single output can be a target; figure panels are `results/figures/Fig_<panel>.pdf`, made by rules named `fig_<panel>`. Outputs go to `results/`. Rule logs are `logs/workflow/<rule>/<sample>.log` (per-sample rules) or `logs/workflow/<rule>.log`, in `logs/workflow_<res_dir>/` instead when `res_dir` is changed; scheduler logs are in `logs/joblogs/`.

Most jobs need under 50 GB of memory and finish in minutes. The largest requests are 256 GB (`fig_ED1G`, all 33 Visium sections at once) and 168 GB (`h04_calculate_vessel_distances`), and the longest job is the stLearn analysis of P29_HD (~4 h).

## Configuration

- `config/workflow/config.yaml`: output folder (`res_dir`), input folder (`h5ad_dir`), sample sheets and the Xenium HEV score thresholds. Any value can be overridden on the command line, e.g. `--config res_dir=results_test`.
- `config/profiles/nci_a56/`: the profile used on NCI Gadi (PBS). Each job's memory, walltime and CPUs come from the rule's `resources` (set from measured peaks); the queue is chosen from the memory request. Other Gadi users need to change `project` and `storage` in its `config.yaml`. On another cluster, adapt the profile; on a single machine, run with `--use-conda --cores <n>` instead of `--profile` (the largest jobs need the memory above).
