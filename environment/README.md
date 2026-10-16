# Conda environments

- **EDF/** - environment definitions, with the main packages and versions: `conda env create -f environment/EDF/<env>.yaml`
- **ESF/** - snapshots of the environments used for the analysis, with every package version: `conda env create -f environment/ESF/<env>.yaml` to reproduce them exactly.

Each script names its environment in its header (`# env: xen-hev`). The workflow activates the environments by name, so create them under the names below; the job submission scripts take the name from `environment` in [config/params.sh](../config/params.sh).

| Environment | Used for |
|---|---|
| `spatial-pilot` | Visium analyses (steps 01, 03), building the scRNA-seq reference, and the Leiden clustering in `05_xenium_mechanisms/08_FDC_recovery` |
| `vis-hd` | Visium HD analyses (steps 01, 04) |
| `xen-hev` | Xenium analyses (steps 01, 05) |
| `spatdat` | Creating and reading the Xenium SpatialData objects (transcript-level analyses and plots) |
| `cell2location` | Cell type deconvolution (step 02) |
| `tangram` | Label transfer to Xenium (step 02) |
| `gseapy` | Gene set enrichment analyses |
| `shapely` | Figure scripts that need `shapely`, `adjustText` or `upsetplot` |
| `snakemake` | Running the [workflow](../workflow/README.md) |
| `data-access` | Pulling the data with DVC from the public archives (see [Getting the data](../README.md#getting-the-data)) |

`spatial-pilot`, `vis-hd` and `xen-hev` were kept separate during development but have similar contents; to run scripts by hand, one merged environment would also work (package versions may then differ slightly from those used), but the workflow expects the three names. Keep `cell2location`, `tangram` and `spatdat` separate: they follow the installation instructions of those packages and pin specific versions.
