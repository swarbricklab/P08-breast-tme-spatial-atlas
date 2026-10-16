# A cellular and vascular atlas of the breast tumor microenvironment with multi-resolution spatial transcriptomics

Analysis repository for spatial transcriptomics profiling of tertiary lymphoid structures (TLS) and high endothelial venules (HEV) in breast cancer using Visium, Visium HD, and Xenium platforms.

**Data:** Zenodo [10.5281/zenodo.18917100](https://doi.org/10.5281/zenodo.18917100) · BioImage Archive [S-BIAD3080](https://doi.org/10.6019/S-BIAD3080) · EGA [EGAS00001008558](https://ega-archive.org/studies/EGAS00001008558) (sequencing reads, controlled access)

## Getting the data

All input data are tracked as DVC imports that point directly at the public archives holding the data (Zenodo, BioImage Archive, GEO), so no access to our own storage is needed. `dvc pull` fetches files from the archive and verifies them against the checksums pinned in this repository. [data/README.md](data/README.md#data-availability) lists each `.dvc` file, what it contains, and its size.

```bash
git clone https://github.com/swarbricklab/P08-breast-tme-spatial-atlas.git
cd P08-breast-tme-spatial-atlas

conda env create -f environment/EDF/data-access.yaml
conda activate data-access
archivefs dvc enable    # once per environment: lets DVC resolve chiro://, zenodo://, bia:// and geo:// URLs
archivefs dvc status    # should report "installed: True"

dvc list -R --size . data/pathology_annotations    # inspect what a target contains, without downloading it
dvc pull data/pathology_annotations.dvc            # fetch it
```

- `dvc pull` with no target fetches everything (~1.4 TB), so pulling selectively is usually more practical. `data/processed.dvc` cannot yet be pulled from outside the lab (see [data/README.md](data/README.md#geo)).
- Run `dvc` from inside the repository: the archive manifests the pins refer to live in `.chirograph/`.
- `Unsupported URL type chiro://` means `archivefs dvc enable` has not been run in the active environment.
- The remote in `.dvc/config` is the lab's private DVC cache; it is not needed to pull the data.

## Running the analyses

A Snakemake workflow runs the analysis (steps 03–05) and figure scripts in dependency order from the published inputs. See [workflow/README.md](workflow/README.md) for requirements, how to run it, targets and configuration.

## Repository structure

| Folder | Contents |
|---|---|
| [scripts/](scripts/README.md) | Analysis scripts in five sequential steps, the figure scripts, and job submission scripts for NCI Gadi |
| [workflow/](workflow/README.md) | Snakemake workflow that runs steps 03–05 and the figure scripts |
| [config/](config/) | Sample sheets (one per platform), workflow configuration and cluster profile, `params.sh` for the job submission scripts |
| [data/](data/README.md) | Input data (AnnData and SpatialData objects, images, pathology annotations), pulled with DVC; layout in [docs/tree.txt](docs/tree.txt) |
| [resources/](resources/README.md) | Gene signatures, cell type hierarchy, colour palettes, regions of interest and tissue metadata |
| [environment/](environment/README.md) | Conda environments |
| results/ | Outputs, one folder per step; `results/02_cell_type_mapping` (cell type mapping) is pulled with DVC, the rest is written by the workflow |
