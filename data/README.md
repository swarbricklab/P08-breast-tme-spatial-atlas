# Data

This directory contains the input data, organized by data type and platform; [docs/tree.txt](../docs/tree.txt) shows the layout after `dvc pull`.

For an overview of which samples are available per modality, see the sample sheets in [../config/](../config/) — one per platform (Visium, Visium HD, lymph nodes, and Xenium).

## Data availability

Each dataset below is tracked by the `.dvc` file shown and is fetched into the location shown; see [Getting the data](../README.md#getting-the-data) for how.

### Zenodo

A curated subset of this data is published on Zenodo: **[10.5281/zenodo.18917100](https://doi.org/10.5281/zenodo.18917100)** — *"Spatial Transcriptomics Data for Breast Cancer TLS and HEV Analysis"*. Each archive unpacks to the same folder layout used in this repository (DVC fetches the unpacked files directly):

| Zenodo archive | Unpacks to | Location in this repo | DVC file | Size |
|---|---|---|---|---|
| `anndata_objects.zip` | `visium/`, `visium_hd/`, `xenium/` | `data/anndata_objects/` | `data/anndata_objects.dvc` | 20 GB |
| `pathology_annotations.zip` | `visium/`, `visium_hd/`, `visium_lymph_nodes/` | `data/pathology_annotations/` | `data/pathology_annotations.dvc` | 102 MB |
| `spatialdata_objects.zip` | `xenium/` (6 SpatialData `.zip` objects) | `data/spatialdata_objects/` | `data/spatialdata_objects.dvc` | 27 GB |
| `cell_type_mapping.zip` | `deconvolution/`, `label_transfer/` | `results/02_cell_type_mapping/` (outside `data/`) | `results/02_cell_type_mapping.dvc` | 2.2 GB |

### BioImage Archive

Imaging data — H&E images, alignment files, CytAssist images, and IHC images (`data/images/`) — is available from the **BioImage Archive**: [S-BIAD3080](https://doi.org/10.6019/S-BIAD3080).

| Contents | Location in this repo | DVC file | Size |
|---|---|---|---|
| H&E images, Visium | `data/images/he/aligned_images/visium/` | `data/images/he/aligned_images/visium.dvc` | 290 GB |
| H&E images, Visium HD | `data/images/he/aligned_images/visium_hd/` | `data/images/he/aligned_images/visium_hd.dvc` | 54 GB |
| H&E images, Xenium | `data/images/he/aligned_images/xenium/` | `data/images/he/aligned_images/xenium.dvc` | 59 GB |
| Alignment files | `data/images/he/alignment_files/` | `data/images/he/alignment_files.dvc` | 34 MB |
| CytAssist images | `data/images/cytassist/` | `data/images/cytassist.dvc` | 1.1 GB |
| IHC images | `data/images/ihc/` | `data/images/ihc.dvc` | 398 GB |


### EGA

Sequencing reads (FASTQ) are deposited at the **European Genome-phenome Archive** under study [EGAS00001008558](https://ega-archive.org/studies/EGAS00001008558). Access is controlled: apply to the study's Data Access Committee through EGA. The reads are not tracked in this repository; the analyses here start from the pipeline outputs in `processed/`.

### GEO

The single-cell reference used for cell type mapping is pulled from its original GEO deposit, [GSE176078](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE176078) ([Wu et al., 2021](https://doi.org/10.1038/s41588-021-00911-1)):

| Contents | Location in this repo | DVC file | Size |
|---|---|---|---|
| scRNA-seq reference (compressed archive) | `resources/sc_ref/` | `resources/sc_ref/GSE176078_Wu_etal_2021_BRCA_scRNASeq.tar.gz.dvc` | 533 MB |

The 10x Genomics pipeline outputs (Space Ranger, Xenium Onboard Analysis), from which the AnnData and SpatialData objects are built, are on **GEO** under SuperSeries [GSE351850](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE351850), one sub-series per platform. Each is pulled as GEO serves it: the per-section files bundled in `<series>_RAW.tar`, plus a table linking each section to its other deposits. The workflow unpacks them into [processed/](#processed) when it needs to build an object.

| GEO series | Contents | Location in this repo | DVC file | Size |
|---|---|---|---|---|
| [GSE351843](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE351843) | Visium, including the lymph node sections | `data/geo/visium/` | `data/geo/visium.dvc` | 4.5 GB |
| [GSE351849](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE351849) | Visium HD | `data/geo/visium_hd/` | `data/geo/visium_hd.dvc` | 5.3 GB |
| [GSE350571](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE350571) | Xenium | `data/geo/xenium/` | `data/geo/xenium.dvc` | 167 GB |

### Recreate from source

- `data/anndata_objects/` and `data/spatialdata_objects/` — objects not pulled from Zenodo are built by the workflow from the GEO deposit (target `objects`; see [../workflow/README.md](../workflow/README.md#building-the-objects-from-geo)).
- `data/anndata_objects/sc_ref/` — built by the workflow (rule `s02_build_sc_ref`, see [../workflow/README.md](../workflow/README.md)) from the [GEO: GSE176078](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE176078) archive in `resources/sc_ref/`.
- `data/anndata_objects/visium_hd/square_002um/` — the 2 µm bins are not in the Zenodo record, to save space. No analysis here uses them (the Visium HD analyses use the 8, 16 and 48 µm bins); [01_preprocessing](../scripts/01_preprocessing/README.md) writes them when it builds the Visium HD objects.

## Directory Structure

### anndata_objects/

AnnData objects (`.h5ad`) with expression, spatial coordinates and metadata:

- **visium/** - Visium (55 µm spots), made by [scripts/01_preprocessing/visium/](../scripts/01_preprocessing/visium/).
- **visium_hd/** - Visium HD at 2, 8, 16 and 48 µm bins (`square_<size>um/`), made by [scripts/01_preprocessing/visium_hd/](../scripts/01_preprocessing/visium_hd/).
- **xenium/** - Xenium (single cells), made by [scripts/01_preprocessing/xenium/](../scripts/01_preprocessing/xenium/).
- **sc_ref/** - the scRNA-seq reference (`miniatlas_GEO.h5ad`), built by the workflow (see [Recreate from source](#recreate-from-source)).

### spatialdata_objects/

Xenium SpatialData objects with transcripts and cell boundaries, made by [01b_make_spatialdata_objects.py](../scripts/01_preprocessing/xenium/01_make_objects/01b_make_spatialdata_objects.py). They are stored as `.zip` and must be unpacked to `.zarr` before reading; the workflow does this, or by hand: `shutil.unpack_archive("P22_Xe5k.zip", "P22_Xe5k.zarr")`, then `spatialdata.read_zarr("P22_Xe5k.zarr")`.

### images/

- **he/aligned_images/** - H&E images (pyramidal OME-TIFF, `{section_id}.ome.tif`) in the same coordinate space as the AnnData objects, so spatial data can be plotted on the full-resolution image; [3A.py](../scripts/figures/figure3/3A.py) shows how to load one and zoom in with `crop_coord`. They open in [QuPath](https://qupath.github.io/) by drag and drop.
- **he/alignment_files/** - how each H&E image was registered:
  - `visium/`: Loupe Browser manual alignment JSON files (H&E to CytAssist image), usable to re-run Space Ranger. Some lymph node sections were aligned automatically and have no file.
  - `xenium/`: Xenium Explorer transformation matrices (`{section_id}_matrix.csv`) from Xenium to H&E coordinates; the Xenium preprocessing scripts read them from here. See the [10x Genomics guide](https://www.10xgenomics.com/support/software/xenium-explorer/latest/tutorials/xe-image-alignment).
- **cytassist/** - CytAssist images (`visium/`, `visium_hd/`), used by Space Ranger with the alignment files to place spots on the H&E image.
- **ihc/** - immunohistochemistry for cGAS, STING and phosphorylated TBK1 (`{patient_id}_{marker}.ome.tif`) for P17, P22, P24, P26, P29 and P32.

### pathology_annotations/

Pathologist annotations per spot or bin (CSV):

- **visium/** - annotated in separate layers (epithelial, immune, stroma, additional features; categories in [resources/pathology_annotations/](../resources/pathology_annotations/pathology_categories_per_layer.json)), plus a `single_layer` column with one label per spot.
- **visium_hd/** - one label per 8 µm bin.
- **visium_lymph_nodes/** - germinal center annotations for the lymph node sections, per spot (CSV) and as polygons (GeoJSON; open the H&E image in QuPath, then drag in the GeoJSON).

### processed/

The 10x Genomics pipeline outputs, the input to [01_preprocessing](../scripts/01_preprocessing/README.md): `spaceranger_outs/visium/`, `spaceranger_outs/visium_hd/` (with `binned_outputs/square_<size>um/` per bin size) and `xenium_outs/`. [00_unpack_geo.py](../scripts/01_preprocessing/00_unpack_geo.py) unpacks them from the [GEO](#geo) deposit in `geo/` (the workflow runs it when needed). See the 10x documentation for the [Space Ranger](https://www.10xgenomics.com/support/software/space-ranger/latest/analysis/outputs/output-overview) and [Xenium](https://www.10xgenomics.com/support/software/xenium-onboard-analysis/latest/analysis/xoa-output-understanding-outputs) output files; the Xenium outputs also open in [Xenium Explorer](https://www.10xgenomics.com/support/software/xenium-explorer).
