# Preprocessing Scripts

This directory contains scripts for preprocessing spatial transcriptomics datasets from Visium, Visium HD, and Xenium platforms. The workflow is organized by data type, with each subdirectory containing scripts for basic quality control and normalization.

## Directory Structure

- `visium/`
  - `01_import_basic_qc.py`: Imports Visium data and performs basic quality control (QC).
  - `02_log_norm.py`: Applies log normalization to the Visium data.

- `visium_hd/`
  - `01_import_basic_qc.py`: Imports Visium HD data and performs basic QC.
  - `02_log_norm.py`: Applies log normalization to the Visium HD data.

- `xenium/`
  - `01a_make_adata_objects.py`: Creates AnnData objects from Xenium raw data. These objects contain cell-by-gene expression and a single spatial coordinate per cell (no transcript-level info or segmentation masks).
  - `01b_make_spatialdata_objects.py`: Creates SpatialData objects for Xenium. These objects include cell-by-gene expression, transcript spatial coordinates, and segmentation masks.

## Workflow

### Visium & Visium HD
1. **Import & QC**
   - Run `01_import_basic_qc.py` to import raw data and perform quality control.
2. **Normalization**
   - Run `02_log_norm.py` to apply log normalization to the QC'd data.

### Xenium
1. **Create AnnData Objects**
   - Run `01a_make_adata_objects.py` to generate AnnData objects from raw Xenium data.
   - Run `01b_make_spatialdata_objects.py` to generate spatial data objects for downstream analysis.

