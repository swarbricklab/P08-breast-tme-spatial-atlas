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
  - `01_make_objects/`
    - `01a_make_adata_objects.py`: Creates AnnData objects from Xenium raw data. These objects contain cell-by-gene expression and a single spatial coordinate per cell (no transcript-level info or segmentation masks).
    - `01b_make_spatialdata_objects.py`: Creates SpatialData objects for Xenium. These objects include cell-by-gene expression, transcript spatial coordinates, and segmentation masks.
  - `02_handle_images/`
    - `01_make_adata_objects_with_images.py`: Transforms Xenium coordinates to image space using an alignment matrix (see [10x Genomics guide](https://www.10xgenomics.com/support/software/xenium-explorer/latest/tutorials/xe-image-alignment) on manual alignment using Xenium Explorer), adds downsampled H&E images to the AnnData object, and saves the result.
    - `02_transform_spatialdata_object.py`: Applies alignment transformations to SpatialData objects so they are registered with the H&E image.

## Workflow

### Visium & Visium HD
1. **Import & QC**
   - Run `01_import_basic_qc.py` to import raw data and perform quality control.
2. **Normalization**
   - Run `02_log_norm.py` to apply log normalization to the QC'd data.

**Note:** All required images for Visium and Visium HD processing are included in the Space Ranger outputs, so no additional image files are needed.


### Xenium
1. **Create AnnData and SpatialData Objects**
  - Run `01_make_objects/01a_make_adata_objects.py` to generate AnnData objects from raw Xenium data.
  - Run `01_make_objects/01b_make_spatialdata_objects.py` to generate SpatialData objects for downstream analysis.
2. **Image Handling and Alignment**
  - Run `02_handle_images/01_make_adata_objects_with_images.py` to align Xenium coordinates to image space and add downsampled H&E images to the AnnData object.
  - Run `02_handle_images/02_transform_spatialdata_object.py` to apply alignment transformations to the SpatialData object for registration with the H&E image.

**Note:** The image handling scripts require full-resolution H&E images to be available in `../../data/images/he/aligned_images/xenium/` (as OME-TIFF pyramidal images). These scripts load and downsample the images to create versions suitable for inclusion in AnnData objects.

