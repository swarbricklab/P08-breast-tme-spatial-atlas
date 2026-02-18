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
    - `01_make_adata_objects_with_images.py`: Transforms Xenium coordinates to image space using an alignment matrix (see [10x Genomics guide](https://www.10xgenomics.com/support/software/xenium-explorer/latest/tutorials/xe-image-alignment) on manual alignment using Xenium Explorer), adds H&E images to the AnnData object, and saves the result.
    - `02_crop_HE_image.py`: Crops large H&E images to the tissue of interest (for convenience) and updates the AnnData object and image files accordingly.
    - `03_transform_spatialdata_object.py`: Applies alignment and cropping transformations to SpatialData objects so they are registered with the cropped H&E image.

## Workflow

### Visium & Visium HD
1. **Import & QC**
   - Run `01_import_basic_qc.py` to import raw data and perform quality control.
2. **Normalization**
   - Run `02_log_norm.py` to apply log normalization to the QC'd data.


### Xenium
1. **Create AnnData and SpatialData Objects**
  - Run `01_make_objects/01a_make_adata_objects.py` to generate AnnData objects from raw Xenium data.
  - Run `01_make_objects/01b_make_spatialdata_objects.py` to generate SpatialData objects for downstream analysis.
2. **Image Handling and Alignment**
  - Run `02_handle_images/01_make_adata_objects_with_images.py` to align Xenium coordinates to image space and add H&E images to the AnnData object.
  - Run `02_handle_images/02_crop_HE_image.py` to crop the H&E image and update the AnnData object and image files.
  - Run `02_handle_images/03_transform_spatialdata_object.py` to apply alignment and cropping transformations to the SpatialData object for registration with the cropped H&E image.

