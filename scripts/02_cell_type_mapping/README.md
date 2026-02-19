# Cell Type Mapping Scripts

This directory contains scripts for cell type mapping in spatial transcriptomics datasets (Visium, Visium HD & Xenium) based on a single-cell RNA-seq reference. The workflow is organized by data type and mapping method (deconvolution or label transfer).

## Directory Structure

- `visium_deconvolution/`
  - `01_prepare_sc_ref/`
    - `01_create_GEO_miniatlas_object.py`: Creates a single-cell reference object from the raw count matrix and metadata obtained from [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE176078) and saves as `.h5ad`.
    - `02_train_model.py`: Trains a cell2location model using the reference object.
  - `02_mapping/`
    - `visium/`
      - `01_cell_type_mapping.py`: Maps cell types for standard Visium data using the trained model.
    - `visium_hd/`
      - `01a_cell_type_mapping_visium_hd_48um.py`: Maps cell types for Visium HD (48um resolution).
      - `01b_cell_type_mapping_visium_hd_16um_chunked.py`: Maps cell types for Visium HD (16um resolution), splitting the object to avoid GPU memory issues.

- `xenium_label_transfer/`
  - `01_prepare_sc_ref_for_tangram.py`: Prepares the single-cell reference for Tangram label transfer.
  - `02_label_transfer_sc_ref.py`: Performs label transfer using Tangram.


## Single-cell RNA-seq Reference Object

The single-cell RNA-seq reference object used in this workflow is constructed from raw data obtained from [GEO: GSE176078](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE176078). This reference includes hierarchical cell type annotations at three levels:

- **major**: Broad cell type categories
- **minor**: Finer subtypes within major categories
- **subset**: Most granular annotation level

All scripts in this directory support running analyses at any of these annotation levels. The desired annotation level (e.g., `major`, `minor`, or `subset`) must be defined within each script, typically as a parameter or variable. This allows flexible mapping and analysis depending on the resolution required for a given study.


---

## Workflow

### Visium/Visium HD Deconvolution (cell2location)
1. **Prepare Reference**
   - Run `01_create_GEO_miniatlas_object.py` to create the reference object.
   - Run `02_train_model.py` to train the cell2location model.
2. **Cell Type Mapping**
   - For standard Visium: run `01_cell_type_mapping.py`.
   - For Visium HD 48um: run `cell_type_mapping_visium_hd_48um.py`.
   - For Visium HD 16um: run `cell_type_mapping_visium_hd_16um_chunked.py` (handles large objects by chunking).

### Xenium Label Transfer (Tangram)
1. **Prepare Reference**
   - Run `01_prepare_sc_ref_for_tangram.py` to subset the reference to genes present in the Xenium panel and find marker genes per cell type.
2. **Label Transfer**
   - Run `02_label_transfer_sc_ref.py` to perform label transfer using Tangram.
