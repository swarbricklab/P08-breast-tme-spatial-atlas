# Resources

This directory contains reference data, metadata, and auxiliary files used throughout the analysis workflow.

## Directory Structure

### gene_lists/

Curated gene signatures and marker gene sets used for analysis:

- **HEV_markers/** - High endothelial venule (HEV) marker genes from [Sawada et al., 2022](https://doi.org/10.1158/2326-6066.CIR-21-0369)
- **STING/** - STING pathway-related gene signatures
- **TLS_niche_signatures/** - Tertiary lymphoid structure (TLS) niche-specific gene signatures (T cell zone, B cell zone, Germinal Center, HEV) extracted from the original analysis. Use these signatures to reproduce manuscript results.
- **andersson_2021/** - Published TLS predictive model gene coefficients ([Andersson et al., 2021](https://doi.org/10.1038/s41467-021-26271-2))
- **enrichment_analysis/** - Gene sets for pathway enrichment analyses
- **marker_genes/** - Endothelial marker genes from [Kumar et al., 2023](https://doi.org/10.1038/s41586-023-06252-9)

### cell_types/

Cell type reference information:

- **miniatlas/** - Cell type hierarchies and metadata from the [single-cell reference atlas](https://doi.org/10.1038/s41588-021-00911-1) (major, minor, subset annotation levels)

### sc_ref/

Single-cell reference data for cell type deconvolution:

**Data retrieval:**
The single-cell reference data is a compressed archive from GEO ([GSE176078](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE176078)), tracked with DVC (see [Getting the data](../README.md#getting-the-data)). Pull and extract it:

```bash
dvc pull resources/sc_ref/GSE176078_Wu_etal_2021_BRCA_scRNASeq.tar.gz.dvc
cd resources/sc_ref
tar -xzf GSE176078_Wu_etal_2021_BRCA_scRNASeq.tar.gz    # creates Wu_etal_2021_BRCA_scRNASeq/
```

Without DVC, the same file can be downloaded directly from GEO (then extract as above):

```bash
cd resources/sc_ref
wget "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE176078&format=file&file=GSE176078%5FWu%5Fetal%5F2021%5FBRCA%5FscRNASeq%2Etar%2Egz" \
    -O GSE176078_Wu_etal_2021_BRCA_scRNASeq.tar.gz
```

**Contents:**
- **Wu_etal_2021_BRCA_scRNASeq/** - Breast cancer single-cell RNA-seq data from [Wu et al., 2021](https://doi.org/10.1038/s41588-021-00911-1)
  - `count_matrix_sparse.mtx` - Raw count matrix (sparse format)
  - `count_matrix_genes.tsv` - Gene identifiers
  - `count_matrix_barcodes.tsv` - Cell barcodes
  - `metadata.csv` - Cell metadata with annotations
- **ensembl_mapping/** - Gene symbol to Ensembl ID mapping files for resolving naming inconsistencies
  - `map.csv` - Curated gene mapping file
  - `genes.tsv` - Gene mapping from the original GEO dataset

**Usage:**
The raw count matrices are processed into an AnnData object using [01_create_GEO_miniatlas_object.py](../scripts/02_cell_type_mapping/visium_deconvolution/01_prepare_sc_ref/01_create_GEO_miniatlas_object.py), which creates `h5ad_dir/sc_ref/miniatlas_GEO.h5ad` for downstream cell type deconvolution and label transfer.

### probe_set/

Xenium gene panel information:

- **xenium_panel.csv** - Complete list of genes included in the Xenium spatial transcriptomics panel
  - `gene_name` - Gene symbol
  - `ensembl_id` - Ensembl gene identifier
  - `num_probesets` - Number of probe sets targeting this gene
  - `panel_source` - Panel origin (Xenium_Prime_Human_5Kpan or Custom_Breast_100g)

This file is used to subset the single-cell reference data to genes measured by Xenium for label transfer analysis (see [01_prepare_sc_ref_for_tangram.py](../scripts/02_cell_type_mapping/xenium_label_transfer/01_prepare_sc_ref_for_tangram.py)).

### plotting_colors/

Standardized color palettes for consistent visualization (CSV format):

- **miniatlas_major_colors.csv**, **miniatlas_minor_colors.csv**, **miniatlas_subset_colors.csv** - Color mappings for cell types at different annotation levels
- **miniatlas_lineage_colors.csv** - Colors for cell lineages
- **endothelial_colors.csv** - Colors for endothelial cell subtypes
- **TLS_HEV_colors.csv** - Colors for TLS and HEV annotations
- **SCSubtype_colors.csv** - Colors for breast cancer subtypes
- **transcript_colors.csv** - Colors for specific transcripts/genes

### areas_of_interest/

Manually defined regions of interest for focused analysis:

- **visium_hd/** - Regions of interest for Visium HD samples
- **xenium/** - Regions of interest for Xenium samples

### image_alignment/

Image alignment parameters and transformations:

- **xenium/** - Alignment parameters for registering Xenium images to spatial coordinates

### pathology_annotations/

Pathology annotation metadata:

- **pathology_categories_per_layer.json** - Defines the annotation categories used in each pathology layer (epithelial, immune, stroma, additional features)

### tissue_metadata.csv

Tissue-level clinical and experimental metadata including:
- Donor IDs
- Tissue IDs  
- Tissue types (primary tumor, lymph node, etc.)
- Clinical subtypes (Luminal A, Luminal B, Triple Negative, etc.)