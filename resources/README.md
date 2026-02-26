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