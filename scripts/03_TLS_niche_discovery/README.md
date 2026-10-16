# TLS Niche Discovery Workflow

This directory contains scripts and workflows for identifying, scoring, and characterizing tertiary lymphoid structure (TLS) niches and related immune aggregates in Visium data.

## Directory Structure

- **01_identify_TLS/**  
  Initial TLS scoring and spatial grouping.
  - `01_run_TLS_model.py`: Apply the published TLS predictive model (Andersson et al., 2021) to Visium spots.
  - `02_group_TLS_high.py`: Visualize and group spatially contiguous regions with high TLS scores.
  - `03_TLS_spatial_organization/`: Characterise discrete TLS across the cohort — number, size, spatial clustering (DBSCAN) and tumour topology (intratumoral / peritumoral / stromal-distant). Reproduces **Extended Data Fig. 3b**. See its [README](01_identify_TLS/03_TLS_spatial_organization/README.md).

- **02_identify_TLS_niches/**  
  Refined clustering and identification of internal TLS architectures.
  - `01_cluster_TLS_high.py`: Integrate, cluster, and visualize TLS-high spots to identify core niches (T cell zone, B cell zone, Germinal Center, HEV). *(Generates data for Supplementary Table 4)*.

- **03_extract_niche_signatures/**  
   Extraction and biological filtering of niche marker gene signatures. *(Generates data for Supplementary Table 5)*.
   - `01_extract_niche_signatures.py`: Derive raw marker gene signatures from TLS clustering results.
   - `02_filter_GC_signature.py`: Filter germinal center (GC) signature genes against an epithelial background to remove shared markers (e.g., proliferation, DNA damage repair).
   - `03a_identify_endothelial_high.py`: Identify a background of pan-endothelial high spots across the tissue.
   - `03b_create_endothelial_object.py`: Integrate endothelial-high spots across samples to create a vascular reference.
   - `03c_filter_HEV_signature.py`: Filter HEV marker genes against the generic endothelial background to isolate HEV-specific markers.
   - `04_organize_signature_files.py`: Copy and organize the filtered gene signature files (top 300 genes per signature) into a single directory for easier downstream access.

- **04_map_niche_signatures/**  
  Scoring and spatial mapping of the refined niche signatures across the tissue.
  - `01_score_signatures.py`: Score all tissue spots for the derived TLS niche gene signatures. This script contains a parameter (`use_resource_signatures`) to either load newly generated signatures from the results directory (set to `False`) or use the gene signatures extracted by the original analysis as reported in the manuscript (set to `True`, loads from [../../resources/gene_lists/TLS_niche_signatures/niche_signatures.csv](../../resources/gene_lists/TLS_niche_signatures/niche_signatures.csv)).
  - `02_group_high_score_areas.py`: Group contiguous regions with high TLS niche scores to enable region quantification (e.g., count GC reactions per tissue).

## Workflow Overview

1. **TLS Scoring:**  
   Apply published TLS predictive models to Visium data to assign a baseline TLS probability score to each spot.

2. **Spatial Grouping:**  
   Identify and visualize regions with high TLS scores, using spatial graphs to group neighboring spots into discrete aggregate structures.

3. **Clustering & Niche Identification:**  
   Integrate and cluster the TLS-high spots to resolve internal architecture. Extract marker genes and apply background filtering to isolate highly specific signatures for GC and HEV niches.

4. **Signature Mapping:**  
   Score the whole tissue for the refined TLS niche signatures, group high-scoring regions, and assess spatial continuity.

5. **Output:**  
   - Spatially annotated `.h5ad` objects with TLS and niche scores.
   - CSV files containing raw differentially expressed genes (**Supplementary Table 4**) and final filtered spatial signatures (**Supplementary Table 5**).
   - A consolidated directory of the top 300 filtered gene signature CSVs for each niche, for convenient downstream analysis.
   - Spatial plots mapping the T cell zone, B cell zone, and GC niches.