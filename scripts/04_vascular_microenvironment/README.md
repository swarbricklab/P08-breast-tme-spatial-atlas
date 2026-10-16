# Vascular Microenvironment Analysis Workflow

This directory contains scripts for identifying, classifying, and analyzing vascular niches in high-resolution spatial transcriptomics data (**Visium HD, 16 µm resolution**). The workflow resolves specific endothelial lineages (arterial, capillary, venous, lymphatic), segments individual vessels, and investigates ligand-receptor cross-talk within the high endothelial venule (HEV) niche.

## Directory Structure

- **01_identify_HEV_niches/** 
  Identification and spatial grouping of HEV structures.
  - `01_score_signature.py`: Score Visium HD bins for HEV and TLS zone signatures. This script contains a parameter (`use_resource_signatures`) to either load newly generated signatures from the results directory (set to `False`) or use the gene signatures extracted by the original analysis as reported in the manuscript (set to `True`, loads from [../../resources/gene_lists/TLS_niche_signatures/niche_signatures.csv](../../resources/gene_lists/TLS_niche_signatures/niche_signatures.csv)).
  - `02a_group_HEV_bins.py`: Identify high-scoring HEV "seed" bins and group contiguous regions.
  - `02b_expand_HEV_groups.py`: Apply the "seed-and-expand" algorithm to include neighboring bins with lower HEV scores into discrete vascular structures.
  - `03a_identify_TLS_areas.py`: Segment and group continuous TLS regions.
  - `03b_calculate_HEV_TLS_distance.py`: Calculate distances from segmented HEVs to TLS areas to classify them as "TLS-associated" or "stromal HEV-like".
  - `04_assess_HEV_pathology_annotations.py`: Map HEV locations against pathologist-annotated tissue compartments.
  - `05_calculate_all_HEV_distances.py`: For each HEV, calculate the distance to surrounding spots to analyze cellular composition across distance. Outputs per-spot distance metrics and summary plots for spatial analysis.

- **02_identify_all_vasculature/** 
  Classification of the broader tumor vascular network.
  - `01_classify_vessels.py`: Classify segmented vessels based on the dominant deconvolved endothelial cell type (arterial, capillary, venous, lymphatic).
  - `02_investigate_vessels.py`: Analyze spatial overlap and associations between vessel types and specific perivascular-like (PVL) subsets.
  - `03_calculate_distances.py`: Calculate radial distances from segmented vessels to surrounding tissue bins.
  - `04_assess_PVL_association.py`: Generate spatial enrichment profiles of PVL composition across distance bins from vessels.

- **03_compare_HEV_to_other_venous/** 
  Targeted comparison between specialized HEVs and canonical venous structures.
  - `01_compare_microenvironment.py`: Compare the pathology-defined tissue composition surrounding HEV vs. canonical venous vessels.
  - `02_compare_gene_expression.py`: Perform differential gene expression analysis between HEV and non-HEV venous bins. *(Generates data for **Supplementary Table 8** and Fig. 3i)*.
  - `03_compare_PVL_association.py`: Quantify the enrichment of inflammatory PVLs in the immediate HEV niche vs. the canonical venous niche.

- **04_ligand_receptor_analysis/** 
  Ligand-receptor (L-R) interaction analysis within the vascular microenvironment. *(Generates data for **Supplementary Table 11** and Fig. 4)*.
  - `01_run_LR_analysis_stlearn.py`: Execute `stLearn` to identify significant L-R co-expression in HEV-localized (0 µm) and HEV-perivascular (25 µm) neighborhoods (one run each, `--bin_size 0` / `--bin_size 25`; default 0).
  - `02_analyze_LR_all.py`: Aggregate interaction results across samples, calculating directional signaling vectors (e.g., Perivascular → HEV) and identifying TLS-enriched pathways (`--bin_size 0` / `--bin_size 25`, matching the stLearn run; default 25).
  - `03_plot_LR_per_sample.py`: Generate per-sample visualizations of significant ligand-receptor interactions, including spatial plots and dotplots stratified by distance to HEV (`--bin_size`, default 25).
  - `04_check_composition_LTB_signaling_bins.py`: Deconstruct the cellular composition of spatial bins sending LTB signals to endothelial LTBR (e.g., cDC2 enrichment).

## Workflow Overview

1. **HEV & TLS Segmentation:** 
   Map HEV signatures at 16 µm resolution, apply a graph-based seed-and-expand algorithm to segment individual vessels, and classify their spatial context relative to TLSs.

2. **Vessel Classification & PVL Coupling:** 
   Segment the entire vascular network using deconvolved endothelial abundances and map the specific spatial coupling of distinct PVL subsets to each vessel lineage.

3. **Venous Specialization:** 
   Isolate the HEV phenotype by directly comparing the gene expression and perivascular remodeling of HEVs against generic veins.

4. **Signaling Mechanisms:** 
   Use permutation-based spatial L-R analysis to identify interactions driving HEV maturation (e.g., LTB-LTBR) and immune recruitment (e.g., chemokines, integrins).

5. **Outputs:** 
   - Spatially annotated `.h5ad` objects containing vessel segmentations and distance metrics.
   - CSV files containing HEV vs. Venous DEGs (**Supplementary Table 8**) and significant L-R interactions (**Supplementary Table 11**).
   - Spatial enrichment profiles and L-R dot plots.
