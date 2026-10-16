# Resources

Reference files, metadata and settings used by the scripts. Throughout the code, "miniatlas" refers to the scRNA-seq breast cancer atlas of [Wu et al. (2021)](https://doi.org/10.1038/s41588-021-00911-1), used as the cell type reference.

## gene_lists/

- **TLS_niche_signatures/** - the TLS niche signatures (T cell zone, B cell zone, germinal center, HEV) published with the manuscript (Supplementary Table 5). The scoring scripts use these by default; see [03_TLS_niche_discovery](../scripts/03_TLS_niche_discovery/README.md).
- **andersson_2021/** - gene coefficients of the TLS model of [Andersson et al. (2021)](https://doi.org/10.1038/s41467-021-26271-2).
- **HEV_markers/** - HEV marker genes from [Sawada et al. (2022)](https://doi.org/10.1158/2326-6066.CIR-21-0369).
- **marker_genes/** - endothelial marker genes from [Kumar et al. (2023)](https://doi.org/10.1038/s41586-023-06252-9).
- **STING/** - STING pathway and type I interferon signatures from MSigDB (`signature_links.csv` links each to its MSigDB page and gives the number of genes on the Xenium panel).
- **enrichment_analysis/** - Gene Ontology Biological Process and Molecular Function gene sets (2025, `.gmt`) for GSEA.

## sc_ref/

The scRNA-seq reference: the GEO archive of [GSE176078](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE176078) (pulled with `dvc pull resources/sc_ref/GSE176078_Wu_etal_2021_BRCA_scRNASeq.tar.gz.dvc`) and gene symbol to Ensembl ID mappings (`ensembl_mapping/`: `genes.tsv` from the GEO dataset, `map.csv` curated). The workflow (rule `s02_build_sc_ref`) extracts the archive into `Wu_etal_2021_BRCA_scRNASeq/` and builds `data/anndata_objects/sc_ref/miniatlas_GEO.h5ad` with [01_create_GEO_miniatlas_object.py](../scripts/02_cell_type_mapping/visium_deconvolution/01_prepare_sc_ref/01_create_GEO_miniatlas_object.py).

Without DVC, download the archive from GEO into this folder:

```bash
cd resources/sc_ref
wget "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE176078&format=file&file=GSE176078%5FWu%5Fetal%5F2021%5FBRCA%5FscRNASeq%2Etar%2Egz" \
    -O GSE176078_Wu_etal_2021_BRCA_scRNASeq.tar.gz
```

## cell_types/

- **miniatlas/miniatlas_subset_classifications.csv** - the reference's cell type hierarchy: each `subset` cell type with its `minor` and `major` type and lineage (Epithelial, Immune, Stromal); used to order and color cell types.

## probe_set/

- **xenium_panel.csv** - the genes on the Xenium panel (`gene_name`, `ensembl_id`, `num_probesets`, `panel_source`: Xenium_Prime_Human_5Kpan or Custom_Breast_100g), used to subset the reference for label transfer.

## plotting_colors/

Color palettes (CSV): cell types per annotation level (`miniatlas_{lineage,major,minor,subset}_colors.csv`), endothelial subtypes (`endothelial_colors.csv`), TLS and HEV annotations (`TLS_HEV_colors.csv`), breast cancer subtypes (`SCSubtype_colors.csv`), pathology annotation layers (`path_colors.csv`) and transcripts (`transcript_colors.csv`).

## areas_of_interest/

Manually drawn regions (GeoJSON) for the zoomed panels of single representative sections, in `visium_hd/` and `xenium/`.

## pathology_annotations/

- **pathology_categories_per_layer.json** - the annotation categories of each pathology layer (epithelial, immune, stroma, additional features).

## tissue_metadata.csv

One row per tissue: `donor_id`, `tissue_id`, `tissue_type` (Primary tumor or Lymph node) and `clinical_subtype` (Luminal A, Luminal B, Luminal B (HER2+), HER2+ or TNBC; `NA` for lymph nodes).
