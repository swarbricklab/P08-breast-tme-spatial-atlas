# 02_cell_type_mapping

Map cell types from the scRNA-seq breast cancer atlas of [Wu et al. (2021)](https://doi.org/10.1038/s41588-021-00911-1) (GEO [GSE176078](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE176078); called "miniatlas" in the code) onto the spatial data: deconvolution with cell2location for Visium and Visium HD, label transfer with Tangram for Xenium. The results (`results/02_cell_type_mapping/`) are published, so this step only needs to be run to regenerate them. Scripts run in the order listed; the cell2location and Tangram scripts need a GPU (`submit_array_jobs_nci_gpu.sh`, see [scripts/README.md](../README.md#job-submission-on-nci-gadi)).

| Script | Env | What it does |
|---|---|---|
| `visium_deconvolution/01_prepare_sc_ref/01_create_GEO_miniatlas_object.py` | spatial-pilot | Build the reference AnnData object (`data/anndata_objects/sc_ref/miniatlas_GEO.h5ad`) from the GEO archive; also run by the workflow |
| `visium_deconvolution/01_prepare_sc_ref/02_train_model.py` | cell2location | Train the cell2location reference model |
| `visium_deconvolution/02_mapping/visium/01_cell_type_mapping.py` | cell2location | Deconvolve Visium spots (per sample) |
| `visium_deconvolution/02_mapping/visium_hd/01a_cell_type_mapping_visium_hd_48um.py` | cell2location | Deconvolve Visium HD 48 µm bins (per sample) |
| `visium_deconvolution/02_mapping/visium_hd/01b_cell_type_mapping_visium_hd_16um_chunked.py` | cell2location | Deconvolve Visium HD 16 µm bins, in chunks to fit in GPU memory (per sample) |
| `xenium_label_transfer/01_prepare_sc_ref_for_tangram.py` | tangram | Subset the reference to the Xenium panel genes and find marker genes per cell type |
| `xenium_label_transfer/02_label_transfer_sc_ref.py` | tangram | Transfer cell type labels to Xenium cells (per sample) |

The reference has three annotation levels: `major`, `minor` and `subset`. Each mapping script sets the level in a `level = ...` line near the top. The published results include all three levels for both methods (one run per level).
