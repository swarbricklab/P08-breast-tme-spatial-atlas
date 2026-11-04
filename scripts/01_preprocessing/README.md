# 01_preprocessing

Build the AnnData and SpatialData objects in `data/anndata_objects/` and `data/spatialdata_objects/` from the 10x pipeline outputs on GEO ([GSE351850](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE351850); see [data/README.md](../../data/README.md#geo)). The objects are published, so this step only needs to be run to regenerate them. The workflow builds any object that has not been pulled from Zenodo (see [workflow/README.md](../../workflow/README.md#building-the-objects-from-geo)); to run the scripts by hand, run `00_unpack_geo.py` once per platform, then the others per sample (`--sample_id`), in the order listed.

| Script | Env | What it does |
|---|---|---|
| `00_unpack_geo.py` | vis-hd | Unpack the GEO deposit of one platform (`--modality`; `data/geo/<platform>/`) into the Space Ranger / Xenium Onboard Analysis layout in `data/processed/` |
| `visium/01_import_basic_qc.py` | spatial-pilot | Import the Space Ranger outputs into AnnData and run basic QC |
| `visium/02_log_norm.py` | spatial-pilot | Normalize (total counts) and log-transform |
| `visium_hd/01_import_basic_qc.py` | vis-hd | As above, for each Visium HD bin size |
| `visium_hd/02_log_norm.py` | vis-hd | As above, for each Visium HD bin size |
| `xenium/01_make_objects/01a_make_adata_objects.py` | xen-hev | AnnData object: cell-by-gene counts and one coordinate per cell |
| `xenium/01_make_objects/01b_make_spatialdata_objects.py` | spatdat | SpatialData object: also transcript locations and segmentation masks |
| `xenium/02_handle_images/01_make_adata_objects_with_images.py` | xen-hev | Transform the cell coordinates to H&E image space and add the downsampled H&E image |
| `xenium/02_handle_images/02_transform_spatialdata_object.py` | spatdat | Register the SpatialData object to the H&E image |

The Visium and Visium HD images come with the Space Ranger outputs. The Xenium image scripts need the aligned H&E images (`data/images/he/aligned_images/xenium/`) and the Xenium Explorer alignment matrices (`data/images/he/alignment_files/xenium/`; see the [10x Genomics guide](https://www.10xgenomics.com/support/software/xenium-explorer/latest/tutorials/xe-image-alignment)).

## Objects built from GEO

- **Visium HD 16 and 48 µm**: GEO holds the 2 and 8 µm bins only, so `00_unpack_geo.py` rebuilds these bins from `hd_feature_slice.h5`.