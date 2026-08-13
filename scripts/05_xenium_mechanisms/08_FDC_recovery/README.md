# 08_FDC_recovery

Recover **follicular dendritic cells (FDCs)** — absent from the single-cell reference — from the label-transferred Tfh pool, by per-sample subclustering and a complement-receptor–based call. Provides the FDC/Tfh cell labels used by the CXCL13 source attribution (09/07).

| Script | Env | Mode | Figure |
|---|---|---|---|
| `01_identify_FDC.py` | spatial-pilot | `--sample_id` | — (upstream) |
| `02_FDC_overview.py` | spatial-pilot | `--samples` | S9J |

**`01_identify_FDC.py`** — subclusters the Tfh pool (scale clip 10 → PCA 30 → kNN k=15 → Leiden res 1.0, seed 0) and calls a cluster FDC if CD3E & CD3G are both below the top-1000 markers, CR1+CR2 combined rank ≤ 250, and ≥ 5 % of cells express CR2. Writes `{s}_FDC_Tfh_barcodes.csv` (cols `cell_id, group ∈ {FDC, Tfh}`) + the subcluster DEG table + a per-sample diagnostic PDF.

**`02_FDC_overview.py`** — cohort FDC-vs-Tfh subcluster marker dot plot (S9J): FDC markers CR1/CR2/FCER2, the CXCL13–CXCR5 axis, pan-T CD3E/CD3G/CD4, and Tfh PDCD1/ICOS/BCL6/IL21/TIGIT. Samples with no recovered FDC are handled gracefully.

Note: Leiden clustering needs the `spatial-pilot` environment (it ships igraph/leidenalg), not `xen-hev`.

**Inputs:** `data/anndata_objects/xenium/{s}.h5ad`; cell labels `results/02_cell_type_mapping/label_transfer/xenium/subset_level/{s}_subset_pred_cell_types.csv` (Tfh label `T_cells_c3_CD4+_Tfh_CXCL13`).

**Outputs:** `results/05_xenium_mechanisms/08_FDC_recovery/{per_sample,all_samples}/`.

**Run order:** `01` (per sample) → `02` (cohort). Run before `09/07`. FDC is recovered in 5 of 6 samples; the non-FDC sample is auto-skipped downstream (no sample ID is hard-coded).
