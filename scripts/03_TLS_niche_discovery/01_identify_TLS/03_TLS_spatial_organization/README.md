# 03_TLS_spatial_organization

Per-sample and cohort characterisation of **discrete tertiary lymphoid structures (TLS)** across the Visium cohort: TLS number, TLS size, spatial clustering (DBSCAN) and topology relative to the tumour. **Reproduces Extended Data Fig. 3b** (and the reported cohort total of 433 TLS across the 33 Visium sections; median 8/section, range 0–37).

## Discrete-TLS definition

These scripts operate on the *stringent* Extended Data Fig. 3b definition of a discrete TLS, derived from the permissive per-spot grouping written by `../02_group_TLS_high.py` (`group_singleton_removed` = contiguous groups of ≥2 directly adjacent TLS-high spots, score > 0.01):

1. keep only groups of **≥3 directly adjacent** TLS-high spots, then
2. **merge surviving groups separated by a single spot**.

`discrete_tls.py` implements this once for `01`, `02` and `04`, so all panels use the same TLS set; it is imported by those scripts, not run on its own.

## Scripts

| Script | Env | Mode | Output |
|---|---|---|---|
| `01_TLS_distances_topology.py` | spatial-pilot | cohort (`-s` / all) | `distances_per_TLS.csv` |
| `02_TLS_clustering.py` | spatial-pilot | cohort (`-s` / all) | `clustering_per_TLS.csv`, `clustering_per_sample.csv` |
| `03_TLS_number_size_overview.py` | spatial-pilot | cohort (`-s` / all) | overview PDF (**Extended Data Fig. 3b**) + summary CSVs |
| `04_TLS_per_sample_plots.py` | spatial-pilot | cohort (`-s` / all) | per-sample spatial PDF |
| `discrete_tls.py` | — | helper (imported) | — |

**`01_TLS_distances_topology.py`** — for each discrete TLS, classifies topology relative to cancer-annotated spots (intratumoral if ≥30% of TLS spots overlap cancer or >50% of the two-spot perimeter ring is cancer; peritumoral if any cancer spot lies in the one-spot ring; else stromal/distant; cancer = invasive + in situ), and records distances to the nearest other TLS and nearest cancer spot (µm + hop distance). Writes one row per TLS to `distances_per_TLS.csv` (`tls_topology_2hop`).

**`02_TLS_clustering.py`** — DBSCAN (eps = 500 µm, minimum 2 members) on the pairwise edge-to-edge (minimum spot-to-spot) distances between discrete TLS, labelling each TLS `clustered` or `lone`. Writes `clustering_per_TLS.csv` (per-TLS labels) and `clustering_per_sample.csv` (per-sample summary).

**`03_TLS_number_size_overview.py`** — the cohort overview figure (**Extended Data Fig. 3b**): a per-sample sorted bar of TLS number with stacked panels for TLS size, clustering composition (clustered vs lone) and topology, plus subtype/capture-area strips. Writes the single-page overview PDF and `summary_per_sample.csv` / `sizes_per_aggregate.csv`.

**`04_TLS_per_sample_plots.py`** — per-section 3-panel spatial overview (epithelial annotation, TLS topology, DBSCAN clustering) over the H&E, all sections in one multi-page PDF.

## Inputs

- Permissive TLS groupings: `results/03_TLS_niche_discovery/01_identify_TLS/02_group_TLS_high/threshold_01/csv_files/{s}_TLS_groups.csv`
- TLS-model scores: `results/03_TLS_niche_discovery/01_identify_TLS/01_run_TLS_model/{s}_TLS_model.csv`
- Visium objects: `data/anndata_objects/visium/{s}.h5ad`
- Pathology annotations: `data/pathology_annotations/visium/{s}_pathology_annotations.csv` (`epithelial` column; `invasive`/`in_situ` = cancer)

## Outputs

`results/03_TLS_niche_discovery/01_identify_TLS/03_TLS_spatial_organization/{01_TLS_distances_topology, 02_TLS_clustering, 03_TLS_number_size_overview, 04_TLS_per_sample_plots}/threshold_01/…`

## Run order

Requires `../01_run_TLS_model.py` and `../02_group_TLS_high.py` to have been run for the whole Visium cohort. Then, within this directory: run `01` and `02` first (they default to all sections; a `-s`/`--samples` subset is also accepted), then `03` (cohort overview) and `04` (per-sample plots), which read the CSVs written by `01`/`02`.
