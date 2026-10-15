# Rules for scripts/05_xenium_mechanisms (Xenium Prime 5K breast tumours):
#   01 identify TLS regions     - TLS polygons + per-cell in_TLS (per sample); composition (S7B-C)
#                                 and CD40/CD40LG expression in TLS (S7D)
#   02 identify HEV regions     - HEV niches per sample (thresholds from config)
#   03 HEV microenvironment     - distance of every cell to the nearest HEV niche (Fig 4E)
#   04 score signatures         - STING / type I IFN signature scores (Table S9, Fig 4C)
#   05 gene-gene correlations   - TMEM173 / SPP1 pseudobulk correlations + GSEA (Table S10, Fig 4D,
#                                 ED Fig 10d/f)
#   06 cellular enrichment      - HEV SPP1 vs cell type abundance <200 um (Table S12, Fig 5H)
#   07 TLS zonation             - B/T zones + GCs on a 50 um grid (S9C-D, S9F-G, S10G-H)
#   08 FDC recovery             - FDC/Tfh subclustering (S9J)
#   09 chemokine sources        - CXCL13/CCL19/CCL21 sources and gradients (S9A-B, S9E, S9H-I, S9K-L)
#   10 HEV-GC proximity         - close vs far from GC DEG (S10I-J)
#   11 SPP1 venous specificity  - SPP1 across venous EC groups (S10C), macrophage SPP1 zoom (S10E)
#
# No f-strings in this file: Snakemake 7.32 on Python 3.12 mangles them, so paths are built by
# concatenation (expand(..., allow_missing=True) keeps the {sample} wildcard).
# Inputs from outside the workflow: {H5AD}/xenium/{sample}.h5ad, {H5AD}/sc_ref/miniatlas_GEO.h5ad,
# data/spatialdata_objects/xenium/{sample}.zip, data/images/he/aligned_images/xenium/{sample}.ome.tif,
# {RES}/02_cell_type_mapping/label_transfer/xenium/{level}_level, resources/.
#
# Notes:
# - The spatialdata .zip is unpacked once per sample by x05_unzip_spatialdata into
#   data/spatialdata_objects/xenium/{sample}.zarr (same call as the scripts:
#   shutil.unpack_archive(zip, zarr)); every rule whose script reads the .zarr (07/03, 07/04, 09/02,
#   09/05, 09/07, 11/02, fig_5B, fig_5E) takes it as input, so the scripts always find a complete
#   .zarr and never unpack it themselves (07/03, 07/04, 11/02, 5B, 5E would unpack next to the .zip;
#   09/02, 09/05, 09/07 into a hard-coded tmp dir, /scratch/a56/sl6147/tmp).

X05_SCRIPTS = "scripts/05_xenium_mechanisms"
X05 = RES + "/05_xenium_mechanisms"

X05_ARGS = base_args("xenium")
X05_SAMPLE_ARGS = " --samples " + " ".join(XEN_SAMPLES)
X05_RE = "|".join(XEN_SAMPLES)

# Published inputs
X05_H5AD = H5AD + "/xenium/{sample}.h5ad"
X05_SC_REF = H5AD + "/sc_ref/miniatlas_GEO.h5ad"
X05_ZIP = "data/spatialdata_objects/xenium/{sample}.zip"
X05_ZARR = "data/spatialdata_objects/xenium/{sample}.zarr"     # unpacked by x05_unzip_spatialdata
X05_HE = "data/images/he/aligned_images/xenium/{sample}.ome.tif"
X05_SHEET = SHEETS["xenium"]
X05_METADATA = "resources/tissue_metadata.csv"
X05_CLASSIF = "resources/cell_types/miniatlas/miniatlas_subset_classifications.csv"
X05_COLORS = "resources/plotting_colors"
X05_AOI = "resources/areas_of_interest/xenium"


def x05_labels(level):
    """Label transfer CSV of a sample at an annotation level (published input under res_dir)."""
    return (RES + "/02_cell_type_mapping/label_transfer/xenium/" + level + "_level/{sample}_"
            + level + "_pred_cell_types.csv")


def x05_all(pattern):
    """Expand a {sample} path over all Xenium samples."""
    return expand(pattern, sample=XEN_SAMPLES)


# Upstream outputs read by several rules
X05_TLS_DIR = X05 + "/01_identify_TLS_regions/csv_files"
X05_TLS_ANN = X05_TLS_DIR + "/{sample}_cell_TLS_annotations.csv"
X05_TLS_GEOJSON = X05_TLS_DIR + "/{sample}_TLS_regions.geojson"
X05_HEV_DIR = X05 + "/02_identify_HEV_regions/csv_files"
X05_HEV_C2N = X05_HEV_DIR + "/{sample}_HEV_niches_cell_to_niche.csv"
X05_HEV_INFO = X05_HEV_DIR + "/{sample}_HEV_niches_info.csv"
X05_ZONES = X05 + "/07_TLS_zonation/01_classify_zones/csv_files"
X05_BT_BINS = X05_ZONES + "/{sample}_BT_zone_bins.csv"
X05_BT_CELLS = X05_ZONES + "/{sample}_BT_zone_classification.csv"
X05_GC_BINS = X05_ZONES + "/{sample}_GC_zone_bins.csv"
X05_GC_CELLS = X05_ZONES + "/{sample}_GC_zone_classification.csv"
X05_GC_GENES = X05_ZONES + "/{sample}_GC_genes_used.csv"
X05_FDC_BARCODES = X05 + "/08_FDC_recovery/per_sample/csv_files/{sample}_FDC_Tfh_barcodes.csv"
X05_CHEMO = X05 + "/09_chemokine_sources"
X05_PRODUCERS = X05_CHEMO + "/01_reference_producers/csv_files"
X05_TX_GENES = ["CXCL13", "CCL19", "CCL21", "TSPAN7"]     # 02_extract_transcripts.py --genes default
X05_TX = X05_CHEMO + "/02_extract_transcripts/csv_files/{sample}_{gene}_per_cell.csv"
X05_HEV_GC = X05 + "/10_HEV_GC_proximity"

# Single representative samples: only samples with an annotated areas-of-interest geojson
X05_ZONE_ZOOM_SAMPLES = [s for s in glob_wildcards(X05_AOI + "/{sample}_zone_zoom.geojson").sample
                         if s in XEN_SAMPLES]
X05_SPP1_MAC_SAMPLES = [s for s in glob_wildcards(X05_AOI + "/{sample}_SPP1_macrophages.geojson").sample
                        if s in XEN_SAMPLES]


# ── spatialdata objects ──────────────────────────────────────────────────────

rule x05_unzip_spatialdata:
    """Unpack the published SpatialData .zip into a .zarr next to it, exactly as the scripts do
    (shutil.unpack_archive(zip, zarr)), so no script unpacks it itself (07/03, 07/04, 09/02,
    09/05, 09/07, 11/02, fig_5B, fig_5E)."""
    input:
        archive = X05_ZIP,
    output:
        zarr = directory(X05_ZARR),
    wildcard_constraints:
        sample = X05_RE,
    conda:
        "spatdat"
    resources:
        mem_mb  = 8000,        # peak 3.9 GB (1 jobs)
        runtime = 15,          # peak 2 min
    log:
        LOGW + "/x05_unzip_spatialdata/{sample}.log",
    shell:
        "python -c \"import shutil; shutil.unpack_archive('{input.archive}', '{output.zarr}')\" > {log} 2>&1"


# ── 01 · identify TLS regions ────────────────────────────────────────────────

rule x05_identify_TLS:
    """01/01: TLS polygons from smoothed T & B cell densities; per-cell in_TLS / TLS_id."""
    input:
        h5ad     = X05_H5AD,
        labels   = x05_labels("major"),
        sheet    = X05_SHEET,
        metadata = X05_METADATA,
        colors   = X05_COLORS + "/miniatlas_major_colors.csv",
    output:
        geojson = X05_TLS_GEOJSON,
        cells   = X05_TLS_ANN,
        pdf     = X05 + "/01_identify_TLS_regions/plots/{sample}_TLS_regions_spatial.pdf",
    wildcard_constraints:
        sample = X05_RE,
    conda:
        "xen-hev"
    params:
        args = X05_ARGS,
    resources:
        mem_mb  = 16000,       # peak 11.2 GB (12 jobs)
        runtime = 15,          # peak 6 min
    log:
        LOGW + "/x05_identify_TLS/{sample}.log",
    shell:
        "python {X05_SCRIPTS}/01_identify_TLS_regions/01_identify_TLS.py {params.args} "
        "--sample_id {wildcards.sample} > {log} 2>&1"


rule x05_composition_TLS:
    """01/02: cohort TLS composition, inside vs outside TLS (S7B) and average composition (S7C)."""
    input:
        cells   = x05_all(X05_TLS_ANN),
        labels  = [x05_all(x05_labels(level)) for level in ["major", "minor", "subset"]],
        classif = X05_CLASSIF,
        colors  = X05_COLORS + "/miniatlas_major_colors.csv",
        sheet   = X05_SHEET,
    output:
        pdf = X05 + "/01_identify_TLS_regions/all_samples/cellular_composition/plots/All_samples_TLS_cellular_composition.pdf",
        csv = expand(X05 + "/01_identify_TLS_regions/all_samples/cellular_composition/csv_files/All_samples_{level}_{table}.csv",
                     level=["major", "minor", "subset"],
                     table=["inside_vs_outside_TLS", "inside_vs_outside_TLS_per_sample",
                            "TLS_composition", "TLS_composition_per_sample"]),
    conda:
        "xen-hev"
    params:
        args = X05_ARGS + X05_SAMPLE_ARGS,
    resources:
        mem_mb  = 4000,        # peak 2.0 GB (2 jobs)
        runtime = 15,          # peak 2 min
    log:
        LOGW + "/x05_composition_TLS.log",
    shell:
        "python {X05_SCRIPTS}/01_identify_TLS_regions/02_composition_TLS.py {params.args} > {log} 2>&1"


rule x05_gene_expression_TLS:
    """01/03: cohort CD40 / CD40LG expression across subset-level cell types in TLS (S7D)."""
    input:
        h5ad    = x05_all(X05_H5AD),
        cells   = x05_all(X05_TLS_ANN),
        labels  = x05_all(x05_labels("subset")),
        classif = X05_CLASSIF,
        colors  = X05_COLORS + "/miniatlas_major_colors.csv",
        sheet   = X05_SHEET,
    output:
        pdf = X05 + "/01_identify_TLS_regions/all_samples/gene_expression/plots/All_samples_TLS_gene_expression_CD40LG_CD40.pdf",
    conda:
        "xen-hev"
    params:
        args = X05_ARGS + X05_SAMPLE_ARGS,
    resources:
        mem_mb  = 40000,       # peak 28.4 GB (2 jobs)
        runtime = 15,          # peak 2 min
    log:
        LOGW + "/x05_gene_expression_TLS.log",
    shell:
        "python {X05_SCRIPTS}/01_identify_TLS_regions/03_gene_expression_TLS.py {params.args} > {log} 2>&1"


# ── 02 · identify HEV regions ────────────────────────────────────────────────

rule x05_identify_HEV:
    """02/01: HEV score, HEV niches (core + expansion), distance of each niche to TLS."""
    input:
        h5ad     = X05_H5AD,
        major    = x05_labels("major"),
        subset   = x05_labels("subset"),
        cells    = X05_TLS_ANN,
        sheet    = X05_SHEET,
        metadata = X05_METADATA,
        colors   = X05_COLORS + "/miniatlas_major_colors.csv",
        classif  = X05_CLASSIF,
    output:
        cell_to_niche = X05_HEV_C2N,
        info          = X05_HEV_INFO,
        parameters    = X05_HEV_DIR + "/{sample}_HEV_niches_parameters.json",
        pdf           = X05 + "/02_identify_HEV_regions/plots/{sample}_identify_HEV_regions.pdf",
    wildcard_constraints:
        sample = X05_RE,
    conda:
        "xen-hev"
    params:
        args = (X05_ARGS + " --hev_high_thr " + str(config["hev_high_thr"])
                + " --hev_low_thr " + str(config["hev_low_thr"])),
    resources:
        mem_mb  = 24000,       # peak 16.0 GB (18 jobs)
        runtime = 15,          # peak 8 min
    log:
        LOGW + "/x05_identify_HEV/{sample}.log",
    shell:
        "python {X05_SCRIPTS}/02_identify_HEV_regions/01_identify_HEV.py {params.args} "
        "--sample_id {wildcards.sample} > {log} 2>&1"


# ── 03 · HEV microenvironment ────────────────────────────────────────────────

rule x05_calculate_distance_HEV:
    """03/01: distance of every cell to the nearest HEV niche (all / close to TLS / far from TLS)."""
    input:
        h5ad          = X05_H5AD,
        cell_to_niche = X05_HEV_C2N,
        info          = X05_HEV_INFO,
        cells         = X05_TLS_ANN,
        sheet         = X05_SHEET,
        metadata      = X05_METADATA,
    output:
        closest = X05 + "/03_HEV_microenvironment/01_calculate_distance/csv_files/{sample}_closest_HEV_distances.csv",
        close   = X05 + "/03_HEV_microenvironment/01_calculate_distance/csv_files/{sample}_HEV_distances_close_TLS.csv",
        far     = X05 + "/03_HEV_microenvironment/01_calculate_distance/csv_files/{sample}_HEV_distances_far_from_TLS.csv",
        pdf     = X05 + "/03_HEV_microenvironment/01_calculate_distance/plots/{sample}_HEV_distances.pdf",
    wildcard_constraints:
        sample = X05_RE,
    conda:
        "xen-hev"
    params:
        args = X05_ARGS,
    resources:
        mem_mb  = 16000,       # peak 10.7 GB (18 jobs)
        runtime = 15,          # peak 4 min
    log:
        LOGW + "/x05_calculate_distance_HEV/{sample}.log",
    shell:
        "python {X05_SCRIPTS}/03_HEV_microenvironment/01_calculate_distance_HEV.py {params.args} "
        "--sample_id {wildcards.sample} > {log} 2>&1"


# ── 04 · score signatures ────────────────────────────────────────────────────

rule x05_score_STING_signatures:
    """04/01: STING and type I IFN signature scores per cell (Table S9)."""
    input:
        h5ad          = X05_H5AD,
        major         = x05_labels("major"),
        subset        = x05_labels("subset"),
        cell_to_niche = X05_HEV_C2N,
        info          = X05_HEV_INFO,
        signatures    = ["resources/gene_lists/STING/STING_signatures.csv",
                         "resources/gene_lists/STING/type_I_interferon_signatures.csv"],
        colors        = [X05_COLORS + "/miniatlas_major_colors.csv",
                         X05_COLORS + "/endothelial_colors.csv",
                         X05_COLORS + "/TLS_HEV_colors.csv"],
        classif       = X05_CLASSIF,
        sheet         = X05_SHEET,
        metadata      = X05_METADATA,
    output:
        scores = X05 + "/04_score_signatures/01_score_STING_signatures/csv_files/{sample}_STING_signature_scores.csv",
        genes  = X05 + "/04_score_signatures/01_score_STING_signatures/csv_files/{sample}_STING_signatures_genes_used.csv",
        pdf    = X05 + "/04_score_signatures/01_score_STING_signatures/plots/{sample}_STING_signatures.pdf",
    wildcard_constraints:
        sample = X05_RE,
    conda:
        "xen-hev"
    params:
        args = X05_ARGS,
    resources:
        mem_mb  = 12000,       # peak 9.0 GB (18 jobs)
        runtime = 15,          # peak 6 min
    log:
        LOGW + "/x05_score_STING_signatures/{sample}.log",
    shell:
        "python {X05_SCRIPTS}/04_score_signatures/01_score_STING_signatures.py {params.args} "
        "--sample_id {wildcards.sample} > {log} 2>&1"


# ── 05 · gene-gene correlations ──────────────────────────────────────────────

# Manuscript runs (Supplementary Table 10): TMEM173 x GO MF (Fig 4d), SPP1 x GO MF (ED Fig 10d),
# SPP1 x GO BP (ED Fig 10f)
X05_CORR_GENES = ["TMEM173", "SPP1"]
X05_GSEA_RUNS = {"gene":     ["TMEM173", "SPP1", "SPP1"],
                 "gene_set": ["GO_Molecular_Function_2025", "GO_Molecular_Function_2025",
                              "GO_Biological_Process_2025"]}
X05_GENE_RE = "|".join(X05_CORR_GENES)
X05_GENE_SET_RE = "GO_Molecular_Function_2025|GO_Biological_Process_2025"
X05_CORR = X05 + "/05_gene_gene_correlations/01_investigate_gene_correlations/{gene}_correlations"
X05_GSEA = X05 + "/05_gene_gene_correlations/02_GSEA_analysis/{gene}_correlations/{gene_set}"

rule x05_gene_correlations:
    """05/01: per-HEV-niche pseudobulk correlation of all genes with {gene} (--gene_of_interest)."""
    input:
        h5ad          = x05_all(X05_H5AD),
        cell_to_niche = x05_all(X05_HEV_C2N),
        info          = x05_all(X05_HEV_INFO),
        sheet         = X05_SHEET,
        metadata      = X05_METADATA,
    output:
        csv = X05_CORR + "/csv_files/All_samples_HEV_{gene}_gene_correlations.csv",
        pdf = X05_CORR + "/plots/All_samples_HEV_{gene}_gene_correlations.pdf",
    wildcard_constraints:
        gene = X05_GENE_RE,
    conda:
        "xen-hev"
    params:
        args = X05_ARGS + X05_SAMPLE_ARGS,
    resources:
        mem_mb  = 32000,       # peak 24.6 GB (6 jobs)
        runtime = 15,          # peak 2 min
    log:
        LOGW + "/x05_gene_correlations/{gene}.log",
    shell:
        "python {X05_SCRIPTS}/05_gene_gene_correlations/01_investigate_gene_correlations.py {params.args} "
        "--gene_of_interest {wildcards.gene} > {log} 2>&1"


rule x05_GSEA_correlated_genes:
    """05/02: preranked GSEA ({gene_set}, --gene_set) on the {gene} correlations (--gene_of_interest).
    The top-terms barplot is only written when significant terms exist, so it is not an output."""
    input:
        corr     = rules.x05_gene_correlations.output.csv,
        gene_set = "resources/gene_lists/enrichment_analysis/{gene_set}.gmt",
    output:
        csv = X05_GSEA + "/csv_files/All_samples_HEV_{gene}_preranked_GSEA_results.csv",
    wildcard_constraints:
        gene     = X05_GENE_RE,
        gene_set = X05_GENE_SET_RE,
    conda:
        "gseapy"
    params:
        args = X05_ARGS + X05_SAMPLE_ARGS,
    resources:
        mem_mb  = 4000,        # peak 1.6 GB (12 jobs)
        runtime = 15,          # peak 3 min
    log:
        LOGW + "/x05_GSEA_correlated_genes/{gene}_{gene_set}.log",
    shell:
        "python {X05_SCRIPTS}/05_gene_gene_correlations/02_GSEA_analysis_on_correlated_genes.py {params.args} "
        "--gene_of_interest {wildcards.gene} --gene_set {wildcards.gene_set} > {log} 2>&1"


# ── 06 · cellular enrichment ─────────────────────────────────────────────────

X05_ENRICH = X05 + "/06_cellular_enrichment/SPP1_vs_cell_abundance/maxdist_200um"

rule x05_SPP1_vs_proximal_enrichment:
    """06/01: cell type proportions <200 um of each HEV niche vs niche endothelial SPP1 (Table S12)."""
    input:
        h5ad          = X05_H5AD,
        cell_to_niche = X05_HEV_C2N,
        info          = X05_HEV_INFO,
        cells         = X05_TLS_ANN,
        labels        = x05_labels("subset"),
        colors        = [X05_COLORS + "/miniatlas_major_colors.csv", X05_COLORS + "/TLS_HEV_colors.csv"],
        classif       = X05_CLASSIF,
        sheet         = X05_SHEET,
        metadata      = X05_METADATA,
    output:
        props = X05_ENRICH + "/csv_files/{sample}_HEV_SPP1_and_cell_proportions.csv",
        corr  = X05_ENRICH + "/csv_files/{sample}_HEV_SPP1_and_cell_proportions_corr.csv",
        pdf   = X05_ENRICH + "/plots/{sample}_HEV_SPP1_vs_cell_abundance.pdf",
    wildcard_constraints:
        sample = X05_RE,
    conda:
        "xen-hev"
    params:
        args = X05_ARGS,
    resources:
        mem_mb  = 12000,       # peak 8.7 GB (18 jobs)
        runtime = 15,          # peak 3 min
    log:
        LOGW + "/x05_SPP1_vs_proximal_enrichment/{sample}.log",
    shell:
        "python {X05_SCRIPTS}/06_cellular_enrichment/01_SPP1_vs_proximal_enrichment.py {params.args} "
        "--sample_id {wildcards.sample} > {log} 2>&1"


# ── 07 · TLS zonation ────────────────────────────────────────────────────────

rule x05_classify_zones:
    """07/01: B/T cell zones and GCs on a 50 um grid inside the TLS regions.
    The zone PDF and *_genes_used.csv are only written when TLS bins exist (true for all six
    samples); GC_genes_used is declared because 07/03 reads it."""
    input:
        h5ad       = X05_H5AD,
        geojson    = X05_TLS_GEOJSON,
        cells      = X05_TLS_ANN,
        signatures = "resources/gene_lists/TLS_niche_signatures/niche_signatures.csv",
        sheet      = X05_SHEET,
        metadata   = X05_METADATA,
    output:
        bt_bins  = X05_BT_BINS,
        bt_cells = X05_BT_CELLS,
        gc_bins  = X05_GC_BINS,
        gc_cells = X05_GC_CELLS,
        gc_genes = X05_GC_GENES,
    wildcard_constraints:
        sample = X05_RE,
    conda:
        "xen-hev"
    params:
        args = X05_ARGS,
    resources:
        mem_mb  = 12000,       # peak 8.9 GB (12 jobs)
        runtime = 15,          # peak 2 min
    log:
        LOGW + "/x05_classify_zones/{sample}.log",
    shell:
        "python {X05_SCRIPTS}/07_TLS_zonation/01_classify_zones.py {params.args} "
        "--sample_id {wildcards.sample} > {log} 2>&1"


rule x05_zonation_overview:
    """07/02: cohort B/T zone map (S9C), T:B ratio per bin (S9D), GC overview (S10G), GCs per sample (S10H)."""
    input:
        h5ad     = x05_all(X05_H5AD),
        labels   = x05_all(x05_labels("major")),
        bt_bins  = x05_all(X05_BT_BINS),
        bt_cells = x05_all(X05_BT_CELLS),
        gc_bins  = x05_all(X05_GC_BINS),
        geojson  = x05_all(X05_TLS_GEOJSON),
        cells    = x05_all(X05_TLS_ANN),
        colors   = X05_COLORS + "/SCSubtype_colors.csv",
        sheet    = X05_SHEET,
        metadata = X05_METADATA,
    output:
        pdf    = X05 + "/07_TLS_zonation/02_zonation_overview/plots/cohort_zonation_overview.pdf",
        n_gcs  = X05 + "/07_TLS_zonation/02_zonation_overview/csv_files/all_samples_n_discrete_GCs.csv",
    conda:
        "xen-hev"
    params:
        args = X05_ARGS + X05_SAMPLE_ARGS,
    resources:
        mem_mb  = 32000,       # peak 22.7 GB (2 jobs)
        runtime = 15,          # peak 3 min
    log:
        LOGW + "/x05_zonation_overview.log",
    shell:
        "python {X05_SCRIPTS}/07_TLS_zonation/02_zonation_overview.py {params.args} > {log} 2>&1"


rule x05_zoom_GCs:
    """07/03: per discrete GC, aligned-H&E zoom with the top GC-signature transcripts (S10G).
    Reads the full-resolution H&E into memory (P24 was OOM-killed at 124.7 GB)."""
    input:
        gc_bins  = X05_GC_BINS,
        gc_genes = X05_GC_GENES,
        h5ad     = X05_H5AD,
        image    = X05_HE,
        sdata    = X05_ZARR,
        colors   = X05_COLORS + "/transcript_colors.csv",
        sheet    = X05_SHEET,
        metadata = X05_METADATA,
    output:
        pdf = X05 + "/07_TLS_zonation/03_zoom_GCs/plots/{sample}_GC_zoom.pdf",
    wildcard_constraints:
        sample = X05_RE,
    conda:
        "spatdat"
    params:
        args = X05_ARGS,
    resources:
        mem_mb  = 40000,       # peak 28.4 GB (P32; region reads)
        runtime = 30,          # peak 9 min
    log:
        LOGW + "/x05_zoom_GCs/{sample}.log",
    shell:
        "python {X05_SCRIPTS}/07_TLS_zonation/03_zoom_GCs.py {params.args} "
        "--sample_id {wildcards.sample} > {log} 2>&1"


rule x05_zoom_zones:
    """07/04: adjacent B / T cell zone single-transcript zoom (S9F-G); samples with a
    resources/areas_of_interest/xenium/{sample}_zone_zoom.geojson (P22_Xe5k)."""
    input:
        geojson  = X05_AOI + "/{sample}_zone_zoom.geojson",
        bt_bins  = X05_BT_BINS,
        h5ad     = X05_H5AD,
        image    = X05_HE,
        sdata    = X05_ZARR,
        colors   = X05_COLORS + "/transcript_colors.csv",
        sheet    = X05_SHEET,
        metadata = X05_METADATA,
    output:
        pdf = X05 + "/07_TLS_zonation/04_zoom_zones/plots/{sample}_zone_zoom.pdf",
    wildcard_constraints:
        sample = X05_RE,
    conda:
        "spatdat"
    params:
        args = X05_ARGS,
    resources:
        mem_mb  = 24000,       # peak 16.2 GB (region reads)
        runtime = 15,          # peak 3 min
    log:
        LOGW + "/x05_zoom_zones/{sample}.log",
    shell:
        "python {X05_SCRIPTS}/07_TLS_zonation/04_zoom_zones.py {params.args} "
        "--sample_id {wildcards.sample} > {log} 2>&1"


# ── 08 · FDC recovery ────────────────────────────────────────────────────────

rule x05_identify_FDC:
    """08/01: Tfh-pool subclustering to recover FDC (Leiden needs spatial-pilot).
    The DEG table and PDF are skipped for samples with <20 Tfh cells, so only the barcodes are outputs."""
    input:
        h5ad     = X05_H5AD,
        labels   = x05_labels("subset"),
        sheet    = X05_SHEET,
        metadata = X05_METADATA,
    output:
        barcodes = X05_FDC_BARCODES,
    wildcard_constraints:
        sample = X05_RE,
    conda:
        "spatial-pilot"
    params:
        args = X05_ARGS,
    resources:
        mem_mb  = 16000,       # peak 9.5 GB (12 jobs)
        runtime = 15,          # peak 7 min
    log:
        LOGW + "/x05_identify_FDC/{sample}.log",
    shell:
        "python {X05_SCRIPTS}/08_FDC_recovery/01_identify_FDC.py {params.args} "
        "--sample_id {wildcards.sample} > {log} 2>&1"


rule x05_FDC_overview:
    """08/02: cohort FDC vs Tfh subcluster marker dot plot (S9J)."""
    input:
        h5ad     = x05_all(X05_H5AD),
        barcodes = x05_all(X05_FDC_BARCODES),
        sheet    = X05_SHEET,
        metadata = X05_METADATA,
    output:
        counts = X05 + "/08_FDC_recovery/all_samples/csv_files/FDC_reannotation_counts.csv",
        pdf    = X05 + "/08_FDC_recovery/all_samples/plots/FDC_cohort_summary.pdf",
    conda:
        "spatial-pilot"
    params:
        args = X05_ARGS + X05_SAMPLE_ARGS,
    resources:
        mem_mb  = 28000,       # peak 19.5 GB (2 jobs)
        runtime = 15,          # peak 2 min
    log:
        LOGW + "/x05_FDC_overview.log",
    shell:
        "python {X05_SCRIPTS}/08_FDC_recovery/02_FDC_overview.py {params.args} > {log} 2>&1"


# ── 09 · chemokine sources ───────────────────────────────────────────────────

rule x05_reference_producers:
    """09/01: chemokine producer subsets in the scRNA-seq reference + niche-gene dot plot (S9A)."""
    input:
        sc_ref  = X05_SC_REF,
        classif = X05_CLASSIF,
        colors  = [X05_COLORS + "/miniatlas_lineage_colors.csv", X05_COLORS + "/miniatlas_major_colors.csv"],
    output:
        producers = expand(X05_PRODUCERS + "/{kind}_producers/over_{pct}/{kind}_producers.csv",
                           kind=["ccl", "cxcl13"], pct=["10", "15"]),
        niche_csv = X05_PRODUCERS + "/niche_genes/niche_gene_pct_expressing.csv",
        pdf       = X05_CHEMO + "/01_reference_producers/plots/niche_genes/S9A_miniatlas_niche_genes_dotplots.pdf",
    conda:
        "xen-hev"
    params:
        args = X05_ARGS,
    resources:
        mem_mb  = 12000,       # peak 8.2 GB (2 jobs)
        runtime = 15,          # peak 3 min
    log:
        LOGW + "/x05_reference_producers.log",
    shell:
        "python {X05_SCRIPTS}/09_chemokine_sources/01_reference_producers.py {params.args} > {log} 2>&1"


rule x05_extract_transcripts:
    """09/02: per-cell integer transcript counts (CXCL13, CCL19, CCL21, TSPAN7) from SpatialData."""
    input:
        sdata = X05_ZARR,
    output:
        tx = expand(X05_TX, gene=X05_TX_GENES, allow_missing=True),
    wildcard_constraints:
        sample = X05_RE,
    conda:
        "spatdat"
    params:
        args = X05_ARGS,
    resources:
        mem_mb  = 32000,       # peak 23.6 GB (8 jobs)
        runtime = 15,          # peak 4 min
    log:
        LOGW + "/x05_extract_transcripts/{sample}.log",
    shell:
        "python {X05_SCRIPTS}/09_chemokine_sources/02_extract_transcripts.py {params.args} "
        "--sample_id {wildcards.sample} > {log} 2>&1"


rule x05_impvl_vs_icaf:
    """09/03: imPVL s1 vs iCAF s2 discrimination dot plot in TLS (S9B)."""
    input:
        h5ad     = x05_all(X05_H5AD),
        labels   = x05_all(x05_labels("subset")),
        geojson  = x05_all(X05_TLS_GEOJSON),
        classif  = X05_CLASSIF,
        colors   = X05_COLORS + "/miniatlas_subset_colors.csv",
        sheet    = X05_SHEET,
        metadata = X05_METADATA,
    output:
        counts = X05_CHEMO + "/03_impvl_vs_icaf/csv_files/focus_cell_counts_per_sample.csv",
        pdf    = X05_CHEMO + "/03_impvl_vs_icaf/plots/S9B_impvl_s1_vs_icaf_s2_TLS.pdf",
    conda:
        "xen-hev"
    params:
        args = X05_ARGS + X05_SAMPLE_ARGS,
    resources:
        mem_mb  = 20000,       # peak 12.7 GB (2 jobs)
        runtime = 15,          # peak 4 min
    log:
        LOGW + "/x05_impvl_vs_icaf.log",
    shell:
        "python {X05_SCRIPTS}/09_chemokine_sources/03_impvl_vs_icaf.py {params.args} > {log} 2>&1"


rule x05_zone_expression_density:
    """09/04: CXCL13 / CCL19 CP10K density per zone, paired Wilcoxon (S9E)."""
    input:
        h5ad     = x05_all(X05_H5AD),
        bt_bins  = x05_all(X05_BT_BINS),
        sheet    = X05_SHEET,
        metadata = X05_METADATA,
    output:
        csv = X05_CHEMO + "/04_zone_expression_density/csv_files/cohort_expression_density_avg_BT_seed42.csv",
        pdf = X05_CHEMO + "/04_zone_expression_density/plots/S9E_cohort_expression_density_libnorm_avg_BT_seed42.pdf",
    conda:
        "xen-hev"
    params:
        args = X05_ARGS + X05_SAMPLE_ARGS,
    resources:
        mem_mb  = 24000,       # peak 16.9 GB (2 jobs)
        runtime = 15,          # peak 3 min
    log:
        LOGW + "/x05_zone_expression_density.log",
    shell:
        "python {X05_SCRIPTS}/09_chemokine_sources/04_zone_expression_density.py {params.args} > {log} 2>&1"


X05_GRADIENTS = X05_CHEMO + "/05_chemokine_gradients"

rule x05_chemokine_gradients:
    """09/05: CCL19/CCL21 transcript and imPVL s1 / iCAF s2 density vs distance from HEV (S9H)."""
    input:
        h5ad          = x05_all(X05_H5AD),
        labels        = x05_all(x05_labels("subset")),
        bt_bins       = x05_all(X05_BT_BINS),
        cell_to_niche = x05_all(X05_HEV_C2N),
        sdata         = x05_all(X05_ZARR),
        colors        = [X05_COLORS + "/transcript_colors.csv", X05_COLORS + "/miniatlas_subset_colors.csv"],
        sheet         = X05_SHEET,
        metadata      = X05_METADATA,
    output:
        per_sample = expand(X05_GRADIENTS + "/csv_files/{sample}_{table}.csv", sample=XEN_SAMPLES,
                            table=["transcript_density_vs_dist", "celltype_density_vs_dist",
                                   "coexpr_source_composition", "coexpr_source_transcripts",
                                   "distance_profile_meta"]),
        pdf        = X05_GRADIENTS + "/plots/S9H_CCL21_CCL19_density_vs_dist_bin_15um_excl0.pdf",   # name written by the script (0 µm bin excluded)
    conda:
        "spatdat"
    params:
        args = X05_ARGS + X05_SAMPLE_ARGS,
    resources:
        mem_mb  = 56000,       # peak 41.0 GB (3 jobs)
        runtime = 30,          # peak 11 min
    log:
        LOGW + "/x05_chemokine_gradients.log",
    shell:
        "python {X05_SCRIPTS}/09_chemokine_sources/05_chemokine_gradients.py {params.args} > {log} 2>&1"


rule x05_reference_producer_corr:
    """09/06: reference-producer density vs CCL19/CCL21 transcript density per 50 um bin (S9I)."""
    input:
        producers = X05_PRODUCERS + "/ccl_producers/over_10/ccl_producers.csv",
        tx        = [x05_all(X05_TX.replace("{gene}", gene)) for gene in ["CCL19", "CCL21"]],
        h5ad      = x05_all(X05_H5AD),
        labels    = x05_all(x05_labels("subset")),
        bt_bins   = x05_all(X05_BT_BINS),
        classif   = X05_CLASSIF,
        colors    = [X05_COLORS + "/miniatlas_lineage_colors.csv", X05_COLORS + "/miniatlas_major_colors.csv",
                     X05_COLORS + "/transcript_colors.csv"],
        sheet     = X05_SHEET,
        metadata  = X05_METADATA,
    output:
        csv = X05_CHEMO + "/06_reference_producer_corr/csv_files/reference_producer_corr.csv",
        pdf = X05_CHEMO + "/06_reference_producer_corr/plots/S9I_reference_producer_corr.pdf",
    conda:
        "xen-hev"
    params:
        args = X05_ARGS + X05_SAMPLE_ARGS,
    resources:
        mem_mb  = 16000,       # peak 9.9 GB (2 jobs)
        runtime = 15,          # peak 2 min
    log:
        LOGW + "/x05_reference_producer_corr.log",
    shell:
        "python {X05_SCRIPTS}/09_chemokine_sources/06_reference_producer_corr.py {params.args} > {log} 2>&1"


X05_ATTRIBUTION = X05_CHEMO + "/07_cxcl13_proximity_attribution"

rule x05_cxcl13_proximity_attribution:
    """09/07: Tfh/FDC density per zone (S9K) and CXCL13 transcript source attribution (S9L)."""
    input:
        producers = X05_PRODUCERS + "/cxcl13_producers/over_15/cxcl13_producers.csv",
        barcodes  = x05_all(X05_FDC_BARCODES),
        bt_cells  = x05_all(X05_BT_CELLS),
        bt_bins   = x05_all(X05_BT_BINS),
        labels    = x05_all(x05_labels("subset")),
        h5ad      = x05_all(X05_H5AD),
        sdata     = x05_all(X05_ZARR),
        classif   = X05_CLASSIF,
        colors    = [X05_COLORS + "/miniatlas_lineage_colors.csv", X05_COLORS + "/miniatlas_subset_colors.csv"],
        sheet     = X05_SHEET,
        metadata  = X05_METADATA,
    output:
        csv = expand(X05_ATTRIBUTION + "/csv_files/{table}.csv",
                     table=["cxcl13_tfh_fdc_density_per_zone", "cxcl13_producer_attribution_counts_long",
                            "cxcl13_producer_attribution_per_sample", "cxcl13_producer_attribution_cohort"]),
        pdf = X05_ATTRIBUTION + "/plots/S9K_S9L_cxcl13_proximity_attribution.pdf",
    conda:
        "spatdat"
    params:
        args = X05_ARGS + X05_SAMPLE_ARGS,
    resources:
        mem_mb  = 44000,       # peak 32.8 GB (2 jobs)
        runtime = 30,          # peak 12 min
    log:
        LOGW + "/x05_cxcl13_proximity_attribution.log",
    shell:
        "python {X05_SCRIPTS}/09_chemokine_sources/07_cxcl13_proximity_attribution.py {params.args} "
        "> {log} 2>&1"


# ── 10 · HEV-GC proximity ────────────────────────────────────────────────────

X05_HEV_GC_DIST = X05_HEV_GC + "/01_classify_HEV_by_GC_dist/csv_files/{sample}_intra_TLS_HEV_GC_dist.csv"

rule x05_classify_HEV_by_GC_dist:
    """10/01: intra-TLS HEV niches close (<=500 um) / far from the nearest GC."""
    input:
        h5ad          = X05_H5AD,
        cell_to_niche = X05_HEV_C2N,
        info          = X05_HEV_INFO,
        gc_bins       = X05_GC_BINS,
        sheet         = X05_SHEET,
        metadata      = X05_METADATA,
    output:
        csv = X05_HEV_GC_DIST,
        pdf = X05_HEV_GC + "/01_classify_HEV_by_GC_dist/plots/{sample}_intra_TLS_HEV_GC_dist.pdf",
    wildcard_constraints:
        sample = X05_RE,
    conda:
        "xen-hev"
    params:
        args = X05_ARGS,
    resources:
        mem_mb  = 12000,       # peak 7.5 GB (18 jobs)
        runtime = 15,          # peak 3 min
    log:
        LOGW + "/x05_classify_HEV_by_GC_dist/{sample}.log",
    shell:
        "python {X05_SCRIPTS}/10_HEV_GC_proximity/01_classify_HEV_by_GC_dist.py {params.args} "
        "--sample_id {wildcards.sample} > {log} 2>&1"


rule x05_DEG_close_vs_far:
    """10/02: close vs far from GC Mann-Whitney DEG for ACKR1+ HEV endothelial cells."""
    input:
        dist          = x05_all(X05_HEV_GC_DIST),
        cell_to_niche = x05_all(X05_HEV_C2N),
        labels        = x05_all(x05_labels("subset")),
        h5ad          = x05_all(X05_H5AD),
    output:
        csv = X05_HEV_GC + "/02_DEG_close_vs_far/csv_files/DEG_close_vs_far.csv",
    conda:
        "xen-hev"
    params:
        args = X05_ARGS + X05_SAMPLE_ARGS,
    resources:
        mem_mb  = 24000,       # peak 16.3 GB (4 jobs)
        runtime = 15,          # peak 3 min
    log:
        LOGW + "/x05_DEG_close_vs_far.log",
    shell:
        "python {X05_SCRIPTS}/10_HEV_GC_proximity/02_DEG_close_vs_far.py {params.args} > {log} 2>&1"


rule x05_volcano_dotplot:
    """10/03: close vs far DEG volcano (S10I) and dot plot (S10J; only written when significant
    genes exist, so not an output)."""
    input:
        deg = rules.x05_DEG_close_vs_far.output.csv,
    output:
        volcano = X05_HEV_GC + "/03_volcano_dotplot/plots/S10I_DEG_close_vs_far_volcano.pdf",
    conda:
        "shapely"
    params:
        args = X05_ARGS + X05_SAMPLE_ARGS,
    resources:
        mem_mb  = 4000,        # peak 0.8 GB (4 jobs)
        runtime = 15,          # peak 1 min
    log:
        LOGW + "/x05_volcano_dotplot.log",
    shell:
        "python {X05_SCRIPTS}/10_HEV_GC_proximity/03_volcano_dotplot.py {params.args} > {log} 2>&1"


# ── 11 · SPP1 venous specificity ─────────────────────────────────────────────

rule x05_SPP1_across_venous_groups:
    """11/01: SPP1 + venous/HEV markers across four ACKR1+ venous EC groups (S10C)."""
    input:
        h5ad          = x05_all(X05_H5AD),
        labels        = x05_all(x05_labels("subset")),
        cells         = x05_all(X05_TLS_ANN),
        cell_to_niche = x05_all(X05_HEV_C2N),
        info          = x05_all(X05_HEV_INFO),
        colors        = X05_COLORS + "/TLS_HEV_colors.csv",
        sheet         = X05_SHEET,
        metadata      = X05_METADATA,
    output:
        pdf = X05 + "/11_SPP1_venous_specificity/plots/S10C_SPP1_across_venous_groups.pdf",
    conda:
        "xen-hev"
    params:
        args = X05_ARGS + X05_SAMPLE_ARGS,
    resources:
        mem_mb  = 28000,       # peak 21.0 GB (3 jobs)
        runtime = 15,          # peak 3 min
    log:
        LOGW + "/x05_SPP1_across_venous_groups.log",
    shell:
        "python {X05_SCRIPTS}/11_SPP1_venous_specificity/01_SPP1_across_venous_groups.py {params.args} "
        "> {log} 2>&1"


rule x05_SPP1_macrophage_regions:
    """11/02: macrophage SPP1 region zoom (S10E); samples with a
    resources/areas_of_interest/xenium/{sample}_SPP1_macrophages.geojson (P29_Xe5k)."""
    input:
        geojson  = X05_AOI + "/{sample}_SPP1_macrophages.geojson",
        h5ad     = X05_H5AD,
        labels   = x05_labels("subset"),
        image    = X05_HE,
        sdata    = X05_ZARR,
        colors   = X05_COLORS + "/transcript_colors.csv",
        classif  = X05_CLASSIF,
        sheet    = X05_SHEET,
        metadata = X05_METADATA,
    output:
        pdf = X05 + "/11_SPP1_venous_specificity/plots/S10E_{sample}_SPP1_macrophages.pdf",
    wildcard_constraints:
        sample = X05_RE,
    conda:
        "spatdat"
    params:
        args = X05_ARGS,
    resources:
        mem_mb  = 52000,       # peak 40.0 GB (2 jobs)
        runtime = 15,          # peak 5 min
    log:
        LOGW + "/x05_SPP1_macrophage_regions/{sample}.log",
    shell:
        "python {X05_SCRIPTS}/11_SPP1_venous_specificity/02_SPP1_macrophage_regions.py {params.args} "
        "--sample_id {wildcards.sample} > {log} 2>&1"


# ── targets ──────────────────────────────────────────────────────────────────

for x05_target in (
    rules.x05_composition_TLS.output.pdf,
    rules.x05_gene_expression_TLS.output.pdf,
    x05_all(rules.x05_identify_HEV.output.pdf),
    x05_all(rules.x05_calculate_distance_HEV.output.closest),
    x05_all(rules.x05_score_STING_signatures.output.scores),
    expand(rules.x05_GSEA_correlated_genes.output.csv, zip, **X05_GSEA_RUNS),
    x05_all(rules.x05_SPP1_vs_proximal_enrichment.output.corr),
    rules.x05_zonation_overview.output.pdf,
    x05_all(rules.x05_zoom_GCs.output.pdf),
    expand(rules.x05_zoom_zones.output.pdf, sample=X05_ZONE_ZOOM_SAMPLES),
    rules.x05_FDC_overview.output.pdf,
    rules.x05_reference_producers.output.pdf,
    rules.x05_impvl_vs_icaf.output.pdf,
    rules.x05_zone_expression_density.output.pdf,
    rules.x05_chemokine_gradients.output.pdf,
    rules.x05_reference_producer_corr.output.pdf,
    rules.x05_cxcl13_proximity_attribution.output.pdf,
    rules.x05_volcano_dotplot.output.volcano,
    rules.x05_SPP1_across_venous_groups.output.pdf,
    expand(rules.x05_SPP1_macrophage_regions.output.pdf, sample=X05_SPP1_MAC_SAMPLES),
):
    TARGETS.extend([x05_target] if isinstance(x05_target, str) else x05_target)
