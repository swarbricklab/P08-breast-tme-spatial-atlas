# Rules for scripts/04_vascular_microenvironment (Visium HD, 16 µm bins unless stated):
#   01 identify HEV niches     - score signatures, group/expand HEV bins, TLS areas, HEV-TLS distance,
#                                pathology annotations, distance of all bins to HEVs  (per sample)
#   02 identify all vasculature - classify vessels, HEV overlap, vessel distances, PVL association (per sample)
#   03 compare HEV vs venous    - microenvironment, PVL association (per sample); DEGs (all samples)
#   04 ligand-receptor (stLearn) - run per sample at 0 µm (HEV bins only) and 25 µm (HEV + perivascular),
#                                  aggregate per bin size, per-sample plots and LTB composition (25 µm).
#                                  stLearn takes up to ~24 h per sample, so step 04 is in TARGETS_LONG.
#
# All scripts read and write results under --res_dir (config res_dir). Niche signatures come from the published
# list in resources/gene_lists/TLS_niche_signatures/niche_signatures.csv (01_score_signature.py:
# use_resource_signatures = True), not from a rerun of 03_TLS_niche_discovery.

V04 = RES + "/04_vascular_microenvironment"
S04 = "scripts/04_vascular_microenvironment"
HD_H5AD = H5AD + "/visium_hd/{res}/{{sample}}.h5ad"
HD_RE = "|".join(HD_SAMPLES)
HD_SHEET = SHEETS["visium_hd"]
METADATA = "resources/tissue_metadata.csv"
CLASSES = "resources/cell_types/miniatlas/miniatlas_subset_classifications.csv"
ENDO_COLORS = "resources/plotting_colors/endothelial_colors.csv"
PATH_ANNO = "data/pathology_annotations/visium_hd/{sample}_pathology_annotations.csv"
DECON = RES + "/02_cell_type_mapping/deconvolution/visium_hd/{level}_level/square_016um/{sample}_{level}_square_016um_q05_cell_abundance_w_sf.csv"

# Endothelial cell types (file-name form) from the classifications file, as the scripts derive them
_cls = pd.read_csv(CLASSES)
ENDO = sorted(_cls.loc[_cls["miniatlas_major"] == "Endothelial", "miniatlas_minor"].unique())
ENDO = [e.replace(" ", "_") for e in ENDO]

# 01_score_signature: scores x resolutions
SCORES = ["GC", "B_cell_zone", "T_cell_zone", "HEVs"]
SCORE_RES = ["square_048um", "square_016um", "square_008um"]
SCORE_DIR = V04 + "/01_identify_HEV_niches/01_score_signature"

def score_csv(score, res):
    return SCORE_DIR + "/" + score + "_score/" + res + "/csv_files/{sample}_" + score + "_score.csv"

# Per-sample outputs used by several steps
HEV_GROUPS = V04 + "/01_identify_HEV_niches/02b_expand_HEV_groups/square_016um/threshold_5/csv_files/{sample}_HEVs_groups.csv"
HEV_TLS_DIR = V04 + "/01_identify_HEV_niches/03b_calculate_HEV_TLS_distance/square_016um/threshold_5"
HEV_TLS_CLASS = HEV_TLS_DIR + "/csv_files/{sample}_HEV_closest_TLS_distance_classifications.csv"
HEV_DIST_DIR = V04 + "/01_identify_HEV_niches/05_calculate_all_HEV_distances/square_016um"
VESSEL_DIR = V04 + "/02_identify_all_vasculature/01_classify_vessels"
VESSEL_CLASS = VESSEL_DIR + "/csv_files/{sample}_vessel_classification.csv"
OVERLAP_DIR = V04 + "/02_identify_all_vasculature/02_investigate_vessels"
OVERLAP = [OVERLAP_DIR + "/csv_files/{sample}_vessel_HEV_overlap.csv",
           OVERLAP_DIR + "/csv_files/{sample}_vessel_HEV_overlap_spot_level.csv"]
VDIST_DIR = V04 + "/02_identify_all_vasculature/03_calculate_distances"

def min_dist(res):
    return expand(VDIST_DIR + "/" + res + "/csv_files/min_dist/{{sample}}_distance_closest_{e}_per_spot.csv", e=ENDO)

LR_DIR = V04 + "/04_ligand_receptor_analysis"
LR_BINS = ["0", "25"]       # µm; 0 = HEV bins only (Fig. 5d, ED10), 25 = HEV + perivascular (Fig. 4, ED8)
STLEARN_H5AD = LR_DIR + "/01_run_LR_analysis_stlearn/{lr_bin}um_bin/csv_files/{sample}_HEV_LR_analysis_stlearn.h5ad"


# ── 01 · identify HEV niches ─────────────────────────────────────────────────────

rule h04_score_signature:
    """01: Score HEV/TLS zone signatures (GC, B cell zone, T cell zone, HEV) at 48, 16 and 8 µm."""
    input:
        h5ad     = expand(HD_H5AD, res=SCORE_RES),
        sigs     = "resources/gene_lists/TLS_niche_signatures/niche_signatures.csv",   # published signatures (use_resource_signatures = True)
        sheet    = HD_SHEET,
        metadata = METADATA,
    output:
        csv = [score_csv(s, r) for s in SCORES for r in SCORE_RES],
        pdf = expand(SCORE_DIR + "/{s}_score/{r}/plots/{{sample}}_{s}_score_{r}.pdf", s=SCORES, r=SCORE_RES),
    wildcard_constraints:
        sample = HD_RE,
    conda:
        "vis-hd"
    params:
        args = base_args("visium_hd"),
    resources:
        mem_mb  = 16000,       # peak 12.4 GB (4 jobs)
        runtime = 120,         # peak 54 min
    log:
        LOGW + "/h04_score_signature/{sample}.log",
    shell:
        "python {S04}/01_identify_HEV_niches/01_score_signature.py {params.args} "
        "--sample_id {wildcards.sample} > {log} 2>&1"


rule h04_group_HEV_bins:
    """02a: Group contiguous high-HEV-score 16 µm bins (score > 0.5)."""
    input:
        h5ad     = expand(HD_H5AD, res="square_016um"),
        scores   = [score_csv(s, "square_016um") for s in SCORES],
        sheet    = HD_SHEET,
        metadata = METADATA,
    output:
        csv  = V04 + "/01_identify_HEV_niches/02a_group_HEV_bins/square_016um/threshold_5/csv_files/{sample}_HEVs_groups.csv",
        json = V04 + "/01_identify_HEV_niches/02a_group_HEV_bins/square_016um/threshold_5/csv_files/{sample}_grouping_settings.json",
        pdf  = V04 + "/01_identify_HEV_niches/02a_group_HEV_bins/square_016um/threshold_5/plots/groups/{sample}_HEVs_groups.pdf",
    wildcard_constraints:
        sample = HD_RE,
    conda:
        "vis-hd"
    params:
        args = base_args("visium_hd"),
    resources:
        mem_mb  = 8000,        # peak 5.4 GB (4 jobs)
        runtime = 15,          # peak 4 min
    log:
        LOGW + "/h04_group_HEV_bins/{sample}.log",
    shell:
        "python {S04}/01_identify_HEV_niches/02a_group_HEV_bins.py {params.args} "
        "--sample_id {wildcards.sample} > {log} 2>&1"


rule h04_expand_HEV_groups:
    """02b: Seed-and-expand HEV groups to neighbouring bins with HEV score >= 0.3."""
    input:
        h5ad     = expand(HD_H5AD, res=["square_048um", "square_016um"]),
        groups   = rules.h04_group_HEV_bins.output.csv,
        scores   = [score_csv("B_cell_zone", "square_048um"), score_csv("T_cell_zone", "square_048um"),
                    score_csv("HEVs", "square_016um")],
        sheet    = HD_SHEET,
        metadata = METADATA,
    output:
        csv  = HEV_GROUPS,
        json = V04 + "/01_identify_HEV_niches/02b_expand_HEV_groups/square_016um/threshold_5/csv_files/{sample}_grouping_settings.json",
        pdf  = V04 + "/01_identify_HEV_niches/02b_expand_HEV_groups/square_016um/threshold_5/plots/groups/{sample}_HEVs_expanded_groups.pdf",
    wildcard_constraints:
        sample = HD_RE,
    conda:
        "vis-hd"
    params:
        args = base_args("visium_hd"),
    resources:
        mem_mb  = 12000,       # peak 8.4 GB (4 jobs)
        runtime = 15,          # peak 4 min
    log:
        LOGW + "/h04_expand_HEV_groups/{sample}.log",
    shell:
        "python {S04}/01_identify_HEV_niches/02b_expand_HEV_groups.py {params.args} "
        "--sample_id {wildcards.sample} > {log} 2>&1"


rule h04_identify_TLS_areas:
    """03a: Group contiguous 48 µm bins with high B or T cell zone score into TLS areas."""
    input:
        h5ad     = expand(HD_H5AD, res="square_048um"),
        scores   = [score_csv("B_cell_zone", "square_048um"), score_csv("T_cell_zone", "square_048um")],
        sheet    = HD_SHEET,
        metadata = METADATA,
    output:
        csv  = V04 + "/01_identify_HEV_niches/03a_identify_TLS_areas/square_048um/csv_files/{sample}_TLS_zones_groups.csv",
        json = V04 + "/01_identify_HEV_niches/03a_identify_TLS_areas/square_048um/csv_files/{sample}_grouping_settings.json",
        pdf  = V04 + "/01_identify_HEV_niches/03a_identify_TLS_areas/square_048um/plots/groups/{sample}_TLS_zones_groups.pdf",
    wildcard_constraints:
        sample = HD_RE,
    conda:
        "spatial-pilot"
    params:
        args = base_args("visium_hd"),
    resources:
        mem_mb  = 8000,        # peak 4.8 GB (4 jobs)
        runtime = 15,          # peak 2 min
    log:
        LOGW + "/h04_identify_TLS_areas/{sample}.log",
    shell:
        "python {S04}/01_identify_HEV_niches/03a_identify_TLS_areas.py {params.args} "
        "--sample_id {wildcards.sample} > {log} 2>&1"


rule h04_calculate_HEV_TLS_distance:
    """03b: Distance of each HEV to the closest TLS area; intra-/peri-/non-TLS classification (Fig. 3k)."""
    input:
        h5ad     = expand(HD_H5AD, res=["square_048um", "square_016um"]),
        tls      = rules.h04_identify_TLS_areas.output.csv,
        hev      = HEV_GROUPS,
        scores   = [score_csv("B_cell_zone", "square_048um"), score_csv("T_cell_zone", "square_048um"),
                    score_csv("HEVs", "square_016um")],
        sheet    = HD_SHEET,
        metadata = METADATA,
    output:
        classes = HEV_TLS_CLASS,
        per_tls = HEV_TLS_DIR + "/csv_files/{sample}_HEV_per_TLS_group.csv",
        pdf     = HEV_TLS_DIR + "/plots/{sample}_HEVs_TLS_classification.pdf",
    wildcard_constraints:
        sample = HD_RE,
    conda:
        "vis-hd"
    params:
        args = base_args("visium_hd"),
    resources:
        mem_mb  = 12000,       # peak 7.3 GB (4 jobs)
        runtime = 15,          # peak 2 min
    log:
        LOGW + "/h04_calculate_HEV_TLS_distance/{sample}.log",
    shell:
        "python {S04}/01_identify_HEV_niches/03b_calculate_HEV_TLS_distance.py {params.args} "
        "--sample_id {wildcards.sample} > {log} 2>&1"


rule h04_assess_HEV_pathology_annotations:
    """04: Pathology annotation (8 µm) under each HEV."""
    input:
        h5ad     = expand(HD_H5AD, res=["square_016um", "square_008um"]),
        hev      = HEV_GROUPS,
        classes  = HEV_TLS_CLASS,
        anno     = PATH_ANNO,
        sheet    = HD_SHEET,
        metadata = METADATA,
    output:
        dominant = V04 + "/01_identify_HEV_niches/04_assess_HEV_pathology_annotations/csv_files/{sample}_dominant_pathology_per_HEV_group.csv",
        counts   = V04 + "/01_identify_HEV_niches/04_assess_HEV_pathology_annotations/csv_files/{sample}_pathology_counts_per_HEV_group.csv",
        pdf      = V04 + "/01_identify_HEV_niches/04_assess_HEV_pathology_annotations/plots/{sample}_HEVs_per_pathology_annotation.pdf",
    wildcard_constraints:
        sample = HD_RE,
    conda:
        "vis-hd"
    params:
        args = base_args("visium_hd"),
    resources:
        mem_mb  = 16000,       # peak 10.6 GB (4 jobs)
        runtime = 15,          # peak 6 min
    log:
        LOGW + "/h04_assess_HEV_pathology_annotations/{sample}.log",
    shell:
        "python {S04}/01_identify_HEV_niches/04_assess_HEV_pathology_annotations.py {params.args} "
        "--sample_id {wildcards.sample} > {log} 2>&1"


rule h04_calculate_all_HEV_distances:
    """05: Distance of every 16 µm bin to each HEV (input for the stLearn neighbourhoods)."""
    input:
        h5ad     = expand(HD_H5AD, res="square_016um"),
        hev      = HEV_GROUPS,
        sheet    = HD_SHEET,
        metadata = METADATA,
    output:
        all_spots = HEV_DIST_DIR + "/csv_files/{sample}_HEV_distances_all_spots.csv",
        closest   = HEV_DIST_DIR + "/csv_files/{sample}_distance_closest_HEV_per_spot.csv",
        pdf       = HEV_DIST_DIR + "/plots/{sample}_HEV_distances.pdf",
    wildcard_constraints:
        sample = HD_RE,
    conda:
        "vis-hd"
    params:
        args = base_args("visium_hd"),
    resources:
        mem_mb  = 28000,       # peak 19.1 GB (4 jobs)
        runtime = 45,          # peak 18 min
    log:
        LOGW + "/h04_calculate_all_HEV_distances/{sample}.log",
    shell:
        "python {S04}/01_identify_HEV_niches/05_calculate_all_HEV_distances.py {params.args} "
        "--sample_id {wildcards.sample} > {log} 2>&1"


# ── 02 · identify all vasculature ────────────────────────────────────────────────

rule h04_classify_vessels:
    """01: Segment vessels by dominant deconvolved endothelial type (Fig. 3b-d)."""
    input:
        h5ad     = expand(HD_H5AD, res="square_016um"),
        decon    = expand(DECON, level="minor", allow_missing=True),
        classes  = CLASSES,
        colors   = ENDO_COLORS,
        sheet    = HD_SHEET,
        metadata = METADATA,
    output:
        classification = VESSEL_CLASS,
        size           = VESSEL_DIR + "/csv_files/{sample}_vessel_size.csv",
        numbers        = VESSEL_DIR + "/csv_files/{sample}_vessel_numbers.csv",
        conflict_spot  = VESSEL_DIR + "/csv_files/{sample}_vessel_conflict_spot_level.csv",
        conflict_ves   = VESSEL_DIR + "/csv_files/{sample}_vessel_conflict_vessel_level.csv",
        endo_numbers   = VESSEL_DIR + "/csv_files/{sample}_endothelial_numbers_normalized.csv",
        pdf            = VESSEL_DIR + "/plots/{sample}_vessel_classification.pdf",
    wildcard_constraints:
        sample = HD_RE,
    conda:
        "vis-hd"
    params:
        args = base_args("visium_hd"),
    resources:
        mem_mb  = 12000,       # peak 7.9 GB (4 jobs)
        runtime = 30,          # peak 9 min
    log:
        LOGW + "/h04_classify_vessels/{sample}.log",
    shell:
        "python {S04}/02_identify_all_vasculature/01_classify_vessels.py {params.args} "
        "--sample_id {wildcards.sample} > {log} 2>&1"


rule h04_investigate_vessels:
    """02: Overlap of each vessel with HEVs (Fig. 3g)."""
    input:
        h5ad     = expand(HD_H5AD, res="square_016um"),
        decon    = expand(DECON, level="minor", allow_missing=True),
        vessels  = VESSEL_CLASS,
        hev      = HEV_GROUPS,
        classes  = CLASSES,
        colors   = ENDO_COLORS,
        sheet    = HD_SHEET,
        metadata = METADATA,
    output:
        overlap      = OVERLAP[0],
        overlap_spot = OVERLAP[1],
        pdf          = OVERLAP_DIR + "/plots/{sample}_vessel_HEV_overlap.pdf",
    wildcard_constraints:
        sample = HD_RE,
    conda:
        "vis-hd"
    params:
        args = base_args("visium_hd"),
    resources:
        mem_mb  = 8000,        # peak 3.8 GB (4 jobs)
        runtime = 15,          # peak 1 min
    log:
        LOGW + "/h04_investigate_vessels/{sample}.log",
    shell:
        "python {S04}/02_identify_all_vasculature/02_investigate_vessels.py {params.args} "
        "--sample_id {wildcards.sample} > {log} 2>&1"


rule h04_calculate_vessel_distances:
    """03: Distance of 8 and 16 µm bins to each vessel type (all_dist only written for 16 µm)."""
    input:
        h5ad     = expand(HD_H5AD, res=["square_016um", "square_008um"]),
        vessels  = VESSEL_CLASS,
        classes  = CLASSES,
        sheet    = HD_SHEET,
        metadata = METADATA,
    output:
        min_dist_016 = min_dist("square_016um"),
        min_dist_008 = min_dist("square_008um"),
        all_dist_016 = expand(VDIST_DIR + "/square_016um/csv_files/all_dist/{{sample}}_{e}_distances_all_spots.csv", e=ENDO),
        pdf          = expand(VDIST_DIR + "/{res}/plots/{{sample}}_vessel_distances.pdf", res=["square_008um", "square_016um"]),
    wildcard_constraints:
        sample = HD_RE,
    conda:
        "vis-hd"
    params:
        args = base_args("visium_hd"),
    resources:
        mem_mb  = 168000,      # peak 130.0 GB (4 jobs)
        runtime = 105,         # peak 51 min
    log:
        LOGW + "/h04_calculate_vessel_distances/{sample}.log",
    shell:
        "python {S04}/02_identify_all_vasculature/03_calculate_distances.py {params.args} "
        "--sample_id {wildcards.sample} > {log} 2>&1"


rule h04_assess_PVL_association:
    """04: PVL abundance across distance bins from each vessel type (Fig. 3f, ED6e)."""
    input:
        h5ad     = expand(HD_H5AD, res="square_016um"),
        decon    = expand(DECON, level="minor", allow_missing=True),
        vessels  = VESSEL_CLASS,
        min_dist = min_dist("square_016um"),
        classes  = CLASSES,
        colors   = ENDO_COLORS,
        sheet    = HD_SHEET,
        metadata = METADATA,
    output:
        csv = V04 + "/02_identify_all_vasculature/04_assess_PVL_association/csv_files/{sample}_PVL_composition_vs_distance.csv.gz",
        pdf = V04 + "/02_identify_all_vasculature/04_assess_PVL_association/plots/{sample}_all_vessels_PVL_composition_vs_distance.pdf",
    wildcard_constraints:
        sample = HD_RE,
    conda:
        "vis-hd"
    params:
        args = base_args("visium_hd"),
    resources:
        mem_mb  = 16000,       # peak 11.2 GB (4 jobs)
        runtime = 15,          # peak 6 min
    log:
        LOGW + "/h04_assess_PVL_association/{sample}.log",
    shell:
        "python {S04}/02_identify_all_vasculature/04_assess_PVL_association.py {params.args} "
        "--sample_id {wildcards.sample} > {log} 2>&1"


# ── 03 · compare HEVs with other venous vessels (ACKR1) ──────────────────────────

rule h04_compare_microenvironment:
    """01: Pathology annotations across distance from HEV vs non-HEV ACKR1 vessels (Fig. 3h)."""
    input:
        h5ad     = expand(HD_H5AD, res=["square_016um", "square_008um"]),
        vessels  = VESSEL_CLASS,
        min_dist = min_dist("square_008um"),
        anno     = PATH_ANNO,
        overlap  = OVERLAP,
        classes  = CLASSES,
        colors   = ENDO_COLORS,
        sheet    = HD_SHEET,
        metadata = METADATA,
    output:
        csv = expand(V04 + "/03_compare_HEV_to_other_venous/01_compare_microenvironment/csv_files/{{sample}}_ACKR1_{grp}_composition_vs_distance.csv",
                     grp=["no_overlap", "overlap"]),
        pdf = V04 + "/03_compare_HEV_to_other_venous/01_compare_microenvironment/plots/{sample}_ACKR1_composition_vs_distance.pdf",
    wildcard_constraints:
        sample = HD_RE,
    conda:
        "vis-hd"
    params:
        args = base_args("visium_hd"),
    resources:
        mem_mb  = 96000,       # peak 74.8 GB (4 jobs)
        runtime = 45,          # peak 20 min
    log:
        LOGW + "/h04_compare_microenvironment/{sample}.log",
    shell:
        "python {S04}/03_compare_HEV_to_other_venous/01_compare_microenvironment.py {params.args} "
        "--sample_id {wildcards.sample} > {log} 2>&1"


rule h04_compare_gene_expression:
    """02: DEGs HEV vs non-HEV ACKR1 vessel bins, all samples (Fig. 3i, Table S8; Fig. 5d)."""
    input:
        h5ad     = expand(H5AD + "/visium_hd/square_016um/{sample}.h5ad", sample=HD_SAMPLES),
        vessels  = expand(VESSEL_CLASS, sample=HD_SAMPLES),
        overlap  = expand(OVERLAP, sample=HD_SAMPLES),
        classes  = CLASSES,
        colors   = ENDO_COLORS,
        sheet    = HD_SHEET,
        metadata = METADATA,
    output:
        degs  = V04 + "/03_compare_HEV_to_other_venous/02_compare_gene_expression/csv_files/All_samples_rank_genes_HEV_overlap.csv",
        h5ad  = V04 + "/03_compare_HEV_to_other_venous/02_compare_gene_expression/csv_files/All_samples_HEV_ACKR1.h5ad",
        pdf   = V04 + "/03_compare_HEV_to_other_venous/02_compare_gene_expression/plots/All_samples_ACKR1_HEV_overlap_DEGs.pdf",
        pdf_s = V04 + "/03_compare_HEV_to_other_venous/02_compare_gene_expression/plots/All_samples_ACKR1_HEV_overlap_spatial.pdf",
    conda:
        "vis-hd"
    params:
        args    = base_args("visium_hd"),
        samples = " ".join(HD_SAMPLES),
    resources:
        mem_mb  = 12000,       # peak 7.7 GB (1 jobs)
        runtime = 15,          # peak 1 min
    log:
        LOGW + "/h04_compare_gene_expression.log",
    shell:
        "python {S04}/03_compare_HEV_to_other_venous/02_compare_gene_expression.py {params.args} "
        "--samples {params.samples} > {log} 2>&1"


rule h04_compare_PVL_association:
    """03: PVL abundance near HEV vs non-HEV ACKR1 vessels (Fig. 3j)."""
    input:
        h5ad     = expand(HD_H5AD, res="square_016um"),
        decon    = expand(DECON, level="minor", allow_missing=True),
        vessels  = VESSEL_CLASS,
        overlap  = OVERLAP,
        classes  = CLASSES,
        colors   = ENDO_COLORS,
        sheet    = HD_SHEET,
        metadata = METADATA,
    output:
        csv = V04 + "/03_compare_HEV_to_other_venous/03_compare_PVL_association/csv_files/{sample}_PVL_composition_vs_distance.csv",
        pdf = V04 + "/03_compare_HEV_to_other_venous/03_compare_PVL_association/plots/{sample}_HEV_vs_non_HEV_PVL_composition_vs_distance.pdf",
    wildcard_constraints:
        sample = HD_RE,
    conda:
        "vis-hd"
    params:
        args = base_args("visium_hd"),
    resources:
        mem_mb  = 28000,       # peak 21.6 GB (4 jobs)
        runtime = 15,          # peak 4 min
    log:
        LOGW + "/h04_compare_PVL_association/{sample}.log",
    shell:
        "python {S04}/03_compare_HEV_to_other_venous/03_compare_PVL_association.py {params.args} "
        "--sample_id {wildcards.sample} > {log} 2>&1"


# ── 04 · ligand-receptor analysis (stLearn; long) ────────────────────────────────

rule h04_run_LR_analysis_stlearn:
    """01: stLearn LR analysis per sample in HEV neighbourhoods; lr_bin 0 = HEV bins only,
    25 = HEV + perivascular bins (~24 h for the largest sample)."""
    input:
        h5ad      = expand(HD_H5AD, res="square_016um"),
        hev       = HEV_GROUPS,
        distances = rules.h04_calculate_all_HEV_distances.output.all_spots,
        classes   = HEV_TLS_CLASS,
        celltypes = CLASSES,
        sheet     = HD_SHEET,
        metadata  = METADATA,
    output:
        h5ad       = STLEARN_H5AD,
        summary    = LR_DIR + "/01_run_LR_analysis_stlearn/{lr_bin}um_bin/csv_files/{sample}_LR_summary.csv",
        lrfeatures = LR_DIR + "/01_run_LR_analysis_stlearn/{lr_bin}um_bin/csv_files/{sample}_HEV_LR_analysis_uns_lrfeatures.csv",
        pdf        = LR_DIR + "/01_run_LR_analysis_stlearn/{lr_bin}um_bin/plots/{sample}_HEV_LR_analysis.pdf",
    wildcard_constraints:
        sample = HD_RE,
        lr_bin = "|".join(LR_BINS),
    conda:
        "vis-hd"
    params:
        args = base_args("visium_hd"),
    threads: 16
    resources:
        mem_mb  = 128000,
        runtime = 2880,
    log:
        LOGW + "/h04_run_LR_analysis_stlearn/{sample}_{lr_bin}um_bin.log",
    shell:
        "python {S04}/04_ligand_receptor_analysis/01_run_LR_analysis_stlearn.py {params.args} "
        "--sample_id {wildcards.sample} --bin_size {wildcards.lr_bin} > {log} 2>&1"


rule h04_analyze_LR_all:
    """02: Aggregate stLearn results across samples per bin size (Fig. 4g-h, ED8d,g at 25 µm;
    Fig. 5d, ED10a-b at 0 µm)."""
    input:
        h5ad  = expand(STLEARN_H5AD, sample=HD_SAMPLES, allow_missing=True),
        hev   = expand(HEV_GROUPS, sample=HD_SAMPLES),
        sheet = HD_SHEET,
    output:
        summary      = LR_DIR + "/02_analyze_LR_all/{lr_bin}um_bin/csv_files/All_samples_LR_summary.csv",
        group_counts = LR_DIR + "/02_analyze_LR_all/{lr_bin}um_bin/csv_files/All_samples_group_counts.csv",
        spot_counts  = LR_DIR + "/02_analyze_LR_all/{lr_bin}um_bin/csv_files/All_samples_spot_counts.csv",
        pdf          = LR_DIR + "/02_analyze_LR_all/{lr_bin}um_bin/plots/All_samples_HEV_LR_analysis.pdf",
    wildcard_constraints:
        lr_bin = "|".join(LR_BINS),
    conda:
        "shapely"
    params:
        args    = base_args("visium_hd"),
        samples = " ".join(HD_SAMPLES),
    resources:
        mem_mb  = 32000,
        runtime = 60,
    log:
        LOGW + "/h04_analyze_LR_all/{lr_bin}um_bin.log",
    shell:
        "python {S04}/04_ligand_receptor_analysis/02_analyze_LR_all.py {params.args} "
        "--samples {params.samples} --bin_size {wildcards.lr_bin} > {log} 2>&1"


rule h04_plot_LR_per_sample:
    """03: Per-sample LR plots and HEV spot grouping (only the 25 µm run is needed; Fig. 4i-j)."""
    input:
        h5ad         = expand(HD_H5AD, res="square_016um"),
        stlearn      = STLEARN_H5AD,
        group_counts = rules.h04_analyze_LR_all.output.group_counts,
        classes      = CLASSES,
        sheet        = HD_SHEET,
        metadata     = METADATA,
    output:
        grouping = LR_DIR + "/03_plot_LR_per_sample/{lr_bin}um_bin/csv_files/{sample}_HEV_spot_grouping.csv",
        h5ad     = LR_DIR + "/03_plot_LR_per_sample/{lr_bin}um_bin/csv_files/{sample}_HEV_LR_analysis_stlearn_spatial_groups.h5ad",
        pdf      = LR_DIR + "/03_plot_LR_per_sample/{lr_bin}um_bin/plots/{sample}_HEV_LR_gene_expression.pdf",
        pdf_hres = LR_DIR + "/03_plot_LR_per_sample/{lr_bin}um_bin/plots/{sample}_spot_labels_grouping_hires.pdf",
    wildcard_constraints:
        sample = HD_RE,
        lr_bin = "25",
    conda:
        "vis-hd"
    params:
        args = base_args("visium_hd"),
    resources:
        mem_mb  = 64000,
        runtime = 120,
    log:
        LOGW + "/h04_plot_LR_per_sample/{sample}_{lr_bin}um_bin.log",
    shell:
        "python {S04}/04_ligand_receptor_analysis/03_plot_LR_per_sample.py {params.args} "
        "--sample_id {wildcards.sample} --bin_size {wildcards.lr_bin} > {log} 2>&1"


rule h04_check_composition_LTB_signaling_bins:
    """04: Cell type composition of bins sending LTB to HEV LTBR (25 µm run; Fig. 4j)."""
    input:
        stlearn  = expand(STLEARN_H5AD, sample=HD_SAMPLES, lr_bin="25"),
        decon    = expand(DECON, sample=HD_SAMPLES, level="subset"),
        grouping = expand(rules.h04_plot_LR_per_sample.output.grouping, sample=HD_SAMPLES, lr_bin="25"),
        classes  = CLASSES,
        colors   = "resources/plotting_colors/miniatlas_major_colors.csv",
        sheet    = HD_SHEET,
        metadata = METADATA,
    output:
        csv   = LR_DIR + "/04_check_composition_LTB_signaling_bins/csv_files/All_samples_HEV_LTB_LTBR_cellular_enrichment.csv",
        stats = LR_DIR + "/04_check_composition_LTB_signaling_bins/csv_files/All_samples_HEV_LTB_LTBR_cellular_enrichment_statistics.csv",
        pdf   = LR_DIR + "/04_check_composition_LTB_signaling_bins/plots/All_samples_HEV_LTB_LTBR_cellular_enrichment.pdf",
    conda:
        "vis-hd"
    params:
        args    = base_args("visium_hd"),
        samples = " ".join(HD_SAMPLES),
    resources:
        mem_mb  = 32000,
        runtime = 60,
    log:
        LOGW + "/h04_check_composition_LTB_signaling_bins.log",
    shell:
        "python {S04}/04_ligand_receptor_analysis/04_check_composition_LTB_signaling_bins.py {params.args} "
        "--samples {params.samples} > {log} 2>&1"


# ── targets ─────────────────────────────────────────────────────────────────────

TARGETS += expand(rules.h04_score_signature.output.csv, sample=HD_SAMPLES)
TARGETS += expand(rules.h04_calculate_HEV_TLS_distance.output, sample=HD_SAMPLES)
TARGETS += expand(rules.h04_assess_HEV_pathology_annotations.output, sample=HD_SAMPLES)
TARGETS += expand(rules.h04_calculate_all_HEV_distances.output, sample=HD_SAMPLES)
TARGETS += expand(rules.h04_classify_vessels.output, sample=HD_SAMPLES)
TARGETS += expand(rules.h04_investigate_vessels.output, sample=HD_SAMPLES)
TARGETS += expand(rules.h04_assess_PVL_association.output, sample=HD_SAMPLES)
TARGETS += expand(rules.h04_compare_microenvironment.output, sample=HD_SAMPLES)
TARGETS += rules.h04_compare_gene_expression.output
TARGETS += expand(rules.h04_compare_PVL_association.output, sample=HD_SAMPLES)

TARGETS_LONG += expand(rules.h04_run_LR_analysis_stlearn.output, sample=HD_SAMPLES, lr_bin=LR_BINS)
TARGETS_LONG += expand(rules.h04_analyze_LR_all.output, lr_bin=LR_BINS)
TARGETS_LONG += expand(rules.h04_plot_LR_per_sample.output, sample=HD_SAMPLES, lr_bin="25")
TARGETS_LONG += rules.h04_check_composition_LTB_signaling_bins.output
