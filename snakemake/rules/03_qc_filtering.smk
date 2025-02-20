# Libraries
from os.path import join

# Includes
include: "load_samples.smk"

# Variables
DBL_FILTERED_DIR = join(config["doublets_dir"],"filtered")
MAD_DIR = join(config["qc_dir"],"mad_qc")
MAD_DIR_FILTERED = join(MAD_DIR, "filtered")
CELLTYPES_DIR = join(config["qc_dir"],"celltypes_qc")

FILTERED_DIR = join(config["qc_dir"],"filtered")
PROCESSED_DIR = join(config["qc_dir"],"processed")
STATS_DIR = config["stats_dir"]
PLOTS_DIR = join(config["plots_dir"],"qc")

plots_list = ["umi_vs_genes", "mt_vs_nuclei", "mt_vs_umi", 
              "umap_log_umi", "umap_umi_pass", "umap_ngenes", 
              "umap_ngenes_pass", "umap_garnett", "qc_pass"]

all_output_list = [join(MAD_DIR,"{sample}.rds"),
                   join(CELLTYPES_DIR,"{sample}.rds"),
                   join(FILTERED_DIR,"{sample}.rds"),
                   join(STATS_DIR,"qc_stats.tsv"),
                   *[join(PLOTS_DIR,e,"{sample}.png") for e in plots_list]
                  ]

if ("qc_integrated" in config.keys() and config["qc_integrated"]):
    all_output_list.append(join(PROCESSED_DIR,"integrated.qs"))


# Rules
rule qc_all:
    input:
        expand(all_output_list,
                sample = SAMPLES_DF.index)
    default_target: True


rule filter_samples_by_QC:
    conda: config["conda_env"]
    input:
        join(DBL_FILTERED_DIR,"{sample}.rds")
    output:
        join(MAD_DIR,"{sample}.rds")
    log:
        err="logs/{sample}/filter_samples_by_QC.err",
        log="logs/{sample}/filter_samples_by_QC.log"
    script:
        "../../scripts/get_cells_QC.R"

rule qc_by_cell_types:
    conda: 
        config["conda_env"]
    input: 
        join(DBL_FILTERED_DIR,"{sample}.rds"),
        join(MAD_DIR,"{sample}.rds")
    output: 
        qc_step_1 = join(MAD_DIR_FILTERED,"{sample}.rds"),
        qc_step_2 = join(CELLTYPES_DIR,"{sample}.rds")
    log:
        err="logs/{sample}/cell_types_qc.err",
        log="logs/{sample}/cell_types_qc.log"
    script:
        "../../scripts/cell_types_qc.R"  

rule qc_process:
    conda: 
        config["conda_env"]
    input: 
        join(DBL_FILTERED_DIR,"{sample}.rds"),
        join(MAD_DIR,"{sample}.rds"),
        join(CELLTYPES_DIR,"{sample}.rds")
    output: 
        processed = join(PROCESSED_DIR,"{sample}.rds"),
        filtered = join(FILTERED_DIR,"{sample}.rds")
    log:
        err="logs/{sample}/processed_qc.err",
        log="logs/{sample}/processed_qc.log"
    script:
        "../../scripts/processed_qc.R"

rule qc_make_plots:
    conda: 
        config["conda_env"]
    input: 
        qc_step_1 = rules.qc_by_cell_types.output.qc_step_1
    output: 
        [join(PLOTS_DIR,e,"{sample}.png") for e in plots_list]
    params:
        plots_list = plots_list
    log:
        err="logs/{sample}/make_plots_qc.err",
        log="logs/{sample}/make_plots_qc.log"
    script:
        "../../scripts/make_plots_qc.R"


rule qc_collect_stats:
    conda: 
        config["conda_env"]
    input: 
        step_1 = expand(rules.qc_by_cell_types.output.qc_step_1,
               sample = SAMPLES_DF.index),
        step_2 = expand(rules.qc_process.output.processed,
               sample = SAMPLES_DF.index)
    output: 
        join(STATS_DIR,"qc_stats.tsv")
    params:
        plots_list = plots_list
    log:
        err="logs/collect_stats_qc.err",
        log="logs/collect_stats_qc.log"
    script:
        "../../scripts/collect_stats_qc.R"


rule qc_integrate:
    conda: 
        config["conda_env"]
    input: 
        expand(rules.qc_process.output.processed,
               sample = SAMPLES_DF.index)
    output: 
        integrated = join(PROCESSED_DIR,"integrated.qs"),
        features = join(PROCESSED_DIR,"integrated_features.qs"),
        anchors = join(PROCESSED_DIR,"integrated_anchors.qs")
    log:
        err="logs/integrate_qc.err",
        log="logs/integrate_qc.log"
    script:
        "../../scripts/integrate_qc.R"
