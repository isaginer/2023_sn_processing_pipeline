# Libraries
from os.path import join

# Includes
include: "load_samples.smk"

# Variables
ALIGN_DIR = config["decontx_dir"]
CR_LOCATION = "outs/filtered_feature_bc_matrix"
if ("skip_soup" in config.keys() and config["skip_soup"]):
    ALIGN_DIR = config["align_dir"]
STATS_DIR = config["stats_dir"]
SCDF_RESULTS_DIR = join(config["doublets_dir"],"scDblFinder")
DF_RESULTS_DIR = join(config["doublets_dir"],"DoubletFinder")
DF_PARAMS_DIR = join(DF_RESULTS_DIR,"params")
SCTYPE_DIR = join(config["doublets_dir"],"sctype")
GARNETT_DIR = join(config["doublets_dir"],"garnett")
PROCESSED_DIR = join(config["doublets_dir"],"processed")
FILTERED_DIR = join(config["doublets_dir"],"filtered")
PLOTS_DIR = join(config["plots_dir"],"doublets")

plots_list = ["scDblfinder_density", "DoubletFinder_density",
              "sctype_unknown_density", "doublets_overlap", "doublets_pass",
              "umap_sctype", "umap_sctype_unknown"]

all_output_list = [join(PROCESSED_DIR,"{sample}.rds"),
                   join(FILTERED_DIR,"{sample}.rds"),
                   join(STATS_DIR,"doublets_stats.tsv"),
                   *[join(PLOTS_DIR,e,"{sample}.png") for e in plots_list]
                  ]

if ("doublets_integrated" in config.keys() and config["doublets_integrated"]):
    all_output_list.append(join(PROCESSED_DIR,"integrated.qs"))


# Rules
rule doublets_all:
    input:
        expand(all_output_list,
                sample = SAMPLES_DF.index)
    default_target: True


rule doublets_optimize_pk_df:
    conda: 
        config["conda_env"]
    input:
        join(ALIGN_DIR, "{sample}", "decontx.done")
    output:
        join(DF_PARAMS_DIR,"{sample}.rds")
    params:
        mtx_location = CR_LOCATION
    log:
        err="logs/{sample}/optimize_pk_df.err",
        log="logs/{sample}/optimize_pk_df.log"
    script:
        "../../scripts/optimize_pk_df.R"


rule doublets_find_doublets_scdf:
    conda: 
        config["conda_env"]
    input: 
        join(ALIGN_DIR, "{sample}", "decontx.done")
    output: 
        join(SCDF_RESULTS_DIR,"{sample}.rds")
    params:
        mtx_location = CR_LOCATION
    log:
        err="logs/{sample}/find_doublets_scdf.err",
        log="logs/{sample}/find_doublets_scdf.log"
    script:
        "../../scripts/find_doublets_scdf.R"


rule doublets_find_doublets_df:
    conda: 
        config["conda_env"]
    input:
        join(ALIGN_DIR, "{sample}", "decontx.done"),
        rules.doublets_optimize_pk_df.output,
        rules.doublets_find_doublets_scdf.output
    output:
        join(DF_RESULTS_DIR,"{sample}.rds")
    params:
        mtx_location = CR_LOCATION
    log:
        err="logs/{sample}/find_doublets_df.err",
        log="logs/{sample}/find_doublets_df.log"
    script:
        "../../scripts/find_doublets_df.R"


rule doublets_garnett_predict:
    conda: 
        config["conda_env"]
    input: 
        join(ALIGN_DIR, "{sample}", "decontx.done"),
        rules.doublets_find_doublets_df.output,
        rules.doublets_find_doublets_scdf.output
    output: 
        join(GARNETT_DIR,"{sample}.rds")
    params:
        mtx_location = CR_LOCATION
    log:
        err="logs/{sample}/garnett_predict.err",
        log="logs/{sample}/garnett_predict.log"
    script:
        "../../scripts/garnett_predict.R"


rule doublets_sctype_predict:
    conda: 
        config["conda_env"]
    input: 
        join(ALIGN_DIR, "{sample}", "decontx.done")
    output: 
        join(SCTYPE_DIR,"{sample}.rds")
    params:
        mtx_location = CR_LOCATION
    log:
        err="logs/{sample}/sctype_predict.err",
        log="logs/{sample}/sctype_predict.log"
    script:
        "../../scripts/sctype_predict.R"


rule doublets_process:
    conda: 
        config["conda_env"]
    input: 
        join(ALIGN_DIR, "{sample}", "decontx.done"),
        rules.doublets_sctype_predict.output,
        rules.doublets_find_doublets_scdf.output,
        rules.doublets_find_doublets_df.output
    output: 
        processed = join(PROCESSED_DIR,"{sample}.rds"),
        filtered = join(FILTERED_DIR,"{sample}.rds")
    params:
        mtx_location = CR_LOCATION
    log:
        err="logs/{sample}/processed_doublets.err",
        log="logs/{sample}/processed_doublets.log"
    script:
        "../../scripts/processed_doublets.R"


rule doublets_make_plots:
    conda: 
        config["conda_env"]
    input: 
        rules.doublets_process.output.processed
    output: 
        [join(PLOTS_DIR,e,"{sample}.png") for e in plots_list]
    params:
        plots_list = plots_list
    log:
        err="logs/{sample}/make_plots_doublets.log",
        log="logs/{sample}/make_plots_doublets.err"
    script:
        "../../scripts/make_plots_doublets.R"


rule doublets_collect_stats:
    conda: 
        config["conda_env"]
    input: 
        expand(rules.doublets_process.output.processed,
               sample = SAMPLES_DF.index)
    output: 
        join(STATS_DIR,"doublets_stats.tsv")
    params:
        plots_list = plots_list
    log:
        err="logs/collect_stats_doublets.err",
        log="logs/collect_stats_doublets.log"
    script:
        "../../scripts/collect_stats_doublets.R"


rule doublets_integrate:
    conda: 
        config["conda_env"]
    input: 
        expand(rules.doublets_process.output.processed,
               sample = SAMPLES_DF.index)
    output: 
        integrated = join(PROCESSED_DIR,"integrated.qs"),
        features = join(PROCESSED_DIR,"integrated_features.qs"),
        anchors = join(PROCESSED_DIR,"integrated_anchors.qs")
    params:
        plots_list = plots_list
    log:
        err="logs/integrate_doublets.err",
        log="logs/integrate_doublets.log"
    script:
        "../../scripts/integrate_doublets.R"