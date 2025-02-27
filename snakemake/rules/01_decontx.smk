# Libraries
from os.path import join

# Includes
include: "load_samples.smk"

# Variables
ALIGN_DIR = config["align_dir"]
FILTERED_LOCATION = "outs/filtered_feature_bc_matrix"
RAW_LOCATION = "outs/raw_feature_bc_matrix"
DECONTX_DIR = config["decontx_dir"]

#Rules
rule decontx_all:
    input: 
        expand(join(DECONTX_DIR, "{sample}", "decontx.done"),
                sample = SAMPLES_DF.index)
    default_target: True

rule decontx_run_decontx:
    conda:  
        config["conda_env"]
    input:
        join(ALIGN_DIR, "{sample}", FILTERED_LOCATION),
        join(ALIGN_DIR, "{sample}", RAW_LOCATION)
    output:
        join(DECONTX_DIR, "{sample}", "decontx.done"),
        join(DECONTX_DIR, "{sample}", "unfiltered", "{sample}.rds"),
    log:
        err="logs/{sample}/run_decontx.err",
        log="logs/{sample}/run_decontx.log"
    script:
        "../../scripts/run_decontx.R"
