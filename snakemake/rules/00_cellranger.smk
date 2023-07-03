include: "load_samples.smk"

rule cellranger_all:
    input:
        os.path.join(config["stats_dir"],"cellranger_stats.csv")
    default_target: True

rule cellranger_collect_statistics:
    input:
        html_reports = expand(os.path.join(config["align_dir"],"{sample}","web_summary.html"),
        sample = SAMPLES_DF.index),
        csv_metrics = expand(os.path.join(config["align_dir"],"{sample}","outs","metrics_summary.csv"),
        sample = SAMPLES_DF.index)
    output:
        cellranger_stats = os.path.join(config["stats_dir"],"cellranger_stats.csv")
    params:
        sample = SAMPLES_DF.index
    log:
        "logs/cellranger_collect_statistics.log"
    script:
        "../../scripts/prepare_cellranger_report.py"