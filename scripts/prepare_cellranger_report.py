import os
import json
import csv

with open(snakemake.output["cellranger_stats"], 'w+') as csvfile:
    fieldnames = ["Sample", "Errors", "Estimated Number of Cells",
                  "Mean Reads per Cell","Median Genes per Cell",
                  "Number of Reads","Valid Barcodes",
                  "Sequencing Saturation","Q30 Bases in Barcode",
                  "Q30 Bases in RNA Read","Q30 Bases in UMI",
                  "Reads Mapped to Genome","Reads Mapped Confidently to Genome",
                  "Reads Mapped Confidently to Intergenic Regions",
                  "Reads Mapped Confidently to Intronic Regions",
                  "Reads Mapped Confidently to Exonic Regions",
                  "Reads Mapped Confidently to Transcriptome",
                  "Reads Mapped Antisense to Gene",
                  "Fraction Reads in Cells","Total Genes Detected",
                  "Median UMI Counts per Cell"]
    csvwriter = csv.DictWriter(csvfile, fieldnames=fieldnames)
    csvwriter.writeheader()
    samples_list = zip(snakemake.params["sample"],
                       snakemake.input["html_reports"],
                       snakemake.input["csv_metrics"])
    
    for sample_info in samples_list:
        sample, html_report_path, csv_metrics_path = sample_info
        # Extracting warnings/errors from cellranger summary file
        with open(html_report_path,"r") as report_file:
            alarms_text = ""
            report_content = report_file.readlines()
            for line in report_content:
                if line.lstrip(" ").startswith('const data = {"summary"'):
                    clean_json = line.strip(" ").replace("const data = ","")
                    summary = json.loads(clean_json)["summary"]
                    if "alarms" in summary.keys():
                        alarms = []
                        for alarm in summary["alarms"]["alarms"]:
                            if alarm["level"] != 'INFO':
                                alarms.append( "{0}, {1}: {2}, {3}".format(
                                        alarm["level"],
                                        alarm["title"],
                                        alarm["formatted_value"],
                                        alarm["message"]))
                    break
            alarms_text = "\n".join(alarms)

        # Parse cell ranger csv with metrics
        stats = dict()
        with open(csv_metrics_path, 'r') as metrics_file:
            csvreader = csv.DictReader(metrics_file)
            for row in csvreader:
                stats = row
        csvwriter.writerow({"Sample": sample,
                            "Errors": alarms_text, 
                            **stats})
        