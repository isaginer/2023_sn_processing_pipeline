import os
import json
import csv

with open(snakemake.output["cellranger_stats"], 'w+') as csvfile:
    with open(snakemake.input["csv_metrics"][0], 'r') as metrics_file:
        csvreader = csv.DictReader(metrics_file)
        for row in csvreader:
            fieldnames = ["Sample","Errors"] + list(row.keys())
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
                if line.lstrip(" ").startswith('const data = {'):
                    clean_json = line.strip(" ").replace("const data = ","")
                    summary = json.loads(clean_json)
                    if "alarms" in summary.keys():
                        alarms = []
                        for alarm in summary["alarms"]["alarms"]:
                            if alarm["level"] != 'INFO':
                                alarms.append( "{0}, {1}: {2}, {3}".format(
                                        alarm["level"],
                                        alarm["title"],
                                        alarm["formatted_value"],
                                        alarm["message"]))
                        alarms_text = "\n".join(alarms)
                    break

        # Parse cell ranger csv with metrics
        stats = dict()
        with open(csv_metrics_path, 'r') as metrics_file:
            csvreader = csv.DictReader(metrics_file)
            for row in csvreader:
                stats = row
        csvwriter.writerow({"Sample": sample,
                            "Errors": alarms_text, 
                            **stats})
        
