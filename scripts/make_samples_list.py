import os
import csv

with open("../data/samples.csv", 'w+') as csvfile:
  fieldnames = ['Sample']
  csvwriter = csv.DictWriter(csvfile, fieldnames=fieldnames)
  csvwriter.writeheader()
  for sample_file in os.listdir("../data/cellranger/"):
    csvwriter.writerow({"Sample":sample_file})