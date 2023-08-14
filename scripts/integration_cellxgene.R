suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(sceasy))
suppressPackageStartupMessages(library(reticulate))
suppressPackageStartupMessages(library(qs))
set.seed(123)

merged.integrated <- qread(snakemake@input[[1]])
sceasy::convertFormat(merged.integrated, from = "seurat", to = "anndata",
                       outFile = snakemake@output[[1]])