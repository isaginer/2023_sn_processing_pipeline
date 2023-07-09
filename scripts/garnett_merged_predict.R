suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(garnett))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(org.Hs.eg.db))
suppressPackageStartupMessages(library(qs))
source("scripts/utils/custom_predict.R")
set.seed(123)

model_path <- snakemake@config[["model_path"]]
garnett_cluster <- snakemake@config[["garnett_cluster"]]

merged.integrated <- qread(snakemake@input[[1]])
load(model_path)

counts <- GetAssayData(merged.integrated, assay = "SCT", slot = "counts")
meta <- merged.integrated@meta.data
fdata <- (rownames(counts))
fdata <- as.data.frame(fdata)
colnames(fdata) <- "gene_short_name"
rownames(fdata) <- fdata$gene_short_name
row.names(counts) <- row.names(fdata)
colnames(counts) <- row.names(meta)
pdata <- meta

pd <- new("AnnotatedDataFrame", data = pdata)
fd <- new("AnnotatedDataFrame", data = fdata)
cds_seu <- newCellDataSetProb(counts, phenoData = pd, featureData = fd)
cds_seu <- estimateSizeFactors(cds_seu)

pData(cds_seu)$garnett_cluster <- pData(cds_seu)[, garnett_cluster]

predicted_cds <- classify_cells(cds_seu, trained_classifier,
                                db = org.Hs.eg.db,
                                cluster_extend = TRUE,
                                cds_gene_id_type = "SYMBOL",
                                verbose = TRUE)

merged.integrated <- AddMetaData(merged.integrated,
                   pData(predicted_cds)[c("Size_Factor", "cell_type", "cluster_ext_type")])
qsave(merged.integrated, file = snakemake@output[[1]])