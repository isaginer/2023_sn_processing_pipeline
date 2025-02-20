# processed_qc.R
output_log <- file(snakemake@log[[2]], open="wt")
error_log <- file(snakemake@log[[1]], open="wt")
sink(output_log, type = "output")
sink(error_log, type = "message")

suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(dplyr))
set.seed(123)

ndims <- snakemake@config[["processing_ndims"]]

seu <- readRDS(snakemake@input[[1]])
metadata <- readRDS(snakemake@input[[2]])
metadata_2 <- readRDS(snakemake@input[[3]])

seu <- AddMetaData(seu, metadata)
selected_cells <- rownames(seu@meta.data[seu$keep, ])
seu$pass_QC_1 <- "FAIL"
seu@meta.data[selected_cells, "pass_QC_1"] <- "PASS"
seu <- subset(seu, pass_QC_1 == "PASS")

seu <- AddMetaData(seu, metadata_2)

DefaultAssay(seu) <- "RNA"
seu <- SCTransform(seu, method = "glmGamPoi",
                    vars.to.regress = "subsets_Mito_percent",
                    verbose = FALSE) %>%
    RunPCA(verbose = FALSE) %>%
    RunUMAP(dims = 1:ndims, verbose = FALSE) %>%
    FindNeighbors(dims = 1:ndims, verbose = FALSE) %>%
    FindClusters(verbose = FALSE)

seu$qc_cluster_ext_type <- "Not_predicted"
saveRDS(seu, file = snakemake@output[[1]])

seu_filtered <- subset(seu, pass_QC_2 == "PASS")
saveRDS(seu_filtered, file = snakemake@output[[2]])

sink()
sink()