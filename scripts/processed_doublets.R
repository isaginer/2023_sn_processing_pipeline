output_log <- file(snakemake@log[[2]], open="wt")
error_log <- file(snakemake@log[[1]], open="wt")
sink(output_log, type = "output")
sink(error_log, type = "message")

suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(dplyr))
set.seed(123)

min_cells <- snakemake@config[["min_cells"]]
ndims <- snakemake@config[["processing_ndims"]]

sample_path <- snakemake@input[[1]]
metadata <- readRDS(snakemake@input[[2]])
metadata_2 <- as.data.frame(readRDS(snakemake@input[[3]]))
metadata_3 <- as.data.frame(readRDS(snakemake@input[[4]]))

counts <- Read10X(file.path(dirname(snakemake@input[[1]]),
                            snakemake@params[["mtx_location"]]))

if (length(names(counts)) > 1) {
  counts <- counts$`Gene Expression`
}
seu <- CreateSeuratObject(counts,
                          min.cells = min_cells,
                          project = snakemake@wildcards$sample)

seu <- AddMetaData(seu, metadata)
seu <- AddMetaData(seu, metadata_2)
seu <- AddMetaData(seu, metadata_3)

seu <- PercentageFeatureSet(seu,
                            pattern = "^(MT|mt)-",
                            col.name = "percent.mt") %>%
    SCTransform(method = "glmGamPoi",
                vars.to.regress = "percent.mt", verbose = FALSE) %>%
    RunPCA(verbose = FALSE) %>%
    RunUMAP(dims = 1:ndims, verbose = FALSE) %>%
    FindNeighbors(dims = 1:ndims, verbose = FALSE) %>%
    FindClusters(verbose = FALSE)

seu$garnett_prediction <- "Not_predicted"
seu$merged_doublets <- paste0(seu$doublet_prc, "_", seu$scDblFinder.class)
selected_cells <- rownames(seu@meta.data[(seu$merged_doublets == "Singlet_singlet"), ])
seu$pass_doublets_QC <- "FAIL"
seu@meta.data[selected_cells, "pass_doublets_QC"] <- "PASS"
saveRDS(seu, file = snakemake@output[[1]])

seu_filtered <- subset(seu, pass_doublets_QC == "PASS")
saveRDS(seu_filtered, file = snakemake@output[[2]])

sink()
