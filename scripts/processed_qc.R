# processed_qc.R

suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(dplyr))
set.seed(123)

model_path <- snakemake@config[["model_path"]]
load(model_path)

min_cells <- snakemake@config[["min_cells"]]
ndims <- snakemake@config[["processing_ndims"]]

sample_path <- snakemake@input[[1]]
metadata <- readRDS(snakemake@input[[2]])

counts <- Read10X(sample_path)
if (length(names(counts)) > 1) {
  counts <- counts$`Gene Expression`
}
seu <- CreateSeuratObject(counts,
                          min.cells = min_cells,
                          project = snakemake@wildcards$sample)

seu <- AddMetaData(seu, metadata)

seu <- PercentageFeatureSet(seu,
                            pattern = "^[MT|mt]-",
                            col.name = "percent.mt") %>%
    SCTransform(method = "glmGamPoi",
                vars.to.regress = "percent.mt", verbose = FALSE) %>%
    RunPCA(verbose = FALSE) %>%
    RunUMAP(dims = 1:ndims, verbose = FALSE) %>%
    FindNeighbors(dims = 1:ndims, verbose = FALSE) %>%
    FindClusters(verbose = FALSE)

counts <- GetAssayData(seu, assay = "SCT", slot = "counts")
meta <- seu@meta.data
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
pData(cds_seu)$garnett_cluster <- pData(cds_seu)[, "seurat_clusters"]

predicted_cds <- classify_cells(cds_seu, trained_classifier,
                                db = org.Hs.eg.db,
                                cluster_extend = TRUE,
                                cds_gene_id_type = "SYMBOL",
                                verbose = TRUE)
predictions <- pData(predicted_cds)[,c("cluster_ext_type")]
colnames(predictions) <- paste0("qc_",colnames(predictions))
seu <- AddMetaData(seu, predictions)

selected_cells <- rownames(seu@meta.data[seu$keep, ])
seu$pass_QC <- "PASS"
seu@meta.data[selected_cells,"pass_QC"] <- "FAIL"
saveRDS(seu, file = snakemake@output[[1]])

seu_filtered <- subset(seu, pass_QC == "PASS")
saveRDS(seu_filtered, file = snakemake@output[[2]])