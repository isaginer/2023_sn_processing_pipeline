output_log <- file(snakemake@log[[2]], open="wt")
error_log <- file(snakemake@log[[1]], open="wt")
sink(output_log, type = "output")
sink(error_log, type = "message")

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(garnett))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(org.Hs.eg.db))
source("scripts/utils/custom_predict.R")
set.seed(123)

min_cells <- snakemake@config[["min_cells"]]
model_path <- snakemake@config[["model_path"]]

sample_path <- snakemake@input[[1]]
metadata <- readRDS(snakemake@input[[2]])
metadata_2 <- as.data.frame(readRDS(snakemake@input[[3]]))

load(model_path)

counts <- Read10X(sample_path)
if (length(names(counts)) > 1) {
  counts <- counts$`Gene Expression`
}

seu <- CreateSeuratObject(counts,
                          min.cells = min_cells)

seu <- AddMetaData(seu, metadata)
seu <- AddMetaData(seu, metadata_2)

seu <- PercentageFeatureSet(seu,
                                pattern = "^(MT|mt)-", col.name = "percent.mt") %>%
        SCTransform(method = "glmGamPoi",
                    vars.to.regress = "percent.mt", verbose = TRUE) %>%
        RunPCA(verbose = TRUE) %>%
        RunUMAP(dims = 1:30, verbose = TRUE) %>%
        FindNeighbors(dims = 1:30, verbose = TRUE) %>%
        FindClusters(verbose = TRUE)

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
predicted_cds@probabilities <- predicted_cds@probabilities[, -ncol(predicted_cds@probabilities)]
colnames(predicted_cds@probabilities) <- gsub("[.0-9]+$","",colnames(predicted_cds@probabilities))
seu <- AddMetaData(seu, pData(predicted_cds)[c("Size_Factor", "cell_type", "cluster_ext_type", "garnett_cluster")])

cluster_cell_types <- matrix(0, 
                             nrow = length(unique(seu@meta.data$seurat_cluster)), 
                             ncol = ncol(predicted_cds@probabilities))
idx <- 1
for (cluster_name in unique(seu@meta.data$seurat_cluster)) {
    doublet_cells <- rownames(seu@meta.data %>% filter(garnett_cluster == cluster_name))
    selected <- predicted_cds@probabilities[doublet_cells, ]
    cluster_cell_types[idx, ] <- round(apply(selected, 2, median),6)
    idx <- idx + 1
}
cluster_cell_types <- cluster_cell_types / Biobase::rowMax(cluster_cell_types)
rownames(cluster_cell_types) <- unique(seu@meta.data$seurat_cluster)
colnames(cluster_cell_types) <- colnames(predicted_cds@probabilities)
cluster_cell_types <- cluster_cell_types[order(as.numeric(row.names(cluster_cell_types))), ]

clusters_annotation <- apply(cluster_cell_types, 1, function(x) {
    names(x[which.max(x)])
})

seu@meta.data$garnett_prediction <- sapply(as.numeric(seu@meta.data$seurat_cluster),
function(x){
    clusters_annotation[x]
})
seu <- AddMetaData(seu, predicted_cds@probabilities)

garnett_metadata <- seu@meta.data[,c("seurat_clusters", "Size_Factor", "garnett_cluster", "cell_type", 
                                     "cluster_ext_type", "garnett_prediction", colnames(predicted_cds@probabilities))]
saveRDS(garnett_metadata, snakemake@output[[1]])

sink()