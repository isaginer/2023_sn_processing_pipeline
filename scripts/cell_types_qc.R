suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(HGNChelper))
suppressPackageStartupMessages(library(tidyr))
source("scripts/utils/sc-type-master/R/gene_sets_prepare.R")
source("scripts/utils/sc-type-master/R/sctype_score_.R")
source("scripts/utils/peakfinder.R")
set.seed(123)

min_cells <- snakemake@config[["min_cells"]]
ndims <- snakemake@config[["processing_ndims"]]

seu <- readRDS(snakemake@input[[1]])
metadata <- readRDS(snakemake@input[[2]])

seu <- AddMetaData(seu, metadata)

selected_cells <- rownames(seu@meta.data[seu$keep, ])
seu$pass_QC_1 <- "FAIL"
seu@meta.data[selected_cells, "pass_QC_1"] <- "PASS"
seu <- SCTransform(seu, method = "glmGamPoi",
                    vars.to.regress = "subsets_Mito_percent",
                    verbose = FALSE) %>%
        RunPCA(verbose = FALSE) %>%
        RunUMAP(dims = 1:ndims, verbose = FALSE) %>%
        FindNeighbors(dims = 1:ndims, verbose = FALSE) %>%
        FindClusters(verbose = FALSE)
saveRDS(seu, file = snakemake@output[[1]])

seu_filtered <- subset(seu, pass_QC_1 == "PASS")


DefaultAssay(seu_filtered) <- "RNA"
seu_filtered <- SCTransform(seu_filtered, method = "glmGamPoi",
                    vars.to.regress = "subsets_Mito_percent",
                    verbose = FALSE) %>%
        RunPCA(verbose = FALSE) %>%
        RunUMAP(dims = 1:ndims, verbose = FALSE) %>%
        FindNeighbors(dims = 1:ndims, verbose = FALSE) %>%
        FindClusters(verbose = FALSE)

db_ <- snakemake@config[["markers_file"]]
tissue <- snakemake@config[["tissue"]]
gs_list <- gene_sets_prepare(db_, tissue)

if ("cell_types" %in% names(snakemake@config)) {
  cell_types <- snakemake@config[["cell_types"]]
  gs_list$gs_positive <- gs_list$gs_positive[cell_types]
  gs_list$gs_negative <- gs_list$gs_negative[cell_types]
}

# assign cell types
es.max <- sctype_score(scRNAseqData = seu_filtered[["SCT"]]@scale.data, 
                      scaled = TRUE,
                      gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)
cluster_predictions <- as.data.frame(apply(t(es.max), 1, function(x) {
  rownames(es.max)[which.max(x)]
}))
colnames(cluster_predictions) <- "qc_sctype_prediction"

seu_filtered <- AddMetaData(seu_filtered, cluster_predictions)
DefaultAssay(seu_filtered) <- "RNA"
seu_filtered$nCount_RNA_log10 <- log10(seu_filtered$nCount_RNA)
seu_filtered$nFeature_RNA_log10 <- log10(seu_filtered$nFeature_RNA)

filtered_cells <- c()
for (cell_type in unique(seu_filtered$qc_sctype_prediction)) {
  print(cell_type)
  seu_cell_type <- subset(seu_filtered, qc_sctype_prediction == cell_type)
  peaks <- peakfinder(seu_cell_type$nCount_RNA_log10)
  umi_dist <- seu_cell_type@meta.data[,c("nCount_RNA_log10")]
  # average and merge picks that are too close to each other
  if (length(peaks) > 1) {
    min_range <- median(umi_dist) - sd(umi_dist)
    max_range <- median(umi_dist) + sd(umi_dist)
    in_interval <- between(peaks,min_range,max_range)
    if (any(in_interval)) {
      peaks[which(in_interval)] <- mean(peaks[which(in_interval)])
    }
    #reduce peaks
    peaks <- unique(peaks)
  }
  # filter everything less than first peak
  if (length(peaks) > 1) {
    seu_cell_type <- subset(seu_cell_type, nCount_RNA_log10 <= peaks[1])
    filtered_cells <- c(filtered_cells,colnames(seu_cell_type))
  }
}
seu_filtered$pass_QC_2 <- "PASS"
seu_filtered@meta.data[filtered_cells,"pass_QC_2"] <- "FAIL"
print(table(seu_filtered$pass_QC_2))
metadata <- seu_filtered@meta.data[,c("qc_sctype_prediction",
                                      "nCount_RNA_log10",
                                      "nFeature_RNA_log10",
                                      "pass_QC_2")]

saveRDS(metadata, snakemake@output[[2]])