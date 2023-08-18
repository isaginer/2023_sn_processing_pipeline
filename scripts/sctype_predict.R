suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(HGNChelper))
suppressPackageStartupMessages(library(tidyr))
source("scripts/utils/sc-type-master/R/gene_sets_prepare.R")
source("scripts/utils/sc-type-master/R/sctype_score_.R")
set.seed(123)

min_cells <- snakemake@config[["min_cells"]]
if ("sctype_resolution" %in% names(snakemake@config)) {
  res <- snakemake@config[["sctype_resolution"]]
} else {
  res <- 20
}
if ("stable_clusters_prc" %in% names(snakemake@config)) {
  stable_cluster_prc <- snakemake@config[["stable_clusters_prc"]]
} else {
  stable_cluster_prc <- 0.9
}
ndims <- snakemake@config[["processing_ndims"]]

counts <- Read10X(file.path(dirname(snakemake@input[[1]]),
                            snakemake@params[["mtx_location"]]))

if (length(names(counts)) > 1) {
  counts <- counts$`Gene Expression`
}
seu <- CreateSeuratObject(counts,
                          min.cells = min_cells)
seu <- PercentageFeatureSet(seu, pattern = "^(MT|mt)-",
                            col.name = "percent.mt") %>%
        SCTransform(method = "glmGamPoi",
                    vars.to.regress = "percent.mt", verbose = TRUE) %>%
        RunPCA(verbose = FALSE) %>%
        RunUMAP(dims = 1:ndims, verbose = FALSE) %>%
        FindNeighbors(dims = 1:ndims, verbose = FALSE)

db_ <- snakemake@config[["markers_file"]]
tissue <- snakemake@config[["tissue"]]
gs_list <- gene_sets_prepare(db_, tissue)

if ("cell_types" %in% names(snakemake@config)) {
  cell_types <- snakemake@config[["cell_types"]]
  gs_list$gs_positive <- gs_list$gs_positive[cell_types]
  gs_list$gs_negative <- gs_list$gs_negative[cell_types]
}

# assign cell types
es.max <- sctype_score(scRNAseqData = seu[["SCT"]]@scale.data, scaled = TRUE,
                      gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)
cluster_predictions <- as.data.frame(apply(t(es.max), 1, function(x){
  rownames(es.max)[which.max(x)]
}))
colnames(cluster_predictions) <- "sctype_prediction"

seu <- FindClusters(seu, resolution = res)
cell_types <- gsub(" ", ".",names(gs_list$gs_positive))
seu <- AddMetaData(seu, as.data.frame(t(es.max)))
seu <- AddMetaData(seu, cluster_predictions)

agg_data <- as.data.frame(seu@meta.data %>% group_by(across(paste0("SCT_snn_res.",res))) %>% 
  mutate(across(all_of(cell_types), ~ sum(.x)), .keep = "used") %>% distinct())
agg_data[,paste0("SCT_snn_res.",res)] <- as.character(agg_data[,paste0("SCT_snn_res.",res)])

cluster_stats <- seu@meta.data %>% group_by_at(paste0("SCT_snn_res.",res)) %>% mutate(total_n = n()) %>%
  group_by_at(c(paste0("SCT_snn_res.",res), "sctype_prediction")) %>% reframe(cnt = n(), prc = cnt/total_n) %>% unique()
top_1_cluster_stats <- as.data.frame(cluster_stats %>% group_by_at(paste0("SCT_snn_res.",res)) %>% top_n(1,prc))

toPlot <- merge(agg_data,top_1_cluster_stats, by = c(paste0("SCT_snn_res.",res))) %>% 
  arrange(prc)
#toPlot[toPlot$value < 0, "value"] <- 0

stableClusters <- toPlot %>% filter(prc > stable_cluster_prc) %>%
  group_by(sctype_prediction) %>%
  mutate(across(all_of(cell_types), ~ min(.x)), .keep = "used") %>%
  distinct() %>% pivot_longer(cell_types)
stableClusters$sctype_prediction <- gsub(" ",".",stableClusters$sctype_prediction)

thresholds <- as.data.frame(stableClusters[stableClusters$sctype_prediction==stableClusters$name,])
non_stable_cluster <- thresholds[is.na(thresholds$value),"name"]
if (length(non_stable_cluster) > 0) {
  warning(paste0("No stable clusters for: ",
  paste0(non_stable_cluster,collapse=",")))
}
thresholds <- thresholds[!is.na(thresholds$value), ]

seu$low_quality <- FALSE
for (prediction_ct in thresholds$sctype_prediction) {
  toPlot_sub <- toPlot %>% filter(prc <= stable_cluster_prc) %>%
    filter(sctype_prediction == prediction_ct)
  ct_minimum <- as.numeric(thresholds[thresholds$sctype_prediction == prediction_ct,"value"])
  low_quality_clusters <- as.numeric(toPlot_sub[toPlot_sub[,prediction_ct] < ct_minimum, paste0("SCT_snn_res.",res)])
  seu@meta.data[seu@meta.data[,paste0("SCT_snn_res.",res)] %in% low_quality_clusters, "low_quality"] <- T
}

sctype_predictions <- cbind(as.data.frame(t(es.max)),
seu@meta.data[,c("sctype_prediction", "low_quality")])
colnames(sctype_predictions) <- c(rownames(es.max), "sctype_prediction", "Unknown_celltype")
sctype_predictions <- sctype_predictions[rownames(seu@meta.data), ]

saveRDS(sctype_predictions, snakemake@output[[1]])