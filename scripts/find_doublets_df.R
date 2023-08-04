
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(DoubletFinder))
set.seed(123)

rename_columns <- function(pattern, columns_, new_name) {
    columns_[grepl(pattern, columns_)] <- new_name
    return(columns_)
}

doublets_data <- readRDS(snakemake@input[[3]])

stats_ <- as.list(table(doublets_data$scDblFinder.class))
doublets_prc <- round(stats_$doublet / nrow(doublets_data), 4)

min_cells <- snakemake@config[["min_cells"]]
dims_ <- snakemake@config[["doublet_dims"]]

counts <- Read10X(file.path(dirname(snakemake@input[[1]]),
                            snakemake@params[["mtx_location"]]))

if (length(names(counts)) > 1) {
  counts <- counts$`Gene Expression`
}

seu <- CreateSeuratObject(counts,
                          min.cells = min_cells) %>%
    SCTransform() %>%
    RunPCA() %>%
    RunUMAP(dims = 1:dims_)

# Load optimal parameters
bcmvn <- readRDS(snakemake@input[[2]])

best_pk <- as.numeric(as.character(bcmvn[which.max(bcmvn$BCmetric), "pK"]))

nexp_poi <- round(doublets_prc * nrow(seu@meta.data))

seu <- doubletFinder_v3(seu, PCs = 1:dims_, pN = 0.25,
                        pK = best_pk, nExp = nexp_poi,
                        reuse.pANN = FALSE, sct = TRUE)
columns_ <- colnames(seu@meta.data)
new_column <- paste0("doublet_prob_", doublets_prc)
colnames(seu@meta.data) <- rename_columns("^pANN_",
                                            columns_,
                                            new_column)
columns_ <- colnames(seu@meta.data)
new_column <- paste0("doublet_", doublets_prc)
colnames(seu@meta.data) <- rename_columns("^DF\\.",
                                            columns_,
                                            new_column)

to_save <- seu@meta.data[, grepl("^doublet_",
                              colnames(seu@meta.data))]
colnames(to_save) <- gsub("[0-9.]+$", "prc", colnames(to_save))

saveRDS(to_save,
        file = snakemake@output[[1]])