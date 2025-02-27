output_log <- file(snakemake@log[[2]], open="wt")
error_log <- file(snakemake@log[[1]], open="wt")
sink(output_log, type = "output")
sink(error_log, type = "message")

suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(scater))
set.seed(123)

sum_log_nmad <- snakemake@config[["sum_log_nmad"]]
detected_log_nmad <- snakemake@config[["detected_log_nmad"]]
mt_thresh <- snakemake@config[["mt_thresh"]]
ribo_thresh <- snakemake@config[["ribo_thresh"]]
nuclei_thresh <- snakemake@config[["nuclei_thresh"]]
keep_columns <- snakemake@config[["keep_statement"]]

seu_sub <- readRDS(snakemake@input[[1]])
seu_sub_sce <- as.SingleCellExperiment(seu_sub, assay = "RNA")
is.mito <- grep("^(MT-|mt-)", rownames(seu_sub_sce))
is.nuclei <- grep(paste0("^",paste0(c("MALAT1", "NEAT1", "FTX",
                  "FOXP1", "RBMS3", "ZBTB20", "LRMDA", "PBX1", "ITPR2", "AUTS2", "TTC28", "BNC2", "EXOC4", "RORA",
                  "PRKG1", "ARID1B", "PARD3B", "GPHN", "N4BP2L2", "PKHD1L1", "EXOC6B", "FBXL7", "MED13L",
                  "TBC1D5", "IMMP2L", "SYNE1", "RERE", "MBD5", "EXT1", "WWOX",
                  "Malat1","Neat1","Ftx","Foxp1","Rbms3","Zbtb20","Lrmda",
                  "Pbx1","Itpr2","Auts2","Ttc28","Bnc2","Exoc4","Rora",
                  "Prkg1", "Arid1b", "Pard3b", "Gphn", "N4bp2l2","Pkhd1l1",
                  "Exoc6b", "Fbxl7", "Med13l", "Tbc1d5", "Immp2l", "Syne1",
                  "Rere", "Mbd5", "Ext1", "Wwox"),collapse="|"),"$"), rownames(seu_sub_sce))
is.ribo <- grep("^(RPL|RPS|Rpl|Rps)", rownames(seu_sub_sce))
per_cell.stats <- as.data.frame(perCellQCMetrics(seu_sub_sce,
                                                 percent.top = c(20, 50),
                                                 subsets = list(Mito = is.mito,
                                                                Nuclei = is.nuclei,
                                                                Ribo = is.ribo)))

for (nmad in sum_log_nmad) {
    for (qc_type in c("lower", "both")) {
        column_name <- paste0("outlier_sum_log_", nmad, "_", qc_type)
        low_thresh_name <- paste0("low_thresh_sum_log_", nmad, "_", qc_type)
        high_thresh_name <- paste0("high_thresh_sum_log_", nmad, "_", qc_type)
        outliers_ <- isOutlier(per_cell.stats$sum,
                                                type = qc_type,
                                                nmads = nmad,
                                                log = TRUE)
        per_cell.stats[, low_thresh_name] <- attr(outliers_, "threshold")[1]
        per_cell.stats[, high_thresh_name] <- attr(outliers_, "threshold")[2]
        per_cell.stats[, column_name] <- outliers_
    }
}

for (nmad in detected_log_nmad) {
    for (qc_type in c("lower", "both")) {
        column_name <- paste0("outlier_detected_log_", nmad, "_", qc_type)
        low_thresh_name <- paste0("low_thresh_detected_log_", nmad, "_", qc_type)
        high_thresh_name <- paste0("high_thresh_detected_log_", nmad, "_", qc_type)
        outliers_ <- isOutlier(per_cell.stats$detected,
                                                type = qc_type,
                                                nmads = nmad,
                                                log = TRUE)
        per_cell.stats[, low_thresh_name] <- attr(outliers_, "threshold")[1]
        per_cell.stats[, high_thresh_name] <- attr(outliers_, "threshold")[2]
        per_cell.stats[, column_name] <- outliers_
    }
}

per_cell.stats[,paste0("outlier_mt_thresh_", mt_thresh)] <- per_cell.stats$subsets_Mito_percent > mt_thresh
per_cell.stats[,paste0("outlier_nuclei_thresh_", nuclei_thresh)] <- per_cell.stats$subsets_Nuclei_percent < nuclei_thresh
per_cell.stats[,paste0("outlier_ribo_thresh_", ribo_thresh)] <- per_cell.stats$subsets_Ribo_percent > ribo_thresh

per_cell.stats$keep <- !apply(per_cell.stats[, keep_columns],
                                            1, any)

#seu_sub <- AddMetaData(seu_sub, per_cell.stats)
saveRDS(per_cell.stats, snakemake@output[[1]])

sink()
