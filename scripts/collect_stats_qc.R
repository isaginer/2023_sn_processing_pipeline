# collect_stats_qc.R

suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(dplyr))
set.seed(123)

mt_thresh <- snakemake@config[["mt_thresh"]]
nuclei_thresh <- snakemake@config[["nuclei_thresh"]]
keep_columns <- snakemake@config[["keep_statement"]]

qc_stats <- matrix(0, nrow = length(snakemake@input), ncol = 13)
colnames(qc_stats) <- c("Total", "Fraction_QC_pass", "Removed_cells",
                        "MT_fraction", "MT_fraction_prc",
                        "Nuclei_fraction", "Nuclei_fraction_prc",
                        "nUMI_outliers", "nUMI_outliers_prc",
                        "nGenes_outliers", "nGenes_outliers_prc",
                        "Phase2_Median_UMI", "Phase2_Median_nGenes")

idx <- 1
SAMPLES <- c()

for (SAMPLE_RDS in snakemake@input) {
    seu <- readRDS(SAMPLE_RDS)
    SAMPLE <- seu@project.name
    SAMPLES <- c(SAMPLES, SAMPLE)

    stat_1 <- table(seu@meta.data$keep)
    stat_2 <- table(seu@meta.data[,paste0("outlier_mt_thresh_",mt_thresh)])
    stat_3 <- table(seu@meta.data[,paste0("outlier_nuclei_thresh_",nuclei_thresh)])
    stat_4 <- seu@meta.data %>% filter(keep) %>% 
        mutate(med_UMI = median(sum), med_nGenes = median(detected)) %>% 
        select(med_UMI,med_nGenes) %>% unique()
    stat_5 <- table(seu@meta.data[, keep_columns[1]])
    stat_6 <- table(seu@meta.data[, keep_columns[2]])
  
    total <- ncol(seu)
    qc_stats[idx, ] <- c(total,
                         stat_1["TRUE"],
                         round((stat_1["FALSE"]/total) * 100,2),
                         stat_2["TRUE"],
                         round((stat_2["TRUE"]/total) * 100,2),
                         stat_3["TRUE"],
                         round((stat_3["TRUE"]/total) * 100,2),
                         stat_5["TRUE"],
                         round((stat_5["TRUE"]/total) * 100,2),
                         stat_6["TRUE"],
                         round((stat_6["TRUE"]/total) * 100,2),
                         stat_4[1,1],
                         stat_4[1,2])
    
    idx <- idx + 1
}
rownames(qc_stats) <- SAMPLES

write.table(qc_stats,
            file = snakemake@output[[1]],
            sep = "\t", quote = FALSE,
            row.names = TRUE, col.names = NA)