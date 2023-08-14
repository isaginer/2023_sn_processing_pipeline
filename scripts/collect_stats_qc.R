# collect_stats_qc.R

suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(dplyr))
set.seed(123)

mt_thresh <- snakemake@config[["mt_thresh"]]
nuclei_thresh <- snakemake@config[["nuclei_thresh"]]
ribo_thresh <- snakemake@config[["ribo_thresh"]]
keep_columns <- snakemake@config[["keep_statement"]]

qc_stats <- matrix(0, nrow = length(snakemake@input[["step_1"]]), ncol = 21)
colnames(qc_stats) <- c("Total", "Fraction_QC_pass_step1", "Removed_cells_step1",
                        "Fraction_QC_pass_step2", "Removed_cells_step2",
                        "MT_fraction_step1", "MT_fraction_prc_step1",
                        "Nuclei_fraction_step1", "Nuclei_fraction_prc_step1",
                        "Ribo_fraction_step1", "Ribo_fraction_prc_step1",
                        "nUMI_outliers_step1", "nUMI_outliers_prc_step1",
                        "nGenes_outliers_step1", "nGenes_outliers_prc_step1",
                        "Phase2_Median_UMI_step1", "Phase2_Median_nGenes_step1",
                        "nUMI_outliers_step2", "nUMI_outliers_prc_step2",
                        "Phase2_Median_UMI_step2", "Phase2_Median_nGenes_step2")

idx <- 1
SAMPLES <- c()

for (SAMPLE_RDS in snakemake@input[["step_1"]]) {
    seu_1 <- readRDS(SAMPLE_RDS)
    SAMPLE <- seu_1@project.name
    seu_2 <- readRDS(snakemake@input[["step_2"]][grepl(paste0(SAMPLE,"[.]rds$"), snakemake@input[["step_2"]])])
    
    SAMPLES <- c(SAMPLES, SAMPLE)

    stat_1 <- table(seu_1@meta.data$keep)
    stat_2 <- table(seu_1@meta.data[,paste0("outlier_mt_thresh_", mt_thresh)])
    stat_3 <- table(seu_1@meta.data[,paste0("outlier_nuclei_thresh_", nuclei_thresh)])
    stat_4 <- table(seu_1@meta.data[,paste0("outlier_ribo_thresh_", ribo_thresh)])
    stat_5 <- seu_1@meta.data %>% filter(keep) %>% 
        mutate(med_UMI = median(sum), med_nGenes = median(detected)) %>% 
        select(med_UMI,med_nGenes) %>% unique()
    stat_6 <- table(seu_1@meta.data[, keep_columns[1]])
    stat_7 <- table(seu_1@meta.data[, keep_columns[2]])
  
    total <- ncol(seu_1)

    stat_8 <- table(seu_2@meta.data$pass_QC_2)
    stat_9 <- seu_2@meta.data %>% filter(keep) %>% 
        mutate(med_UMI = median(sum), med_nGenes = median(detected)) %>% 
        select(med_UMI,med_nGenes) %>% unique()

    qc_stats[idx, ] <- c(total,
                         stat_1["TRUE"],
                         round((stat_1["FALSE"]/total) * 100,2),
                         stat_8["PASS"],
                         round((stat_8["FAIL"]/total) * 100,2),
                         stat_2["TRUE"],
                         round((stat_2["TRUE"]/total) * 100,2),
                         stat_3["TRUE"],
                         round((stat_3["TRUE"]/total) * 100,2),
                         stat_4["TRUE"],
                         round((stat_4["TRUE"]/total) * 100,2),
                         stat_6["TRUE"],
                         round((stat_6["TRUE"]/total) * 100,2),
                         stat_7["TRUE"],
                         round((stat_7["TRUE"]/total) * 100,2),
                         stat_5[1,1],
                         stat_5[1,2],
                         stat_8["FAIL"],
                         round((stat_8["FAIL"]/total) * 100,2),
                         stat_9[1,1],
                         stat_9[1,2])
    
    idx <- idx + 1
}
rownames(qc_stats) <- SAMPLES
qc_stats[is.na(qc_stats)] <- 0

write.table(qc_stats,
            file = snakemake@output[[1]],
            sep = "\t", quote = FALSE,
            row.names = TRUE, col.names = NA)