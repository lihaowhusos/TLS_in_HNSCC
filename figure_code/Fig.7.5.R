library(pheatmap)

# Read cluster data
folder_path <- "XXX"
file_path <- paste0(folder_path, "GSE136961_RNA_cell2location.csv")
cluster <- read.csv(file_path, header = TRUE, row.names = 1)

# Reorder columns
new_order <- c(
  "B_IGHM_04", "B_CD74_1_00", "B_CD69_1_01", "B_MHC_II_05", "B_GC_11", "B_ISG_08",
  "Plasma_cell_01_02", "Plasma_cell_02_03", "CD8_T_naive_01", "CD8_T_effector_memory_01",
  "CD8_T_cytotoxic_02", "CD8_T_late_exhausted_01", "CD8_T_ISG", "CD4_T_fh_01_00",
  "CD4_T_naive_1_01", "CD4_T_naive_2_03", "CD4_T_early_exhausted_06", "CD4_T_late_exhausted_08",
  "CD4_T_ISG", "cDC_1_15", "cDC_2_06", "DC_LAMP3_11", "pDC_10", "CD4_T_reg_naive_07",
  "CD4_T_reg_1_04", "CD4_T_reg_2_05", "CD4_T_reg_3_10", "CD4_T_reg_ISG", "CD4_T_reg_exhausted_01_02",
  "CD4_T_reg_exhausted_02_09", "Cleaning_macrophage_13", "Cycling_B_plamsa_cell_09", "Cycling_CD4_T",
  "Cycling_CD8_T", "Cycling_NK", "Cycling_myeloid_cell_14", "Endothelial_17", "Fibroblast_13",
  "ILC_06", "Lymphatic_endothelial_cell_21", "M1_S100A8_07", "M2_COL1A1_04", "M2_CXCL10_12",
  "M2_MARCO_05", "M2_MMP9_08", "M2_SELENOP_02", "M2_STAB1_09", "Mast_cell_00", "Neutrophil_1_01",
  "Neutrophil_1_03"
) 

cluster <- cluster[, new_order]

# Normalize columns
normalize_column <- function(x) {
  (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}

plot_df_normalized <- as.data.frame(lapply(cluster, normalize_column))
plot_df_normalized <- t(plot_df_normalized)

colors33 <- colorRampPalette(c("white", "#ff9517", "#b83835", '#6d4490', "#000000"))(100)

# Plot heatmap
plot <- pheatmap(plot_df_normalized, 
                 cluster_rows = FALSE, 
                 cluster_cols = FALSE,
                 border_color = NA, 
                 color = colors33,
                 show_colnames = FALSE
)

pdf(paste0(folder_path, "GEO_cluster.pdf"))  # Create a new pdf device
print(plot)
dev.off()  # Close the pdf device

# Add group information and perform Mann-Whitney tests
merged_df <- cluster
cell_types <- colnames(merged_df)
new_column <- c(rep(0, 12), rep(1, 9))
merged_df$cluster_type <- new_column

man_whitney_results <- lapply(cell_types, function(i) {
  formula <- as.formula(paste(i, "~ cluster_type"))
  wilcox.test(formula, data = merged_df, alternative = "less")
})

names(man_whitney_results) <- cell_types
p_values <- sapply(man_whitney_results, function(x) x$p.value)

results_df <- data.frame(
  Cell_Type = names(p_values),
  P_Value = p_values
)

write.csv(results_df, file = paste0(folder_path, "p_values.csv"), row.names = FALSE)

# Repeat the process for other datasets
datasets <- list(
  list(file = "GSE93157_to_plot.csv", new_column = c(rep(0, 6), rep(1, 12))),
  list(file = "Cancer_Dis_cell2location_plot.csv", new_column = c(rep(0, 12), rep(1, 5)))
)

for (dataset in datasets) {
  file_path <- paste0(folder_path, dataset$file)
  cluster <- read.csv(file_path, header = TRUE, row.names = 1)
  cluster <- cluster[, new_order]
  
  plot_df_normalized <- as.data.frame(lapply(cluster, normalize_column))
  plot_df_normalized <- t(plot_df_normalized)
  
  plot <- pheatmap(plot_df_normalized, 
                   cluster_rows = FALSE, 
                   cluster_cols = FALSE,
                   border_color = NA, 
                   color = colors33,
                   show_colnames = FALSE
  )
  
  pdf(paste0(folder_path, gsub(".csv", ".pdf", dataset$file)))  # Create a new pdf device
  print(plot)
  dev.off()  # Close the pdf device
  
  merged_df <- cluster
  cell_types <- colnames(merged_df)
  merged_df$cluster_type <- dataset$new_column
  
  man_whitney_results <- lapply(cell_types, function(i) {
    formula <- as.formula(paste(i, "~ cluster_type"))
    wilcox.test(formula, data = merged_df, alternative = "less")
  })
  
  names(man_whitney_results) <- cell_types
  p_values <- sapply(man_whitney_results, function(x) x$p.value)
  
  results_df <- data.frame(
    Cell_Type = names(p_values),
    P_Value = p_values
  )
  
  write.csv(results_df, file = paste0(folder_path, "p_values.csv"), row.names = FALSE)
}