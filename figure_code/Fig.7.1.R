library(ggplot2)
library(RColorBrewer)
library(dplyr)

# Function to read and preprocess data
read_and_preprocess <- function(file_path, custom_group_order, custom_name_order) {
  df <- read.csv(file_path, header = TRUE, stringsAsFactors = FALSE, row.names = 1)
  
  # Scale the data
  df_scaled <- df %>%
    group_by(names) %>%
    mutate(scaled_meanchange = scale(meanchange)) %>%
    mutate(
      pvals_adj = ifelse(pvals_adj < 1E-300, 1E-300, pvals_adj),
      scaled_meanchange = ifelse(scaled_meanchange > 1, 1, scaled_meanchange),
      scaled_meanchange = ifelse(scaled_meanchange < -1, -1, scaled_meanchange)
    )
  
  df_scaled$group <- factor(df_scaled$group, levels = custom_group_order)
  df_scaled <- df_scaled[order(df_scaled$group),]
  
  df_scaled$names <- factor(df_scaled$names, levels = custom_name_order)
  df_scaled <- df_scaled[order(df_scaled$names),]
  
  return(df_scaled)
}

# Function to plot data
plot_data <- function(df_scaled, result_path) {
  plot <- ggplot(df_scaled, aes(x = names, y = factor(group, levels = rev(unique(group))), color = scaled_meanchange, size = -log10(pvals_adj))) +
    geom_point() +
    scale_color_gradientn(
      colors = rev(brewer.pal(9, "RdBu")),
      limits = c(-1, 1),
      labels = c("≤-1", 0, "≥1"),
      breaks = c(-1, 0, 1)
    ) +
    scale_size_continuous(
      limits = c(0, 300),
      breaks = c(0, 150, 300),
      labels = c(0, 150, 300)
    ) +
    scale_x_discrete(drop = FALSE) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  
  ggsave(result_path, plot, width = 10, height = 7)
}

# CD8_T data
cd8_file_path <- "XXX/AUCell_LH_annotation_to_plot.csv"
cd8_group_order <- c("CD8_T_naive_01", "CD8_T_effector_memory_01", "CD8_T_cytotoxic_02", "CD8_T_late_exhausted_01", "CD8_T_ISG")
cd8_name_order <- c("CD8.T.naive.centr.mem", "CD8.T.effector.memory", "CD8.T.cytotoxic", "CD8.T.dysfunc.early", "CD8.T.dysfunc.late", "CD8.T.ISG.early", "CD8.T.ISG.late", "CD8.T.dysfunc.ISG")
cd8_result_path <- "XXX/CD8_T_AUCell.pdf"

df_scaled_cd8 <- read_and_preprocess(cd8_file_path, cd8_group_order, cd8_name_order)
plot_data(df_scaled_cd8, cd8_result_path)

# CD4_T data
cd4_file_path <- "XXX/CD4_T_AUCell.csv"
cd4_group_order <- c("CD4_T_fh_01_00", "CD4_T_naive_1_01", "CD4_T_naive_2_03", "CD4_T_early_exhausted_06", "CD4_T_late_exhausted_08", "CD4_T_ISG", "CD4_T_reg_ISG", "CD4_T_reg_naive_07", "CD4_T_reg_1_04", "CD4_T_reg_2_05", "CD4_T_reg_3_10", "CD4_T_reg_exhausted_01_02", "CD4_T_reg_exhausted_02_09")
cd4_name_order <- c("CD4.T.naive.centr.mem.1", "CD4.T.naive.centr.mem.2", "CD4.T.effector.memory", "CD4.T.dysfunc.early", "CD4.T.dysfunc.late.1", "CD4.T.dysfunc.late.2", "CD4.T.ISG", "CD4.T.reg.ISG", "CD4.T.reg.1", "CD4.T.reg.2", "CD4.T.reg.3", "CD4.Th17.1", "CD4.Th17.2", "CD4.Th17.3")
cd4_result_path <- "XXX/CD4_T_AUCell.pdf"

df_scaled_cd4 <- read_and_preprocess(cd4_file_path, cd4_group_order, cd4_name_order)
plot_data(df_scaled_cd4, cd4_result_path)

# Myeloid data
myeloid_file_path <- "XXX/Myeloid_AUCell.csv"
myeloid_group_order <- c("cDC_1_15", "cDC_2_06", "DC_LAMP3_11", "pDC_10", "M1_S100A8_07", "M2_CXCL10_12", "M2_MARCO_05", "M2_STAB1_09", "M2_SELENOP_02", "M2_MMP9_08", "M2_COL1A1_04", "Cleaning_macrophage_13", "Cycling_myeloid_cell_14", "Neutrophil_1_01", "Neutrophil_1_03", "Mast_cell_00")
myeloid_name_order <- c("cDC1", "cDC2", "mDC", "pDC", "M1.S100A8", "M2.CXCL10", "M2.MARCO", "M2.SELENOP", "M2.MMP9", "M2.COL1A1", "Clearing.M", "Cycling.M", "Mast.cell", "dissociated_1", "dissociated", "doublet.Ovarian.cancer.cell")
myeloid_result_path <- "XXX/Myeloid_AUCell.pdf"

df_scaled_myeloid <- read_and_preprocess(myeloid_file_path, myeloid_group_order, myeloid_name_order)
plot_data(df_scaled_myeloid, myeloid_result_path)