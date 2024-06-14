# Load necessary libraries
library(ggplot2)
library(RColorBrewer)
library(dplyr)

# Function to process and plot data
process_and_plot <- function(cell_type, custom_order_group, custom_order_names, result_width, result_height) {
  cell_type_path <- paste0("XXX/", cell_type, "_AUCell.csv")
  df <- read.csv(cell_type_path, header = TRUE, stringsAsFactors = FALSE, row.names = 1)
  
  # Scale data and adjust values
  df_scaled <- df %>%
    group_by(names) %>%
    mutate(scaled_meanchange = scale(meanchange)) %>%
    mutate(pvals_adj = ifelse(pvals_adj < 1E-300, 1E-300, pvals_adj)) %>%
    mutate(scaled_meanchange = ifelse(scaled_meanchange > 1, 1, scaled_meanchange)) %>%
    mutate(scaled_meanchange = ifelse(scaled_meanchange < -1, -1, scaled_meanchange))
  
  # Custom order for groups and names
  df_scaled$group <- factor(df_scaled$group, levels = custom_order_group)
  df_scaled <- df_scaled[order(df_scaled$group), ]
  
  df_scaled$names <- factor(df_scaled$names, levels = custom_order_names)
  df_scaled <- df_scaled[order(df_scaled$names), ]
  
  # Plot data
  plot_data <- ggplot(df_scaled, aes(x = names, y = factor(group, levels = rev(unique(group))),
                                     color = scaled_meanchange, size = -log10(pvals_adj))) +
    geom_point() +
    scale_color_gradientn(colors = rev(brewer.pal(9, "RdBu")), limits = c(-1, 1),
                          labels = c("≤-1", 0, "≥1"), breaks = c(-1, 0, 1)) +
    scale_size_continuous(limits = c(0, 300), breaks = c(0, 150, 300), labels = c(0, 150, 300)) +
    scale_x_discrete(drop = FALSE) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  
  plot_data
  
  result_path <- paste0("XXX/", cell_type, "_AUCell.pdf")
  ggsave(result_path, plot_data, width = result_width, height = result_height)
}

# CD8_T cell type
cd8_custom_order_group <- c("CD8_T_naive_01", "CD8_T_effector_memory_01", "CD8_T_cytotoxic_02",
                            "CD8_T_late_exhausted_01", "CD8_T_ISG")
cd8_custom_order_names <- c("CD8.T.naive.centr.mem", "CD8.T.effector.memory", "CD8.T.cytotoxic",
                            "CD8.T.dysfunc.early", "CD8.T.dysfunc.late", "CD8.T.ISG.early",
                            "CD8.T.ISG.late", "CD8.T.dysfunc.ISG")
process_and_plot("CD8_T", cd8_custom_order_group, cd8_custom_order_names, 6.5, 4)

# CD4_T cell type
cd4_custom_order_group <- c("CD4_T_fh_01_00", "CD4_T_naive_1_01", "CD4_T_naive_2_03",
                            "CD4_T_early_exhausted_06", "CD4_T_late_exhausted_08", "CD4_T_ISG",
                            "CD4_T_reg_ISG", "CD4_T_reg_naive_07", "CD4_T_reg_1_04", "CD4_T_reg_2_05",
                            "CD4_T_reg_3_10", "CD4_T_reg_exhausted_01_02", "CD4_T_reg_exhausted_02_09")
cd4_custom_order_names <- c("CD4.T.naive.centr.mem.1", "CD4.T.naive.centr.mem.2", "CD4.T.effector.memory",
                            "CD4.T.dysfunc.early", "CD4.T.dysfunc.late.1", "CD4.T.dysfunc.late.2",
                            "CD4.T.ISG", "CD4.T.reg.ISG", "CD4.T.reg.1", "CD4.T.reg.2", "CD4.T.reg.3",
                            "CD4.Th17.1", "CD4.Th17.2", "CD4.Th17.3")
process_and_plot("CD4_T", cd4_custom_order_group, cd4_custom_order_names, 10, 7)

# Myeloid cell type
myeloid_custom_order_group <- c("cDC_1_15", "cDC_2_06", "DC_LAMP3_11", "pDC_10", "M1_S100A8_07",
                                "M2_CXCL10_12", "M2_MARCO_05", "M2_STAB1_09", "M2_SELENOP_02",
                                "M2_MMP9_08", "M2_COL1A1_04", "Cleaning_macrophage_13", "Cycling_myeloid_cell_14",
                                "Neutrophil_1_01", "Neutrophil_1_03", "Mast_cell_00")
myeloid_custom_order_names <- c("cDC1", "cDC2", "mDC", "pDC", "M1.S100A8", "M2.CXCL10", "M2.MARCO",
                                "M2.SELENOP", "M2.MMP9", "M2.COL1A1", "Clearing.M", "Cycling.M",
                                "Mast.cell", "dissociated_1", "dissociated", "doublet.Ovarian.cancer.cell")
process_and_plot("Myeloid", myeloid_custom_order_group, myeloid_custom_order_names, 10, 7)