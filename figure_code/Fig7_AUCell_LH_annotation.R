# remove data
rm(list=ls())

library(ggplot2)
library(RColorBrewer)

################## CD8_T_AUCell #############################

df <- read.csv("G:/TLSdata/TLS_total_data/ICB_varify/GSE200996_scRNA/AUCell_LH_annotation_to_plot.csv", header = TRUE, stringsAsFactor = FALSE, row.names = 1)


library(dplyr)
# Assuming the data is in a data frame called df
df_scaled <- df %>%
  group_by(names) %>%
  mutate(scaled_meanchange = scale(meanchange))

df_scaled <- df_scaled %>%
  mutate(pvals_adj = ifelse(pvals_adj < 1E-300, 1E-300, pvals_adj)) %>%
  mutate(scaled_meanchange = ifelse(scaled_meanchange > 1, 1, scaled_meanchange)) %>%
  mutate(scaled_meanchange = ifelse(scaled_meanchange < -1, -1, scaled_meanchange))



custom_order <- c("CD8_T_naive_01",
                  "CD8_T_effector_memory_01",
                  "CD8_T_cytotoxic_02",
                  "CD8_T_late_exhausted_01",
                  "CD8_T_ISG"
)

df_scaled$group <- factor(df_scaled$group, levels = custom_order)
df_scaled <- df_scaled[order(df_scaled$group),]




custom_order <- c("CD8.T.naive.centr.mem", 
                  "CD8.T.effector.memory", 
                  "CD8.T.cytotoxic",
                  "CD8.T.dysfunc.early",
                  "CD8.T.dysfunc.late",   
                  "CD8.T.ISG.early",
                  "CD8.T.ISG.late",
                  "CD8.T.dysfunc.ISG"
)

df_scaled$names <- factor(df_scaled$names, levels = custom_order)
df_scaled <- df_scaled[order(df_scaled$names),]




plot_data = ggplot(df_scaled, aes(x = names, y = factor(group, levels = rev(unique(group))), color = scaled_meanchange, size = -log10(pvals_adj))) +
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
  scale_x_discrete(drop=FALSE) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

plot_data

result_path = paste0("G:\\TLSdata\\TLS_total_data\\figures\\annotation\\ovarian_cancer_ref\\",cell_type,"_AUCell.pdf")
ggsave(result_path,plot_data,width = 6.5,height = 4)



################## CD4_T_AUCell #############################
cell_type = "CD4_T"
cell_type_path = paste0("G:/TLSdata/TLS_total_data/AUCell/table/ovarian_cancer_ref/",cell_type,"_AUCell.csv")

df <- read.csv(cell_type_path, header = TRUE, stringsAsFactor = FALSE, row.names = 1)


library(dplyr)
# Assuming the data is in a data frame called df
df_scaled <- df %>%
  group_by(names) %>%
  mutate(scaled_meanchange = scale(meanchange))

df_scaled <- df_scaled %>%
  mutate(pvals_adj = ifelse(pvals_adj < 1E-300, 1E-300, pvals_adj)) %>%
  mutate(scaled_meanchange = ifelse(scaled_meanchange > 1, 1, scaled_meanchange)) %>%
  mutate(scaled_meanchange = ifelse(scaled_meanchange < -1, -1, scaled_meanchange))


c(unique(df_scaled$group))
custom_order <- c("CD4_T_fh_01_00",
                  "CD4_T_naive_1_01",
                  "CD4_T_naive_2_03",
                  
                  "CD4_T_early_exhausted_06",
                  "CD4_T_late_exhausted_08",
                  "CD4_T_ISG",
                  "CD4_T_reg_ISG",
                  
                  
                  "CD4_T_reg_naive_07",
                  "CD4_T_reg_1_04",
                  "CD4_T_reg_2_05",
                  "CD4_T_reg_3_10",
                  "CD4_T_reg_exhausted_01_02",
                  "CD4_T_reg_exhausted_02_09"
)

df_scaled$group <- factor(df_scaled$group, levels = custom_order)
df_scaled <- df_scaled[order(df_scaled$group),]



unique(df_scaled$names)
custom_order <- c("CD4.T.naive.centr.mem.1", 
                  "CD4.T.naive.centr.mem.2", 
                  "CD4.T.effector.memory",
                  "CD4.T.dysfunc.early",
                  "CD4.T.dysfunc.late.1",   
                  "CD4.T.dysfunc.late.2",
                  "CD4.T.ISG",
                  "CD4.T.reg.ISG",
                  "CD4.T.reg.1",   
                  "CD4.T.reg.2",
                  "CD4.T.reg.3",
                  "CD4.Th17.1",
                  "CD4.Th17.2",
                  "CD4.Th17.3"
)

df_scaled$names <- factor(df_scaled$names, levels = custom_order)
df_scaled <- df_scaled[order(df_scaled$names),]




plot_data = ggplot(df_scaled, aes(x = names, y = factor(group, levels = rev(unique(group))), color = scaled_meanchange, size = -log10(pvals_adj))) +
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
  scale_x_discrete(drop=FALSE) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

plot_data


result_path = paste0("G:\\TLSdata\\TLS_total_data\\figures\\annotation\\ovarian_cancer_ref\\",cell_type,"_AUCell.pdf")
ggsave(result_path,plot_data,width = 10,height = 7)



################## Myeloid_AUCell #############################
cell_type = "Myeloid"
cell_type_path = paste0("G:/TLSdata/TLS_total_data/AUCell/table/ovarian_cancer_ref/",cell_type,"_AUCell.csv")

df <- read.csv(cell_type_path, header = TRUE, stringsAsFactor = FALSE, row.names = 1)


library(dplyr)
# Assuming the data is in a data frame called df
df_scaled <- df %>%
  group_by(names) %>%
  mutate(scaled_meanchange = scale(meanchange))

df_scaled <- df_scaled %>%
  mutate(pvals_adj = ifelse(pvals_adj < 1E-300, 1E-300, pvals_adj)) %>%
  mutate(scaled_meanchange = ifelse(scaled_meanchange > 1, 1, scaled_meanchange)) %>%
  mutate(scaled_meanchange = ifelse(scaled_meanchange < -1, -1, scaled_meanchange))


c(unique(df_scaled$group))
custom_order <- c("cDC_1_15", 
                  "cDC_2_06", 
                  "DC_LAMP3_11",
                  "pDC_10",
                  "M1_S100A8_07",   
                  "M2_CXCL10_12",
                  "M2_MARCO_05",
                  "M2_STAB1_09",
                  "M2_SELENOP_02",   
                  "M2_MMP9_08",
                  "M2_COL1A1_04",
                  "Cleaning_macrophage_13",
                  "Cycling_myeloid_cell_14",   
                  "Neutrophil_1_01",
                  "Neutrophil_1_03",
                  "Mast_cell_00"
)

df_scaled$group <- factor(df_scaled$group, levels = custom_order)
df_scaled <- df_scaled[order(df_scaled$group),]




unique(df_scaled$names)
custom_order <- c("cDC1", 
                  "cDC2", 
                  "mDC",
                  "pDC",
                  "M1.S100A8",   
                  "M2.CXCL10",
                  "M2.MARCO",
                  "M2.SELENOP",
                  "M2.MMP9",   
                  "M2.COL1A1",
                  "Clearing.M",
                  "Cycling.M",   
                  "Mast.cell",
                  "dissociated_1",
                  "dissociated",
                  "doublet.Ovarian.cancer.cell"
)

df_scaled$names <- factor(df_scaled$names, levels = custom_order)
df_scaled <- df_scaled[order(df_scaled$names),]


plot_data = ggplot(df_scaled, aes(x = names, y = factor(group, levels = rev(unique(group))), color = scaled_meanchange, size = -log10(pvals_adj))) +
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
  scale_x_discrete(drop=FALSE) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

plot_data

result_path = paste0("G:\\TLSdata\\TLS_total_data\\figures\\annotation\\ovarian_cancer_ref\\",cell_type,"_AUCell.pdf")
ggsave(result_path,plot_data,width = 10,height = 7)












