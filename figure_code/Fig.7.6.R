library(pheatmap)
library(survival)
library(survminer)
library(readr)
library(dplyr)

# Read cluster data
folder_path = "XXX"
file_path = paste0(folder_path, "GSE93157_to_plot.csv")
cluster = read.csv(file_path, header = TRUE, row.names = 1)

# Re-order columns
new_order <- c(
  "TLS_imprint",
  "B_CD74_1_00", "B_CD69_1_01",
  "B_MHC_II_05", "B_GC_11",
  "B_ISG_08", "Plasma_cell_01_02",
  "Plasma_cell_02_03", "CD8_T_naive_01",
  "CD8_T_effector_memory_01", "CD8_T_cytotoxic_02",
  "CD8_T_late_exhausted_01", "CD8_T_ISG",
  "CD4_T_fh_01_00", "CD4_T_naive_1_01",
  "CD4_T_naive_2_03", "CD4_T_early_exhausted_06",
  "CD4_T_late_exhausted_08", "CD4_T_ISG",
  "cDC_1_15", "cDC_2_06",
  "DC_LAMP3_11", "pDC_10",
  "CD4_T_reg_naive_07", "CD4_T_reg_1_04",
  "CD4_T_reg_2_05", "CD4_T_reg_3_10",
  "CD4_T_reg_ISG", "CD4_T_reg_exhausted_01_02",
  "CD4_T_reg_exhausted_02_09", "Cleaning_macrophage_13",
  "Cycling_B_plamsa_cell_09", "Cycling_CD4_T",
  "Cycling_CD8_T", "Cycling_NK",
  "Cycling_myeloid_cell_14", "Endothelial_17",
  "Fibroblast_13", "ILC_06",
  "Lymphatic_endothelial_cell_21", "M1_S100A8_07",
  "M2_COL1A1_04", "M2_CXCL10_12",
  "M2_MARCO_05", "M2_MMP9_08",
  "M2_SELENOP_02", "M2_STAB1_09",
  "Mast_cell_00", "Neutrophil_1_01",
  "Neutrophil_1_03"
)

cluster = cluster[, new_order]

# Normalize data
normalize_column <- function(x) {
  (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}

plot_df_normalized <- as.data.frame(lapply(cluster, normalize_column))
plot_df_normalized = t(plot_df_normalized)

# Define color palette
colors33 <- colorRampPalette(c("white", "#ff9517", "#b83835", '#6d4490', "#000000"))(100)

# Plot heatmap
plot = pheatmap(
  plot_df_normalized, 
  cluster_rows = FALSE, 
  cluster_cols = FALSE,
  border_color = NA, 
  color = colors33,
  show_colnames = FALSE
)

# Save heatmap to PDF
pdf(paste0(folder_path, "GEO_cluster.pdf"))
print(plot)
dev.off()

# Add group information
merged_df = cluster
new_column <- c(rep(0, 14), rep(1, 4))
merged_df$cluster_type <- as.numeric(new_column)

# Perform Mann-Whitney U test
man_whitney_results <- lapply(colnames(cluster), function(i) {
  formula = as.formula(paste(i, "~ cluster_type"))
  wilcox.test(formula, data = merged_df, alternative = "less")
})

p_values <- sapply(man_whitney_results, function(x) x$p.value)
adjusted_p_values <- p.adjust(p_values, method = "BH")

results_df <- data.frame(
  Cell_Type = colnames(cluster),
  P_Value = p_values,
  Adjusted_P_Value = adjusted_p_values
)

write.csv(results_df, file = paste0(folder_path, "p_values.csv"), row.names = FALSE)

# Combine clinical and sequencing data
setwd('XXX')
rt <- read.csv("survival.csv")
rt$PFSE <- as.numeric(rt$PFSE)
rt$PFS <- as.numeric(rt$PFS)

cell_choose = 'TLS'
a <- rt[, cell_choose] <= median(rt[, cell_choose], na.rm = TRUE)

diff = survdiff(Surv(PFS, PFSE) ~ a, data = rt)
pValue = round(1 - pchisq(diff$chisq, df = 1), 5)
fit <- survfit(Surv(PFS, PFSE) ~ a, data = rt)

p = ggsurvplot(
  fit,
  conf.int = FALSE,
  pval = TRUE,
  surv.median.line = "hv",
  risk.table = TRUE,
  xlab = "Time (month)",
  legend = c(0.8, 0.75),
  legend.title = "",
  legend.labs = c("high", "low"),
  palette = c("#ff0000", "#0000c0")
)

file_path = paste0("survival/", cell_choose, ".pdf")
pdf(file_path)
print(p, newpage = FALSE)
dev.off()

# Loop through cell types and plot survival analysis
cell_types = colnames(rt)[-c(1:19)]

for (cell_choose in cell_types) {
  a <- rt[, cell_choose] <= median(rt[, cell_choose], na.rm = TRUE)
  
  diff = survdiff(Surv(PFS, PFSE) ~ a, data = rt)
  pValue = round(1 - pchisq(diff$chisq, df = 1), 5)
  fit <- survfit(Surv(PFS, PFSE) ~ a, data = rt)
  
  p = ggsurvplot(
    fit,
    conf.int = FALSE,
    pval = TRUE,
    surv.median.line = "hv",
    risk.table = TRUE,
    xlab = "Time (month)",
    legend = c(0.8, 0.75),
    title = cell_choose,
    legend.title = "",
    legend.labs = c("high", "low"),
    palette = c("#ff0000", "#0000c0")
  )
  
  file_path = paste0("survival/", cell_choose, ".pdf")
  pdf(file_path)
  print(p, newpage = FALSE)
  dev.off()
}