####################### GSE93157 ############################
library(pheatmap)

###read cluster
folder_path = "G:\\TLSdata\\TLS_total_data\\ICB_varify\\GSE93157\\"
file_path = paste0(folder_path,"GSE93157_to_plot.csv")
cluster = read.csv(file_path,header = TRUE,row.names=1)


## re-order
new_order <- c("TLS_imprint",
  "q05cell_abundance_w_sf_B_IGHM_04", "q05cell_abundance_w_sf_B_CD74_1_00",
  "q05cell_abundance_w_sf_B_CD69_1_01", "q05cell_abundance_w_sf_B_MHC_II_05",
  "q05cell_abundance_w_sf_B_GC_11", "q05cell_abundance_w_sf_B_ISG_08",
  "q05cell_abundance_w_sf_Plasma_cell_01_02", "q05cell_abundance_w_sf_Plasma_cell_02_03",
  "q05cell_abundance_w_sf_CD8_T_naive_01", "q05cell_abundance_w_sf_CD8_T_effector_memory_01",
  "q05cell_abundance_w_sf_CD8_T_cytotoxic_02", "q05cell_abundance_w_sf_CD8_T_late_exhausted_01",
  "q05cell_abundance_w_sf_CD8_T_ISG",
  "q05cell_abundance_w_sf_CD4_T_fh_01_00", "q05cell_abundance_w_sf_CD4_T_naive_1_01",
  "q05cell_abundance_w_sf_CD4_T_naive_2_03",
  "q05cell_abundance_w_sf_CD4_T_early_exhausted_06", "q05cell_abundance_w_sf_CD4_T_late_exhausted_08",
  "q05cell_abundance_w_sf_CD4_T_ISG", "q05cell_abundance_w_sf_cDC_1_15",
  "q05cell_abundance_w_sf_cDC_2_06", "q05cell_abundance_w_sf_DC_LAMP3_11", "q05cell_abundance_w_sf_pDC_10",
  "q05cell_abundance_w_sf_CD4_T_reg_naive_07",
  "q05cell_abundance_w_sf_CD4_T_reg_1_04", "q05cell_abundance_w_sf_CD4_T_reg_2_05",
  "q05cell_abundance_w_sf_CD4_T_reg_3_10", "q05cell_abundance_w_sf_CD4_T_reg_ISG",
  "q05cell_abundance_w_sf_CD4_T_reg_exhausted_01_02", "q05cell_abundance_w_sf_CD4_T_reg_exhausted_02_09",
  "q05cell_abundance_w_sf_Cleaning_macrophage_13",
  "q05cell_abundance_w_sf_Cycling_B_plamsa_cell_09", "q05cell_abundance_w_sf_Cycling_CD4_T",
  "q05cell_abundance_w_sf_Cycling_CD8_T", "q05cell_abundance_w_sf_Cycling_NK",
  "q05cell_abundance_w_sf_Cycling_myeloid_cell_14", 
  "q05cell_abundance_w_sf_Endothelial_17", "q05cell_abundance_w_sf_Fibroblast_13",
  "q05cell_abundance_w_sf_ILC_06", "q05cell_abundance_w_sf_Lymphatic_endothelial_cell_21",
  "q05cell_abundance_w_sf_M1_S100A8_07", "q05cell_abundance_w_sf_M2_COL1A1_04",
  "q05cell_abundance_w_sf_M2_CXCL10_12", "q05cell_abundance_w_sf_M2_MARCO_05",
  "q05cell_abundance_w_sf_M2_MMP9_08", "q05cell_abundance_w_sf_M2_SELENOP_02",
  "q05cell_abundance_w_sf_M2_STAB1_09", "q05cell_abundance_w_sf_Mast_cell_00",
  "q05cell_abundance_w_sf_Neutrophil_1_01",
  "q05cell_abundance_w_sf_Neutrophil_1_03"
) 

cluster = cluster[, new_order]





plot_df = as.data.frame(cluster)

normalize_column <- function(x) {
  (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}

plot_df_normalized <- plot_df  # Start with a copy of the original dataframe

for(i in colnames(plot_df)) {
  plot_df_normalized[[i]] <- normalize_column(plot_df[[i]])
}

max(plot_df_normalized)
plot_df_normalized = t(plot_df_normalized)


colors33 <- colorRampPalette(c("white","#ff9517","#b83835",'#6d4490',"#000000"))(100) ###color

#breaks <- seq(0, 5, length.out = length(colors33) + 1)


plot = pheatmap(plot_df_normalized, 
                #scale = "row",
                cluster_rows = FALSE, 
                cluster_cols = FALSE,
                border_color=NA, 
                color=colors33,
                #breaks = breaks,
                show_colnames = FALSE,
                #display_numbers = T,
                #cellheight= 10,
                #gaps_col = c(334,105,96,11),
                #annotation_col = annotation_col,
                #cutree_rows = 4
)


pdf(paste0(folder_path,"GEO_cluster.pdf")) # Create a new pdf device
print(plot)
dev.off() # Close the pdf device


##################################### statistic analysis ##########################
merged_df = cluster
cell_types = colnames(merged_df)

######## add group infomation
new_column <- c(rep(0, 14), rep(1, 4))
merged_df$cluster_type <- new_column



# If you don't want to keep the RowNames column, you can remove it after merging
merged_df$RowNames <- NULL

merged_df$cluster_type = as.numeric(merged_df$cluster_type)

man_whitney_results <- list()
for(i in cell_types) {
  # Create the formula with the current column name
  formula = as.formula(paste(i, "~ cluster_type"))
  
  # Run the Mann-Whitney U test with the constructed formula
  # Note: in wilcox.test(), by default, alternative="two.sided", we set it to "less" for a one-tailed test
  man_whitney_results[[i]] <- wilcox.test(formula, data = merged_df, alternative="less")
  
  # Print the results with the current cell type name
  print(paste("Results for", i))
  print(man_whitney_results[[i]])
}

# Assuming that you have already run your tests and the results are in man_whitney_results,
# Extracting p-values from each man_whitney_result
p_values <- sapply(man_whitney_results, \(x) x$p.value)

# Adjusting the p-values for multiple testing
adjusted_p_values <- p.adjust(p_values, method = "BH")  # Or your method of choice

# If you need to keep track of which p-value corresponds to which test:
#names(p_values) <- cell_types

# Check the adjusted p-values
#print(adjusted_p_values)
#adjusted_p_values


results_df <- data.frame(
  Cell_Type = names(adjusted_p_values),
  P_Value = p_values,
  Adjuest_P_Value = adjusted_p_values
)

# Save the dataframe to a CSV file
write.csv(results_df, file = paste0(folder_path,"p_values.csv"), row.names = FALSE)



########### survival analysis ###########
####combine clinical and sequencing data
library("survival")
library("survminer")
library("readr")
library(dplyr)
library(stringr)

setwd('G:/TLSdata/TLS_total_data/ICB_varify/GSE93157/')

rt <- read.csv("survival.csv")

rt$PFSE <- as.numeric(rt$PFSE)
rt$PFS <- as.numeric(rt$PFS)

#write.csv(rt, 'cell2location_survival_data.csv', row.names = FALSE)

print(colnames(rt))

cell_choose = 'q05cell_abundance_w_sf_CD4_T_early_exhausted_06'
cell_choose = 'TLS'

a <- rt[,cell_choose] <= median(rt[,cell_choose], na.rm = TRUE)

diff=survdiff(Surv(PFS, PFSE) ~a,data = rt)
pValue=1-pchisq(diff$chisq,df=1)
pValue=round(pValue,5)
fit <- survfit(Surv(PFS, PFSE) ~ a, data = rt)




#summary(fit)    #查看五年生存率

p = ggsurvplot(fit, # 创建的拟合对象
               conf.int = FALSE, # 显示置信区间
               pval = TRUE, # 添加P值
               surv.median.line = "hv",  # 添加中位生存时间线
               risk.table = TRUE, # 添加风险表
               xlab = "Time (month)", # 指定x轴标签
               legend = c(0.8,0.75), # 指定图例位置
               legend.title = "", # 设置图例标题
               legend.labs = c("high", "low"), # 指定图例分组标签
               palette = c("#ff0000", "#0000c0")
               ) # 设置x轴刻度间距

file_path = paste0("survival/",cell_choose,".pdf")
pdf(file_path)
print(p, newpage = FALSE)
dev.off()


############### 循环输出 OS ##############
cell_types = colnames(rt)
cell_types <- cell_types[-c(1:19)]
cell_types

for (cell_choose in cell_types){
  cell_choose
  a <- rt[,cell_choose] <= median(rt[,cell_choose], na.rm = TRUE)
  
  diff=survdiff(Surv(PFS, PFSE) ~a,data = rt)
  pValue=1-pchisq(diff$chisq,df=1)
  pValue=round(pValue,5)
  fit <- survfit(Surv(PFS, PFSE) ~ a, data = rt)
  p = ggsurvplot(fit, # 创建的拟合对象
                 conf.int = FALSE, # 显示置信区间
                 pval = TRUE, # 添加P值
                 surv.median.line = "hv",  # 添加中位生存时间线
                 risk.table = TRUE, # 添加风险表
                 xlab = "Time (month)", # 指定x轴标签
                 legend = c(0.8,0.75), # 指定图例位置
                 title = cell_choose,
                 #tables.height = ,# numeric value (in [0 - 1])
                 legend.title = "", # 设置图例标题
                 legend.labs = c("high", "low"), # 指定图例分组标签
                 palette = c("#ff0000", "#0000c0")
  ) # 设置x轴刻度间距
  
  file_path = paste0("survival/",cell_choose,".pdf")
  pdf(file_path)
  print(p, newpage = FALSE)
  dev.off()
}


