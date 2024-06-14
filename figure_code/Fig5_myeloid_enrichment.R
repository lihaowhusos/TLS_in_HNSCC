# 删除无关数据
rm(list=ls())


library(dplyr)
library(ggplot2)
library(tidyverse)
library(stringr)
library(tidyr)
library(reshape2)
library(forcats)
library(RColorBrewer)

df = read.csv("G:\\TLSdata\\TLS_total_data\\table\\decoupleR_go_biological_process_Myeloid.csv")
# 假设您的数据已经作为一个名为df的数据框导入了R环境中
# 现在按照group和pvals_adj进行排序，并筛选出每个组pvals_adj最小的前五项




unique(df$group)
custom_order <- c('cDC_1_15',
                  'cDC_2_06',
                  'DC_LAMP3_11',
                  'pDC_10',
                  'M1_S100A8_07',
                  'M2_CXCL10_12',
                  'M2_MARCO_05',
                  'M2_STAB1_09',
                  'M2_SELENOP_02',
                  'M2_MMP9_08',
                  'M2_COL1A1_04',
                  'Cleaning_macrophage_13',
                  'Cycling_myeloid_cell_14',
                  'Neutrophil_1_01',
                  'Neutrophil_1_03',
                  'Mast_cell_00'
)

df$group <- factor(df$group, levels = custom_order)
df <- df[order(df$group),]


# 1. 对每个组筛选statistic为正且pvals_adj最小的前五项
df_top <- df %>%
  filter(statistic > 0) %>%
  arrange(group, pvals_adj) %>%
  group_by(group) %>%
  slice_head(n = 4) %>%
  ungroup()



y_names = df_top$names

# 使用data.frame来创建一个排序的映射表，确保即使有重复的y_name也能保持顺序
y_names_df <- data.frame(name = y_names, order = seq_along(y_names))

# 选择那些df$names出现在y_names中的行，然后依据y_names的排序进行整理
filtered_sorted_df <- df %>%
  inner_join(y_names_df, by = c("names" = "name")) %>%
  arrange(order) %>%
  select(-order) # 将辅助排序列'order'移除





filtered_sorted_df <- filtered_sorted_df %>%
  group_by(names) %>%
  mutate(z_statistic = (statistic - mean(statistic, na.rm = TRUE)) / 
           sd(statistic, na.rm = TRUE))



filtered_sorted_df <- filtered_sorted_df %>%
  mutate(pvals_adj = ifelse(pvals_adj < 1E-250, 1E-250, pvals_adj)) %>%
  mutate(z_statistic = ifelse(z_statistic > 3, 3, z_statistic))






plot_data <- ggplot(filtered_sorted_df, aes(x = group, y = factor(names, levels = unique(names)), color = z_statistic, size = -log10(pvals_adj))) +
  geom_point() +
  scale_color_gradientn(
    colors = rev(brewer.pal(9, "PuOr")),
    limits = c(-3, 3),
    labels = c("≤-3", 0, "≥3"),
    breaks = c(-3, 0, 3)
  ) +
  scale_size_continuous(
    limits = c(0, 250),
    breaks = c(0, 125, 250),
    labels = c(0, 125, 250)
  ) +
  scale_x_discrete(drop=FALSE) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5))

plot_data


ggsave("G:\\TLSdata\\TLS_total_data\\figures\\Macrophage\\Macrophage_enrichment_top5.pdf",plot_data,height = 17,width = 14,limitsize = FALSE)


ggsave("G:\\TLSdata\\TLS_total_data\\figures\\Macrophage\\Macrophage_enrichment_top4.pdf",plot_data,height = 15,width = 13.6,limitsize = FALSE)





