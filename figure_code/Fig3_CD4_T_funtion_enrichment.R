library(dplyr)
library(ggplot2)
library(tidyverse)
library(stringr)
library(tidyr)
library(reshape2)
library(forcats)
library(RColorBrewer)

df = read.csv("G:\\TLSdata\\TLS_total_data\\table\\decoupleR_go_biological_process_T_NK.csv")
# 假设您的数据已经作为一个名为df的数据框导入了R环境中
# 现在按照group和pvals_adj进行排序，并筛选出每个组pvals_adj最小的前五项


df <- df %>%  filter(str_starts(group, "CD4"))
df <- df[grepl("^CD4", df$group), ]



# 1. 对每个组筛选statistic为正且pvals_adj最小的前五项
df_top <- df %>%
  filter(statistic > 0) %>%
  arrange(group, pvals_adj) %>%
  group_by(group) %>%
  slice_head(n = 10) %>%
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
  mutate(group_names = paste0(row_number(), "_", group, "_", names))


filtered_sorted_df <- filtered_sorted_df %>%
  mutate(pvals_adj = ifelse(pvals_adj < 1E-250, 1E-250, pvals_adj)) %>%
  mutate(z_statistic = ifelse(z_statistic > 3, 3, z_statistic))


filtered_sorted_df$group_names <- factor(filtered_sorted_df$group_names, levels = unique(filtered_sorted_df$group_names))

p <- ggplot(filtered_sorted_df, aes(x = group, y = group_names, color = z_statistic, size = -log10(pvals_adj))) +
  geom_point() +
  scale_size_continuous(
    limits = c(0, 250),
    breaks = c(0, 125, 250),
    labels = c(0, 125, 250)
  ) +
  scale_color_gradientn(
    colors = rev(brewer.pal(9, "RdBu")),
    limits = c(-3, 3),
    labels = c(-3, 0, 3),
    breaks = c(-3, 0, 3)
  ) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_discrete(limits = levels(filtered_sorted_df$group_names))

print(p)


ggsave("G:\\TLSdata\\TLS_total_data\\figures\\TCR_CD8\\CD4_enrichment_top3.pdf",p,height = 17,width = 12,limitsize = FALSE)


y_names





