library(dplyr)
library(ggplot2)
library(tidyr)
library(RColorBrewer)

# Read the data
df <- read.csv("XXX/decoupleR_go_biological_process_Myeloid.csv")

# Custom order for groups
custom_order <- c('cDC_1_15', 'cDC_2_06', 'DC_LAMP3_11', 'pDC_10', 'M1_S100A8_07', 
                  'M2_CXCL10_12', 'M2_MARCO_05', 'M2_STAB1_09', 'M2_SELENOP_02', 
                  'M2_MMP9_08', 'M2_COL1A1_04', 'Cleaning_macrophage_13', 
                  'Cycling_myeloid_cell_14', 'Neutrophil_1_01', 'Neutrophil_1_03', 
                  'Mast_cell_00')

df$group <- factor(df$group, levels = custom_order)
df <- df[order(df$group), ]

# Select top 5 entries with positive statistic and smallest pvals_adj for each group
df_top <- df %>%
  filter(statistic > 0) %>%
  arrange(group, pvals_adj) %>%
  group_by(group) %>%
  slice_head(n = 4) %>%
  ungroup()

y_names <- df_top$names

# Create a sorted mapping table to maintain order even with duplicate y_names
y_names_df <- data.frame(name = y_names, order = seq_along(y_names))

# Select and sort rows where df$names appear in y_names
filtered_sorted_df <- df %>%
  inner_join(y_names_df, by = c("names" = "name")) %>%
  arrange(order) %>%
  select(-order) # Remove auxiliary sorting column

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
  scale_x_discrete(drop = FALSE) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5))

plot_data

ggsave("XXX/Macrophage_enrichment_top5.pdf", plot_data, height = 17, width = 14, limitsize = FALSE)

ggsave("XXX/Macrophage_enrichment_top4.pdf", plot_data, height = 15, width = 13.6, limitsize = FALSE)