library(dplyr)
library(ggplot2)
library(tidyr)
library(stringr)
library(RColorBrewer)

# Read the CSV file
df <- read.csv("XXX/decoupleR_go_biological_process_T_NK.csv")

# Filter and sort data by group
df <- df %>% filter(str_starts(group, "CD8"))

# Custom order for the groups
custom_order <- c(
  "CD8_T_naive_01",
  "CD8_T_effector_memory_01",
  "CD8_T_cytotoxic_02",
  "CD8_T_late_exhausted_01",
  "CD8_T_ISG"
)

df$group <- factor(df$group, levels = custom_order)
df <- df[order(df$group),]

# Select top 5 items with positive statistics and lowest pvals_adj for each group
df_top <- df %>%
  filter(statistic > 0) %>%
  arrange(group, pvals_adj) %>%
  group_by(group) %>%
  slice_head(n = 5) %>%
  ungroup()

y_names <- df_top$names

# Create a sorted mapping table to maintain order even with duplicate y_names
y_names_df <- data.frame(name = y_names, order = seq_along(y_names))

# Filter and sort the original dataframe based on the sorted y_names
filtered_sorted_df <- df %>%
  inner_join(y_names_df, by = c("names" = "name")) %>%
  arrange(order) %>%
  select(-order)  # Remove the auxiliary sorting column 'order'

# Z-score normalization of the statistics
filtered_sorted_df <- filtered_sorted_df %>%
  group_by(names) %>%
  mutate(z_statistic = (statistic - mean(statistic, na.rm = TRUE)) / sd(statistic, na.rm = TRUE))

# Adjust pvals_adj and z_statistic values
filtered_sorted_df <- filtered_sorted_df %>%
  mutate(pvals_adj = ifelse(pvals_adj < 1E-250, 1E-250, pvals_adj)) %>%
  mutate(z_statistic = ifelse(z_statistic > 2, 2, z_statistic))

# Plotting the data
plot_data <- ggplot(filtered_sorted_df, aes(x = group, y = factor(names, levels = rev(unique(names))), color = z_statistic, size = -log10(pvals_adj))) +
  geom_point() +
  scale_color_gradientn(
    colors = rev(brewer.pal(9, "RdBu")),
    limits = c(-2, 2),
    labels = c("≤-2", 0, "≥2"),
    breaks = c(-2, 0, 2)
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

# Save the plot
ggsave("XXX/CD8_enrichment_top5_new.pdf", plot_data, height = 9, width = 8.5, limitsize = FALSE)