library(dplyr)
library(ggplot2)
library(stringr)
library(tidyr)
library(forcats)
library(RColorBrewer)

# Read the data
df <- read.csv("XXX")

# Filter and sort data
df <- df %>%
  filter(str_starts(group, "CD4")) %>%
  filter(grepl("^CD4", group))

df$group <- factor(df$group, levels = custom_order)
df <- df[order(df$group),]

# Select top 5 entries for each group with positive statistic and minimum pvals_adj
df_top <- df %>%
  filter(statistic > 0) %>%
  arrange(group, pvals_adj) %>%
  group_by(group) %>%
  slice_head(n = 5) %>%
  ungroup()

y_names <- df_top$names

# Create a sorted mapping table
y_names_df <- data.frame(name = y_names, order = seq_along(y_names))

# Filter and sort df based on y_names
filtered_sorted_df <- df %>%
  inner_join(y_names_df, by = c("names" = "name")) %>%
  arrange(order) %>%
  select(-order) # Remove the auxiliary sorting column 'order'

# Normalize statistics
filtered_sorted_df <- filtered_sorted_df %>%
  group_by(names) %>%
  mutate(z_statistic = (statistic - mean(statistic, na.rm = TRUE)) / 
           sd(statistic, na.rm = TRUE))

# Adjust pvals_adj and z_statistic
filtered_sorted_df <- filtered_sorted_df %>%
  mutate(pvals_adj = ifelse(pvals_adj < 1E-250, 1E-250, pvals_adj)) %>%
  mutate(z_statistic = ifelse(z_statistic > 3, 3, z_statistic))

# Plot data
plot_data <- ggplot(filtered_sorted_df, aes(x = group, y = factor(names, levels = rev(unique(names))), color = z_statistic, size = -log10(pvals_adj))) +
  geom_point() +
  scale_color_gradientn(
    colors = rev(brewer.pal(9, "RdBu")),
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
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5))

plot_data

ggsave("XXX", plot_data, height = 16, width = 9, limitsize = FALSE)
