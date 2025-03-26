library(tidyverse)
library(ggpubr)

# Read data
data <- read.csv('XXX/CD8_clone_id_group.csv', header = TRUE)

# Define colors
colors <- c("#E53238", "#FDCA30", "#0064D2")

# Create boxplots
plot_boxplot <- function(df, y_var, colors) {
  ggplot(df, aes(x = group, y = !!sym(y_var))) + 
    geom_boxplot(outlier.shape = NA, width = 0.25) +
    geom_jitter(aes(color = group), width = 0.25, size = 1) +
    coord_flip() +
    scale_color_manual(values = colors) +
    theme_classic() +
    theme(legend.position = "bottom")
}

# Generate plots
p_total <- plot_boxplot(data, "total_num", colors)

# Filter non-zero values and generate plots
filter_and_plot <- function(data, y_var, colors) {
  df_filtered <- data %>% filter(!!sym(y_var) != 0)
  plot_boxplot(df_filtered, y_var, colors)
}

p_CD8_ISG <- filter_and_plot(data, "CD8_T_ISG", colors)
p_CD8_T_cytotoxic <- filter_and_plot(data, "CD8_T_cytotoxic_02", colors)
p_CD8_T_effector_memory <- filter_and_plot(data, "CD8_T_effector_memory_01", colors)
p_CD8_T_late_exhausted <- filter_and_plot(data, "CD8_T_late_exhausted_01", colors)
p_CD8_T_naive <- filter_and_plot(data, "CD8_T_naive_01", colors)

# Arrange and save plots
p_all <- ggarrange(p_total, p_CD8_T_effector_memory, p_CD8_T_naive, p_CD8_T_cytotoxic, p_CD8_T_late_exhausted, nrow = 1)
ggsave("XXX/CD8_T_clone_size.pdf", plot = p_all, width = 10, height = 4)

# Statistical analysis
group_pairs <- combn(unique(data$group), 2)
test_results <- list()

# Function to perform Mann-Whitney U test
perform_test <- function(data, cluster_name, group1_name, group2_name) {
  group1 <- data %>% filter(group == group1_name) %>% select(cluster_name) %>% filter(!!sym(cluster_name) != 0)
  group2 <- data %>% filter(group == group2_name) %>% select(cluster_name) %>% filter(!!sym(cluster_name) != 0)
  wilcox.test(group1[[cluster_name]], group2[[cluster_name]], alternative = "less")
}

# List of cluster names
cluster_names <- c("total_num", "CD8_T_ISG", "CD8_T_cytotoxic_02", "CD8_T_effector_memory_01", "CD8_T_naive_01", "CD8_T_late_exhausted_01")

# Perform tests and store results
for (cluster_name in cluster_names) {
  for (i in seq(ncol(group_pairs))) {
    group1_name <- group_pairs[1, i]
    group2_name <- group_pairs[2, i]
    test_result <- perform_test(data, cluster_name, group1_name, group2_name)
    test_results[[paste(cluster_name, group1_name, group2_name, sep = "-")]] <- test_result
  }
}

# Print test results
print(test_results)
