library(tidyverse)
library(ggpubr)
library(dplyr)

# Read data
data <- read.csv('XXX/CD4_T_and_Treg_clone_id_group.csv', header = TRUE)

# Recode and factorize the group column
data$group <- recode(data$sample, 
                     'szj105988' = 'immature',
                     'szj106005' = 'none',
                     'szj106495' = 'none',
                     'szj106560' = 'immature',                
                     'szj106562' = 'mature',
                     'szj107010' = 'mature',
                     'szj107145' = 'mature',
                     'szj107734' = 'mature',
                     'szj107849' = 'none',
                     'szj108352' = 'immature',
                     'szj106121' = 'mature',                  
                     'szj106138' = 'immature',                  
                     'szj106759' = 'none',
                     'szj106771' = 'mature')

data$group <- factor(data$group, ordered = TRUE, levels = c("mature", "immature", "none"))

# Define colors
colors <- c("#E53238", "#FDCA30", "#0064D2")

# Function to create plots
create_plot <- function(data, column_name) {
  df_filtered <- data %>% filter(!!sym(column_name) != 0)
  ggplot(df_filtered, aes(x = group, y = !!sym(column_name))) + 
    geom_boxplot(outlier.shape = NA, width = 0.25) +
    geom_jitter(aes(color = group), width = 0.25, size = 1) +
    coord_flip() +
    scale_color_manual(values = colors) +
    theme_classic() +
    theme(legend.position = "bottom")
}

# Create plots
p_total <- create_plot(data, "total_num")
p_CD4_T_ISG <- create_plot(data, "CD4_T_ISG")
p_CD4_T_early_exhausted_06 <- create_plot(data, "CD4_T_early_exhausted_06")
p_CD4_T_fh_01_00 <- create_plot(data, "CD4_T_fh_01_00")
p_CD4_T_late_exhausted_08 <- create_plot(data, "CD4_T_late_exhausted_08")
p_CD4_T_naive_1_01 <- create_plot(data, "CD4_T_naive_1_01")
p_CD4_T_naive_2_03 <- create_plot(data, "CD4_T_naive_2_03")

# Arrange plots
p_all <- ggarrange(p_total, p_CD4_T_fh_01_00, p_CD4_T_naive_1_01, p_CD4_T_naive_2_03, p_CD4_T_early_exhausted_06, p_CD4_T_late_exhausted_08, p_CD4_T_ISG, nrow = 1)

# Save plot
ggsave("XXX/CD4_T_clone_size.pdf", plot = p_all, width = 10, height = 4)

# Treg Analysis
data <- read.csv('XXX/CD4_T_and_Treg_clone_id_group.csv', header = TRUE)

data$group <- recode(data$sample, 
                     'szj105988' = 'immature',
                     'szj106005' = 'none',
                     'szj106495' = 'none',
                     'szj106560' = 'immature',                
                     'szj106562' = 'mature',
                     'szj107010' = 'mature',
                     'szj107145' = 'mature',
                     'szj107734' = 'mature',
                     'szj107849' = 'none',
                     'szj108352' = 'immature',
                     'szj106121' = 'mature',                  
                     'szj106138' = 'immature',                  
                     'szj106759' = 'none',
                     'szj106771' = 'mature')

data$group <- factor(data$group, ordered = TRUE, levels = c("mature", "immature", "none"))

# Create Treg plots
p_CD4_T_reg_naive_07 <- create_plot(data, "CD4_T_reg_naive_07")
p_CD4_T_reg_1_04 <- create_plot(data, "CD4_T_reg_1_04")
p_CD4_T_reg_2_05 <- create_plot(data, "CD4_T_reg_2_05")
p_CD4_T_reg_3_10 <- create_plot(data, "CD4_T_reg_3_10")
p_CD4_T_reg_exhausted_01_02 <- create_plot(data, "CD4_T_reg_exhausted_01_02")
p_CD4_T_reg_exhausted_02_09 <- create_plot(data, "CD4_T_reg_exhausted_02_09")
p_CD4_T_reg_ISG <- create_plot(data, "CD4_T_reg_ISG")

# Arrange Treg plots
p_all <- ggarrange(p_CD4_T_reg_naive_07, p_CD4_T_reg_1_04, p_CD4_T_reg_2_05, p_CD4_T_reg_3_10, p_CD4_T_reg_exhausted_01_02, p_CD4_T_reg_exhausted_02_09, p_CD4_T_reg_ISG, nrow = 1)

# Save Treg plot
ggsave("XXX/CD4_Treg_clone_size.pdf", plot = p_all, width = 10, height = 4)

# Statistical Analysis
data <- read.csv('XXX/CD4_T_and_Treg_clone_id_group.csv', header = TRUE)

data$group <- recode(data$sample, 
                     'szj105988' = 'immature',
                     'szj106005' = 'none',
                     'szj106495' = 'none',
                     'szj106560' = 'immature',                
                     'szj106562' = 'mature',
                     'szj107010' = 'mature',
                     'szj107145' = 'mature',
                     'szj107734' = 'mature',
                     'szj107849' = 'none',
                     'szj108352' = 'immature',
                     'szj106121' = 'mature',                  
                     'szj106138' = 'immature',                  
                     'szj106759' = 'none',
                     'szj106771' = 'mature')

data$group <- factor(data$group, ordered = TRUE, levels = c("mature", "immature", "none"))

# Create group pairs
group_pairs <- combn(unique(data$group), 2)

# Perform Mann-Whitney U test
perform_test <- function(data, group_pairs, column_name) {
  test_results <- list()
  for (i in seq(ncol(group_pairs))) {
    group1 <- data %>% filter(group == group_pairs[1, i]) %>% select(column_name) %>% filter(!!sym(column_name) != 0)
    group2 <- data %>% filter(group == group_pairs[2, i]) %>% select(column_name) %>% filter(!!sym(column_name) != 0)
    test_result <- wilcox.test(group1[, column_name], group2[, column_name], alternative = "less")
    test_results[[paste(group_pairs[1, i], group_pairs[2, i], sep = "-")]] <- test_result
  }
  return(test_results)
}

# Perform tests
test_results_total_num <- perform_test(data, group_pairs, "total_num")
test_results_CD4_T_reg_1_04 <- perform_test(data, group_pairs, "CD4_T_reg_1_04")

print(test_results_total_num)
print(test_results_CD4_T_reg_1_04)

colnames(data)