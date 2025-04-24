library(tidyverse)
library(ggpubr)

# Load data
data <- read.csv('XXX/BCR_clone_id_group.csv', header = TRUE)


data$group <- factor(data$group, ordered = TRUE, levels = c("mature", "immature", "none"))

# Define colors
colors <- c("#E53238", "#FDCA30", "#0064D2")

# Function to create plots
create_plot <- function(data, y_var) {
  df_filtered <- data %>% filter(!!sym(y_var) != 0)
  
  ggplot(df_filtered, aes(x = group, y = !!sym(y_var))) + 
    geom_boxplot(outlier.shape = NA, width = 0.25) +
    geom_jitter(aes(color = group), width = 0.25, size = 1) +
    coord_flip() +
    scale_color_manual(values = colors) +
    theme_classic() +
    theme(
      legend.position = "bottom",
      axis.text.y = element_blank()
    )
}

# Create individual plots
p_total <- create_plot(data, "total_num")
p_B_IGHM_04 <- create_plot(data, "B_IGHM_04")
p_B_CD74_1_00 <- create_plot(data, "B_CD74_1_00")
p_B_CD69_1_01 <- create_plot(data, "B_CD69_1_01")
p_B_MHC_II_05 <- create_plot(data, "B_MHC_II_05")
p_B_GC_11 <- create_plot(data, "B_GC_11")
p_B_ISG_08 <- create_plot(data, "B_ISG_08")
p_Cycling_B_plamsa_cell_09 <- create_plot(data, "Cycling_B_plamsa_cell_09")
p_Plasma_cell_02_03 <- create_plot(data, "Plasma_cell_02_03")
p_Plasma_cell_01_02 <- create_plot(data, "Plasma_cell_01_02")

# Arrange all plots into one figure
p_all <- ggarrange(p_total, p_B_IGHM_04, p_B_CD74_1_00, p_B_CD69_1_01, p_B_MHC_II_05, p_B_GC_11, p_B_ISG_08, p_Cycling_B_plamsa_cell_09, p_Plasma_cell_02_03, p_Plasma_cell_01_02, nrow = 1) 

# Save the combined plot
ggsave("XXX/B_clone_size.pdf", plot = p_all, width = 15, height = 4)
