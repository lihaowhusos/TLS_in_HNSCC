library(ggplot2)
library(ggrepel)

# Load data
CD4_TCR_table <- read.csv('XXX/CD4_T_and_Treg_plot_clone_size_UMAP_for_R.csv', header = TRUE)

# Subset data by group
CD4_TCR_table_immature <- subset(CD4_TCR_table, group == 'immature')
CD4_TCR_table_mature <- subset(CD4_TCR_table, group == 'mature')
CD4_TCR_table_none <- subset(CD4_TCR_table, group == 'none')

# Define color palette
color <- c('#ea615dff', '#f18988ff', '#6a120c', '#f2cece', '#d8433cff', '#b12921ff', 
                      '#f18da7', '#ec6595', '#c72e82', '#f2ced4', '#580944', '#e0478c', '#9c196e')
                      
# Function to create plot
create_plot <- function(data, fill = FALSE) {
  aes_params <- aes(x = X, y = Y, size = clone_id_size, color = cluster_2)
  if (fill) {
    aes_params <- aes(x = X, y = Y, size = clone_id_size, fill = cluster_2)
  }
  ggplot(data, aes_params) +
    geom_point(color = ifelse(fill, "black", NA), shape = ifelse(fill, 21, 16), alpha = ifelse(fill, 1, 0.25)) +
    scale_color_manual(values = color) +
    scale_fill_manual(values = color) +
    scale_size_continuous(
      limits = c(0, ifelse(fill, 60, 250)),
      breaks = c(0, ifelse(fill, 30, 125), ifelse(fill, 60, 250)),
      labels = c(0, ifelse(fill, 30, 125), ifelse(fill, 60, 250))
    ) +
    theme_bw()
}

# Create plots
p_im_fill <- create_plot(CD4_TCR_table_immature, fill = TRUE)
p_m_fill <- create_plot(CD4_TCR_table_mature, fill = TRUE)
p_n_fill <- create_plot(CD4_TCR_table_none, fill = TRUE)

# Save plots
ggsave("XXX/CD4_T_clone_type_immature.pdf", plot = p_im_fill, width = 7, height = 5)
ggsave("XXX/CD4_T_clone_type_mature.pdf", plot = p_m_fill, width = 7, height = 5)
ggsave("XXX/CD4_T_clone_type_none.pdf", plot = p_n_fill, width = 7, height = 5)