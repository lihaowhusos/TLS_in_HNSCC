library(ggplot2)

# Read data
CD8_TCR_table <- read.csv('XXX/CD8_T_plot_clone_size_UMAP_for_R.csv', header = TRUE)
CD8_TCR_table$clone_id_size10 <- CD8_TCR_table$clone_id_size / 10

# Subset data based on group
CD8_TCR_table_immature <- subset(CD8_TCR_table, group == 'immature')
CD8_TCR_table_mature <- subset(CD8_TCR_table, group == 'mature')
CD8_TCR_table_none <- subset(CD8_TCR_table, group == 'none')

# Define color palette
color <- c('#d4b44eff', '#e9dd7bff', '#f1f1c1ff', '#a87d2bff', '#603b10ff')

# Function to generate plots
generate_plot <- function(data, fill = FALSE) {
  aes_params <- aes(x = X, y = Y, size = clone_id_size, color = cluster_2)
  if (fill) {
    aes_params <- aes(x = X, y = Y, size = clone_id_size, fill = cluster_2)
  }
  
  p <- ggplot(data, aes_params) +
    geom_point(alpha = ifelse(fill, 1, 0.25), colour = ifelse(fill, "black", NA), shape = ifelse(fill, 21, 16)) +
    scale_color_manual(values = color) +
    scale_fill_manual(values = color) +
    scale_size_continuous(
      range = ifelse(fill, c(1, 9), NULL),
      limits = c(0, 250),
      breaks = c(0, 125, 250),
      labels = c(0, 125, 250)
    ) +
    theme_bw()
  return(p)
}

# Generate and save plots
p_im_fill <- generate_plot(CD8_TCR_table_immature, fill = TRUE)
p_m_fill <- generate_plot(CD8_TCR_table_mature, fill = TRUE)
p_n_fill <- generate_plot(CD8_TCR_table_none, fill = TRUE)

ggsave("XXX/CD8_T_clone_type_immature.pdf", plot = p_im_fill, width = 7, height = 5)
ggsave("XXX/CD8_T_clone_type_mature.pdf", plot = p_m_fill, width = 7, height = 5)
ggsave("XXX/CD8_T_clone_type_none.pdf", plot = p_n_fill, width = 7, height = 5)