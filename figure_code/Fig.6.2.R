library(dplyr)
library(Seurat)
library(ggplot2)
library(ggrastr)
library(MASS)
library(scales)
library(tidyr)
library(stringr)

# Function to compute KDE and tidy the result
kde2d_tidy <- function(data_tbl, kernel_group, kernel_subset, n = 200, h = 0.4) {
  data_tbl <- filter(data_tbl, kernel_group %in% kernel_subset)
  x <- as.numeric(data_tbl$umap_1)
  y <- as.numeric(data_tbl$umap_2)
  kde_data <- kde2d(x, y, n = n, h = h)
  result_tbl <- kde_data$z %>%
    as.data.frame() %>%
    as_tibble(rownames = "x") %>%
    gather(y, value, -x) %>%
    mutate(y = str_remove_all(y, "V") %>% as.numeric,
           x = as.numeric(x),
           var = paste0(kernel_subset, collapse = "_"),
           x = rescale(x),
           y = rescale(y))
  return(result_tbl)
}

# Function to compute KDE contrast
kde2d_contrast <- function(data_tbl, kernel_group, kernel_subsets_a, kernel_subsets_b, 
                           n = 200, h = 0.2, quench = 0.015) {
  kernel_a <- kde2d_tidy(data_tbl, kernel_group, kernel_subsets_a, n = n, h = h)
  kernel_b <- kde2d_tidy(data_tbl, kernel_group, kernel_subsets_b, n = n, h = h)
  
  kernel_tbl <- kernel_a %>%
    left_join(kernel_b, by = c("x", "y")) %>%
    mutate(value = value.x - value.y,
           delta_group = value > 0,
           value_quenched = ifelse(value > quench, quench, ifelse(value < -quench, -quench, value)),
           value_scaled = ifelse(delta_group, rescale(value_quenched, c(0, 2)), rescale(value_quenched, c(-2, 0)))) %>%
    ungroup()
  
  return(kernel_tbl)
}

# Load data
cluster_type = 'Myeloid'
file <- paste0('XXX/', cluster_type, '_after_seurat.Rdata')
load(file)

# Regroup data
combined.sct@meta.data$group <- recode(combined.sct@meta.data$sample, 
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

# Create table with UMAP embeddings and group
table <- combined.sct@reductions$umap@cell.embeddings %>%
  as.data.frame() %>%
  cbind(group = combined.sct@meta.data$group)

# Plot data preparation
plot_data_sub <- table %>%
  mutate(umapscaled_1 = rescale(umap_1),
         umapscaled_2 = rescale(umap_2))

# Compute kernel tables
kernel_tbl_m <- kde2d_contrast(
  data_tbl = table, kernel_group = table$group,
  kernel_subsets_a = c("mature"), kernel_subsets_b = c("immature", "none")
) %>% filter(value_quenched > 0.0001 | value_quenched < -0.0001)

kernel_tbl_im <- kde2d_contrast(
  data_tbl = table, kernel_group = table$group,
  kernel_subsets_a = c("immature"), kernel_subsets_b = c("mature", "none")
) %>% filter(value_quenched > 0.0001 | value_quenched < -0.0001)

kernel_tbl_n <- kde2d_contrast(
  data_tbl = table, kernel_group = table$group,
  kernel_subsets_a = c("none"), kernel_subsets_b = c("mature", "immature")
) %>% filter(value_quenched > 0.0001 | value_quenched < -0.0001)

# Function to create UMAP plot
create_umap_plot <- function(data_sub, kernel_tbl, title) {
  ggplot() +
    geom_point_rast(aes(umapscaled_1, umapscaled_2), color = "grey80", size = 0.01, alpha = 0.02, data = data_sub) +
    geom_point_rast(aes(x, y, color = value_scaled), data = kernel_tbl, size = 0.01) +
    scale_color_gradientn(colours = c("#2166AC", "#67A9CF", "#D1E5F0", "#FFFFFF", "#FDDBC7", "#EF8A62", "#B2182B"),
                          values = rescale(c(min(kernel_tbl$value_scaled), 0, max(kernel_tbl$value_scaled)))) +
    guides(color = guide_colorbar(label.position = "right", title.position = "top", title.hjust = 0, title.vjust = 1, direction = "vertical")) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          aspect.ratio = 1,
          legend.key.height = unit(0.05, "npc"),
          legend.key.width = unit(0.03, "npc"),
          legend.position = c(-0.1, 0.95),
          legend.justification = c("left", "top"),
          plot.title = element_text(hjust = 0.5, vjust = 0.5, face = "plain", size = 18)) +
    labs(color = title)
}

# Generate UMAP plots
kernel_umap_m <- create_umap_plot(plot_data_sub, kernel_tbl_m, "Enrichment\nin mature")
kernel_umap_im <- create_umap_plot(plot_data_sub, kernel_tbl_im, "Enrichment\nin immature")
kernel_umap_n <- create_umap_plot(plot_data_sub, kernel_tbl_n, "Enrichment\nin none")

# Save figures
ggsave(paste0('XXX/', cluster_type, '_mature.pdf'), plot = kernel_umap_m, width = 6, height = 5)
ggsave(paste0('XXX/', cluster_type, '_immature.pdf'), plot = kernel_umap_im, width = 6, height = 5)
ggsave(paste0('XXX/', cluster_type, '_none.pdf'), plot = kernel_umap_n, width = 6, height = 5)