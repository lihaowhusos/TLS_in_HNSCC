library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(MASS)
library(Seurat)
library(stringr)
library(ggrastr)

# 自定义函数:核密度估计
kernel_density_estimation <- function(data, group_var, subset_var, grid_points = 200, bandwidth = 0.4) {
  filtered_data <- filter(data, {{group_var}} %in% subset_var)
  x <- as.numeric(filtered_data$umap_1)
  y <- as.numeric(filtered_data$umap_2)
  density_data <- MASS::kde2d(x, y, n = grid_points, h = bandwidth)
  result <- density_data$z %>%
    as.data.frame() %>%
    as_tibble(rownames = "x") %>%
    pivot_longer(cols = -x, names_to = "y", values_to = "value") %>%
    mutate(y = str_remove_all(y, "V") %>% as.numeric,
           x = as.numeric(x)) %>%
    mutate(subset = paste0(subset_var, collapse = "_"),
           x = scales::rescale(x),
           y = scales::rescale(y))
  return(result)
}

# 自定义函数:核密度对比
kernel_density_contrast <- function(data, group_var, subset_a, subset_b, grid_points = 200, bandwidth = 0.2, quench_factor = 0.015) {
  density_a <- kernel_density_estimation(data, {{group_var}}, subset_a, grid_points, bandwidth)
  density_b <- kernel_density_estimation(data, {{group_var}}, subset_b, grid_points, bandwidth)
  
  contrast_data <- density_a %>%
    left_join(density_b, by = c("x", "y")) %>%
    mutate(value = value.x - value.y) %>%
    mutate(delta_group = value > 0) %>%
    group_by(delta_group) %>%
    mutate(value_quenched = ifelse(value > quench_factor, quench_factor, ifelse(value < -quench_factor, -quench_factor, value))) %>%
    mutate(value_scaled = ifelse(delta_group, scales::rescale(value_quenched, c(0, 2)), scales::rescale(value_quenched, c(-2, 0)))) %>%
    ungroup()
  
  return(contrast_data)
}

# 读取数据
load("G:/TLSdata/TLS_total_data/Processing data/subset/T_NK/T_NK_after_seurat.Rdata")

# 重新分组
combined.sct@meta.data$group <- combined.sct@meta.data$sample
combined.sct@meta.data$group <- recode(combined.sct@meta.data$group,
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

# 创建数据表
data_table <- combined.sct@reductions$umap@cell.embeddings %>%
  as.data.frame() %>%
  cbind(group = combined.sct@meta.data$group)

# 缩放UMAP坐标
data_table$umapscaled_1 <- scales::rescale(data_table$umap_1)
data_table$umapscaled_2 <- scales::rescale(data_table$umap_2)

# 计算核密度对比
contrast_mature <- kernel_density_contrast(data_table, group, c("mature"), c("immature", "none")) %>%
  filter(value_quenched > 0.0001 | value_quenched < -0.0001)

contrast_immature <- kernel_density_contrast(data_table, group, c("immature"), c("mature", "none")) %>%
  filter(value_quenched > 0.0001 | value_quenched < -0.0001)

contrast_none <- kernel_density_contrast(data_table, group, c("none"), c("mature", "immature")) %>%
  filter(value_quenched > 0.0001 | value_quenched < -0.0001)

# 绘制UMAP图
plot_umap <- function(contrast_data, title) {
  ggplot() +
    ggrastr::geom_point_rast(aes(umapscaled_1, umapscaled_2), color = "grey80", size = 0.01, alpha = 0.02, data = data_table) +
    ggrastr::geom_point_rast(aes(x, y, color = value_scaled), data = contrast_data, size = 0.01) +
    scale_color_gradientn(colors = c("#2166AC", "#67A9CF", "#D1E5F0", "#FFFFFF", "#FDDBC7", "#EF8A62", "#B2182B"),
                          values = scales::rescale(c(min(contrast_data$value_scaled), 0, max(contrast_data$value_scaled)))) +
    guides(color = guide_colorbar(label.position = "right", title.position = "top", title.hjust = 0, title.vjust = 1, direction = "vertical")) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          aspect.ratio = 1,
          legend.key.height = unit(0.05, "npc"),
          legend.key.width = unit(0.03, "npc"),
          legend.position = c(-0.1, 0.95),
          legend.justification = c("left", "top"),
          plot.title = element_text(hjust = 0.5, vjust = 0.5, face = "plain", size = 18)) +
    labs(color = paste0("Enrichment\nin ", title))
}

# 保存图像
ggsave("G:/TLSdata/TLS_total_data/figures/enrichment_UMAP/T_NK_mature_new.pdf", plot = plot_umap(contrast_mature, "mature"), width = 6, height = 5)
ggsave("G:/TLSdata/TLS_total_data/figures/enrichment_UMAP/T_NK_immature_new.pdf", plot = plot_umap(contrast_immature, "immature"), width = 6, height = 5)
ggsave("G:/TLSdata/TLS_total_data/figures/enrichment_UMAP/T_NK_none_new.pdf", plot = plot_umap(contrast_none, "none"), width = 6, height = 5)