library(dplyr)
library(tidyverse) 
library(magrittr)
library(Seurat)
library(readxl)
library(cowplot)
library(viridis)
library(grid)
library(rlang)
library(ggplot2)
library(ggrastr)



################### function definition####################
kde2d_tidy <- function(data_tbl, kernel_group, kernel_subset,  
                       n = 200, h = 0.4) {
  
  data_tbl <- filter(data_tbl, kernel_group %in% kernel_subset)
  x <- as.numeric(data_tbl$umap_1)
  y <- as.numeric(data_tbl$umap_2)
  kde_data <- MASS::kde2d(x, y, n = n, h = h) 
  var <- as.name(kernel_group)
  result_tbl <- kde_data$z %>% 
    as.data.frame() %>% 
    as_tibble(rownames = "x") %>% 
    gather(y, value, -x) %>% 
    mutate(y = str_remove_all(y, "V") %>% as.numeric,
           x = as.numeric(x)) %>%
    mutate(var := paste0(kernel_subset, collapse = "_"),
           x = scales::rescale(x),
           y = scales::rescale(y))
  return(result_tbl)
  
}



kde2d_contrast <- function(data_tbl, kernel_group, kernel_subsets_a, kernel_subsets_b, 
                           n = 200, h = 0.2, quench = 0.015) { #################change quench
  
  kernel_a <- kde2d_tidy(data_tbl, kernel_group, kernel_subsets_a, n = n, h = h)
  kernel_b <- kde2d_tidy(data_tbl, kernel_group, kernel_subsets_b, n = n, h = h)
  
  kernel_tbl <- kernel_a %>% 
    left_join(kernel_b, by = c("x", "y")) %>% 
    mutate(value = value.x - value.y) %>% 
    mutate(delta_group = value > 0) %>% 
    group_by(delta_group) %>% 
    mutate(value_quenched = ifelse(value > quench, quench, ifelse(value < -quench, -quench, value))) %>% 
    mutate(value_scaled = ifelse(delta_group, scales::rescale(value_quenched, c(0, 2)), scales::rescale(value_quenched, c(-2, 0)))) %>% 
    ungroup()
  
  return(kernel_tbl)
}




########################## Read data ###################
cluster_type = 'Myeloid'

file <- paste0('g:/TLSdata/TLS_total_data/Processing data/subset/',cluster_type,'/',cluster_type,'_after_seurat.Rdata') 
load(file)

###########regroup############
combined.sct@meta.data$group <- combined.sct@meta.data$sample
combined.sct@meta.data$group <- recode(combined.sct@meta.data$group, 
                                       'szj105988'= 'immature',
                                       'szj106005'= 'none',
                                       'szj106495'= 'none',
                                       'szj106560'= 'immature',                
                                       'szj106562'= 'mature',
                                       'szj107010'= 'mature',
                                       'szj107145'= 'mature',
                                       'szj107734'= 'mature',
                                       'szj107849'= 'none',
                                       'szj108352'= 'immature',
                                       'szj106121'= 'mature',                  
                                       'szj106138'= 'immature',                  
                                       'szj106759'= 'none',
                                       'szj106771'= 'mature'
)



#### table ###########
table = combined.sct@reductions$umap@cell.embeddings %>% 
  as.data.frame() %>%
  cbind(group = combined.sct@meta.data$group)




##### plot_data_sub ########
plot_data_sub = table
plot_data_sub$umapscaled_1 = scales::rescale(plot_data_sub$umap_1)
plot_data_sub$umapscaled_2 = scales::rescale(plot_data_sub$umap_2)





###########kernel_table#############
kernel_tbl_m <- kde2d_contrast(
  data_tbl = table, kernel_group = table$group,
  kernel_subsets_a = c("mature"), kernel_subsets_b = c("immature", "none")
) %>%   filter(value_quenched > 0.0001 | value_quenched < -0.0001)

kernel_tbl_im <- kde2d_contrast(
  data_tbl = table, kernel_group = table$group,
  kernel_subsets_a = c("immature"), kernel_subsets_b = c("mature", "none")
) %>%   filter(value_quenched > 0.0001 | value_quenched < -0.0001)

kernel_tbl_n <- kde2d_contrast(
  data_tbl = table, kernel_group = table$group,
  kernel_subsets_a = c("none"), kernel_subsets_b = c("mature", "immature")
) %>%   filter(value_quenched > 0.0001 | value_quenched < -0.0001)




###########kernel_umap_manual#############
kernel_umap_m <- ggplot() +
  ggrastr::geom_point_rast(aes(umapscaled_1, umapscaled_2), color = "grey80", size = 0.01, alpha = 0.02,
                           data = plot_data_sub) +
  ggrastr::geom_point_rast(aes(x, y, color = value_scaled), data = kernel_tbl_m, size = 0.01) +
  
  scale_color_gradientn(colours = c("#2166AC", "#67A9CF", "#D1E5F0", "#FFFFFF", "#FDDBC7", "#EF8A62", "#B2182B"),
                        values = scales::rescale(c(min(kernel_tbl_m$value_scaled), 0,
                                                   max(kernel_tbl_m$value_scaled)))) +
  guides(color = guide_colorbar(label.position = "right",
                                title.position = "top",
                                title.hjust = 0,
                                title.vjust = 1,
                                direction = "vertical")) +
  theme_bw()+
  theme(panel.grid=element_blank())+
  theme(aspect.ratio = 1,
        legend.key.height = unit(0.05, "npc"),
        legend.key.width = unit(0.03, "npc"),
        legend.position = c(-0.1, 0.95),
        legend.justification = c("left", "top"),
        plot.title = element_text(hjust = 0.5, vjust = 0.5, face = "plain", size = 18)) + 
  labs(color = "Enrichment\nin mature"
  ) 

kernel_umap_im <- ggplot() +
  ggrastr::geom_point_rast(aes(umapscaled_1, umapscaled_2), color = "grey80", size = 0.01, alpha = 0.02,
                           data = plot_data_sub) +
  ggrastr::geom_point_rast(aes(x, y, color = value_scaled), data = kernel_tbl_im, size = 0.01) +
  
  scale_color_gradientn(colours = c("#2166AC", "#67A9CF", "#D1E5F0", "#FFFFFF", "#FDDBC7", "#EF8A62", "#B2182B"),
                        values = scales::rescale(c(min(kernel_tbl_im$value_scaled), 0,
                                                   max(kernel_tbl_im$value_scaled)))) +
  guides(color = guide_colorbar(label.position = "right",
                                title.position = "top",
                                title.hjust = 0,
                                title.vjust = 1,
                                direction = "vertical")) +
  theme_bw()+
  theme(panel.grid=element_blank())+
  theme(aspect.ratio = 1,
        legend.key.height = unit(0.05, "npc"),
        legend.key.width = unit(0.03, "npc"),
        legend.position = c(-0.1, 0.95),
        legend.justification = c("left", "top"),
        plot.title = element_text(hjust = 0.5, vjust = 0.5, face = "plain", size = 18)) + 
  labs(color = "Enrichment\nin mature"
  ) 



kernel_umap_n <- ggplot() +
  ggrastr::geom_point_rast(aes(umapscaled_1, umapscaled_2), color = "grey80", size = 0.01, alpha = 0.02,
                           data = plot_data_sub) +
  ggrastr::geom_point_rast(aes(x, y, color = value_scaled), data = kernel_tbl_n, size = 0.01) +
  
  scale_color_gradientn(colours = c("#2166AC", "#67A9CF", "#D1E5F0", "#FFFFFF", "#FDDBC7", "#EF8A62", "#B2182B"),
                        values = scales::rescale(c(min(kernel_tbl_n$value_scaled), 0,
                                                   max(kernel_tbl_n$value_scaled)))) +
  guides(color = guide_colorbar(label.position = "right",
                                title.position = "top",
                                title.hjust = 0,
                                title.vjust = 1,
                                direction = "vertical")) +
  theme_bw()+
  theme(panel.grid=element_blank())+
  theme(aspect.ratio = 1,
        legend.key.height = unit(0.05, "npc"),
        legend.key.width = unit(0.03, "npc"),
        legend.position = c(-0.1, 0.95),
        legend.justification = c("left", "top"),
        plot.title = element_text(hjust = 0.5, vjust = 0.5, face = "plain", size = 18)) + 
  labs(color = "Enrichment\nin mature"
  ) 


########### save figure ################
getwd()


ggsave(paste0('G:/TLSdata/TLS_total_data/figures/enrichment_UMAP/',cluster_type,'_mature.pdf') ,plot = kernel_umap_m,width = 6, height = 5)
ggsave(paste0('G:/TLSdata/TLS_total_data/figures/enrichment_UMAP/',cluster_type,'_immature.pdf') ,plot = kernel_umap_im,width = 6, height = 5)
ggsave(paste0('G:/TLSdata/TLS_total_data/figures/enrichment_UMAP/',cluster_type,'_none.pdf') ,plot = kernel_umap_n,width = 6, height = 5)




