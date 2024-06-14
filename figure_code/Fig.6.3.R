library(Seurat)
library(dplyr)
library(ggplot2)
library(tibble)
library(grid)
library(png)
library(cowplot)
library(tidyverse)
library(pheatmap)

# Define sample name
sample = "XXX"

# Load spatial data
spatial = Load10X_Spatial(paste0("XXX", sample, "-V_outs"))

# Add metadata
metadata_add = read.csv(paste0("XXX", sample, "_cell2location_cell_type.csv"))
metadata_add[,'spot_id'] = sapply(metadata_add[,'spot_id'], function(str) substring(str, 11, nchar(str)-10))
rownames(metadata_add) = metadata_add[,'spot_id']
spatial <- AddMetaData(object = spatial, metadata = metadata_add)

# Extract metadata and coordinates
metadata_ds <- data.frame(spatial@meta.data)
spatial_coord <- spatial@images[["slice1"]]@coordinates %>%
  rownames_to_column("barcodeID") %>%
  mutate(imagerow_scaled = imagerow * spatial@images[["slice1"]]@scale.factors$lowres,
         imagecol_scaled = imagecol * spatial@images[["slice1"]]@scale.factors$lowres) %>%
  inner_join(metadata_ds %>% rownames_to_column("barcodeID"), by = "barcodeID")

# Load tissue image
img <- readPNG(paste0("XXX", sample, "-V_outs", "/spatial/tissue_lowres_image.png"))
img_grob <- rasterGrob(img, interpolate = FALSE, width = unit(1, "npc"), height = unit(1, "npc"))

# Normalize cell data
cell = "B_GC_11"
cell1 = "CD4_T_fh_01_00"
spatial_coord[,cell] = spatial_coord[,cell] / max(spatial_coord[,cell])
spatial_coord[,cell1] = spatial_coord[,cell1] / max(spatial_coord[,cell1])

# Create scatter plot
scatterpie_plt <- ggplot() +   
  annotation_custom(grob = img_grob, xmin = 0, xmax = ncol(img), ymin = 0, ymax = -nrow(img)) +
  geom_point(data = spatial_coord, aes(x = imagecol_scaled, y = imagerow_scaled, size = B_GC_11, alpha = B_GC_11), color = 'yellow') + 
  geom_point(data = spatial_coord, aes(x = imagecol_scaled, y = imagerow_scaled, size = CD4_T_fh_01_00, alpha = CD4_T_fh_01_00), color = 'blue') +
  scale_y_reverse() + ylim(nrow(img), 0) + xlim(0, ncol(img)) + 
  theme_half_open(11, rel_small = 1) + theme_void() + 
  coord_fixed(ratio = 1) +
  scale_size_continuous(range = c(0, 2)) +
  scale_alpha_continuous(range = c(0, 1)) +
  labs(size = cell) + labs(size = cell1) + guides(alpha = "none")

# Save scatter plot
ggsave(paste0("XXX", sample, "_cell2location_cell_type.pdf"), scatterpie_plt, height = 15, width = 20)

# Read metadata for dotplot
metadata_add = read.csv(paste0("XXX/cell2location_cell_type.csv"))

# Extract unique sample names
unique_values <- unique(metadata_add$sample)

# Normalize metadata
metadata_add[, 16:ncol(metadata_add)] <- apply(metadata_add[, 16:ncol(metadata_add)], 2, function(x) x / max(x))
metadata_add[,'spot_id'] = sapply(metadata_add[,'spot_id'], function(str) substring(str, 11, nchar(str)-10))

# Generate dotplots for each sample
for (sample_name in unique_values) {
  folder_name = paste0("XXX/single_per_sample/", sample_name)
  if (!dir.exists(folder_name)){
    dir.create(folder_name)
  } 
  selected_rows <- subset(metadata_add, sample == sample_name)
  rownames(selected_rows) = selected_rows[,'spot_id']
  cell_names <- colnames(selected_rows)[16:ncol(selected_rows)]
  spatial = Load10X_Spatial(paste0("XXX", sample_name, "-V_outs"))
  spatial <- AddMetaData(object = spatial, metadata = selected_rows)
  metadata_ds <- data.frame(spatial@meta.data)
  
  spatial_coord <- spatial@images[["slice1"]]@coordinates %>%
    rownames_to_column("barcodeID") %>%
    mutate(imagerow_scaled = imagerow * spatial@images[["slice1"]]@scale.factors$lowres,
           imagecol_scaled = imagecol * spatial@images[["slice1"]]@scale.factors$lowres) %>%
    inner_join(metadata_ds %>% rownames_to_column("barcodeID"), by = "barcodeID")
  
  img <- readPNG(paste0("XXX", sample_name, "-V_outs", "/spatial/tissue_lowres_image.png"))
  img_grob <- rasterGrob(img, interpolate = FALSE, width = unit(1, "npc"), height = unit(1, "npc"))
  
  for (cell_name in cell_names){
    scatterpie_plt <- ggplot() + 
      annotation_custom(grob = img_grob, xmin = 0, xmax = ncol(img), ymin = 0, ymax = -nrow(img)) +
      geom_point(data = spatial_coord, aes(x = imagecol_scaled, y = imagerow_scaled, alpha = !!sym(cell_name), color = !!sym(cell_name)), size = 3.1) + 
      scale_y_reverse() + ylim(nrow(img), 0) + xlim(0, ncol(img)) + 
      theme_half_open(11, rel_small = 1) + theme_void() + 
      coord_fixed(ratio = 1) +
      scale_color_gradientn(colors = c("white", "white", "red", "red"), values = scales::rescale(c(0, 0.35, 0.6, 1)), limits = c(0, 1)) +
      scale_alpha_continuous(range = c(0, 0.7), limits = c(0.3, 1)) +
      guides(alpha = "none")
    
    ggsave(paste0(folder_name, "/", sample_name, "_cell2location_", cell_name, ".pdf"), scatterpie_plt, height = 12, width = 15)
  }
}

# Co-occurrence analysis
metadata = read.csv(paste0("XXX/cell2location_cell_type.csv"))

# Select columns
metadata <- subset(metadata, select = c("spot_id", "in_tissue", "array_row", "array_col", "sample",
                                        "n_genes_by_counts", "total_counts", "total_counts_mt", "pct_counts_mt",
                                        "total_counts_hb", "pct_counts_hb", "library_id", "X_indices", "X_scvi_batch",
                                        "X_scvi_labels", "B_IGHM_04", "B_CD74_1_00", "B_CD69_1_01", "B_MHC_II_05",
                                        "B_GC_11", "B_ISG_08", "Cycling_B_plamsa_cell_09", "Plasma_cell_02_03",
                                        "Plasma_cell_01_02", "CD8_T_naive_01", "CD8_T_effector_memory_01",
                                        "CD8_T_cytotoxic_02", "CD8_T_late_exhausted_01", "CD8_T_ISG", "CD4_T_fh_01_00",
                                        "CD4_T_naive_1_01", "CD4_T_naive_2_03", "CD4_T_early_exhausted_06",
                                        "CD4_T_late_exhausted_08", "CD4_T_ISG", "CD4_T_reg_ISG", "CD4_T_reg_naive_07",
                                        "CD4_T_reg_1_04", "CD4_T_reg_2_05", "CD4_T_reg_3_10", "CD4_T_reg_exhausted_01_02",
                                        "CD4_T_reg_exhausted_02_09", "cDC_1_15", "cDC_2_06", "DC_LAMP3_11",
                                        "pDC_10", "M1_S100A8_07", "M2_CXCL10_12", "M2_MARCO_05", "M2_STAB1_09",
                                        "M2_SELENOP_02", "M2_MMP9_08", "M2_COL1A1_04", "Cleaning_macrophage_13",
                                        "Cycling_myeloid_cell_14", "Neutrophil_1_01", "Neutrophil_1_03", "Mast_cell_00"))

# Extract unique sample names
unique_values <- unique(metadata$sample)

# Normalize values
metadata[, 16:ncol(metadata)] <- apply(metadata[, 16:ncol(metadata)], 2, function(x) x / max(x))
metadata[, 16:ncol(metadata)] <- scale(metadata[, 16:ncol(metadata)])

colours = colorRampPalette(c('white', '#ed1c24', '#ce181e', '#7f181b'))(500)

# Generate heatmaps for co-occurrence
for (sample_name in unique_values) {
  selected_rows <- subset(metadata, sample == sample_name)
  selected_rows <- selected_rows[,-(1:15)]
  names <- colnames(selected_rows)
  dataframe_hm = data.frame(matrix(ncol = 0, nrow = ncol(selected_rows)))
  for (j in 1:ncol(selected_rows)){
    n=c()
    for (k in 1:ncol(selected_rows)){
      c_name = names[j]
      d_name = names[k]
      n_sub = (sum(selected_rows[[c_name]] * selected_rows[[d_name]])) / nrow(selected_rows)
      n = c(n, n_sub)
    }
    dataframe_hm = cbind(dataframe_hm, n)
    dataframe_hm[j, j] = 0
  }
  
  hm = pheatmap(dataframe_hm, color = colours, scale = "none", cluster_rows = FALSE,
                cluster_cols = FALSE, border_color = NA, display_numbers = F,	
                number_color = "black", labels_row = names, labels_col = names, angle_col = 45,
                breaks = seq(0, 3, length.out = 500))
  ggsave(paste0("XXX/co-occurrence-T_B_Myeloid/", sample_name, "co-occurrence.pdf"), hm, height = 14, width = 14.5)
  ggsave(paste0("XXX/co-occurrence-T_B_Myeloid/", sample_name, "co-occurrence.jpg"), hm, height = 14, width = 14.5)
}