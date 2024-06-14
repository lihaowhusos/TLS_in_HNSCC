library(Seurat)
library(dplyr)
library(ggplot2)
library(png)
library(grid)
library(cowplot)

sample <- "XXX"

# Load spatial data
spatial <- Load10X_Spatial(paste0("XXX", sample, "-V_outs"))

# Add metadata
metadata_add <- read.csv(paste0("XXX", sample, "_cell2location_cell_type.csv"))
metadata_add$spot_id <- sapply(metadata_add$spot_id, function(str) substring(str, 11, nchar(str) - 10))
rownames(metadata_add) <- metadata_add$spot_id

spatial <- AddMetaData(object = spatial, metadata = metadata_add)
spatial <- SCTransform(spatial, assay = "Spatial", verbose = FALSE)

# Plot spatial feature
p1 <- SpatialFeaturePlot(spatial, features = "TCF7", image.alpha = 0) +
  scale_color_gradientn(colors = c("white", "white", "red", "red"),
                        values = scales::rescale(c(0, 0.35, 0.6, 1)),
                        limits = c(0, 1))

# Extract information
metadata_ds <- spatial@meta.data %>% as.data.frame()
colnames(metadata_ds) <- colnames(spatial@meta.data)

spatial_coord <- spatial@images[["slice1"]]@coordinates %>%
  tibble::rownames_to_column("barcodeID") %>%
  mutate(imagerow_scaled = imagerow * spatial@images[["slice1"]]@scale.factors$lowres,
         imagecol_scaled = imagecol * spatial@images[["slice1"]]@scale.factors$lowres) %>%
  inner_join(metadata_ds %>% tibble::rownames_to_column("barcodeID"), by = "barcodeID")

img <- readPNG(paste0("XXX", sample, "-V_outs/spatial/tissue_lowres_image.png"))
img_grob <- rasterGrob(img, interpolate = FALSE, width = unit(1, "npc"), height = unit(1, "npc"))

cell <- "B_GC_11"
cell1 <- "CD4_T_fh_01_00"
Max <- max(spatial_coord[[cell]])
spatial_coord[[cell]] <- spatial_coord[[cell]] / Max
Max1 <- max(spatial_coord[[cell1]])
spatial_coord[[cell1]] <- spatial_coord[[cell1]] / Max1

# Plot scatter pie
scatterpie_plt <- ggplot() +
  annotation_custom(grob = img_grob, xmin = 0, xmax = ncol(img), ymin = 0, ymax = -nrow(img)) +
  geom_point(data = spatial_coord, aes(x = imagecol_scaled, y = imagerow_scaled, size = B_GC_11, alpha = B_GC_11), color = 'yellow') +
  geom_point(data = spatial_coord, aes(x = imagecol_scaled, y = imagerow_scaled, size = CD4_T_fh_01_00, alpha = CD4_T_fh_01_00), color = 'blue') +
  scale_y_reverse() +
  ylim(nrow(img), 0) +
  xlim(0, ncol(img)) +
  theme_void() +
  coord_fixed(ratio = 1) +
  scale_size_continuous(range = c(0, 2)) +
  scale_alpha_continuous(range = c(0, 1)) +
  labs(size = cell) +
  labs(size = cell1) +
  guides(alpha = "none")

# Save plot
ggsave(paste0("XXX", sample, "_cell2location_cell_type.pdf"), scatterpie_plt, height = 15, width = 20)