##Fig6_plot_spatial_data
library(Seurat)
library(argparse)
library(dplyr)
library(ggplot2)

sample = "szj105988"



spatial = Load10X_Spatial(paste0("G:\\TLSdata\\spatial\\",sample,"-V_outs"))

#add metadata
metadata_add = read.csv(paste0("G:\\spatial\\result\\cell2location\\T_NK_B_Myeloid_other\\table\\",sample,"_cell2location_cell_type.csv"))

metadata_add[,'spot_id'] = sapply(metadata_add[,'spot_id'], function(str) substring(str, 11, nchar(str)-10))
rownames(metadata_add) = metadata_add[,'spot_id']

spatial <- AddMetaData(
  object = spatial,
  metadata = metadata_add,
)


##extract information

metadata_ds <- data.frame(spatial@meta.data)
colnames(metadata_ds) <- colnames(spatial@meta.data)


spatial_coord <- spatial@images[["slice1"]]@coordinates %>%
  tibble::rownames_to_column("barcodeID") %>% dplyr::mutate(imagerow_scaled = imagerow *
  spatial@images[["slice1"]]@scale.factors$lowres, imagecol_scaled = imagecol *
  spatial@images[["slice1"]]@scale.factors$lowres) %>% dplyr::inner_join(metadata_ds %>%
  tibble::rownames_to_column("barcodeID"), by = "barcodeID")


img <- png::readPNG(paste0("G:\\TLSdata\\spatial\\",sample,"-V_outs","\\spatial\\tissue_lowres_image.png"))
img_grob <- grid::rasterGrob(img, interpolate = FALSE, width = grid::unit(1,    "npc"), height = grid::unit(1, "npc"))

cell = "B_GC_11"
cell1 = "CD4_T_fh_01_00"
cell2 = 
Max = max(spatial_coord[,cell])
spatial_coord[,cell] = spatial_coord[,cell]/Max
Max1 = max(spatial_coord[,cell1])
spatial_coord[,cell1] = spatial_coord[,cell1]/Max1



scatterpie_plt <- ggplot() +   annotation_custom(grob = img_grob,xmin = 0, xmax = ncol(img), ymin = 0, ymax = -nrow(img))+
  geom_point(data = spatial_coord, size = 2,aes(x = imagecol_scaled, y = imagerow_scaled,size = B_GC_11,alpha = B_GC_11), color = 'yellow') + 
  geom_point(data = spatial_coord, size = 2,aes(x = imagecol_scaled, y = imagerow_scaled,size = CD4_T_fh_01_00,alpha = CD4_T_fh_01_00), color = 'blue') +
  scale_y_reverse() + ylim(nrow(img),0) + 
  xlim(0, ncol(img)) + cowplot::theme_half_open(11,rel_small = 1) + theme_void() + 
  coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE, clip = "on") +
  scale_size_continuous(range=c(0,2)) +
  scale_alpha_continuous(range=c(0,1)) +
  labs(size = cell) + 
  labs(size = cell1) + 
  guides(alpha = "none")


ggsave(paste0("G:\\spatial\\result\\cell2location\\T_NK_B_Myeloid_other\\multiple\\",sample,"_cell2location_cell_type.pdf"),scatterpie_plt,height =15, width=20)

                  
metadata_add = read.csv(paste0("G:\\spatial\\result\\cell2location\\T_NK_B_Myeloid_other\\table\\",sample,"_cell2location_cell_type.csv"))





############################################### dotplot in HE slice
#read all metadata
metadata_add = read.csv(paste0("G:\\spatial\\result\\cell2location\\T_NK_B_Myeloid_other\\table\\cell2location_cell_type.csv"))

#extract sample names
unique_values <- unique(metadata_add$sample)

# max to min = 1 to 0
metadata_add[,16:ncol(metadata_add)] <- apply(metadata_add[, 16:ncol(metadata_add)], 2, function(x) x / max(x))

metadata_add[,'spot_id'] = sapply(metadata_add[,'spot_id'], function(str) substring(str, 11, nchar(str)-10))


# export images pdf
for (sample_name in unique_values) {
  ###creat folder
  folder_name = paste0("G:\\spatial\\result\\cell2location\\T_NK_B_Myeloid_other\\single_per_sample\\",sample_name)
  if (!dir.exists(folder_name)){
    dir.create(folder_name)
  } 
  ###select metadata
  selected_rows <- subset(metadata_add, sample == sample_name)
  rownames(selected_rows) = selected_rows[,'spot_id']
  
  ###select cell type
  cell_names <- colnames(selected_rows)[16:ncol(selected_rows)]

  ### read spatial data
  spatial = Load10X_Spatial(paste0("G:\\TLSdata\\spatial\\",sample_name,"-V_outs"))
  ### add metadata
  spatial <- AddMetaData(object = spatial,metadata = selected_rows)
  
  ### extract metadata
  metadata_ds <- data.frame(spatial@meta.data)
  colnames(metadata_ds) <- colnames(spatial@meta.data)
  
  spatial_coord <- spatial@images[["slice1"]]@coordinates %>%
    tibble::rownames_to_column("barcodeID") %>% dplyr::mutate(imagerow_scaled = imagerow *
    spatial@images[["slice1"]]@scale.factors$lowres, imagecol_scaled = imagecol *
    spatial@images[["slice1"]]@scale.factors$lowres) %>% dplyr::inner_join(metadata_ds %>%
    tibble::rownames_to_column("barcodeID"), by = "barcodeID")
  
  ### read image of HE
  img <- png::readPNG(paste0("G:\\TLSdata\\spatial\\",sample_name,"-V_outs","\\spatial\\tissue_lowres_image.png"))
  img_grob <- grid::rasterGrob(img, interpolate = FALSE, width = grid::unit(1,    "npc"), height = grid::unit(1, "npc"))
  
  for (cell_name in cell_names){
    scatterpie_plt <- ggplot() + annotation_custom(grob = img_grob,xmin = 0, xmax = ncol(img), ymin = 0, ymax = -nrow(img)) +
      geom_point(data = spatial_coord,aes(x = imagecol_scaled, y = imagerow_scaled,alpha = !!sym(cell_name), color = !!sym(cell_name)),size = 3.1) + 
      scale_y_reverse() + ylim(nrow(img),0) + 
      xlim(0, ncol(img)) + cowplot::theme_half_open(11,rel_small = 1) + theme_void() + 
      coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE, clip = "on") +
      #scale_size_continuous(range=c(0,2)) +
      scale_color_gradientn(colors = c("white","white", "red","red"),
                            values = scales::rescale(c(0,0.35,0.6,1)),
                           limits = c(0, 1)) +
      scale_alpha_continuous(range = c(0, 0.7),limits = c(0.3, 1)) 
      guides(alpha = "none")
      
    scatterpie_plt = scatterpie_plt + geom_point(data = spatial_coord,aes(x = imagecol_scaled, y = imagerow_scaled,alpha = !!sym(cell_name), color = !!sym(cell_name)),size = 3.1) + 
      scale_y_reverse() + ylim(nrow(img),0) + 
      xlim(0, ncol(img)) + cowplot::theme_half_open(11,rel_small = 1) + theme_void() + 
      coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE, clip = "on") +
      #scale_size_continuous(range=c(0,2)) +
      scale_color_gradientn(colors = c("white","white", "red","red"),
                            values = scales::rescale(c(0,0.35,0.6,1)),
                            limits = c(0, 1))
    scale_alpha_continuous(range = c(0, 0.7),limits = c(0.3, 1)) 
    
    
     ggsave(paste0(folder_name,"\\",sample_name,"_cell2location_",cell_name,".pdf"),scatterpie_plt,height =12, width=15)
  }
}




############################# co-occurrence ############
library(tidyverse)
library(stringr)
library(tidyr)
library(reshape2)
library(forcats)
library(pheatmap)
library(dplyr)


#read all metadata
metadata = read.csv(paste0("G:\\spatial\\result\\cell2location\\T_NK_B_Myeloid_other\\table\\cell2location_cell_type.csv"))
metadata = read.csv(paste0("G:\\spatial\\result\\cell2location\\T_NK_B_Myeloid_other\\table\\cell2location_cell_type_T_B.csv"))


# 选择列
metadata <- subset(metadata, select = c("spot_id",
                                    "in_tissue",
                                    "array_row",
                                    "array_col",
                                    "sample",
                                    "n_genes_by_counts",
                                    "total_counts",
                                    "total_counts_mt",
                                    "pct_counts_mt",
                                    "total_counts_hb",
                                    "pct_counts_hb",
                                    "library_id",
                                    "X_indices",
                                    "X_scvi_batch",
                                    "X_scvi_labels",
                                    "B_IGHM_04",
                                    "B_CD74_1_00",
                                    "B_CD69_1_01",
                                    "B_MHC_II_05",
                                    "B_GC_11",
                                    "B_ISG_08",
                                    "Cycling_B_plamsa_cell_09",
                                    "Plasma_cell_02_03",
                                    "Plasma_cell_01_02",
                                    "CD8_T_naive_01",
                                    "CD8_T_effector_memory_01",
                                    "CD8_T_cytotoxic_02",
                                    "CD8_T_late_exhausted_01",
                                    "CD8_T_ISG",
                                    "CD4_T_fh_01_00",
                                    "CD4_T_naive_1_01",
                                    "CD4_T_naive_2_03",
                                    "CD4_T_early_exhausted_06",
                                    "CD4_T_late_exhausted_08",
                                    "CD4_T_ISG",
                                    "CD4_T_reg_ISG",
                                    "CD4_T_reg_naive_07",
                                    "CD4_T_reg_1_04",
                                    "CD4_T_reg_2_05",
                                    "CD4_T_reg_3_10",
                                    "CD4_T_reg_exhausted_01_02",
                                    "CD4_T_reg_exhausted_02_09",
                                    "cDC_1_15", 
                                    "cDC_2_06", 
                                    "DC_LAMP3_11",
                                    "pDC_10",
                                    "M1_S100A8_07",   
                                    "M2_CXCL10_12",
                                    "M2_MARCO_05",
                                    "M2_STAB1_09",
                                    "M2_SELENOP_02",   
                                    "M2_MMP9_08",
                                    "M2_COL1A1_04",
                                    "Cleaning_macrophage_13",
                                    "Cycling_myeloid_cell_14",   
                                    "Neutrophil_1_01",
                                    "Neutrophil_1_03",
                                    "Mast_cell_00"
                                    ))




#extract sample names
unique_values <- unique(metadata$sample)

# max to min = 1 to 0
metadata[,16:ncol(metadata)] <- apply(metadata[, 16:ncol(metadata)], 2, function(x) x / max(x))
metadata[, 16:ncol(metadata)] <- scale(metadata[, 16:ncol(metadata)])

names <- colnames(metadata)

colours=colorRampPalette(c("white", "#1874CD","#FFC125", "#EE4000","#CD0000"))(500)

colours=colorRampPalette(c("white", "#1874CD", "#FFC125","#FFC125","#EE4000", "#EE4000","#CD0000","#CD0000"))(500)

colours=colorRampPalette(c("white", "#FFC125","#EE4000","#CD0000"))(500)

colours=colorRampPalette(c("white", "#FFC125","#FFC125","#EE4000", "#EE4000","#CD0000","#CD0000"))(500)

colours = colorRampPalette(c('white', '#ed1c24','#ce181e','#7f181b'))(500)



for (sample_name in unique_values) {
  ###select metadata
  selected_rows <- subset(metadata, sample == sample_name)
  selected_rows <- selected_rows[,-(1:15)]
  names <- colnames(selected_rows)
  dataframe_hm = data.frame(matrix(ncol = 0, nrow = ncol(selected_rows)))
  for (j in 1:ncol(selected_rows)){
    n=c()
    for (k in 1:ncol(selected_rows)){
      c_name = names[j]
      d_name = names[k]
      n_sub = (sum(selected_rows[[c_name]] * selected_rows[[d_name]]))/nrow(selected_rows)*1000
      n = c(n,n_sub)
    }
    dataframe_hm = cbind(dataframe_hm, n)
    dataframe_hm[j,j] = 0
  }
  any(is.na(dataframe_hm))
  
  hm = pheatmap(dataframe_hm, color = colours,
                              scale = "none", cluster_rows = FALSE,
                              cluster_cols = FALSE,border_color = NA,display_numbers=F,	
                              number_color = "black",
                              labels_row = names,labels_col = names,angle_col = 45,  # 添加这一行来旋转x轴标签
                breaks = seq(0, 10, length.out=500,)
  )
  ggsave(paste0("G:\\spatial\\result\\cell2location\\T_NK_B_Myeloid_other\\co-occurrence-T_B_Myeloid\\",sample_name,"co-occurrence.pdf"), hm, height =14, width=14.5)
  ggsave(paste0("G:\\spatial\\result\\cell2location\\T_NK_B_Myeloid_other\\co-occurrence-T_B_Myeloid\\",sample_name,"co-occurrence.jpg"), hm, height =14, width=14.5)
}


for (sample_name in unique_values) {
  ###select metadata
  selected_rows <- subset(metadata, sample == sample_name)
  selected_rows <- selected_rows[,-(1:15)]
  names <- colnames(selected_rows)
  dataframe_hm = data.frame(matrix(ncol = 0, nrow = ncol(selected_rows)))
  for (j in 1:ncol(selected_rows)){
    n=c()
    for (k in 1:ncol(selected_rows)){
      c_name = names[j]
      d_name = names[k]
      n_sub = (sum(selected_rows[[c_name]] * selected_rows[[d_name]]))/nrow(selected_rows)
      n = c(n,n_sub)
    }
    dataframe_hm = cbind(dataframe_hm, n)
    dataframe_hm[j,j] = 0
  }
  any(is.na(dataframe_hm))
  
  hm = pheatmap(dataframe_hm, color = colours,
                scale = "none", cluster_rows = FALSE,
                cluster_cols = FALSE,border_color = NA,display_numbers=F,	
                number_color = "black",
                labels_row = names,labels_col = names,angle_col = 45,  # 添加这一行来旋转x轴标签
                breaks = seq(0, 3, length.out=500)
  )
  ggsave(paste0("G:\\spatial\\result\\cell2location\\T_NK_B_Myeloid_other\\co-occurrence-T_B_Myeloid\\",sample_name,"co-occurrence.pdf"), hm, height =14, width=14.5)
  ggsave(paste0("G:\\spatial\\result\\cell2location\\T_NK_B_Myeloid_other\\co-occurrence-T_B_Myeloid\\",sample_name,"co-occurrence.jpg"), hm, height =14, width=14.5)
}








