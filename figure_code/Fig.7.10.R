library(NMF)
library(pheatmap)
library(dplyr)
library(stringr)

# Remove irrelevant data
rm(list=ls())

# Set working directory
setwd("XXX")

# Read data
rt <- read.csv("XXX/Bulk_RNA_cell2location_20240202.csv", sep=",", header=TRUE, check.names=FALSE)
names(rt)[1] <- "Id"

# Filter data
rt <- rt %>% filter(str_detect(Id, "01.{1}$"))
length(unique(rt$Id))

# Collate matrix format
rt <- as.matrix(rt)
rownames(rt) <- rt[,1]
exp <- rt[,2:ncol(rt)]
exp <- t(exp)

dimnames <- list(rownames(exp), colnames(exp))
data <- matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)

# Standardize function
standardize <- function(x) {
  rowmean <- apply(x, 1, mean)
  rowsd <- apply(x, 1, sd)
  rv <- sweep(x, 1, rowmean, "-")
  rv <- sweep(rv, 1, rowsd, "/")
  return(rv)
}

# Normalize function
normalize <- function(x) {
  rowmin <- apply(x, 1, min)
  rowmax <- apply(x, 1, max)
  rowmax.min <- rowmax - rowmin
  rv <- sweep(x, 1, rowmin, "-")
  rv <- sweep(rv, 1, rowmax.min, "/")
  return(rv)
}

data <- normalize(data)
data.standardize <- standardize(data)

# Remove rows containing "Cycling"
rows_to_remove <- grep("Cycling", rownames(data))
data <- data[-rows_to_remove, ]

# NMF analysis
ranks <- 2:10
estim.coad <- nmf(data, ranks, nrun=50)
plot(estim.coad)

# Save plot
pdf("XXX/NMF/estim.coad_20240203.pdf")
plot(estim.coad)
dev.off()

# NMF with rank = 3
seed <- 20231120
rank <- 3
nmf.rank <- nmf(data, rank=rank, nrun=50, seed=seed, method="brunet")

# Set colors
jco <- c("#2874C5", "#EABF00", "#C6524A", "#868686")
index <- extractFeatures(nmf.rank, "max")
sig.order <- unlist(index)
nmf.Exp.rank <- data[sig.order, ]
nmf.Exp.rank <- na.omit(nmf.Exp.rank)
group <- predict(nmf.rank)
table(group)

# Plot consensus map
consensusmap(nmf.rank, labRow=NA, labCol=NA, annCol=data.frame("cluster"=group[colnames(nmf.rank)]), annColors=list(cluster=c("1"=jco[1], "2"=jco[2], "3"=jco[3], "4"=jco[4])))
consensusmap(nmf.rank)
coefmap(nmf.rank)

# Save group information
group <- predict(nmf.rank)
group <- as.data.frame(group)
group$group <- paste0('Cluster', group$group)
group$sample <- rownames(group)
group <- group[order(group$group), ]
table(group$group)

# Save group data to CSV
write.csv(group, paste0("XXX/NMF/group_", rank, ".csv"))

# Plot with pheatmap
file_path <- paste0("XXX/NMF/group_", rank, ".csv")
cluster <- read.csv(file_path, header=TRUE, row.names=1)
cluster <- select(cluster, -last_col())

data_df <- as.data.frame(t(data))
cluster$ID <- rownames(cluster)
data_df$ID <- rownames(data_df)

merged_df <- merge(cluster, data_df, by="ID", all.x=TRUE)
merged_df <- merged_df[order(merged_df$group), ]

# Create cluster_new
cluster_new <- data.frame(merged_df[,"group"], row.names=merged_df[,"ID"])
names(cluster_new) <- 'group'

# Create matrix
rownames(merged_df) <- merged_df$ID
merged_df <- select(merged_df, -group)
merged_df <- select(merged_df, -ID)

merged_df_matrix <- as.matrix(merged_df)
merged_df_matrix <- t(as.matrix(merged_df))
rownames_merged_df_matrix <- substr(rownames(merged_df_matrix), 24, nchar(rownames(merged_df_matrix)))
colnames_merged_df_matrix <- colnames(merged_df_matrix)

merged_df_matrix <- t(apply(merged_df_matrix, 1, scale))

rownames(merged_df_matrix) <- rownames_merged_df_matrix
colnames(merged_df_matrix) <- colnames_merged_df_matrix

# Re-order cell type
desired_rownames <- c("B_IGHM_04", "B_CD74_1_00", "B_CD69_1_01", "B_MHC_II_05", "B_GC_11", "B_ISG_08", "Plasma_cell_02_03", "Plasma_cell_01_02", "CD8_T_naive_01", "CD8_T_effector_memory_01", "CD8_T_cytotoxic_02", "CD8_T_late_exhausted_01", "CD8_T_ISG", "CD4_T_fh_01_00", "CD4_T_naive_1_01", "CD4_T_naive_2_03", "CD4_T_early_exhausted_06", "CD4_T_late_exhausted_08", "CD4_T_ISG", "CD4_T_reg_naive_07", "CD4_T_reg_1_04", "CD4_T_reg_2_05", "CD4_T_reg_3_10", "CD4_T_reg_exhausted_01_02", "CD4_T_reg_exhausted_02_09", "CD4_T_reg_ISG", "NKT_1_00", "NKT_2_02", "NK_cytotoxic_04", "NK_reg_CRTAM_03", "NK_reg_KRT86_01", "ILC_06", "dg_T_05", "cDC_1_15", "cDC_2_06", "DC_LAMP3_11", "pDC_10", "M1_S100A8_07", "M2_CXCL10_12", "M2_MARCO_05", "M2_STAB1_09", "M2_SELENOP_02", "M2_MMP9_08", "M2_COL1A1_04", "Cleaning_macrophage_13", "Neutrophil_1_01", 'Neutrophil_1_03', 'Mast_cell_00')

sorted_indices <- match(desired_rownames, rownames(merged_df_matrix))
sorted_matrix <- merged_df_matrix[sorted_indices, ]

# Add TLS imprint
TLS_imprint <- read.csv("XXX/TLS_imprint/TCGA_HNSCC_TLS_imprint_aucell.csv")
TLS_imprint <- t(as.matrix(TLS_imprint))
colnames(TLS_imprint) <- TLS_imprint[1,]

matching_colnames <- intersect(colnames(sorted_matrix), colnames(TLS_imprint))
new_row <- TLS_imprint[2, matching_colnames]
sorted_matrix <- rbind(sorted_matrix, as.numeric(new_row))

# Plot heatmap
colors33 <- colorRampPalette(c("white", "#ff9517", "#b83835", '#6d4490', "#000000"))(500)
colorRdBu <- colorRampPalette(c("white", "#f57814", "#e2001a", "#b41428", "black"))(500)

breaks <- seq(0, 5, length.out=length(colors33) + 1)

plot.scale <- pheatmap(sorted_matrix, scale="row", cluster_rows=FALSE, cluster_cols=FALSE, border_color=NA, color=colors33, breaks=breaks, annotation_col=cluster_new, show_colnames=FALSE)

pdf("XXX/NMF/cluster.pdf")
print(plot.scale)
dev.off()

breaks <- seq(0, 5, length.out=length(colorRdBu) + 1)

plot.TLS.imprint <- pheatmap(sorted_matrix, scale="row", cluster_rows=FALSE, cluster_cols=FALSE, border_color=NA, color=colorRdBu, breaks=breaks, annotation_col=cluster_new, show_colnames=FALSE)

pdf("XXX/NMF/TLS_imprint.pdf")
print(plot.TLS.imprint)
dev.off()