####NMF
##based on Cell2location
#options(BioC_mirror="https://mirrors.tuna.tsinghua.edu.cn/bioconductor")

#install.packages("NMF")

# 删除无关数据
rm(list=ls())

setwd("G:\\TLSdata\\TLS_total_data\\Bulk_RNA")      

rt=read.csv("Bulk_RNA_cell2location_20240202.csv",sep=",",header=T,check.names=F)
names(rt)[1] <- "Id"
#cell_type2 <- read.csv("cell_type_xcell.csv")
rt <- rt %>% 
  filter(str_detect(Id, "01.{1}$"))
length(unique(rt$Id))


##Collated matrix format##
rt=as.matrix(rt)
rownames(rt)=rt[,1] 
exp=rt[,2:ncol(rt)]
exp=t(exp)


dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)


#####??׼####
standardize <- function(x) {
  rowmean <- apply(x, 1, mean)
  rowsd <- apply(x, 1, sd)  
  rv <- sweep(x, 1, rowmean,"-")  #????��-??ֵ
  rv <- sweep(rv, 1, rowsd, "/")  #?ٳ??Ա?׼??
  return(rv)
}
######
normalize <- function(x) {
  rowmin <- apply(x, 1, min)
  rowmax <- apply(x, 1, max)
  rowmax.min <- rowmax-rowmin
  rv <- sweep(x, 1, rowmin,"-")  
  rv <- sweep(rv, 1, rowmax.min, "/")  
  return(rv)
}



data = normalize(data)
data.standardize = standardize(data)

rows_to_remove <- grep("Cycling", rownames(data))
data <- data[-rows_to_remove, ]

dim(data)



############ NMF ##################
library(NMF) # 加NMF包

ranks <- 2:10
estim.coad <- nmf(data,ranks, nrun=50)
duplicated(colnames(data))
plot(estim.coad)

#Estimation of the rank: Quality measures computed from 10 runs for each value of r.
pdf("NMF/estim.coad_20240203.pdf")
plot(estim.coad)
dev.off()

#再次NMF,rank=4
seed = 20231120
rank = 3
nmf.rank <- nmf(data, 
                 rank = rank, 
                 nrun=50,
                 seed = seed, 
                 method = "brunet")


#设置颜色
jco <- c("#2874C5","#EABF00","#C6524A","#868686")
index <- extractFeatures(nmf.rank,"max") 
sig.order <- unlist(index)
nmf.Exp.rank <- data[sig.order,]
nmf.Exp.rank <- na.omit(nmf.Exp.rank) #sig.order有时候会有缺失值
group <- predict(nmf.rank) # 提出亚型
table(group)
consensusmap(nmf.rank,
             labRow = NA,
             labCol = NA,
             annCol = data.frame("cluster"=group[colnames(nmf.rank)]),
             annColors = list(cluster=c("1"=jco[1],"2"=jco[2],"3"=jco[3],"4"=jco[4])))

consensusmap(nmf.rank)
coefmap(nmf.rank)


####### Save group information ######
group <- predict(nmf.rank)
group <- as.data.frame(group)
group$group <- paste0('Cluster',group$group)
group$sample <- rownames(group)
group<- group[order(group$group),]
table(group$group)
head(group)
#保存分组情况
write.csv(group,paste0("NMF/group_",rank,".csv")) #也可以以csv格式导出到本地，查看具体分类情况




############ plot with pheatmap ####################
library(pheatmap)
library(tidyr)

### read cluster
rank = 3
file_path = paste0("G:\\TLSdata\\TLS_total_data\\Bulk_RNA\\NMF\\group_",rank,".csv")
cluster = read.csv(file_path,header = TRUE,row.names=1)

cluster <- select(cluster, -last_col())


data_df <- as.data.frame(t(data))

cluster$ID <- rownames(cluster)
data_df$ID <- rownames(data_df)


merged_df <- merge(cluster, data_df, by = "ID", all.x = TRUE)
merged_df <- merged_df[order(merged_df$group), ]

### make cluster_new
cluster_new <- data.frame(merged_df[,"group"],row.names = merged_df[,"ID"])
names(cluster_new) <- 'group'

### make matrix
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




####### Plot ######

### re-order cell type ######
# 定义期望的行名排序
desired_rownames <- c("B_IGHM_04",
                      "B_CD74_1_00",
                      "B_CD69_1_01",
                      "B_MHC_II_05",
                      "B_GC_11",
                      "B_ISG_08",
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
                      "CD4_T_reg_naive_07",
                      "CD4_T_reg_1_04",
                      "CD4_T_reg_2_05",
                      "CD4_T_reg_3_10",
                      "CD4_T_reg_exhausted_01_02",
                      "CD4_T_reg_exhausted_02_09",
                      "CD4_T_reg_ISG",
                      "NKT_1_00",
                      "NKT_2_02",
                      "NK_cytotoxic_04",
                      "NK_reg_CRTAM_03",
                      "NK_reg_KRT86_01",
                      "ILC_06",
                      "dg_T_05",
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
                      "Neutrophil_1_01",
                      'Neutrophil_1_03',
                      'Mast_cell_00'
                      )

# 根据desired_rownames对矩阵行进行排序
sorted_indices <- match(desired_rownames, rownames(merged_df_matrix))
sorted_matrix <- merged_df_matrix[sorted_indices, ]







####### add TLS imprint ############
TLS_imprint = read.csv("G:\\TLSdata\\TLS_total_data\\Bulk_RNA\\TLS_imprint\\TCGA_HNSCC_TLS_imprint_aucell.csv")

TLS_imprint = t(as.matrix(TLS_imprint))
colnames(TLS_imprint) = TLS_imprint[1,]



# 确定匹配的列名
matching_colnames <- intersect(colnames(sorted_matrix), colnames(TLS_imprint))

# 从matrix2中选择一个行向量，只包含与matrix1匹配的列
# 这里我们选择matrix2的第一行作为例子
new_row <- TLS_imprint[2, matching_colnames]

# 将该行向量的顺序调整为与matrix1匹配的列名顺序
# 此步骤是可选的，因为在上一步中已经确保了顺序
# new_row <- new_row[match(colnames(matrix1), names(new_row))]

# 将这个行向量添加到matrix1作为新的一行
sorted_matrix <- rbind(sorted_matrix, as.numeric(new_row))






library(RColorBrewer)
colors33 <- colorRampPalette(c("white","#ff9517","#b83835",'#6d4490',"#000000"))(500) ###color

colorRdBu <- colorRampPalette(c("#f5f8fa","#e1e8ed","#aab8c2","#657786","#14171a"))(500) ###color
colorRdBu <- colorRampPalette(c("white","#f57814","#e2001a","#b41428","black"))(500) ###color
colorRdBu <- colorRampPalette(c("aliceblue","antiquewhite","antiquewhite3","antiquewhite4"))(500) ###color
colorRdBu <- colorRampPalette(c("white","bisque3","bisque3","antiquewhite4","black"))(500) ###color



breaks <- seq(0, 5, length.out = length(colors33) + 1)

plot.scale = pheatmap(sorted_matrix, 
                      scale = "row",
                      cluster_rows = FALSE, 
                      cluster_cols = FALSE,
                      border_color=NA, 
                      color=colors33,
                      breaks = breaks,
                      #cellheight= 10,
                      #gaps_col = c(334,105,96,11),
                      
                      annotation_col = cluster_new,
                      #cutree_rows = 4
                      show_colnames = FALSE
)


pdf("NMF/cluster.pdf") # Create a new pdf device
print(plot.scale)
dev.off() # Close the pdf device



breaks <- seq(0, 5, length.out = length(colorRdBu) + 1)

plot.TLS.imprint = pheatmap(sorted_matrix, 
                      scale = "row",
                      cluster_rows = FALSE, 
                      cluster_cols = FALSE,
                      border_color=NA, 
                      color=colorRdBu,
                      breaks = breaks,
                      #cellheight= 10,
                      #gaps_col = c(334,105,96,11),
                      
                      annotation_col = cluster_new,
                      #cutree_rows = 4
                      show_colnames = FALSE
)


pdf("NMF/TLS_imprint.pdf") # Create a new pdf device
print(plot.TLS.imprint)
dev.off() # Close the pdf device






