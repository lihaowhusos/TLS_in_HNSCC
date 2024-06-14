## try http:// if https:// URLs are not supported

#install.packages("ggpubr")
#install.packages("ggthemes")

library(ggpubr)
library(ggthemes)
library(paletteer)
library(dplyr)


setwd("G:\\TLSdata\\TLS_total_data\\table")        

deg.data=read.csv("decoupleR_GO_BP_B.csv",sep=",",header=T,check.names=F)

deg.data$logP = -log10(deg.data$pvals_adj)
deg.data$rank = deg.data[,1]

indices_positive = which(deg.data$meanchange > 0)
indices_negative = which(deg.data$meanchange < 0)

deg.data$logP[indices_positive] = deg.data$logP[indices_positive]
deg.data$logP[indices_negative] = - deg.data$logP[indices_negative]


#[1] "B_CD69_1_01"              "B_CD74_1_00"              "B_GC_11"                  "B_IGHM_04"                "B_ISG_08"                 "B_MHC_II_05"             
#[7] "Cycling_B_plamsa_cell_09" "Plasma_cell_01_02"        "Plasma_cell_02_03"

celltypes = unique(deg.data$group)


names(deg.data)[which(is.na(names(deg.data)))] <- "new_name_for_NA"
names(deg.data)[which(names(deg.data) == "")] <- "new_name_for_empty"

deg.data <- deg.data %>% 
  mutate(logP = ifelse(logP > 300, 300, logP))
  

for(celltype in celltypes){
  deg.data_plot = subset(deg.data, group == celltype)
  deg.data$Label = ""
  up.genes = head(deg.data$Symbol[which(deg.data$Group == "up-regulated")], 20)
  down.genes = head(deg.data$Symbol[which(deg.data$Group == "down-regulated")], 20)
  deg.top20.genes = c(as.character(up.genes), as.character(down.genes))
  deg.data$Label[match(deg.top20.genes, deg.data$Symbol)] <- deg.top20.genes
  
  ggscatter(deg.data_plot, x = "rank", y = "logP", 
            size = 2,
            font.label = 9,
            repel = T,
            xlab = "Rank",
            #          fill = "rank",
            color = "logP",
            ylab = "-log10(P.adj)"
  ) + 
    scale_color_gradient2(midpoint = 0, low = "dodgerblue3", mid = "gray100", high = "#CD2626")+
    theme_base()
  
  # Save pdf
  file_path = paste0("G:\\TLSdata\\TLS_total_data\\figures\\B\\snake_plot\\",celltype,"_snake_plot.pdf")
  ggsave(file_path, height = 10, width = 4)
}


