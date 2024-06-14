## try http:// if https:// URLs are not supported

#install.packages("ggpubr")
#install.packages("ggthemes")

library(ggpubr)
library(ggthemes)


setwd("G:\\TLSdata\\TLS_total_data\\Results")              #???ù???Ŀ¼

deg.data=read.csv("Find_All_Markers_scanpy_B_between_group.csv",sep=",",header=T,check.names=F)
deg.data$Symbol = deg.data$mature_names

deg.data$logP = deg.data$mature_pvals_adj + 1
deg.data$logP = -log10(deg.data$mature_pvals_adj)

deg.data <- deg.data[!(deg.data$mature_logfoldcha > 10 | deg.data$mature_logfoldcha < -10), ]

deg.data$Group = "not-significant"
deg.data$Group[which((deg.data$mature_pvals_adj < 0.05)& (deg.data$mature_logfoldcha > 0.5))] = "up-regulated"
deg.data$Group[which((deg.data$mature_pvals_adj < 0.05)& (deg.data$mature_logfoldcha < -0.5))] = "down-regulated"



deg.data$Label = ""
deg.data = deg.data[order(deg.data$mature_pvals_adj),]

up.genes = head(deg.data$Symbol[which(deg.data$Group == "up-regulated")], 20)
down.genes = head(deg.data$Symbol[which(deg.data$Group == "down-regulated")], 20)

deg.top20.genes = c(as.character(up.genes), as.character(down.genes))
deg.data$Label[match(deg.top20.genes, deg.data$Symbol)] <- deg.top20.genes




ggscatter(deg.data, x = "mature_logfoldcha", y = "logP", 
          color = "Group", 
          palette = c("#00a4e4","#BBBBBB", "#ff0000"),
          size = 2,
          font.label = 9,
          
          repel = T,
          xlab = "log2FC",
          ylab = "-log10(FDR)",) + theme_base() +
  geom_hline(yintercept = 1.30, linetype = "dashed") +
  geom_vline(xintercept = c(-0.5,0.5), linetype = "dashed")+
  geom_text(aes(label = Label),size = 3, vjust = 1,hjust = -0.1)






ggsave("G:\\TLSdata\\TLS_total_data\\figures\\B\\B_vocanol_plot.pdf", height = 9, width = 10)


#############????################
library(ggplot2)

# ??ȡ???ݣ?
data <- read.table("edgerOut.xls",sep="\t",header=T,check.names=F)
data$label <- c(rownames(data)[1:10],rep(NA,(nrow(data)-10)))



ggplot(data,aes(log2FoldChange, -log10(padj)))+
  # ????ˮƽ?ο??ߣ?
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#999999")+
  # ??????ֱ?ο??ߣ?
  geom_vline(xintercept = c(-1.2,1.2), linetype = "dashed", color = "#999999")+
  # ɢ??ͼ:
  geom_point(aes(size=-log10(padj), color= -log10(padj)))+
  # ָ????ɫ????ģʽ??
  scale_color_gradientn(values = seq(0,1,0.2),
                        colors = c("#39489f","#39bbec","#f9ed36","#f38466","#b81f25"))+
  # ָ??ɢ????С????ģʽ??
  scale_size_continuous(range = c(1,3))+
  # ??????????
  theme_bw()+
  # ??????????ͼ??λ?ã?
  theme(panel.grid = element_blank(),
        legend.position = c(0.01,0.7),
        legend.justification = c(0,1)
  )+
  # ???ò???ͼ??????ʾ??
  guides(col = guide_colourbar(title = "-Log10_q-value"),
         size = "none")+
  # ???ӱ?ǩ??
  geom_text(aes(label=label, color = -log10(padj)), size = 3, vjust = 1.5, hjust=1)+
  # ?޸??????᣺
  xlab("Log2FC")+
  ylab("-Log10(FDR q-value)")

# ????ͼƬ??
ggsave("vocanol_plot.pdf", height = 9, width = 10)



