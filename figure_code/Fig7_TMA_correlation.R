library(tidyverse)
library(stringr)
library(tidyr)
library(reshape2)
library(forcats)
library(ggplot2)
library(corrgram)
library(corrplot)
library(GGally)

setwd("G:\\TLSdata\\TLS_total_data\\multiplex_data")
rt=read.csv("cell_type_for_correlation.csv",sep=",",header=T,check.names=F,row.names=1)
cor_cell_type = cor(rt)


library(pheatmap)
library(RColorBrewer)

colours=colorRampPalette(c("white", "#1874CD", "#FFC125","#EE4000","#CD0000"))(1000)



hm = pheatmap(cor_cell_type, color = colours,
              scale = "none", cluster_rows = FALSE,
              cluster_cols = FALSE,border_color = NA,
              number_color = "black",display_numbers = T
              # labels_row = rownames(dataframe_hm)
              )





#hm2 = ggpairs(rt, title="correlogram")+ +theme_classic()

ggsave("hm_total_correlation.pdf",plot = hm,  width = 6.6, height = 6)

#ggsave("hm_total_cell_correlation.pdf",plot = hm2,  width = 6.6, height = 6)



library(ggplot2)
setwd("G:\\TLSdata\\TLS_total_data\\multiplex_data")
rt=read.csv("cell_type_for_correlation_cluster.csv",sep=",",header=T,check.names=F,row.names=1)

colnames(rt)

colors = c( "#FDCA30","#E53238","#0064D2" )

#rt$TLS = as.factor(rt$TLS)


ggplot(rt, aes(x=ratio_CD20, y=ratio_CD4_TCF1,fill=TLS)) +
  geom_point(shape = 21,alpha=1,size = 3) +
  scale_fill_manual(values = colors) +
  scale_color_manual(values = colors)+
  geom_smooth(method=lm , color="black", fill="grey", se=TRUE) +
  theme_classic()

p_im = ggplot(rt, aes(x=X, y=Y, size=clone_id_size, color=cluster_2)) +
  geom_point(alpha=0.25) +
  scale_color_manual(values = colors)


p_im_fill = ggplot(CD4_TCR_table_immature, aes(x=X, y=Y, size=clone_id_size, fill=cluster_2)) +
  geom_point(colour = "black", shape = 21,alpha=1) +
  scale_fill_manual(values = color)+
  scale_size_continuous(
    range =  c(1,9),
    limits = c(0, 40),
    breaks = c(0, 20, 40),
    labels = c(0, 20, 40)
  )+
  theme_bw()


colors = c("#0064D2", "#FDCA30", "#E53238")

scatter_CD4 = ggplot() +
  geom_scatterpie(aes(x = x,y = y,r =r), 
  data = table, 
  cols = c("n","i","m")) +
  geom_scatterpie_legend(table$r, n = 3, x = 2.5, y = 4) + #######x,y is the location of legend
  scale_y_reverse()+
  scale_fill_manual(values = colors)+
  labs(x  = "X=size x 25")+
  theme_classic()








