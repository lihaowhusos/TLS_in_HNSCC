library(tidyverse)
library(stringr)
library(tidyr)
library(reshape2)
library(forcats)

data = read.csv('G:/TLSdata/TLS_total_data/table/CD8_clone_id_group.csv',header = TRUE)

data$group <- data$sample
data$group <- recode(data$group, 
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

####################clone type share in CD4 T clusters#################
print(colnames(data))


# 统计cyto列和effector列中数值同时大于0的行的数量

names = c("CD8_T_naive_01",
          "CD8_T_effector_memory_01" ,
          "CD8_T_cytotoxic_02",
          "CD8_T_late_exhausted_01",
          "CD8_T_ISG"
)

number = as.numeric(length(names))


dataframe_hm = data.frame(matrix(ncol = 0, nrow = number))###得第一次有6行
n=c()
for (j in 1:number){
  for (k in 1:number){
    c_name = names[j]
    d_name = names[k]
    n_sub = nrow(subset(data, data[[c_name]] > 0 & data[[d_name]] > 0))
    n = c(n,n_sub)
    
  }
  dataframe_hm = cbind(dataframe_hm, n)
  #  dataframe_hm[j,j] = 0
  n=c()
}

rownames(dataframe_hm) = names
colnames(dataframe_hm) = names

dataframe_hm = log10(dataframe_hm+1)


disp_num = as.matrix(10**dataframe_hm-1)



library(pheatmap)
library(RColorBrewer)

colours=colorRampPalette(c("white", "#1874CD", "#FFC125","#EE4000","#CD0000"))(500)
colours = colorRampPalette(c("white", "#EBE5ED", "#D2D1E7","#AABDD9","#73ABC9","#4A8FBD","#338186","#29685B","#174233"))(500)

#colours=colorRampPalette(c("white", "#1874CD", "#FFC125","#CD0000"))(500)

my_breaks <- seq(1, 3.1, length.out = 500)

CD8_T_clone_share_hm = pheatmap(dataframe_hm, color = colours,
                                scale = "none", cluster_rows = FALSE,
                                cluster_cols = FALSE,border_color = NA,display_numbers = disp_num,	
                                number_color = "black",
                                legend_breaks = c(0,1,2,3),
                                breaks = my_breaks
                                )




ggsave("G:/TLSdata/TLS_total_data/figures/TCR_CD8/CD8_T_clone_share_hm.pdf",plot = CD8_T_clone_share_hm,  width = 4.5, height = 4)





################ 饼状图 ################
#install.packages("scatterpie")
library(ggplot2)
library(scatterpie)

print(colnames(data))

names = c("CD8_T_naive_01",
          "CD8_T_effector_memory_01" ,
          "CD8_T_cytotoxic_02",
          "CD8_T_late_exhausted_01",
          "CD8_T_ISG"
)

number = as.numeric(length(names))

x = c()
y = c()

for (i in 2:number){
  k = i-1
  x = c(x, rep(i,k))
  y = c(y, c(1:k))
}




n = c()
i = c()
m = c()
total = c()
for (j in 2:number){
  s = j-1
  for (k in 1:s){
    c_name = names[j]
    d_name = names[k]
    n_sub = nrow(subset(data, group == "none" & data[[c_name]] > 0 & data[[d_name]] > 0))
    i_sub = nrow(subset(data, group == "immature" & data[[c_name]] > 0 & data[[d_name]] > 0))
    m_sub = nrow(subset(data, group == "mature" & data[[c_name]] > 0 & data[[d_name]] > 0))
    total_sub = nrow(subset(data, data[[c_name]] > 0 & data[[d_name]] > 0))
    n = c(n,n_sub)
    i = c(i,i_sub)
    m = c(m,m_sub)
    total = c(total,total_sub)
  }
}


table = data.frame(x,y,n,i,m,total)

table$r = sqrt(table$total)/50 #盘子大小


colors = c("#0064D2", "#FDCA30", "#E53238")

scatter_CD8 = ggplot() +
  geom_scatterpie(aes(x = x,y = y,r =r), 
                  data = table, 
                  cols = c("n","i","m")) +
  geom_scatterpie_legend(table$r, n = 3, x = 2.5, y = 4) + #######x,y is the location of legend
  scale_y_reverse()+
  scale_fill_manual(values = colors)+
  labs(x  = "X=size x 50")+
  theme_classic()

scatter_CD8

ggsave("G:/TLSdata/TLS_total_data/figures/TCR_CD8/CD8_clone_share_scatter.pdf",plot = scatter_CD8,  width = 4.5, height = 4)




