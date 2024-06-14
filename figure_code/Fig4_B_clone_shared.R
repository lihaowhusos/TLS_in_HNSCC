library(tidyverse)
library(stringr)
library(tidyr)
library(reshape2)
library(forcats)

data = read.csv('G:/TLSdata/TLS_total_data/table/BCR_clone_id_group.csv',header = TRUE)

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



####################clone type share in B clusters#################
print(colnames(data))

# 统计cyto列和effector列中数值同时大于0的行的数量

names = c("B_IGHM_04",
          "B_CD74_1_00" ,
          "B_CD69_1_01",
          "B_MHC_II_05",
          "B_GC_11",
          "B_ISG_08",
          "Cycling_B_plamsa_cell_09",
          "Plasma_cell_02_03",
          "Plasma_cell_01_02")

dataframe_hm = data.frame(matrix(ncol = 0, nrow = 9))###得第一次有9行
n=c()
for (j in 1:9){
  for (k in 1:9){
    c_name = names[j]
    d_name = names[k]
    n_sub = nrow(subset(data, data[[c_name]] > 0 & data[[d_name]] > 0))
    n = c(n,n_sub)
  }
  dataframe_hm = cbind(dataframe_hm, n)
  #dataframe_hm[j,j] = 0
  n=c()
}

rownames(dataframe_hm) = names
colnames(dataframe_hm) = names

dataframe_hm = log10(dataframe_hm+1)


disp_num = as.matrix(10**dataframe_hm-1)



library(pheatmap)
library(RColorBrewer)

colours=colorRampPalette(c("white", "#1874CD", "#FFC125","#EE4000","#CD0000"))(500)

B_clone_share_hm = pheatmap(dataframe_hm, color = colours,
                            scale = "none", cluster_rows = FALSE,
                            cluster_cols = FALSE,border_color = NA,
                            display_numbers = disp_num,	
                            number_color = "black")




ggsave("G:/TLSdata/TLS_total_data/figures/B/B_clone_share_hm.pdf",plot = B_clone_share_hm,  width = 6.5, height = 6)




################ 饼状图 ################
#install.packages("scatterpie")



library(ggplot2)
library(scatterpie)

x = c()
y = c()


for (i in 2:9){
  k = i-1
  x = c(x, rep(i,k))
  y = c(y, c(1:k))
}


names = colnames(data)
names = c("B_IGHM_04",
          "B_CD74_1_00" ,
          "B_CD69_1_01",
          "B_MHC_II_05",
          "B_ISG_08",
          "B_GC_11",
          "Cycling_B_plamsa_cell_09",
          "Plasma_cell_02_03",
          "Plasma_cell_01_02")

n = c()
i = c()
m = c()
total = c()
for (j in 2:9){
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

table$r = sqrt(table$total)/30


colors = c("#0064D2", "#FDCA30", "#E53238")

scatter_B = ggplot() +
  geom_scatterpie(aes(x = x,y = y,r =r), 
                  data = table, 
                  cols = c("n","i","m")) +
  geom_scatterpie_legend(table$r, n = 3, x = 3, y = 7) +
  scale_y_reverse()+
  scale_fill_manual(values = colors)+
  labs(x  = "X=size x 30")+
  theme_classic()

scatter_B




ggsave("G:/TLSdata/TLS_total_data/figures/B/B_clone_share_scatter.pdf",plot = scatter_B,  width = 4.5, height = 4)





