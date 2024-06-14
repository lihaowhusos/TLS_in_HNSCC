rm(list=ls())


library(tidyverse)
library(stringr)
library(tidyr)
library(reshape2)
library(forcats)
library(ggpubr)



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

data$group <-factor(data$group,ordered=TRUE,levels=c("mature","immature","none")) #修改因子水平 

colors = c("#E53238", "#FDCA30", "#0064D2")

# total clone size box plot
p_total = ggplot(data, aes(x = group, y = total_num)) + 
  geom_boxplot(outlier.shape = NA,width = 0.25) +
  geom_jitter(aes(color=group),width = 0.25,size = 1) +
  coord_flip() +
  scale_color_manual(values = colors)+
  theme_classic()+
  theme(
    legend.position = "bottom",
    axis.text.y = element_blank(), # 不显示y轴标签
  )


colnames(data)

#### p_B_IGHM_04
library(dplyr)
df_filtered <- data %>%
  filter(B_IGHM_04 != 0)

# Now use df_filtered in your ggplot call
p_B_IGHM_04 = ggplot(df_filtered, aes(x = group, y = B_IGHM_04)) + 
  geom_boxplot(outlier.shape = NA,width = 0.25) +
  geom_jitter(aes(color=group),width = 0.25,size = 1) +
  coord_flip() +
  scale_color_manual(values = colors)+
  theme_classic()+
  theme(
    legend.position = "bottom",
    axis.text.y = element_blank(), # 不显示y轴标签
  )


#### p_B_CD74_1_00
library(dplyr)
df_filtered <- data %>%
  filter(B_CD74_1_00 != 0)

# Now use df_filtered in your ggplot call
p_B_CD74_1_00 = ggplot(df_filtered, aes(x = group, y = B_CD74_1_00)) + 
  geom_boxplot(outlier.shape = NA,width = 0.25) +
  geom_jitter(aes(color=group),width = 0.25,size = 1) +
  coord_flip() +
  scale_color_manual(values = colors)+
  theme_classic()+
  theme(
    legend.position = "bottom",
    axis.text.y = element_blank(), # 不显示y轴标签
  )



#### B_CD69_1_01
library(dplyr)
df_filtered <- data %>%
  filter(B_CD69_1_01 != 0)

# Now use df_filtered in your ggplot call
p_B_CD69_1_01 = ggplot(df_filtered, aes(x = group, y = B_CD69_1_01)) + 
  geom_boxplot(outlier.shape = NA,width = 0.25) +
  geom_jitter(aes(color=group),width = 0.25,size = 1) +
  coord_flip() +
  scale_color_manual(values = colors)+
  theme_classic()+
  theme(
    legend.position = "bottom",
    axis.text.y = element_blank(), # 不显示y轴标签
  )


#### B_MHC_II_05
library(dplyr)
df_filtered <- data %>%
  filter(B_MHC_II_05 != 0)

# Now use df_filtered in your ggplot call
p_B_MHC_II_05 = ggplot(df_filtered, aes(x = group, y = B_MHC_II_05)) + 
  geom_boxplot(outlier.shape = NA,width = 0.25) +
  geom_jitter(aes(color=group),width = 0.25,size = 1) +
  coord_flip() +
  scale_color_manual(values = colors)+
  theme_classic()+
  theme(
    legend.position = "bottom",
    axis.text.y = element_blank(), # 不显示y轴标签
  )



#### B_ISG_08
library(dplyr)
df_filtered <- data %>%
  filter(B_ISG_08 != 0)

# Now use df_filtered in your ggplot call
p_B_ISG_08 = ggplot(df_filtered, aes(x = group, y = B_ISG_08)) + 
  geom_boxplot(outlier.shape = NA,width = 0.25) +
  geom_jitter(aes(color=group),width = 0.25,size = 1) +
  coord_flip() +
  scale_color_manual(values = colors)+
  theme_classic()+
  theme(
    legend.position = "bottom",
    axis.text.y = element_blank(), # 不显示y轴标签
  )


####del 0 value
library(dplyr)
df_filtered <- data %>%
  filter(B_GC_11 != 0)

# Now use df_filtered in your ggplot call
p_B_GC_11 = ggplot(df_filtered, aes(x = group, y = B_GC_11)) + 
  geom_boxplot(outlier.shape = NA,width = 0.25) +
  geom_jitter(aes(color=group),width = 0.25,size = 1) +
  coord_flip() +
  scale_color_manual(values = colors)+
  theme_classic()+
  theme(
    legend.position = "bottom",
    axis.text.y = element_blank(), # 不显示y轴标签
  )


#### Cycling_B_plamsa_cell_09
library(dplyr)
df_filtered <- data %>%
  filter(Cycling_B_plamsa_cell_09 != 0)

# Now use df_filtered in your ggplot call
p_Cycling_B_plamsa_cell_09 = ggplot(df_filtered, aes(x = group, y = Cycling_B_plamsa_cell_09)) + 
  geom_boxplot(outlier.shape = NA,width = 0.25) +
  geom_jitter(aes(color=group),width = 0.25,size = 1) +
  coord_flip() +
  scale_color_manual(values = colors)+
  theme_classic()+
  theme(
    legend.position = "bottom",
    axis.text.y = element_blank(), # 不显示y轴标签
  )

#### Plasma_cell_02_03
library(dplyr)
df_filtered <- data %>%
  filter(Plasma_cell_02_03 != 0)

# Now use df_filtered in your ggplot call
p_Plasma_cell_02_03 = ggplot(df_filtered, aes(x = group, y = Plasma_cell_02_03)) + 
  geom_boxplot(outlier.shape = NA,width = 0.25) +
  geom_jitter(aes(color=group),width = 0.25,size = 1) +
  coord_flip() +
  scale_color_manual(values = colors)+
  theme_classic()+
  theme(
    legend.position = "bottom",
    axis.text.y = element_blank(), # 不显示y轴标签
  )


#### Plasma_cell_01_02
library(dplyr)
df_filtered <- data %>%
  filter(Plasma_cell_01_02 != 0)

# Now use df_filtered in your ggplot call
p_Plasma_cell_01_02 = ggplot(df_filtered, aes(x = group, y = Plasma_cell_01_02)) + 
  geom_boxplot(outlier.shape = NA,width = 0.25) +
  geom_jitter(aes(color=group),width = 0.25,size = 1) +
  coord_flip() +
  scale_color_manual(values = colors)+
  theme_classic()+
  theme(
    legend.position = "bottom",
    axis.text.y = element_blank(), # 不显示y轴标签
  )




p_all = ggarrange(p_total, p_B_IGHM_04, p_B_CD74_1_00, p_B_CD69_1_01, p_B_MHC_II_05, p_B_GC_11, p_B_ISG_08, p_Cycling_B_plamsa_cell_09, p_Plasma_cell_02_03, p_Plasma_cell_01_02, nrow = 1) 

ggsave("G:/TLSdata/TLS_total_data/figures/B/B_clone_size.pdf",plot = p_all,  width = 15, height = 4)


