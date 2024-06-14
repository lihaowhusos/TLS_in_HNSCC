library(tidyverse)
library(stringr)
library(tidyr)
library(reshape2)
library(forcats)
library(ggpubr)



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

data$group <-factor(data$group,ordered=TRUE,levels=c("mature","immature","none")) #修改因子水平 


colors = c("#E53238", "#FDCA30", "#0064D2")

# total clone size box plot
p_total = ggplot(data, aes(x = group, y = total_num)) + 
  geom_boxplot(outlier.shape = NA,width = 0.25) +
  geom_jitter(aes(color=group),width = 0.25,size = 1) +
  coord_flip() +
  scale_color_manual(values = colors)+
  theme_classic()+
  theme(legend.position="bottom")


####del 0 value
library(dplyr)
df_filtered <- data %>%
  filter(CD8_T_ISG != 0)

# Now use df_filtered in your ggplot call
p_CD8_ISG = ggplot(df_filtered, aes(x = group, y = CD8_T_ISG)) + 
  geom_boxplot(outlier.shape = NA,width = 0.25) +
  geom_jitter(aes(color=group),width = 0.25,size = 1) +
  coord_flip() +
  scale_color_manual(values = colors)+
  theme_classic()+
  theme(legend.position="bottom")


####del 0 value
library(dplyr)
df_filtered <- data %>%
  filter(CD8_T_cytotoxic_02 != 0)

# Now use df_filtered in your ggplot call
p_CD8_T_cytotoxic = ggplot(df_filtered, aes(x = group, y = CD8_T_cytotoxic_02)) + 
  geom_boxplot(outlier.shape = NA,width = 0.25) +
  geom_jitter(aes(color=group),width = 0.25,size = 1) +
  coord_flip() +
  scale_color_manual(values = colors)+
  theme_classic()+
  theme(legend.position="bottom")



####del 0 value
library(dplyr)
df_filtered <- data %>%
  filter(CD8_T_effector_memory_01 != 0)

# Now use df_filtered in your ggplot call
p_CD8_T_effector_memory_01= ggplot(df_filtered, aes(x = group, y = CD8_T_effector_memory_01)) + 
  geom_boxplot(outlier.shape = NA,width = 0.25) +
  geom_jitter(aes(color=group),width = 0.25,size = 1) +
  coord_flip() +
  scale_color_manual(values = colors)+
  theme_classic()+
  theme(legend.position="bottom")


####del 0 value
library(dplyr)
df_filtered <- data %>%
  filter(CD8_T_late_exhausted_01 != 0)

# Now use df_filtered in your ggplot call
p_CD8_T_late_exhausted_01= ggplot(df_filtered, aes(x = group, y = CD8_T_late_exhausted_01)) + 
  geom_boxplot(outlier.shape = NA,width = 0.25) +
  geom_jitter(aes(color=group),width = 0.25,size = 1) +
  coord_flip() +
  scale_color_manual(values = colors)+
  theme_classic()+
  theme(legend.position="bottom")



####del 0 value
library(dplyr)
df_filtered <- data %>%
  filter(CD8_T_naive_01 != 0)

# Now use df_filtered in your ggplot call
p_CD8_T_naive_01= ggplot(df_filtered, aes(x = group, y = CD8_T_naive_01)) + 
  geom_boxplot(outlier.shape = NA,width = 0.25) +
  geom_jitter(aes(color=group),width = 0.25,size = 1) +
  coord_flip() +
  scale_color_manual(values = colors)+
  theme_classic()+
  theme(legend.position="bottom")




p_all = ggarrange(p_total, p_CD8_T_effector_memory_01, p_CD8_T_naive_01, p_CD8_T_cytotoxic, p_CD8_T_late_exhausted_01
,nrow = 1) 

ggsave("G:/TLSdata/TLS_total_data/figures/TCR_CD8/CD8_T_clone_size.pdf",plot = p_all,  width = 10, height = 4)





################################### statistic ##########################################
library(dplyr)

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

data$group <-factor(data$group,ordered=TRUE,levels=c("mature","immature","none")) #修改因子水平 


# 使用combn创建不同组合
group_pairs <- combn(unique(data$group), 2)

# 初始化一个列表来存储测试结果
test_results <- list()

################################### statistic ##########################################
library(dplyr)

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

data$group <-factor(data$group,ordered=TRUE,levels=c("mature","immature","none")) #修改因子水平 


# 使用combn创建不同组合
group_pairs <- combn(unique(data$group), 2)

# 初始化一个列表来存储测试结果
test_results <- list()

# 对每对组合运行Mann-Whitney U测试
for(i in seq(ncol(group_pairs))) {
  group1 <- data %>% filter(group == group_pairs[1, i])
  group2 <- data %>% filter(group == group_pairs[2, i])
  
  test_result <- wilcox.test(group1$total_num, group2$total_num, alternative = "less")
  
  # 存储结果
  test_results[[paste(group_pairs[1, i], group_pairs[2, i], sep = "-")]] <- test_result
}
print(test_results)


cluster_name = 'CD8_T_ISG'
# 对每对组合运行Mann-Whitney U测试
for(i in seq(ncol(group_pairs))) {
  group1 <- data %>% filter(group == group_pairs[1, i]) %>% select(cluster_name) %>% filter(!!sym(cluster_name) != 0)
  group2 <- data %>% filter(group == group_pairs[2, i]) %>% select(cluster_name) %>% filter(!!sym(cluster_name) != 0)
  
  test_result <- wilcox.test(group1[,cluster_name], group2[,cluster_name], alternative = "less")
  
  # 存储结果
  test_results[[paste(group_pairs[1, i], group_pairs[2, i], sep = "-")]] <- test_result
}
print(test_results)



cluster_name = 'CD8_T_cytotoxic_02'
# 对每对组合运行Mann-Whitney U测试
for(i in seq(ncol(group_pairs))) {
  group1 <- data %>% filter(group == group_pairs[1, i]) %>% select(cluster_name) %>% filter(!!sym(cluster_name) != 0)
  group2 <- data %>% filter(group == group_pairs[2, i]) %>% select(cluster_name) %>% filter(!!sym(cluster_name) != 0)
  
  test_result <- wilcox.test(group1[,cluster_name], group2[,cluster_name], alternative = "less")
  
  # 存储结果
  test_results[[paste(group_pairs[1, i], group_pairs[2, i], sep = "-")]] <- test_result
}
print(test_results)



cluster_name = 'CD8_T_effector_memory_01'
# 对每对组合运行Mann-Whitney U测试
for(i in seq(ncol(group_pairs))) {
  group1 <- data %>% filter(group == group_pairs[1, i]) %>% select(cluster_name) %>% filter(!!sym(cluster_name) != 0)
  group2 <- data %>% filter(group == group_pairs[2, i]) %>% select(cluster_name) %>% filter(!!sym(cluster_name) != 0)
  
  test_result <- wilcox.test(group1[,cluster_name], group2[,cluster_name], alternative = "less")
  
  # 存储结果
  test_results[[paste(group_pairs[1, i], group_pairs[2, i], sep = "-")]] <- test_result
}
print(test_results)



cluster_name = 'CD8_T_naive_01'
# 对每对组合运行Mann-Whitney U测试
for(i in seq(ncol(group_pairs))) {
  group1 <- data %>% filter(group == group_pairs[1, i]) %>% select(cluster_name) %>% filter(!!sym(cluster_name) != 0)
  group2 <- data %>% filter(group == group_pairs[2, i]) %>% select(cluster_name) %>% filter(!!sym(cluster_name) != 0)
  
  test_result <- wilcox.test(group1[,cluster_name], group2[,cluster_name], alternative = "less")
  
  # 存储结果
  test_results[[paste(group_pairs[1, i], group_pairs[2, i], sep = "-")]] <- test_result
}
print(test_results)





cluster_name = 'CD8_T_late_exhausted_01'
# 对每对组合运行Mann-Whitney U测试
for(i in seq(ncol(group_pairs))) {
  group1 <- data %>% filter(group == group_pairs[1, i]) %>% select(cluster_name) %>% filter(!!sym(cluster_name) != 0)
  group2 <- data %>% filter(group == group_pairs[2, i]) %>% select(cluster_name) %>% filter(!!sym(cluster_name) != 0)
  
  test_result <- wilcox.test(group1[,cluster_name], group2[,cluster_name], alternative = "less")
  
  # 存储结果
  test_results[[paste(group_pairs[1, i], group_pairs[2, i], sep = "-")]] <- test_result
}
print(test_results)


