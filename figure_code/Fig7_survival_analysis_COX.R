# 删除无关数据
rm(list=ls())

####combine clinical and sequencing data
library("survival")
library("survminer")
library("readr")
library(dplyr)
library(stringr)

setwd('G:/TLSdata/TLS_total_data/Bulk_RNA/')


#clinical  <- read_table("TCGA临床数据.txt")
#clinical <- clinical[clinical$Id != 1, ]

#处理clinical2
clinical2  <- read.csv("G:/TLSdata/TLS_total_data/Bulk_RNA/TCGA_HNSCC/clinical_data/survival_HNSC_survival.csv")
names(clinical2)[2] <- "Id"
clinical2 <- clinical2 %>% distinct(Id, .keep_all = TRUE)

length(unique(clinical2$Id))


#cell_type <- read.csv("proportions_NuSVR.csv")
#cell_type <- read.csv("Bulk_RNA_cell2location.csv")

cell_type <- read.csv("Bulk_RNA_cell2location_20240202.csv")
names(cell_type)[1] <- "Id"
#cell_type2 <- read.csv("cell_type_xcell.csv")
cell_type <- cell_type %>% 
  filter(str_detect(Id, "01.{1}$"))
length(unique(cell_type$Id))


cell_type$Id <- substr(cell_type$Id, 1, nchar(cell_type$Id)-4)
length(unique(cell_type$Id))

#merged_df <- merge(clinical2, cell_type, by = "Id", all.x = TRUE)
#length(unique(merged_df$Id))


#names(cell_type2)[1] <- "Id"
#cell_type2$Id <- substr(cell_type2$Id, 1, nchar(cell_type2$Id)-4)
#merged_df <- merge(merged_df, cell_type2, by = "Id", all.x = TRUE)

#print(colnames(merged_df))


## plot annotation from cell2location and xcell
#plot(merged_df$q05cell_abundance_w_sf_B_GC_11, merged_df$B.cells, pch=16)



# 假设 "Column_A" 是你要处理的列
# 使用 na.omit() 函数
rt <- inner_join(clinical2, cell_type, by = "Id")
length(unique(rt$Id))

rt$OS <- as.numeric(rt$OS)
rt$OS.time <- as.numeric(rt$OS.time)
rt$OS.time <- rt$OS.time /30    #如果以月为单位，除以30；以年为单位，除以365



rt$DSS <- as.numeric(rt$DSS)
rt$DSS.time <- as.numeric(rt$DSS.time)
rt$DSS.time <- rt$DSS.time /30    #如果以月为单位，除以30；以年为单位，除以365


rt$PFI <- as.numeric(rt$PFI)
rt$PFI.time <- as.numeric(rt$PFI.time)
rt$PFI.time <- rt$PFI.time /30    #如果以月为单位，除以30；以年为单位，除以365





#write.csv(rt, 'cell2location_survival_data.csv', row.names = FALSE)

print(colnames(rt))

cell_choose = 'q05cell_abundance_w_sf_CD4_T_early_exhausted_06'
a <- rt[,cell_choose] <= median(rt[,cell_choose], na.rm = TRUE)

diff=survdiff(Surv(OS.time, OS) ~a,data = rt)
pValue=1-pchisq(diff$chisq,df=1)
pValue=round(pValue,5)
fit <- survfit(Surv(OS.time, OS) ~ a, data = rt)




#summary(fit)    #查看五年生存率

p = ggsurvplot(fit, # 创建的拟合对象
           conf.int = FALSE, # 显示置信区间
           pval = TRUE, # 添加P值
           surv.median.line = "hv",  # 添加中位生存时间线
           risk.table = TRUE, # 添加风险表
           xlab = "Time (month)", # 指定x轴标签
           legend = c(0.8,0.75), # 指定图例位置
           legend.title = "", # 设置图例标题
           legend.labs = c("high", "low"), # 指定图例分组标签
           palette = c("#ff0000", "#0000c0")
           ) # 设置x轴刻度间距

file_path = paste0("G:/TLSdata/TLS_total_data/Bulk_RNA/survival_analysis_figure/",cell_choose,".pdf")
pdf(file_path)
print(p, newpage = FALSE)
dev.off()

############### 循环输出 OS ##############
cell_types = colnames(rt)
cell_types <- cell_types[-c(1:12)]
cell_types

for (cell_choose in cell_types){
  cell_choose
  a <- rt[,cell_choose] <= median(rt[,cell_choose], na.rm = TRUE)
  
  diff=survdiff(Surv(OS.time, OS) ~a,data = rt)
  pValue=1-pchisq(diff$chisq,df=1)
  pValue=round(pValue,5)
  fit <- survfit(Surv(OS.time, OS) ~ a, data = rt)
  p = ggsurvplot(fit, # 创建的拟合对象
                 conf.int = FALSE, # 显示置信区间
                 pval = TRUE, # 添加P值
                 surv.median.line = "hv",  # 添加中位生存时间线
                 risk.table = TRUE, # 添加风险表
                 xlab = "Time (month)", # 指定x轴标签
                 legend = c(0.8,0.75), # 指定图例位置
                 title = cell_choose,
                 #tables.height = ,# numeric value (in [0 - 1])
                 legend.title = "", # 设置图例标题
                 legend.labs = c("high", "low"), # 指定图例分组标签
                 palette = c("#ff0000", "#0000c0")
  ) # 设置x轴刻度间距
  
  file_path = paste0("G:/TLSdata/TLS_total_data/Bulk_RNA/survival_analysis_OS/",cell_choose,".pdf")
  pdf(file_path)
  print(p, newpage = FALSE)
  dev.off()
}


############### 循环输出 DSS ##############
cell_types = colnames(rt)
cell_types <- cell_types[-c(1:12)]
cell_types

for (cell_choose in cell_types){
  cell_choose
  a <- rt[,cell_choose] <= median(rt[,cell_choose], na.rm = TRUE)
  
  diff=survdiff(Surv(DSS.time, DSS) ~a,data = rt)
  pValue=1-pchisq(diff$chisq,df=1)
  pValue=round(pValue,5)
  fit <- survfit(Surv(DSS.time, DSS) ~ a, data = rt)
  p = ggsurvplot(fit, # 创建的拟合对象
                 conf.int = FALSE, # 显示置信区间
                 pval = TRUE, # 添加P值
                 surv.median.line = "hv",  # 添加中位生存时间线
                 risk.table = TRUE, # 添加风险表
                 xlab = "Time (month)", # 指定x轴标签
                 legend = c(0.8,0.75), # 指定图例位置
                 title = cell_choose,
                 #tables.height = ,# numeric value (in [0 - 1])
                 legend.title = "", # 设置图例标题
                 legend.labs = c("high", "low"), # 指定图例分组标签
                 palette = c("#ff0000", "#0000c0")
  ) # 设置x轴刻度间距
  
  file_path = paste0("G:/TLSdata/TLS_total_data/Bulk_RNA/survival_analysis_DSS/",cell_choose,".pdf")
  pdf(file_path)
  print(p, newpage = FALSE)
  dev.off()
}


############### 循环输出 PFI ##############
cell_types = colnames(rt)
cell_types <- cell_types[-c(1:12)]
cell_types

for (cell_choose in cell_types){
  cell_choose
  a <- rt[,cell_choose] <= median(rt[,cell_choose], na.rm = TRUE)
  
  diff=survdiff(Surv(PFI.time, PFI) ~a,data = rt)
  pValue=1-pchisq(diff$chisq,df=1)
  pValue=round(pValue,5)
  fit <- survfit(Surv(PFI.time, PFI) ~ a, data = rt)
  p = ggsurvplot(fit, # 创建的拟合对象
                 conf.int = FALSE, # 显示置信区间
                 pval = TRUE, # 添加P值
                 surv.median.line = "hv",  # 添加中位生存时间线
                 risk.table = TRUE, # 添加风险表
                 xlab = "Time (month)", # 指定x轴标签
                 legend = c(0.8,0.75), # 指定图例位置
                 title = cell_choose,
                 #tables.height = ,# numeric value (in [0 - 1])
                 legend.title = "", # 设置图例标题
                 legend.labs = c("high", "low"), # 指定图例分组标签
                 palette = c("#ff0000", "#0000c0")
  ) # 设置x轴刻度间距
  
  file_path = paste0("G:/TLSdata/TLS_total_data/Bulk_RNA/survival_analysis_PFI/",cell_choose,".pdf")
  pdf(file_path)
  print(p, newpage = FALSE)
  dev.off()
}







##################################### for immune cluster #######################################################
library("survival")
library("survminer")
library("readr")

setwd('G:/TLSdata/TLS_total_data/Bulk_RNA/')

#处理clinical2
clinical2  <- read.csv("G:/TLSdata/TLS_total_data/Bulk_RNA/TCGA_HNSCC/clinical_data/survival_HNSC_survival.csv")
names(clinical2)[2] <- "Id"
clinical2 <- clinical2 %>% distinct(Id, .keep_all = TRUE)


###read cluster
cell_type <- read.csv("NMF/group_3.csv")
names(cell_type)[1] <- "Id"
#cell_type2 <- read.csv("cell_type_xcell.csv")
cell_type <- cell_type %>% 
  filter(str_detect(Id, "01.{1}$"))
length(unique(cell_type$Id))





###read cluster
iterations = 3
file_path = paste0("NMF/group_",iterations,".csv")
cluster_type = read.csv(file_path,header = T,row.names=1)
colnames(cluster_type) = c("cluster_type","Id")

cluster_type$Id <- substr(cluster_type$Id, 1, nchar(cluster_type$Id)-4)

merged_df <- merge(clinical2, cluster_type, by = "Id", all.x = TRUE)

merged_df <- merged_df[!is.na(merged_df$cluster_type), ]


merged_df$cluster_type <- substr(merged_df$cluster_type, 8, nchar(merged_df$cluster_type))




# 假设 "Column_A" 是你要处理的列
# 使用 na.omit() 函数
rt <- merged_df[!is.na(merged_df$cluster_type), ]

# OS  DSS  PFI
survival_type = 'OS'
survival_time = paste0(survival_type, '.time')

rt$futime <- as.numeric(rt[,survival_time])
rt$futime <- rt$futime /30    #如果以月为单位，除以30；以年为单位，除以365
rt$fustat <- as.numeric(rt[,survival_type])

#write.csv(rt, 'cell2location_survival_data.csv', row.names = FALSE)

print(colnames(rt))

fit <- survfit(Surv(futime, fustat) ~ cluster_type, data = rt)

#summary(fit)    #查看五年生存率




p = ggsurvplot(fit, # 创建的拟合对象
               data = rt,
               conf.int = FALSE, # 显示置信区间
               pval = TRUE, # 添加P值
               surv.median.line = "hv",  # 添加中位生存时间线
               risk.table = TRUE, # 添加风险表
               xlab = "Time (month)", # 指定x轴标签
               title = "cluster",
               legend = c(0.8,0.75), # 指定图例位置
               legend.title = "", # 设置图例标题
               #legend.labs = c("n","im","m"), # 指定图例分组标签
               #palette = c("#0064d2","#fdca30", "#e53238")
) 

p = ggsurvplot(fit, # 创建的拟合对象
               data = rt,
               conf.int = FALSE, # 显示置信区间
               pval = TRUE, # 添加P值
               surv.median.line = "hv",  # 添加中位生存时间线
               risk.table = TRUE, # 添加风险表
               xlab = "Time (month)", # 指定x轴标签
               title = "cluster",
               legend = c(0.8,0.75), # 指定图例位置
               legend.title = "", # 设置图例标题
               legend.labs = c("m","n","im"), # 指定图例分组标签
               palette = c("#e53238","#0064d2", "#fdca30")

) 

p

# write.csv(merged_df, file = "NMF/survival_data.csv")


file_path = paste0("NMF/survival_cluster.pdf")
pdf(file_path)
print(p, newpage = FALSE)
dev.off()




##################################### for TMA #######################################################
library("survival")
library("survminer")
library("readr")

setwd('G:/TLSdata/TLS_total_data/multiplex_data/')

rt  <- read.csv("TMA_survival.csv")
colnames(rt)


#rt$death <- as.numeric(rt$fustat)
#rt$suvival_time <- as.numeric(rt$futime)

fit <- survfit(Surv(suvival_time, death) ~ TLS, data = rt)

#summary(fit)    #查看五年生存率




p = ggsurvplot(fit, # 创建的拟合对象
               data = rt,
               conf.int = FALSE, # 显示置信区间
               pval = TRUE, # 添加P值
               surv.median.line = "hv",  # 添加中位生存时间线
               risk.table = TRUE, # 添加风险表
               xlab = "Time (month)", # 指定x轴标签
               title = "cluster",
               legend = c(0.8,0.75), # 指定图例位置
               legend.title = "", # 设置图例标题
               #legend.labs = c("n","im","m"), # 指定图例分组标签
               palette = c("#0064d2","#fdca30", "#e53238")
) 
p


file_path = paste0("G:/TLSdata/TLS_total_data/multiplex_data/survival_TMA_TLS_state.pdf")
pdf(file_path)
print(p, newpage = FALSE)
dev.off()






################# Multivariate Cox regression analysis ###########
library("survival")
library("survminer")
library("forestmodel")
#install.packages("finalfit")
library("finalfit")

setwd('G:/TLSdata/TLS_total_data/Bulk_RNA/')

#处理clinical2
clinical2  <- read.csv("G:/TLSdata/TLS_total_data/Bulk_RNA/TCGA_HNSCC/clinical_data/survival_HNSC_survival.csv")
names(clinical2)[2] <- "Id"
clinical2 <- clinical2 %>% distinct(Id, .keep_all = TRUE)

clinical_matrix <- read.csv("G:/TLSdata/TLS_total_data/Bulk_RNA/TCGA_HNSCC/clinical_data/TCGA.HNSC.sampleMap_HNSC_clinicalMatrix.csv")


###read cluster
cell_type <- read.csv("NMF/group_3.csv")
names(cell_type)[1] <- "Id"
#cell_type2 <- read.csv("cell_type_xcell.csv")
cell_type <- cell_type %>% 
  filter(str_detect(Id, "01.{1}$"))
length(unique(cell_type$Id))


#### read cluster
iterations = 3
file_path = paste0("NMF/group_",iterations,".csv")
cluster_type = read.csv(file_path,header = T,row.names=1)
colnames(cluster_type) = c("cluster_type","Id")

cluster_type$Id <- substr(cluster_type$Id, 1, nchar(cluster_type$Id)-4)

merged_df <- merge(clinical2, cluster_type, by = "Id", all.x = TRUE)

merged_df <- merged_df[!is.na(merged_df$cluster_type), ]


merged_df$cluster_type <- substr(merged_df$cluster_type, 8, nchar(merged_df$cluster_type))




merged_df <- merge(merged_df, clinical_matrix, by = "sample", all.x = TRUE)




# 假设 "Column_A" 是你要处理的列
# 使用 na.omit() 函数
rt <- merged_df[!is.na(merged_df$cluster_type), ]
rt$cluster_type = factor(rt$cluster_type,levels=c("2","3",'1'))
#rt$age_group <- ifelse(rt$age_at_initial_pathologic_diagnosis <= 60, "low", "high")
rt <- rt[!is.na(rt$OS.time), ]
rt$clinical_stage[rt$clinical_stage == ""] <- "None"

rt$cluster_type = factor(rt$cluster_type,levels=c("2","3",'1'))


rt$clinical_stage = factor(rt$clinical_stage,levels=c("Stage I","Stage II",'Stage III', 'Stage IVA', 'Stage IVB', 'Stage IVC','None'))
 
# OS  DSS  PFI
survival_type = 'PFI'
survival_time = paste0(survival_type, '.time')

rt$futime <- as.numeric(rt[,survival_time])
rt$futime <- rt$futime /30    #如果以月为单位，除以30；以年为单位，除以365
rt$fustat <- as.numeric(rt[,survival_type])


valid_rows <- !is.na(rt[[survival_type]])

rt_non_NA <- rt[valid_rows, ]



colnames(rt_non_NA)


res.cox <- coxph(Surv(futime, fustat) ~ age_at_initial_pathologic_diagnosis + gender + clinical_stage + cluster_type, data =  rt_non_NA)
summary(res.cox)


sur.forest = ggforest(res.cox,
          data = rt_non_NA,  #数据集
          main = paste0('Hazard ratio of TCGA_HNSC_',iterations),  #标题
          noDigits = 3
                      
)


pdf(paste0("COX/TCGA_HNSC_COX_",survival_type,".pdf")) # Create a new pdf device
print(sur.forest)
dev.off() # Close the pdf device



############ Universal and multi ############

dependent_os  <- "Surv(OS.time, OS)"
dependent_dss <- "Surv(DSS.time, DSS)"
dependent_pfi <- "Surv(PFI.time, PFI)"
explanatory   <- c("age_at_initial_pathologic_diagnosis", "gender", "clinical_stage", "cluster_type")


sur.univariable.multivariable.OS =   rt %>% 
  finalfit(dependent_os, explanatory, add_dependent_label = FALSE) %>% 
  rename("Overall survival" = label) %>% 
  rename(" " = levels) %>% 
  rename("  " = all)

sur.univariable.multivariable.DSS =   rt %>% 
  finalfit(dependent_dss, explanatory, add_dependent_label = FALSE) %>% 
  rename("DSS" = label) %>% 
  rename(" " = levels) %>% 
  rename("  " = all)

sur.univariable.multivariable.PFI =   rt %>% 
  finalfit(dependent_pfi, explanatory, add_dependent_label = FALSE) %>% 
  rename("PFI" = label) %>% 
  rename(" " = levels) %>% 
  rename("  " = all)



