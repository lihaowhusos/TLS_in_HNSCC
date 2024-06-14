library(tidyverse)
library(stringr)
library(tidyr)
library(reshape2)
library(forcats)
library(ggplot2)
library(corrgram)
library(corrplot)
library(GGally)
library(ggplot2)


setwd("G:\\TLSdata\\TLS_total_data\\multiplex_data")
set.seed(20240207) # 确保示例可重复


rt=read.csv("cell_type_for_correlation_cluster.csv",sep=",",header=T,check.names=F,row.names=1)

rt$TLS = as.factor(rt$TLS)

colnames(rt)

group_pairs <- combn(levels(rt$TLS), 2, simplify = FALSE)

# 准备结果的容器
results <- list()

# 对每一个变量和每对组进行Mann-Whitney U检验
for(variable in setdiff(names(rt), "TLS")) {
  for(pair in group_pairs) {
    group1 <- rt %>% filter(TLS == pair[1]) %>% pull(variable)
    group2 <- rt %>% filter(TLS == pair[2]) %>% pull(variable)
    
    # Mann-Whitney U检验
    test <- wilcox.test(group1, group2, exact = FALSE, alternative = "two.sided")
    
    # 保存结果
    results[[paste(variable, paste(pair, collapse=" vs "), sep=": ")]] <- c(PValue=test$p.value)
  }
}

# 输出结果前不调整
print("Results before adjustment:")
print(results)