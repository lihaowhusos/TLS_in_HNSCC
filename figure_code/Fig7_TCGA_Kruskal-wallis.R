#### Kruskal-Wallis #########
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ggridges)
library(data.table)
library(tidytext)
library(cowplot)
library(ggthemes)
library(ggrepel)
library(corrplot)
library(viridisLite)
library(grid)
library(RColorBrewer)
library(lme4)
library(emmeans)
library(openxlsx)
##Collated matrix format##

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


########
standardize <- function(x) {
  rowmean <- apply(x, 1, mean)
  rowsd <- apply(x, 1, sd)  
  rv <- sweep(x, 1, rowmean,"-")  #????��-??ֵ
  rv <- sweep(rv, 1, rowsd, "/")  #?ٳ??Ա?׼??
  return(rv)
}
########
normalize <- function(x) {
  rowmin <- apply(x, 1, min)
  rowmax <- apply(x, 1, max)
  rowmax.min <- rowmax-rowmin
  rv <- sweep(x, 1, rowmin,"-")  #????��-??ֵ
  rv <- sweep(rv, 1, rowmax.min, "/")  #?ٳ??Ա?׼??
  return(rv)
}


data = normalize(data)
data.standardize = standardize(data)

data.standardize = t(data.standardize)
data.standardize = as.data.frame(data.standardize)

colnames(data.standardize) = substr(colnames(data.standardize), 24, nchar(colnames(data.standardize)))

# Convert row names to a column
data.standardize$sample <- rownames(data.standardize)

library(tidyr)

### read cluster
rank = 3
file_path = paste0("G:\\TLSdata\\TLS_total_data\\Bulk_RNA\\NMF\\group_",rank,".csv")
cluster_type = read.csv(file_path,header = TRUE,row.names=1)






####### add TLS imprint ############
TLS_imprint = read.csv("G:\\TLSdata\\TLS_total_data\\Bulk_RNA\\TLS_imprint\\TCGA_HNSCC_TLS_imprint_aucell.csv")
colnames(TLS_imprint) = c("sample","TLS")
TLS_imprint$TLS = as.numeric(TLS_imprint$TLS)
data.standardize <- merge(data.standardize, TLS_imprint, by = "sample")


# Merge the data frames based on the new RowNames column
merged_df <- merge(data.standardize, cluster_type, by = "sample")
merged_df$sample <- NULL


merged_df$group = as.numeric(gsub("Cluster", "", merged_df$group))



########### calculate kruskal_results ########
#cell_types <- colnames(merged_df)
cell_types <- colnames(merged_df)
cell_types <- cell_types[-length(cell_types)]

kruskal_results <- list()
for(i in cell_types) {
  # Create the formula with the current column name
  formula = as.formula(paste(i, "~ group"))
  
  # Run the Kruskal-Wallis test with the constructed formula
  kruskal_results[[i]] <- kruskal.test(formula, data = merged_df)
  
  # Print the results with the current cell type name
  print(paste("Results for", i))
  print(kruskal_results[[i]])
}

# Assuming that you have already run your tests and the results are in kruskal_results,
# Extracting p-values from each kruskal_result
p_values <- sapply(kruskal_results, function(x) x$p.value)

# Adjusting the p-values for multiple testing
adjusted_p_values <- p.adjust(p_values, method = "BH")  # Or your method of choice

# If you need to keep track of which p-value corresponds to which test:
names(adjusted_p_values) <- cell_types

# Check the adjusted p-values
print(adjusted_p_values)
adjusted_p_values


results_df <- data.frame(
  Cell_Type = names(adjusted_p_values),
  Adjusted_P_Value = adjusted_p_values
)

# Save the dataframe to a CSV file
write.csv(results_df, file = "adjusted_p_values.csv", row.names = FALSE)
