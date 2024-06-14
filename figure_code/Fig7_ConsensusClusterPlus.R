####ConsensusClusterPlus
##based on xCell
####??װ
#BiocManager::install("ConsensusClusterPlus")
#BiocManager::install("ComplexHeatmap")

#options(BioC_mirror="https://mirrors.tuna.tsinghua.edu.cn/bioconductor")


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


#####??׼####
standardize <- function(x) {
  rowmean <- apply(x, 1, mean)
  rowsd <- apply(x, 1, sd)  
  rv <- sweep(x, 1, rowmean,"-")  #????��-??ֵ
  rv <- sweep(rv, 1, rowsd, "/")  #?ٳ??Ա?׼??
  return(rv)
}
######
normalize <- function(x) {
  rowmin <- apply(x, 1, min)
  rowmax <- apply(x, 1, max)
  rowmax.min <- rowmax-rowmin
  rv <- sweep(x, 1, rowmin,"-")  
  rv <- sweep(rv, 1, rowmax.min, "/")  
  return(rv)
}



d=data   #data
data = normalize(data)
data.standardize = standardize(data)





library(ConsensusClusterPlus)

results = ConsensusClusterPlus(data.standardize,maxK=10,reps=400,pItem=0.8,pFeature=1,clusterAlg="hc",distance="pearson",seed=12621213.6666666,plot="png", writeTable = TRUE)

###########PAC for best K###########
Kvec = 2:10
x1 = 0.1; x2 = 0.9        # threshold defining the intermediate sub-interval
PAC = rep(NA,length(Kvec)) 
names(PAC) = paste("K=",Kvec,sep="")  # from 2 to maxK
for(i in Kvec){
  M = results[[i]]$consensusMatrix
  Fn = ecdf(M[lower.tri(M)])          # M 为计算出共识矩阵
  PAC[i-1] = Fn(x2) - Fn(x1)
} 

optK = Kvec[which.min(PAC)]  # 理想的K值





########################## plot with pheatmap ############################
library(pheatmap)

###read cluster
iterations = 3
file_path = paste0("G:\\TLSdata\\TLS_total_data\\Bulk_RNA\\untitled_consensus_cluster\\untitled_consensus_cluster.k=",iterations,".consensusClass.csv")
cluster = read.csv(file_path,header = FALSE,row.names=1)



annotation_col = t(cluster)####annotation????????Ϊfactor
rownames(annotation_col) = "Group"

data <- rbind(data, `group` = annotation_col)


last_row <- tail(data, 1)
ordered_columns <- order(last_row)

data.paixu <- data[, ordered_columns]# 假设df是你的数据框变量
# 假设df是你的DataFrame
rownames(data.paixu) <- substr(rownames(data.paixu), 24, nchar(rownames(data.paixu)))

table(tail(data, 1))



data.clear <- data.paixu[-nrow(data.paixu), ]
data.clear <- t(apply(data.clear, 1, scale))

last_row <- tail(data.paixu, 1)


#annotation_col <- tail(data.paixu, n = 1)
#colnames(annotation_col) = "Group"


#annotation_col = data.frame(factor(tail(data, n = 1)))####annotation????????Ϊfactor
#colnames(annotation_col) = "Group"



#data.normalize = normalize(data)
#data.meannormalize = meannormalize(data)
#data.stan = standardize(data)
#######??һ??????######

colors33 <- colorRampPalette(c("white","#ff9517","#b83835",'#6d4490',"#000000"))(500) ###color

breaks <- seq(0, 5, length.out = length(colors33) + 1)

plot = pheatmap(data.paixu, 
                scale = "row",
                cluster_rows = FALSE, 
                cluster_cols = FALSE,
                border_color=NA, 
                color=colors33,
                breaks = breaks,
                show_colnames = FALSE,
                #cellheight= 10,
                #gaps_col = c(334,105,96,11),
                
                #annotation_col = annotation_col,
                #cutree_rows = 4
)


pdf("cluster3.pdf") # Create a new pdf device
print(plot.scale)
dev.off() # Close the pdf device








s