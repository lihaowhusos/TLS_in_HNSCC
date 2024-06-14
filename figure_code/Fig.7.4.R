library(ConsensusClusterPlus)
library(pheatmap)
library(dplyr)
library(stringr)

# Set working directory
setwd("XXX")

# Read and preprocess data
rt <- read.csv("XXX.csv", sep=",", header=TRUE, check.names=FALSE)
names(rt)[1] <- "Id"
rt <- rt %>% filter(str_detect(Id, "01.{1}$"))

# Collate matrix format
rt <- as.matrix(rt)
rownames(rt) <- rt[,1]
exp <- rt[,2:ncol(rt)]
exp <- t(exp)

dimnames <- list(rownames(exp), colnames(exp))
data <- matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)

# Standardize and normalize functions
standardize <- function(x) {
  rowmean <- apply(x, 1, mean)
  rowsd <- apply(x, 1, sd)
  rv <- sweep(x, 1, rowmean, "-")
  rv <- sweep(rv, 1, rowsd, "/")
  return(rv)
}

normalize <- function(x) {
  rowmin <- apply(x, 1, min)
  rowmax <- apply(x, 1, max)
  rowmax.min <- rowmax - rowmin
  rv <- sweep(x, 1, rowmin, "-")
  rv <- sweep(rv, 1, rowmax.min, "/")
  return(rv)
}

data <- normalize(data)
data.standardize <- standardize(data)

# Consensus clustering
results <- ConsensusClusterPlus(data.standardize, maxK=10, reps=400, pItem=0.8, pFeature=1, clusterAlg="hc", distance="pearson", seed=12621213.6666666, plot="png", writeTable=TRUE)

# PAC for best K
Kvec <- 2:10
x1 <- 0.1
x2 <- 0.9
PAC <- sapply(Kvec, function(i) {
  M <- results[[i]]$consensusMatrix
  Fn <- ecdf(M[lower.tri(M)])
  Fn(x2) - Fn(x1)
})
names(PAC) <- paste("K=", Kvec, sep="")
optK <- Kvec[which.min(PAC)]

# Plot with pheatmap
iterations <- 3
file_path <- paste0("XXX", iterations, ".consensusClass.csv")
cluster <- read.csv(file_path, header=FALSE, row.names=1)

annotation_col <- t(cluster)
rownames(annotation_col) <- "Group"
data <- rbind(data, `group` = annotation_col)

last_row <- tail(data, 1)
ordered_columns <- order(last_row)
data.paixu <- data[, ordered_columns]

rownames(data.paixu) <- substr(rownames(data.paixu), 24, nchar(rownames(data.paixu)))

data.clear <- data.paixu[-nrow(data.paixu), ]
data.clear <- t(apply(data.clear, 1, scale))

colors33 <- colorRampPalette(c("white", "#ff9517", "#b83835", '#6d4490', "#000000"))(500)
breaks <- seq(0, 5, length.out = length(colors33) + 1)

plot <- pheatmap(data.paixu, 
                 scale="row",
                 cluster_rows=FALSE, 
                 cluster_cols=FALSE,
                 border_color=NA, 
                 color=colors33,
                 breaks=breaks,
                 show_colnames=FALSE)

pdf("cluster3.pdf")
print(plot)
dev.off()