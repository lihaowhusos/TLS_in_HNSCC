library(dplyr)
library(ggplot2)
library(ggpubr)
library(data.table)
library(tidyr)

# Clear the workspace
rm(list = ls())

# Set working directory
setwd("XXX")

# Load data
rt <- read.csv("XXX", sep = ",", header = TRUE, check.names = FALSE)
names(rt)[1] <- "Id"

# Filter data
rt <- rt %>% filter(str_detect(Id, "01.{1}$"))
length(unique(rt$Id))

# Convert to matrix format
rt <- as.matrix(rt)
rownames(rt) <- rt[, 1]
exp <- rt[, -1]
exp <- t(exp)

dimnames <- list(rownames(exp), colnames(exp))
data <- matrix(as.numeric(as.matrix(exp)), nrow = nrow(exp), dimnames = dimnames)

# Define standardize function
standardize <- function(x) {
  rowmean <- apply(x, 1, mean)
  rowsd <- apply(x, 1, sd)
  rv <- sweep(x, 1, rowmean, "-")
  rv <- sweep(rv, 1, rowsd, "/")
  return(rv)
}

# Define normalize function
normalize <- function(x) {
  rowmin <- apply(x, 1, min)
  rowmax <- apply(x, 1, max)
  rowmax.min <- rowmax - rowmin
  rv <- sweep(x, 1, rowmin, "-")
  rv <- sweep(rv, 1, rowmax.min, "/")
  return(rv)
}

# Normalize and standardize data
data <- normalize(data)
data.standardize <- standardize(data)

data.standardize <- t(data.standardize)
data.standardize <- as.data.frame(data.standardize)

colnames(data.standardize) <- substr(colnames(data.standardize), 24, nchar(colnames(data.standardize)))

# Convert row names to a column
data.standardize$sample <- rownames(data.standardize)

# Read cluster data
rank <- 3
file_path <- paste0("XXX", rank, ".csv")
cluster_type <- read.csv(file_path, header = TRUE, row.names = 1)

# Add TLS imprint
TLS_imprint <- read.csv("XXX")
colnames(TLS_imprint) <- c("sample", "TLS")
TLS_imprint$TLS <- as.numeric(TLS_imprint$TLS)
data.standardize <- merge(data.standardize, TLS_imprint, by = "sample")

# Merge data frames
merged_df <- merge(data.standardize, cluster_type, by = "sample")
merged_df$sample <- NULL

merged_df$group <- as.numeric(gsub("Cluster", "", merged_df$group))

# Calculate Kruskal-Wallis test results
cell_types <- colnames(merged_df)
cell_types <- cell_types[-length(cell_types)]

kruskal_results <- list()
for (i in cell_types) {
  formula <- as.formula(paste(i, "~ group"))
  kruskal_results[[i]] <- kruskal.test(formula, data = merged_df)
  print(paste("Results for", i))
  print(kruskal_results[[i]])
}

# Extract and adjust p-values
p_values <- sapply(kruskal_results, function(x) x$p.value)
adjusted_p_values <- p.adjust(p_values, method = "BH")
names(adjusted_p_values) <- cell_types

# Check adjusted p-values
print(adjusted_p_values)

# Save results to a CSV file
results_df <- data.frame(
  Cell_Type = names(adjusted_p_values),
  Adjusted_P_Value = adjusted_p_values
)
write.csv(results_df, file = "adjusted_p_values.csv", row.names = FALSE)