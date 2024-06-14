library(stats)

# Set working directory
setwd("XXX")

# Read the data
rt <- read.csv("XXX", sep = ",", header = TRUE, check.names = FALSE)

# Collated matrix format
rt <- as.matrix(rt)
rownames(rt) <- rt[, 1]
exp <- rt[, -1]
exp <- t(exp)

dimnames <- list(rownames(exp), colnames(exp))
data <- matrix(as.numeric(as.matrix(exp)), nrow = nrow(exp), dimnames = dimnames)

# Standardize function
standardize <- function(x) {
  rowmean <- apply(x, 1, mean)
  rowsd <- apply(x, 1, sd)
  rv <- sweep(x, 1, rowmean, "-")  # Subtract mean
  rv <- sweep(rv, 1, rowsd, "/")  # Divide by standard deviation
  return(rv)
}

# Normalize function
normalize <- function(x) {
  rowmin <- apply(x, 1, min)
  rowmax <- apply(x, 1, max)
  rowmax.min <- rowmax - rowmin
  rv <- sweep(x, 1, rowmin, "-")  # Subtract min
  rv <- sweep(rv, 1, rowmax.min, "/")  # Divide by range
  return(rv)
}

data <- normalize(data)
data.standardize <- standardize(data)

data.standardize <- t(data.standardize)
data.standardize <- as.data.frame(data.standardize)

cell_types <- colnames(data.standardize)

iterations <- 3
file_path <- paste0("XXX", iterations, ".consensusClass.csv")
cluster_type <- read.csv(file_path, header = FALSE, row.names = 1)
colnames(cluster_type) <- "cluster_type"

# Convert row names to a column
data.standardize$RowNames <- rownames(data.standardize)
cluster_type$RowNames <- rownames(cluster_type)

# Merge the data frames based on the new RowNames column
merged_df <- merge(data.standardize, cluster_type, by = "RowNames")

# Remove RowNames column after merging
merged_df$RowNames <- NULL
merged_df$cluster_type <- as.numeric(merged_df$cluster_type)

kruskal_results <- list()
for (i in cell_types) {
  # Create the formula with the current column name
  formula <- as.formula(paste(i, "~ cluster_type"))
  
  # Run the Kruskal-Wallis test with the constructed formula
  kruskal_results[[i]] <- kruskal.test(formula, data = merged_df)
  
  # Print the results with the current cell type name
  print(paste("Results for", i))
  print(kruskal_results[[i]])
}

# Extracting p-values from each kruskal_result
p_values <- sapply(kruskal_results, function(x) x$p.value)

# Adjusting the p-values for multiple testing
adjusted_p_values <- p.adjust(p_values, method = "BH")

# Keep track of which p-value corresponds to which test
names(adjusted_p_values) <- cell_types

# Check the adjusted p-values
print(adjusted_p_values)

results_df <- data.frame(
  Cell_Type = names(adjusted_p_values),
  Adjusted_P_Value = adjusted_p_values
)

# Save the dataframe to a CSV file
write.csv(results_df, file = "XXX/adjusted_p_values.csv", row.names = FALSE)