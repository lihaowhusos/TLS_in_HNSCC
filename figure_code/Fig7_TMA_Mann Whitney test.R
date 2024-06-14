library(tidyverse)
library(ggplot2)

# Set working directory and seed
set.seed(20240207) # Ensure reproducibility

# Read data
rt <- read.csv("XXX/cell_type_for_correlation_cluster.csv", sep = ",", header = TRUE, check.names = FALSE, row.names = 1)

# Convert TLS column to factor
rt$TLS <- as.factor(rt$TLS)

# Display column names
print(colnames(rt))

# Generate all pairwise combinations of TLS levels
group_pairs <- combn(levels(rt$TLS), 2, simplify = FALSE)

# Container for results
results <- list()

# Perform Mann-Whitney U test for each variable and each pair of groups
for (variable in setdiff(names(rt), "TLS")) {
  for (pair in group_pairs) {
    group1 <- rt %>% filter(TLS == pair[1]) %>% pull(variable)
    group2 <- rt %>% filter(TLS == pair[2]) %>% pull(variable)
    
    # Mann-Whitney U test
    test <- wilcox.test(group1, group2, exact = FALSE, alternative = "two.sided")
    
    # Save results
    results[[paste(variable, paste(pair, collapse = " vs "), sep = ": ")]] <- c(PValue = test$p.value)
  }
}

# Print results before adjustment
print("Results before adjustment:")
print(results)