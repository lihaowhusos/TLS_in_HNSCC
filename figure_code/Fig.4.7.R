library(ggpubr)
library(ggthemes)
library(dplyr)

# Set working directory
setwd("XXX")  # Set working directory

# Load data
deg.data <- read.csv("XXX", sep = ",", header = TRUE, check.names = FALSE)

# Calculate -log10 of adjusted p-values
deg.data$logP <- -log10(deg.data$pvals_adj)
deg.data$rank <- deg.data[, 1]

# Adjust logP based on meanchange
deg.data$logP[deg.data$meanchange < 0] <- -deg.data$logP[deg.data$meanchange < 0]

# Replace NA and empty column names
names(deg.data)[is.na(names(deg.data))] <- "new_name_for_NA"
names(deg.data)[names(deg.data) == ""] <- "new_name_for_empty"

# Cap logP values at 300
deg.data <- deg.data %>% mutate(logP = ifelse(logP > 300, 300, logP))

# Get unique cell types
celltypes <- unique(deg.data$group)

# Generate plots for each cell type
for (celltype in celltypes) {
  deg.data_plot <- subset(deg.data, group == celltype)
  deg.data$Label <- ""
  
  # Identify top 20 up-regulated and down-regulated genes
  up.genes <- head(deg.data$Symbol[deg.data$Group == "up-regulated"], 20)
  down.genes <- head(deg.data$Symbol[deg.data$Group == "down-regulated"], 20)
  deg.top20.genes <- c(as.character(up.genes), as.character(down.genes))
  deg.data$Label[match(deg.top20.genes, deg.data$Symbol)] <- deg.top20.genes
  
  # Generate scatter plot
  ggscatter(deg.data_plot, x = "rank", y = "logP", 
            size = 2,
            font.label = 9,
            repel = TRUE,
            xlab = "Rank",
            color = "logP",
            ylab = "-log10(P.adj)"
  ) + 
    scale_color_gradient2(midpoint = 0, low = "dodgerblue3", mid = "gray100", high = "#CD2626") +
    theme_base()
  
  # Save plot as PDF
  file_path <- paste0("XXX", celltype, "_snake_plot.pdf")
  ggsave(file_path, height = 10, width = 4)
}