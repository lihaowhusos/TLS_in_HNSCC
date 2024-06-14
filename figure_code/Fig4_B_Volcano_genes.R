library(ggpubr)
library(ggthemes)
library(ggplot2)

# Set working directory
setwd("XXX")

# Read data
deg.data <- read.csv("Find_All_Markers_scanpy_B_between_group.csv", sep = ",", header = TRUE, check.names = FALSE)
deg.data$Symbol <- deg.data$mature_names

# Calculate logP
deg.data$logP <- -log10(deg.data$mature_pvals_adj)

# Filter data
deg.data <- deg.data[!(deg.data$mature_logfoldcha > 10 | deg.data$mature_logfoldcha < -10), ]

# Assign groups based on significance
deg.data$Group <- "not-significant"
deg.data$Group[deg.data$mature_pvals_adj < 0.05 & deg.data$mature_logfoldcha > 0.5] <- "up-regulated"
deg.data$Group[deg.data$mature_pvals_adj < 0.05 & deg.data$mature_logfoldcha < -0.5] <- "down-regulated"

# Add labels for top genes
deg.data$Label <- ""
deg.data <- deg.data[order(deg.data$mature_pvals_adj), ]
up.genes <- head(deg.data$Symbol[deg.data$Group == "up-regulated"], 20)
down.genes <- head(deg.data$Symbol[deg.data$Group == "down-regulated"], 20)
deg.top20.genes <- c(as.character(up.genes), as.character(down.genes))
deg.data$Label[match(deg.top20.genes, deg.data$Symbol)] <- deg.top20.genes

# Create scatter plot
ggscatter(deg.data, x = "mature_logfoldcha", y = "logP",
          color = "Group",
          palette = c("#00a4e4", "#BBBBBB", "#ff0000"),
          size = 2,
          font.label = 9,
          repel = TRUE,
          xlab = "log2FC",
          ylab = "-log10(FDR)") + 
  theme_base() +
  geom_hline(yintercept = 1.30, linetype = "dashed") +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed") +
  geom_text(aes(label = Label), size = 3, vjust = 1, hjust = -0.1)

# Save scatter plot
ggsave("XXX/B_vocanol_plot.pdf", height = 9, width = 10)

# Read data for second plot
data <- read.table("edgerOut.xls", sep = "\t", header = TRUE, check.names = FALSE)
data$label <- c(rownames(data)[1:10], rep(NA, nrow(data) - 10))

# Create ggplot
ggplot(data, aes(log2FoldChange, -log10(padj))) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#999999") +
  geom_vline(xintercept = c(-1.2, 1.2), linetype = "dashed", color = "#999999") +
  geom_point(aes(size = -log10(padj), color = -log10(padj))) +
  scale_color_gradientn(values = seq(0, 1, 0.2),
                        colors = c("#39489f", "#39bbec", "#f9ed36", "#f38466", "#b81f25")) +
  scale_size_continuous(range = c(1, 3)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = c(0.01, 0.7),
        legend.justification = c(0, 1)) +
  guides(col = guide_colourbar(title = "-Log10_q-value"),
         size = "none") +
  geom_text(aes(label = label, color = -log10(padj)), size = 3, vjust = 1.5, hjust = 1) +
  xlab("Log2FC") +
  ylab("-Log10(FDR q-value)")

# Save ggplot
ggsave("XXX/vocanol_plot.pdf", height = 9, width = 10)