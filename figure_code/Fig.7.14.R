library(tidyverse)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)

# Read data
rt <- read.csv("XXX/cell_type_for_correlation.csv", sep = ",", header = TRUE, check.names = FALSE, row.names = 1)

# Calculate correlation
cor_cell_type <- cor(rt)

# Define color palette
colours <- colorRampPalette(c("white", "#1874CD", "#FFC125", "#EE4000", "#CD0000"))(1000)

# Create heatmap
hm <- pheatmap(cor_cell_type, 
               color = colours,
               scale = "none", 
               cluster_rows = FALSE,
               cluster_cols = FALSE,
               border_color = NA,
               number_color = "black",
               display_numbers = TRUE)

# Save heatmap
ggsave("XXX/hm_total_correlation.pdf", plot = hm, width = 6.6, height = 6)

# Read data for scatter plots
rt <- read.csv("XXX/cell_type_for_correlation_cluster.csv", sep = ",", header = TRUE, check.names = FALSE, row.names = 1)

# Define colors
colors <- c("#FDCA30", "#E53238", "#0064D2")

# Create scatter plot
ggplot(rt, aes(x = ratio_CD20, y = ratio_CD4_TCF1, fill = TLS)) +
  geom_point(shape = 21, alpha = 1, size = 3) +
  scale_fill_manual(values = colors) +
  scale_color_manual(values = colors) +
  geom_smooth(method = lm, color = "black", fill = "grey", se = TRUE) +
  theme_classic()

# Create scatter plot with size and color
p_im <- ggplot(rt, aes(x = X, y = Y, size = clone_id_size, color = cluster_2)) +
  geom_point(alpha = 0.25) +
  scale_color_manual(values = colors)

# Create scatter plot with size and fill
p_im_fill <- ggplot(CD4_TCR_table_immature, aes(x = X, y = Y, size = clone_id_size, fill = cluster_2)) +
  geom_point(colour = "black", shape = 21, alpha = 1) +
  scale_fill_manual(values = colors) +
  scale_size_continuous(
    range = c(1, 9),
    limits = c(0, 40),
    breaks = c(0, 20, 40),
    labels = c(0, 20, 40)
  ) +
  theme_bw()

# Define colors for scatter pie
colors <- c("#0064D2", "#FDCA30", "#E53238")

# Create scatter pie plot
scatter_CD4 <- ggplot() +
  geom_scatterpie(aes(x = x, y = y, r = r), 
                  data = table, 
                  cols = c("n", "i", "m")) +
  geom_scatterpie_legend(table$r, n = 3, x = 2.5, y = 4) + # x, y is the location of legend
  scale_y_reverse() +
  scale_fill_manual(values = colors) +
  labs(x = "X=size x 25") +
  theme_classic()