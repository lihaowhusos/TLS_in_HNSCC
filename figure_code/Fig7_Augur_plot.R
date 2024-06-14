library(ggplot2)
library(ggrepel)

# Load data
data_1 <- read.csv("XXX/Augur_NR_single.csv", row.names = 1)
data_2 <- read.csv("XXX/Augur_R_single.csv", row.names = 1)

# Select the first row and transpose
data_1 <- as.data.frame(t(data_1[1,]))
data_2 <- as.data.frame(t(data_2[1,]))

# Rename row names
rownames(data_1) <- "NR"
rownames(data_2) <- "R"

# Add Key column
data_1$Key <- rownames(data_1)
data_2$Key <- rownames(data_2)

# Merge data
data <- merge(data_1, data_2, by = "Key")

# Calculate delta_AUC
data$delta_AUC <- data$R - data$NR

# Assign colors to each point
colors <- colorRampPalette(c("#FF0000", "#FFFF00", "#00FF00", '#00FFFF', "#0000FF", "#FF00FF"))(45)
rownames_colors <- setNames(colors, data$Key)

# Calculate distances
distances <- abs(data$R - data$NR)

# Set color to gray for points with a distance less than 0.015
colors[distances < 0.015] <- "lightgray"
  
# Plot scatter plot
plot(data$NR, data$R, main = "Augur", xlab = "NR", ylab = "R", pch = 19, col = colors)
abline(a = 0, b = 1, col = "red")  # Add y=x line
text(data$NR, data$R, labels = data$Key, pos = 4, cex = 0.6)

# Save plot to PDF
pdf("XXX/Augur_plot.pdf", width = 4.5, height = 8)
plot(data$R, data$delta_AUC, main = "Augur", xlab = "R", ylab = "delta_Augur_score", pch = 19, col = colors)
abline(a = 0, b = 0, col = "red")  # Add y=0 line
text(data$R, data$delta_AUC, labels = data$Key, pos = 4, cex = 0.6)
dev.off()

# Plot scatter plot using ggplot2
ggplot(data, aes(x = NR, y = R)) + 
  geom_point(aes(colour = colors)) +  # Plot points
  geom_text_repel(aes(label = Key), box.padding = 0.35, point.padding = 0.5, max.overlaps = Inf) + 
  geom_abline(intercept = 0, slope = 1, col = "black") +  # Add y=x line
  theme_minimal() +
  ggtitle("Scatter Plot Example")