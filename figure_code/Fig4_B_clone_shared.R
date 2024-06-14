library(tidyverse)
library(pheatmap)
library(RColorBrewer)
library(ggplot2)
library(scatterpie)

# Read data
data <- read.csv('XXX/BCR_clone_id_group.csv', header = TRUE)

# Recode groups
data$group <- recode(data$sample, 
                     'szj105988' = 'immature',
                     'szj106005' = 'none',
                     'szj106495' = 'none',
                     'szj106560' = 'immature',                
                     'szj106562' = 'mature',
                     'szj107010' = 'mature',
                     'szj107145' = 'mature',
                     'szj107734' = 'mature',
                     'szj107849' = 'none',
                     'szj108352' = 'immature',
                     'szj106121' = 'mature',                  
                     'szj106138' = 'immature',                  
                     'szj106759' = 'none',
                     'szj106771' = 'mature'
)

# Define clone names
clone_names <- c("B_IGHM_04", "B_CD74_1_00", "B_CD69_1_01", "B_MHC_II_05", "B_GC_11", 
                 "B_ISG_08", "Cycling_B_plamsa_cell_09", "Plasma_cell_02_03", "Plasma_cell_01_02")

# Initialize heatmap dataframe
dataframe_hm <- data.frame(matrix(ncol = 0, nrow = 9))

# Populate heatmap dataframe
for (j in 1:9) {
  n <- sapply(1:9, function(k) nrow(subset(data, data[[clone_names[j]]] > 0 & data[[clone_names[k]]] > 0)))
  dataframe_hm <- cbind(dataframe_hm, n)
}

rownames(dataframe_hm) <- clone_names
colnames(dataframe_hm) <- clone_names

# Log transform data and prepare display numbers
dataframe_hm <- log10(dataframe_hm + 1)
disp_num <- as.matrix(10^dataframe_hm - 1)

# Define heatmap colors
colours <- colorRampPalette(c("white", "#1874CD", "#FFC125", "#EE4000", "#CD0000"))(500)

# Create heatmap
B_clone_share_hm <- pheatmap(dataframe_hm, color = colours,
                             scale = "none", cluster_rows = FALSE,
                             cluster_cols = FALSE, border_color = NA,
                             display_numbers = disp_num,	
                             number_color = "black")

# Save heatmap
ggsave("XXX/B_clone_share_hm.pdf", plot = B_clone_share_hm, width = 6.5, height = 6)

# Prepare data for scatter pie chart
x <- unlist(lapply(2:9, function(i) rep(i, i-1)))
y <- unlist(lapply(2:9, function(i) 1:(i-1)))

# Initialize vectors
n <- i <- m <- total <- c()

# Populate vectors
for (j in 2:9) {
  for (k in 1:(j-1)) {
    c_name <- clone_names[j]
    d_name <- clone_names[k]
    n <- c(n, nrow(subset(data, group == "none" & data[[c_name]] > 0 & data[[d_name]] > 0)))
    i <- c(i, nrow(subset(data, group == "immature" & data[[c_name]] > 0 & data[[d_name]] > 0)))
    m <- c(m, nrow(subset(data, group == "mature" & data[[c_name]] > 0 & data[[d_name]] > 0)))
    total <- c(total, nrow(subset(data, data[[c_name]] > 0 & data[[d_name]] > 0)))
  }
}

# Create table for scatter pie chart
table <- data.frame(x, y, n, i, m, total)
table$r <- sqrt(table$total) / 30

# Define colors for scatter pie chart
colors <- c("#0064D2", "#FDCA30", "#E53238")

# Create scatter pie chart
scatter_B <- ggplot() +
  geom_scatterpie(aes(x = x, y = y, r = r), data = table, cols = c("n", "i", "m")) +
  geom_scatterpie_legend(table$r, n = 3, x = 3, y = 7) +
  scale_y_reverse() +
  scale_fill_manual(values = colors) +
  labs(x = "X=size x 30") +
  theme_classic()

# Save scatter pie chart
ggsave("XXX/B_clone_share_scatter.pdf", plot = scatter_B, width = 4.5, height = 4)