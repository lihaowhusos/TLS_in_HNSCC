library(tidyverse)
library(pheatmap)
library(RColorBrewer)
library(ggplot2)
library(scatterpie)

# Read data
data <- read.csv('XXX/CD8_clone_id_group.csv', header = TRUE)

# Recode sample groups
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

# Define names of interest
names <- c("CD8_T_naive_01",
           "CD8_T_effector_memory_01",
           "CD8_T_cytotoxic_02",
           "CD8_T_late_exhausted_01",
           "CD8_T_ISG")

number <- length(names)

# Initialize dataframe for heatmap
dataframe_hm <- data.frame(matrix(ncol = number, nrow = number))

# Populate the heatmap dataframe
for (j in 1:number) {
  for (k in 1:number) {
    c_name <- names[j]
    d_name <- names[k]
    n_sub <- nrow(subset(data, data[[c_name]] > 0 & data[[d_name]] > 0))
    dataframe_hm[j, k] <- n_sub
  }
}

rownames(dataframe_hm) <- names
colnames(dataframe_hm) <- names

# Apply log transformation
dataframe_hm <- log10(dataframe_hm + 1)

# Prepare data for display
disp_num <- as.matrix(10^dataframe_hm - 1)

# Define color palette
colours <- colorRampPalette(c("white", "#EBE5ED", "#D2D1E7", "#AABDD9", "#73ABC9", "#4A8FBD", "#338186", "#29685B", "#174233"))(500)

# Define breaks for the heatmap
my_breaks <- seq(1, 3.1, length.out = 500)

# Generate heatmap
CD8_T_clone_share_hm <- pheatmap(dataframe_hm, color = colours,
                                 scale = "none", cluster_rows = FALSE,
                                 cluster_cols = FALSE, border_color = NA, display_numbers = disp_num,	
                                 number_color = "black",
                                 legend_breaks = c(0, 1, 2, 3),
                                 breaks = my_breaks)

# Save heatmap to file
ggsave("XXX/CD8_T_clone_share_hm.pdf", plot = CD8_T_clone_share_hm, width = 4.5, height = 4)

# Initialize x and y coordinates
x <- c()
y <- c()

# Populate coordinates
for (i in 2:number) {
  k <- i - 1
  x <- c(x, rep(i, k))
  y <- c(y, c(1:k))
}

# Initialize counts
n <- c()
i <- c()
m <- c()
total <- c()

# Populate counts based on group and clone presence
for (j in 2:number) {
  s <- j - 1
  for (k in 1:s) {
    c_name <- names[j]
    d_name <- names[k]
    n_sub <- nrow(subset(data, group == "none" & data[[c_name]] > 0 & data[[d_name]] > 0))
    i_sub <- nrow(subset(data, group == "immature" & data[[c_name]] > 0 & data[[d_name]] > 0))
    m_sub <- nrow(subset(data, group == "mature" & data[[c_name]] > 0 & data[[d_name]] > 0))
    total_sub <- nrow(subset(data, data[[c_name]] > 0 & data[[d_name]] > 0))
    n <- c(n, n_sub)
    i <- c(i, i_sub)
    m <- c(m, m_sub)
    total <- c(total, total_sub)
  }
}

# Create dataframe for scatter pie chart
table <- data.frame(x, y, n, i, m, total)

# Define radius for pie chart
table$r <- sqrt(table$total) / 50

# Define colors
colors <- c("#0064D2", "#FDCA30", "#E53238")

# Generate scatter pie chart
scatter_CD8 <- ggplot() +
  geom_scatterpie(aes(x = x, y = y, r = r), 
                  data = table, 
                  cols = c("n", "i", "m")) +
  geom_scatterpie_legend(table$r, n = 3, x = 2.5, y = 4) +  # Legend location
  scale_y_reverse() +
  scale_fill_manual(values = colors) +
  labs(x = "X=size x 50") +
  theme_classic()

# Display scatter pie chart
print(scatter_CD8)

# Save scatter pie chart to file
ggsave("XXX/CD8_clone_share_scatter.pdf", plot = scatter_CD8, width = 4.5, height = 4)