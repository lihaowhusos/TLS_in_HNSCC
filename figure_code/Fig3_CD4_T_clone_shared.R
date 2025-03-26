library(tidyverse)
library(pheatmap)
library(RColorBrewer)
library(ggplot2)
library(scatterpie)

# Load data
data <- read.csv('XXX/CD4_T_and_Treg_clone_id_group.csv', header = TRUE)

# Define cluster names
names <- c("CD4_T_fh_01_00",
           "CD4_T_naive_1_01",
           "CD4_T_naive_2_03",
           "CD4_T_early_exhausted_06",
           "CD4_T_late_exhausted_08",
           "CD4_T_ISG",
           "CD4_T_reg_naive_07",
           "CD4_T_reg_1_04",
           "CD4_T_reg_2_05",
           "CD4_T_reg_3_10",
           "CD4_T_reg_ISG",
           "CD4_T_reg_exhausted_01_02",
           "CD4_T_reg_exhausted_02_09")

# Initialize dataframe for heatmap
dataframe_hm <- data.frame(matrix(ncol = 0, nrow = 13))
for (j in 1:13) {
  n <- sapply(1:13, function(k) nrow(subset(data, data[[names[j]]] > 0 & data[[names[k]]] > 0)))
  dataframe_hm <- cbind(dataframe_hm, n)
}

rownames(dataframe_hm) <- names
colnames(dataframe_hm) <- names

dataframe_hm <- log10(dataframe_hm + 1)
disp_num <- as.matrix(10^dataframe_hm - 1)

# Define colors for heatmap
colours <- colorRampPalette(c("white", "#1874CD", "#FFC125", "#EE4000", "#CD0000"))(500)

# Create heatmap
CD4_T_clone_share_hm <- pheatmap(dataframe_hm, color = colours,
                                 scale = "none", cluster_rows = FALSE,
                                 cluster_cols = FALSE, border_color = NA, display_numbers = disp_num,	
                                 number_color = "black")

# Save heatmap
ggsave("XXX/CD4_T_clone_share_hm.pdf", plot = CD4_T_clone_share_hm, width = 6.6, height = 6)

# Pie chart data preparation
x <- unlist(lapply(2:13, function(i) rep(i, i - 1)))
y <- unlist(lapply(2:13, function(i) 1:(i - 1)))

# Calculate group counts
group_counts <- function(group) {
  sapply(2:13, function(j) {
    sapply(1:(j - 1), function(k) {
      nrow(subset(data, group == group & data[[names[j]]] > 0 & data[[names[k]]] > 0))
    })
  })
}

n <- group_counts("none")
i <- group_counts("immature")
m <- group_counts("mature")
total <- group_counts(NULL)

# Create table for pie chart
table <- data.frame(x, y, n = unlist(n), i = unlist(i), m = unlist(m), total = unlist(total))
table$r <- sqrt(table$total) / 30 # Size of pie charts

# Define colors for pie chart
colors <- c("#0064D2", "#FDCA30", "#E53238")

# Create scatter pie chart
scatter_CD4 <- ggplot() +
  geom_scatterpie(aes(x = x, y = y, r = r), data = table, cols = c("n", "i", "m")) +
  geom_scatterpie_legend(table$r, n = 3, x = 2.5, y = 4) +
  scale_y_reverse() +
  scale_fill_manual(values = colors) +
  labs(x = "X=size x 25") +
  theme_classic()

# Display scatter pie chart
print(scatter_CD4)

# Save scatter pie chart
ggsave("XXX/CD4_clone_share_scatter.pdf", plot = scatter_CD4, width = 4.5, height = 4)
