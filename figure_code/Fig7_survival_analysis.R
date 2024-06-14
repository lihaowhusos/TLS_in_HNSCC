library("survival")
library("survminer")
library("readr")

# Set working directory
setwd('XXX')

# Load clinical data
clinical <- read_table("XXX")
clinical <- clinical[clinical$Id != 1, ]

# Load cell type data
cell_type <- read.csv("XXX")
cell_type2 <- read.csv("XXX")

# Merge data
names(cell_type)[1] <- "Id"
cell_type$Id <- substr(cell_type$Id, 1, nchar(cell_type$Id) - 4)
merged_df <- merge(clinical, cell_type, by = "Id", all.x = TRUE)

names(cell_type2)[1] <- "Id"
cell_type2$Id <- substr(cell_type2$Id, 1, nchar(cell_type2$Id) - 4)
merged_df <- merge(merged_df, cell_type2, by = "Id", all.x = TRUE)

# Remove rows with NA values in specific column
rt <- merged_df[!is.na(merged_df$B.cells), ]
rt$futime <- as.numeric(rt$futime) / 30  # Convert time to months
rt$fustat <- as.numeric(rt$fustat)

cell_choose <- 'CD4_T_early_exhausted_06'
a <- rt[, cell_choose] <= median(rt[, cell_choose], na.rm = TRUE)

diff <- survdiff(Surv(futime, fustat) ~ a, data = rt)
pValue <- round(1 - pchisq(diff$chisq, df = 1), 5)
fit <- survfit(Surv(futime, fustat) ~ a, data = rt)

p <- ggsurvplot(fit,
                conf.int = FALSE,
                pval = TRUE,
                surv.median.line = "hv",
                risk.table = TRUE,
                xlab = "Time (month)",
                legend = c(0.8, 0.75),
                legend.title = "",
                legend.labs = c("high", "low"),
                palette = c("#ff0000", "#0000c0"))

file_path <- paste0("XXX", cell_choose, ".pdf")
pdf(file_path)
print(p, newpage = FALSE)
dev.off()

# Loop through cell types
cell_types <- colnames(rt)[-c(1:12)]
for (cell_choose in cell_types) {
  a <- rt[, cell_choose] <= median(rt[, cell_choose], na.rm = TRUE)
  
  diff <- survdiff(Surv(futime, fustat) ~ a, data = rt)
  pValue <- round(1 - pchisq(diff$chisq, df = 1), 5)
  fit <- survfit(Surv(futime, fustat) ~ a, data = rt)
  
  p <- ggsurvplot(fit,
                  conf.int = FALSE,
                  pval = TRUE,
                  surv.median.line = "hv",
                  risk.table = TRUE,
                  xlab = "Time (month)",
                  title = cell_choose,
                  legend = c(0.8, 0.75),
                  legend.title = "",
                  legend.labs = c("high", "low"),
                  palette = c("#ff0000", "#0000c0"))
  
  file_path <- paste0("XXX", cell_choose, ".pdf")
  pdf(file_path)
  print(p, newpage = FALSE)
  dev.off()
}

# For immune cluster
setwd('XXX')

clinical <- read_table("XXX")
clinical <- clinical[clinical$Id != 1, ]

# Read cluster
iterations <- 3
file_path <- paste0("XXX", iterations, ".consensusClass.csv")
cluster_type <- read.csv(file_path, header = FALSE, row.names = 1)
colnames(cluster_type) <- "cluster_type"

cluster_type$Id <- rownames(cluster_type)
cluster_type$Id <- substr(cluster_type$Id, 1, nchar(cluster_type$Id) - 4)

merged_df <- merge(clinical, cluster_type, by = "Id", all.x = TRUE)

# Remove rows with NA values in specific column
rt <- merged_df[!is.na(merged_df$cluster_type), ]
rt$futime <- as.numeric(rt$futime) / 30  # Convert time to months
rt$fustat <- as.numeric(rt$fustat)

fit <- survfit(Surv(futime, fustat) ~ cluster_type, data = rt)

p <- ggsurvplot(fit,
                data = rt,
                conf.int = FALSE,
                pval = TRUE,
                surv.median.line = "hv",
                risk.table = TRUE,
                xlab = "Time (month)",
                title = "cluster",
                legend = c(0.8, 0.75),
                legend.title = "")

file_path <- paste0("XXX", "1_cluster_conse.pdf")
pdf(file_path)
print(p, newpage = FALSE)
dev.off()

# For TMA
setwd('XXX')

rt <- read.csv("XXX")

fit <- survfit(Surv(suvival_time, death) ~ TLS, data = rt)

p <- ggsurvplot(fit,
                data = rt,
                conf.int = FALSE,
                pval = TRUE,
                surv.median.line = "hv",
                risk.table = TRUE,
                xlab = "Time (month)",
                title = "cluster",
                legend = c(0.8, 0.75),
                legend.title = "",
                palette = c("#0064d2", "#fdca30", "#e53238"))

file_path <- paste0("XXX", "survival_TMA_TLS_state.pdf")
pdf(file_path)
print(p, newpage = FALSE)
dev.off()