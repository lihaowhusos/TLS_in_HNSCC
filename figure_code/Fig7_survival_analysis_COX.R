library(survival)
library(survminer)
library(readr)
library(dplyr)
library(stringr)

# Set working directory
setwd("XXX")

# Load and preprocess clinical data
clinical2 <- read.csv("XXX/clinical_data/survival_HNSC_survival.csv")
names(clinical2)[2] <- "Id"
clinical2 <- clinical2 %>% distinct(Id, .keep_all = TRUE)

# Load and preprocess cell type data
cell_type <- read.csv("Bulk_RNA_cell2location_20240202.csv")
names(cell_type)[1] <- "Id"
cell_type <- cell_type %>% filter(str_detect(Id, "01.{1}$"))
cell_type$Id <- substr(cell_type$Id, 1, nchar(cell_type$Id) - 4)

# Merge clinical and cell type data
rt <- inner_join(clinical2, cell_type, by = "Id")

# Convert and normalize survival time
rt <- rt %>%
  mutate(
    OS = as.numeric(OS),
    OS.time = as.numeric(OS.time) / 30,
    DSS = as.numeric(DSS),
    DSS.time = as.numeric(DSS.time) / 30,
    PFI = as.numeric(PFI),
    PFI.time = as.numeric(PFI.time) / 30
  )

# Function to plot survival analysis
plot_survival <- function(rt, cell_choose, time_col, status_col, output_dir) {
  a <- rt[, cell_choose] <= median(rt[, cell_choose], na.rm = TRUE)
  diff <- survdiff(Surv(rt[[time_col]], rt[[status_col]]) ~ a, data = rt)
  pValue <- round(1 - pchisq(diff$chisq, df = 1), 5)
  fit <- survfit(Surv(rt[[time_col]], rt[[status_col]]) ~ a, data = rt)
  
  p <- ggsurvplot(
    fit,
    conf.int = FALSE,
    pval = TRUE,
    surv.median.line = "hv",
    risk.table = TRUE,
    xlab = "Time (month)",
    legend = c(0.8, 0.75),
    legend.title = "",
    legend.labs = c("high", "low"),
    palette = c("#ff0000", "#0000c0")
  )
  
  file_path <- paste0("XXX/", output_dir, "/", cell_choose, ".pdf")
  pdf(file_path)
  print(p, newpage = FALSE)
  dev.off()
}

# Loop for OS, DSS, and PFI survival analysis
cell_types <- colnames(rt)[-c(1:12)]
for (cell_choose in cell_types) {
  plot_survival(rt, cell_choose, "OS.time", "OS", "survival_analysis_OS")
  plot_survival(rt, cell_choose, "DSS.time", "DSS", "survival_analysis_DSS")
  plot_survival(rt, cell_choose, "PFI.time", "PFI", "survival_analysis_PFI")
}

# For immune cluster analysis
clinical2 <- read.csv("XXX/clinical_data/survival_HNSC_survival.csv")
names(clinical2)[2] <- "Id"
clinical2 <- clinical2 %>% distinct(Id, .keep_all = TRUE)

# Read cluster data
iterations <- 3
file_path <- paste0("NMF/group_", iterations, ".csv")
cluster_type <- read.csv(file_path, header = TRUE, row.names = 1)
colnames(cluster_type) <- c("cluster_type", "Id")
cluster_type$Id <- substr(cluster_type$Id, 1, nchar(cluster_type$Id) - 4)

merged_df <- merge(clinical2, cluster_type, by = "Id", all.x = TRUE)
merged_df <- merged_df[!is.na(merged_df$cluster_type), ]
merged_df$cluster_type <- substr(merged_df$cluster_type, 8, nchar(merged_df$cluster_type))

rt <- merged_df[!is.na(merged_df$cluster_type), ]
rt <- rt %>%
  mutate(futime = as.numeric(OS.time) / 30, fustat = as.numeric(OS))

fit <- survfit(Surv(futime, fustat) ~ cluster_type, data = rt)

p <- ggsurvplot(
  fit,
  data = rt,
  conf.int = FALSE,
  pval = TRUE,
  surv.median.line = "hv",
  risk.table = TRUE,
  xlab = "Time (month)",
  title = "cluster",
  legend = c(0.8, 0.75),
  legend.title = "",
  legend.labs = c("m", "n", "im"),
  palette = c("#e53238", "#0064d2", "#fdca30")
)

file_path <- paste0("NMF/survival_cluster.pdf")
pdf(file_path)
print(p, newpage = FALSE)
dev.off()

# For TMA analysis
setwd("XXX")
rt <- read.csv("TMA_survival.csv")
fit <- survfit(Surv(suvival_time, death) ~ TLS, data = rt)

p <- ggsurvplot(
  fit,
  data = rt,
  conf.int = FALSE,
  pval = TRUE,
  surv.median.line = "hv",
  risk.table = TRUE,
  xlab = "Time (month)",
  title = "cluster",
  palette = c("#0064d2", "#fdca30", "#e53238")
)

file_path <- paste0("XXX/survival_TMA_TLS_state.pdf")
pdf(file_path)
print(p, newpage = FALSE)
dev.off()

# Multivariate Cox regression analysis
clinical2 <- read.csv("XXX/clinical_data/survival_HNSC_survival.csv")
names(clinical2)[2] <- "Id"
clinical2 <- clinical2 %>% distinct(Id, .keep_all = TRUE)
clinical_matrix <- read.csv("XXX/clinical_data/TCGA.HNSC.sampleMap_HNSC_clinicalMatrix.csv")

# Read cluster data
iterations <- 3
file_path <- paste0("NMF/group_", iterations, ".csv")
cluster_type <- read.csv(file_path, header = TRUE, row.names = 1)
colnames(cluster_type) <- c("cluster_type", "Id")
cluster_type$Id <- substr(cluster_type$Id, 1, nchar(cluster_type$Id) - 4)

merged_df <- merge(clinical2, cluster_type, by = "Id", all.x = TRUE)
merged_df <- merged_df[!is.na(merged_df$cluster_type), ]
merged_df$cluster_type <- substr(merged_df$cluster_type, 8, nchar(merged_df$cluster_type))
merged_df <- merge(merged_df, clinical_matrix, by = "sample", all.x = TRUE)

rt <- merged_df[!is.na(merged_df$cluster_type), ]
rt <- rt %>%
  mutate(
    cluster_type = factor(cluster_type, levels = c("2", "3", "1")),
    clinical_stage = factor(clinical_stage, levels = c("Stage I", "Stage II", "Stage III", "Stage IVA", "Stage IVB", "Stage IVC", "None")),
    futime = as.numeric(PFI.time) / 30,
    fustat = as.numeric(PFI)
  )

valid_rows <- !is.na(rt$PFI)
rt_non_NA <- rt[valid_rows, ]

res.cox <- coxph(Surv(futime, fustat) ~ age_at_initial_pathologic_diagnosis + gender + clinical_stage + cluster_type, data = rt_non_NA)
sur.forest <- ggforest(res.cox, data = rt_non_NA, main = paste0('Hazard ratio of TCGA_HNSC_', iterations), noDigits = 3)

pdf(paste0("COX/TCGA_HNSC_COX_PFI.pdf"))
print(sur.forest)
dev.off()

# Universal and multivariate Cox regression analysis
dependent_os <- "Surv(OS.time, OS)"
dependent_dss <- "Surv(DSS.time, DSS)"
dependent_pfi <- "Surv(PFI.time, PFI)"
explanatory <- c("age_at_initial_pathologic_diagnosis", "gender", "clinical_stage", "cluster_type")

sur.univariable.multivariable.OS <- rt %>%
  finalfit(dependent_os, explanatory, add_dependent_label = FALSE) %>%
  rename("Overall survival" = label, " " = levels, "  " = all)

sur.univariable.multivariable.DSS <- rt %>%
  finalfit(dependent_dss, explanatory, add_dependent_label = FALSE) %>%
  rename("DSS" = label, " " = levels, "  " = all)

sur.univariable.multivariable.PFI <- rt %>%
  finalfit(dependent_pfi, explanatory, add_dependent_label = FALSE) %>%
  rename("PFI" = label, " " = levels, "  " = all)