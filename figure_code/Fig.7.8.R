library(pheatmap)
library(broom)
library(forestplot)

# Define the new order of columns
new_order <- c("response",
               "B_IGHM_04", "B_CD74_1_00", "B_CD69_1_01", "B_MHC_II_05", "B_GC_11", "B_ISG_08",
               "Plasma_cell_01_02", "Plasma_cell_02_03", "CD8_T_naive_01", "CD8_T_effector_memory_01",
               "CD8_T_cytotoxic_02", "CD8_T_late_exhausted_01", "CD8_T_ISG", "CD4_T_fh_01_00",
               "CD4_T_naive_1_01", "CD4_T_naive_2_03", "CD4_T_early_exhausted_06", "CD4_T_late_exhausted_08",
               "CD4_T_ISG", "cDC_1_15", "cDC_2_06", "DC_LAMP3_11", "pDC_10", "CD4_T_reg_naive_07",
               "CD4_T_reg_1_04", "CD4_T_reg_2_05", "CD4_T_reg_3_10", "CD4_T_reg_ISG", "CD4_T_reg_exhausted_01_02",
               "CD4_T_reg_exhausted_02_09", "Cleaning_macrophage_13", "Cycling_B_plamsa_cell_09", "Cycling_CD4_T",
               "Cycling_CD8_T", "Cycling_NK", "Cycling_myeloid_cell_14", "Endothelial_17", "Fibroblast_13",
               "ILC_06", "Lymphatic_endothelial_cell_21", "M1_S100A8_07", "M2_COL1A1_04", "M2_CXCL10_12",
               "M2_MARCO_05", "M2_MMP9_08", "M2_SELENOP_02", "M2_STAB1_09", "Mast_cell_00", "Neutrophil_1_01",
               "Neutrophil_1_03")

# Read cluster data
folder_path <- "XXX"
file_path <- paste0(folder_path, "GSE136961_RNA_cell2location_for_GLM.csv")
cluster <- read.csv(file_path, header = TRUE, row.names = 1)

df <- cluster[, new_order]
df$response <- ifelse(df$response == "R", 1, 0)

# Create formula for logistic regression
features <- paste(colnames(df)[-1], collapse = " + ")
formula <- as.formula(paste('response ~', features))

# Fit the logistic regression model
model <- glm(formula, data = df, family = binomial)

# Check for model convergence
if (!model$converged) {
  stop("Model failed to converge")
}

# Get summary of the model
summary_model <- summary(model)
coefficients_summary <- summary_model$coefficients

# Calculate Wald confidence intervals
wald_ci <- coefficients_summary[, 1] + c(-1, 1) * qnorm(0.975) * coefficients_summary[, 2]

# Create a data frame of results
results_df <- data.frame(
  Feature = colnames(df)[-1],
  Estimate = coefficients_summary[-1, 1], 
  Std.Error = coefficients_summary[-1, 2],
  Z.Value = coefficients_summary[-1, 3],
  P.Value = coefficients_summary[-1, 4],
  OR = exp(coefficients_summary[-1, 1]),
  LowerCI95 = exp(wald_ci[-1, 1]),
  UpperCI95 = exp(wald_ci[-1, 2])
)

# Print the results data frame
print(results_df)

# Order results by P.Value
results_df <- results_df[order(results_df$P.Value), ]

# Create the forest plot
forestplot(labeltext = results_df$Feature,
           mean = results_df$OR,
           lower = results_df$LowerCI95,
           upper = results_df$UpperCI95,
           is.summary = c(TRUE, rep(FALSE, nrow(results_df) - 1)),
           xlog = TRUE,
           title = "Forest plot of Odds Ratios",
           clip = c(0.1, 10),
           zero = 1,
           xlab = "Odds Ratio (log scale)",
           col = fpColors(box = "royalblue", line = "darkblue", summary = "royalblue"))

# Read cluster data for Cancer_Dis
file_path <- paste0(folder_path, "Cancer_Dis_cell2location_for_GLM.csv")
cluster <- read.csv(file_path, header = TRUE, row.names = 1)

df <- cluster[, new_order]
df$response <- ifelse(df$response == "R", 1, 0)

# Prepare an empty dataframe to store results
results_df <- data.frame(
  Feature = character(),
  Estimate = numeric(),
  Std.Error = numeric(),
  Statistic = numeric(),
  P.Value = numeric(),
  OR = numeric(),
  LowerCI95 = numeric(),
  UpperCI95 = numeric(),
  stringsAsFactors = FALSE
)

# Loop through the feature columns
for (feature in colnames(df)[-1]) { 
  formula <- as.formula(paste('response ~', feature))
  model <- glm(formula, data = df, family = binomial)
  tidy_model <- broom::tidy(model)
  
  results_df <- rbind(results_df, data.frame(
    Feature = feature,
    Estimate = tidy_model$estimate[2],
    Std.Error = tidy_model$std.error[2],
    Statistic = tidy_model$statistic[2],
    P.Value = tidy_model$p.value[2],
    OR = exp(tidy_model$estimate[2]),
    LowerCI95 = exp(tidy_model$estimate[2] - 1.96 * tidy_model$std.error[2]),
    UpperCI95 = exp(tidy_model$estimate[2] + 1.96 * tidy_model$std.error[2])
  ))
}

# Order results by P.Value
results_df <- results_df[order(results_df$P.Value), ]

# Create the forest plot
forestplot(labeltext = results_df$Feature,
           mean = results_df$OR,
           lower = results_df$LowerCI95,
           upper = results_df$UpperCI95,
           is.summary = c(TRUE, rep(FALSE, nrow(results_df) - 1)),
           xlog = TRUE,
           title = "Forest plot of Odds Ratios",
           clip = c(0.1, 10),
           zero = 1,
           xlab = "Odds Ratio (log scale)",
           col = fpColors(box = "royalblue", line = "darkblue", summary = "royalblue"))

# Read cluster data for GSE136961
file_path <- paste0(folder_path, "GSE136961_RNA_cell2location_for_GLM.csv")
cluster <- read.csv(file_path, header = TRUE, row.names = 1)

data <- cluster[, -1]
data <- t(data)

# Standardization function
standardize <- function(x) {
  rowmean <- apply(x, 1, mean)
  rowsd <- apply(x, 1, sd)
  rv <- sweep(x, 1, rowmean, "-")
  rv <- sweep(rv, 1, rowsd, "/")
  return(rv)
}

# Normalization function
normalize <- function(x) {
  rowmin <- apply(x, 1, min)
  rowmax <- apply(x, 1, max)
  rowmax.min <- rowmax - rowmin
  rv <- sweep(x, 1, rowmin, "-")
  rv <- sweep(rv, 1, rowmax.min, "/")
  return(rv)
}

data <- normalize(data)
data.standardize <- standardize(data)
data.standardize <- t(data.standardize)
data.standardize <- as.data.frame(data.standardize)

data.standardize$response <- cluster$response
data.standardize$response <- ifelse(data.standardize$response == "R", 1, 0)