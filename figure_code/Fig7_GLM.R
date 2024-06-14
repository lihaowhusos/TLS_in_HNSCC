########################## GSE136961 ############################
library(pheatmap)
library(broom)
library(forestplot)

new_order <- c("response",
  "q05cell_abundance_w_sf_B_IGHM_04", "q05cell_abundance_w_sf_B_CD74_1_00",
  "q05cell_abundance_w_sf_B_CD69_1_01", "q05cell_abundance_w_sf_B_MHC_II_05",
  "q05cell_abundance_w_sf_B_GC_11", "q05cell_abundance_w_sf_B_ISG_08",
  "q05cell_abundance_w_sf_Plasma_cell_01_02", "q05cell_abundance_w_sf_Plasma_cell_02_03",
  "q05cell_abundance_w_sf_CD8_T_naive_01", "q05cell_abundance_w_sf_CD8_T_effector_memory_01",
  "q05cell_abundance_w_sf_CD8_T_cytotoxic_02", "q05cell_abundance_w_sf_CD8_T_late_exhausted_01",
  "q05cell_abundance_w_sf_CD8_T_ISG",
  "q05cell_abundance_w_sf_CD4_T_fh_01_00", "q05cell_abundance_w_sf_CD4_T_naive_1_01",
  "q05cell_abundance_w_sf_CD4_T_naive_2_03",
  "q05cell_abundance_w_sf_CD4_T_early_exhausted_06", "q05cell_abundance_w_sf_CD4_T_late_exhausted_08",
  "q05cell_abundance_w_sf_CD4_T_ISG", "q05cell_abundance_w_sf_cDC_1_15",
  "q05cell_abundance_w_sf_cDC_2_06", "q05cell_abundance_w_sf_DC_LAMP3_11", "q05cell_abundance_w_sf_pDC_10",
  "q05cell_abundance_w_sf_CD4_T_reg_naive_07",
  "q05cell_abundance_w_sf_CD4_T_reg_1_04", "q05cell_abundance_w_sf_CD4_T_reg_2_05",
  "q05cell_abundance_w_sf_CD4_T_reg_3_10", "q05cell_abundance_w_sf_CD4_T_reg_ISG",
  "q05cell_abundance_w_sf_CD4_T_reg_exhausted_01_02", "q05cell_abundance_w_sf_CD4_T_reg_exhausted_02_09",
  "q05cell_abundance_w_sf_Cleaning_macrophage_13",
  "q05cell_abundance_w_sf_Cycling_B_plamsa_cell_09", "q05cell_abundance_w_sf_Cycling_CD4_T",
  "q05cell_abundance_w_sf_Cycling_CD8_T", "q05cell_abundance_w_sf_Cycling_NK",
  "q05cell_abundance_w_sf_Cycling_myeloid_cell_14", 
  "q05cell_abundance_w_sf_Endothelial_17", "q05cell_abundance_w_sf_Fibroblast_13",
  "q05cell_abundance_w_sf_ILC_06", "q05cell_abundance_w_sf_Lymphatic_endothelial_cell_21",
  "q05cell_abundance_w_sf_M1_S100A8_07", "q05cell_abundance_w_sf_M2_COL1A1_04",
  "q05cell_abundance_w_sf_M2_CXCL10_12", "q05cell_abundance_w_sf_M2_MARCO_05",
  "q05cell_abundance_w_sf_M2_MMP9_08", "q05cell_abundance_w_sf_M2_SELENOP_02",
  "q05cell_abundance_w_sf_M2_STAB1_09", "q05cell_abundance_w_sf_Mast_cell_00",
  "q05cell_abundance_w_sf_Neutrophil_1_01",
  "q05cell_abundance_w_sf_Neutrophil_1_03"
) 



###read cluster
folder_path = "G:\\TLSdata\\TLS_total_data\\ICB_varify\\GSE136961\\"
file_path = paste0(folder_path,"GSE136961_RNA_cell2location_for_GLM.csv")
cluster = read.csv(file_path,header = TRUE,row.names=1)


df = cluster[, new_order]
df$response <- ifelse(df$response == "R", 1, 0)



# Assuming the first column is the 'response' and all other columns are features
features <- paste(colnames(df)[-1], collapse = " + ")
formula <- as.formula(paste('response ~', features))

# Fit the logistic regression model with all features included
model <- glm(formula, data = df, family = binomial)

# Check for model convergence
if (!model$converged) {
  stop("Model failed to converge")
}

# Use summary to get the Wald confidence intervals
summary_model <- summary(model)

# Directly extract the coefficients, standard errors, z values and p-values
coefficients_summary <- summary_model$coefficients

# Calculate the Wald confidence intervals
wald_ci <- coefficients_summary[, 1] + c(-1, 1) * qnorm(0.975) * coefficients_summary[, 2]

# Create a data frame of results
results_df <- data.frame(
  Feature = colnames(df)[-1],
  Estimate = coefficients_summary[-1, 1], # Excluding Intercept
  Std.Error = coefficients_summary[-1, 2], # Excluding Intercept
  Z.Value = coefficients_summary[-1, 3], # Excluding Intercept
  P.Value = coefficients_summary[-1, 4], # Excluding Intercept
  OR = exp(coefficients_summary[-1, 1]), # Excluding Intercept
  LowerCI95 = exp(wald_ci[-1, 1]), # Excluding Intercept
  UpperCI95 = exp(wald_ci[-1, 2]) # Excluding Intercept
)

# Print the results data frame
print(results_df)


# Order results by P.Value or OR
results_df <- results_df[order(results_df$P.Value), ]


# Create the forest plot
forestplot(labeltext = results_df$Feature,
           mean = results_df$OR,
           lower = results_df$LowerCI95,
           upper = results_df$UpperCI95,
           is.summary=c(TRUE, rep(FALSE, nrow(results_df)-1)),
           xlog = TRUE, # Log scale for Odds Ratio
           title = "Forest plot of Odds Ratios",
           clip = c(0.1, 10), # Adjust these values as needed for your data
           zero = 1, # Zero line for Odds Ratio of 1
           xlab = "Odds Ratio (log scale)",
           #cex.lab=.75,
           col=fpColors(box="royalblue", line="darkblue", summary="royalblue"))




########################################################## Cancer_dis ##############################################
library(pheatmap)

###read cluster
folder_path = "G:\\TLSdata\\TLS_total_data\\ICB_varify\\Cancer_Dis\\"
file_path = paste0(folder_path,"Cancer_Dis_cell2location_for_GLM.csv")
cluster = read.csv(file_path,header = TRUE,row.names=1)


df = cluster[, new_order]
df$response <- ifelse(df$response == "R", 1, 0)


# Assuming 'df' is your dataframe and 'response' is the binary outcome variable

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
for (feature in colnames(df)[-1]) { # assuming the first column is the 'response'
  
  # Fit the logistic regression model
  formula <- as.formula(paste('response ~', feature))
  model <- glm(formula, data = df, family = binomial)
  
  # Use broom to tidy the model and extract necessary information
  tidy_model <- broom::tidy(model)
  
  # Add results to the result dataframe
  results_df <- rbind(results_df, data.frame(
    Feature = feature,
    Estimate = tidy_model$estimate[2],
    Std.Error = tidy_model$std.error[2],
    Statistic = tidy_model$statistic[2],
    P.Value = tidy_model$p.value[2],
    OR = exp(tidy_model$estimate[2]), # Odds Ratio
    LowerCI95 = exp(tidy_model$estimate[2] - 1.96 * tidy_model$std.error[2]), # Lower 95% CI for OR
    UpperCI95 = exp(tidy_model$estimate[2] + 1.96 * tidy_model$std.error[2])  # Upper 95% CI for OR
  ))
}

# Order results by P.Value or OR
results_df <- results_df[order(results_df$P.Value), ]


# Create the forest plot
forestplot(labeltext = results_df$Feature,
           mean = results_df$OR,
           lower = results_df$LowerCI95,
           upper = results_df$UpperCI95,
           is.summary=c(TRUE, rep(FALSE, nrow(results_df)-1)),
           xlog = TRUE, # Log scale for Odds Ratio
           title = "Forest plot of Odds Ratios",
           clip = c(0.1, 10), # Adjust these values as needed for your data
           zero = 1, # Zero line for Odds Ratio of 1
           xlab = "Odds Ratio (log scale)",
           #cex.lab=.75,
           col=fpColors(box="royalblue", line="darkblue", summary="royalblue"))










########################## GSE136961 ############################
library(pheatmap)
library(broom)
library(forestplot)
library(rstanarm)


new_order <- c("response",
               "q05cell_abundance_w_sf_B_IGHM_04", "q05cell_abundance_w_sf_B_CD74_1_00",
               "q05cell_abundance_w_sf_B_CD69_1_01", "q05cell_abundance_w_sf_B_MHC_II_05",
               "q05cell_abundance_w_sf_B_GC_11", "q05cell_abundance_w_sf_B_ISG_08",
               "q05cell_abundance_w_sf_Plasma_cell_01_02", "q05cell_abundance_w_sf_Plasma_cell_02_03",
               "q05cell_abundance_w_sf_CD8_T_naive_01", "q05cell_abundance_w_sf_CD8_T_effector_memory_01",
               "q05cell_abundance_w_sf_CD8_T_cytotoxic_02", "q05cell_abundance_w_sf_CD8_T_late_exhausted_01",
               "q05cell_abundance_w_sf_CD8_T_ISG",
               "q05cell_abundance_w_sf_CD4_T_fh_01_00", "q05cell_abundance_w_sf_CD4_T_naive_1_01",
               "q05cell_abundance_w_sf_CD4_T_naive_2_03",
               "q05cell_abundance_w_sf_CD4_T_early_exhausted_06", "q05cell_abundance_w_sf_CD4_T_late_exhausted_08",
               "q05cell_abundance_w_sf_CD4_T_ISG", "q05cell_abundance_w_sf_cDC_1_15",
               "q05cell_abundance_w_sf_cDC_2_06", "q05cell_abundance_w_sf_DC_LAMP3_11", "q05cell_abundance_w_sf_pDC_10",
               "q05cell_abundance_w_sf_CD4_T_reg_naive_07",
               "q05cell_abundance_w_sf_CD4_T_reg_1_04", "q05cell_abundance_w_sf_CD4_T_reg_2_05",
               "q05cell_abundance_w_sf_CD4_T_reg_3_10", "q05cell_abundance_w_sf_CD4_T_reg_ISG",
               "q05cell_abundance_w_sf_CD4_T_reg_exhausted_01_02", "q05cell_abundance_w_sf_CD4_T_reg_exhausted_02_09",
               "q05cell_abundance_w_sf_Cleaning_macrophage_13",
               "q05cell_abundance_w_sf_Cycling_B_plamsa_cell_09", "q05cell_abundance_w_sf_Cycling_CD4_T",
               "q05cell_abundance_w_sf_Cycling_CD8_T", "q05cell_abundance_w_sf_Cycling_NK",
               "q05cell_abundance_w_sf_Cycling_myeloid_cell_14", 
               "q05cell_abundance_w_sf_Endothelial_17", "q05cell_abundance_w_sf_Fibroblast_13",
               "q05cell_abundance_w_sf_ILC_06", "q05cell_abundance_w_sf_Lymphatic_endothelial_cell_21",
               "q05cell_abundance_w_sf_M1_S100A8_07", "q05cell_abundance_w_sf_M2_COL1A1_04",
               "q05cell_abundance_w_sf_M2_CXCL10_12", "q05cell_abundance_w_sf_M2_MARCO_05",
               "q05cell_abundance_w_sf_M2_MMP9_08", "q05cell_abundance_w_sf_M2_SELENOP_02",
               "q05cell_abundance_w_sf_M2_STAB1_09", "q05cell_abundance_w_sf_Mast_cell_00",
               "q05cell_abundance_w_sf_Neutrophil_1_01",
               "q05cell_abundance_w_sf_Neutrophil_1_03"
) 



###read cluster
folder_path = "G:\\TLSdata\\TLS_total_data\\ICB_varify\\GSE136961\\"
file_path = paste0(folder_path,"GSE136961_RNA_cell2location_for_GLM.csv")
cluster = read.csv(file_path,header = TRUE,row.names=1)

data = cluster[, -1]

data = t(data)

#####??׼####
standardize <- function(x) {
  rowmean <- apply(x, 1, mean)
  rowsd <- apply(x, 1, sd)  
  rv <- sweep(x, 1, rowmean,"-")  #????��-??ֵ
  rv <- sweep(rv, 1, rowsd, "/")  #?ٳ??Ա?׼??
  return(rv)
}
######
normalize <- function(x) {
  rowmin <- apply(x, 1, min)
  rowmax <- apply(x, 1, max)
  rowmax.min <- rowmax-rowmin
  rv <- sweep(x, 1, rowmin,"-")  #????��-??ֵ
  rv <- sweep(rv, 1, rowmax.min, "/")  #?ٳ??Ա?׼??
  return(rv)
}




data = normalize(data)
data.standardize = standardize(data)

data.standardize = t(data.standardize)
data.standardize = as.data.frame(data.standardize)


data.standardize <- as.data.frame(data.standardize)
data.standardize$response <- cluster$response  # Add the response column back
data.standardize$response <- ifelse(data.standardize$response == "R", 1, 0)
