library(pheatmap)
library(broom)
library(forestplot)
library(rstanarm)
library(bayesplot)
library(rstan)

# Define the new order for columns
new_order <- c(
  "response",
  "B_IGHM_04", "B_CD74_1_00", "B_CD69_1_01", "B_MHC_II_05",
  "B_GC_11", "B_ISG_08", "Plasma_cell_01_02", "Plasma_cell_02_03",
  "CD8_T_naive_01", "CD8_T_effector_memory_01", "CD8_T_cytotoxic_02",
  "CD8_T_late_exhausted_01", "CD8_T_ISG", "CD4_T_fh_01_00",
  "CD4_T_naive_1_01", "CD4_T_naive_2_03", "CD4_T_early_exhausted_06",
  "CD4_T_late_exhausted_08", "CD4_T_ISG", "cDC_1_15", "cDC_2_06",
  "DC_LAMP3_11", "pDC_10", "CD4_T_reg_naive_07", "CD4_T_reg_1_04",
  "CD4_T_reg_2_05", "CD4_T_reg_3_10", "CD4_T_reg_ISG",
  "CD4_T_reg_exhausted_01_02", "CD4_T_reg_exhausted_02_09",
  "Cleaning_macrophage_13", "Cycling_B_plamsa_cell_09", "Cycling_CD4_T",
  "Cycling_CD8_T", "Cycling_NK", "Cycling_myeloid_cell_14",
  "Endothelial_17", "Fibroblast_13", "ILC_06", "Lymphatic_endothelial_cell_21",
  "M1_S100A8_07", "M2_COL1A1_04", "M2_CXCL10_12", "M2_MARCO_05",
  "M2_MMP9_08", "M2_SELENOP_02", "M2_STAB1_09", "Mast_cell_00",
  "Neutrophil_1_01", "Neutrophil_1_03"
) 

# Load and process data
folder_path = "XXX"
file_path = paste0(folder_path, "GSE136961_RNA_cell2location_for_GLM.csv")
cluster = read.csv(file_path, header = TRUE, row.names = 1)

# Remove prefix from column names
colnames(cluster) <- gsub(colnames(cluster))

data = cluster[, -1]
data = t(data)

# Standardize function
standardize <- function(x) {
  rowmean <- apply(x, 1, mean)
  rowsd <- apply(x, 1, sd)  
  rv <- sweep(x, 1, rowmean, "-") 
  rv <- sweep(rv, 1, rowsd, "/") 
  return(rv)
}

# Normalize function
normalize <- function(x) {
  rowmin <- apply(x, 1, min)
  rowmax <- apply(x, 1, max)
  rowmax.min <- rowmax - rowmin
  rv <- sweep(x, 1, rowmin, "-") 
  rv <- sweep(rv, 1, rowmax.min, "/") 
  return(rv)
}

data = normalize(data)
data.standardize = standardize(data)
data.standardize = t(data.standardize)
data.standardize = as.data.frame(data.standardize)

data.standardize$response <- cluster$response 
data.standardize$response <- ifelse(data.standardize$response == "R", 1, 0)

features <- paste(colnames(data.standardize)[-ncol(data.standardize)], collapse = " + ")
formula <- as.formula(paste('response ~', features))

model <- stan_glm(formula, data = data.standardize, family = binomial(link = "logit"))

bayesplot::color_scheme_set("viridis")
plot(model)
plot(model, plotfun = "areas", prob = 0.9)

pp_check(model)
summary(model)
monitor(model)

# Load and process second dataset
file_path = paste0(folder_path, "Cancer_Dis_cell2location_for_GLM.csv")
cluster = read.csv(file_path, header = TRUE, row.names = 1)

# Remove prefix from column names
colnames(cluster) <- gsub(colnames(cluster))

df = cluster[, new_order]
df$response <- ifelse(df$response == "R", 1, 0)

features <- paste(colnames(df)[-1], collapse = " + ")
formula <- as.formula(paste('response ~', features))

model <- stan_glm(formula, data = df, family = binomial(link = "logit"))

bayesplot::color_scheme_set("brightblue")
plot(model)
posterior_vs_prior(model)
pp_check(model, plotfun = "error_binned")

pp_check(model)
summary(model)
monitor(model)