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

# 设定宽泛的先验分布
# 以下的参数需要根据你自己的数据进行调整
#prior_for_intercept <- normal(0, 10)   # 截距的宽泛正态先验分布
#prior_for_coefficients <- student_t(5, 0, 10) # 斜率的宽泛学生t先验分布

library(rstanarm)
library(broom.mixed)



features <- paste(colnames(data.standardize)[ -ncol(data.standardize)], collapse = " + ")
formula <- as.formula(paste('response ~', features))


model <- stan_glm(formula, data = data.standardize, family = binomial(link = "logit"),
                  #prior_intercept = normal(0, 1),
                  #prior = student_t(5, 0, 2.5)
                  )


bayesplot::color_scheme_set("viridis")
plot(model)


plot(model, plotfun = "areas", prob = 0.9)





library(bayesplot)
pp_check(model)

summary(model)

library(rstan) 
monitor(model)





########################################################## Cancer_dis ##############################################
library(pheatmap)

###read cluster
folder_path = "G:\\TLSdata\\TLS_total_data\\ICB_varify\\Cancer_Dis\\"
file_path = paste0(folder_path,"Cancer_Dis_cell2location_for_GLM.csv")
cluster = read.csv(file_path,header = TRUE,row.names=1)


df = cluster[, new_order]
df$response <- ifelse(df$response == "R", 1, 0)

# 设定宽泛的先验分布
# 以下的参数需要根据你自己的数据进行调整
#prior_for_intercept <- normal(0, 10)   # 截距的宽泛正态先验分布
#prior_for_coefficients <- student_t(5, 0, 10) # 斜率的宽泛学生t先验分布

library(rstanarm)
library(broom.mixed)

# Model fitting loop

results_list <- vector("list", ncol(df) - 1)
names(results_list) <- colnames(df)[-1]


features <- paste(colnames(df)[-1], collapse = " + ")
formula <- as.formula(paste('response ~', features))


model <- stan_glm(formula, data = df, family = binomial(link = "logit")
)


bayesplot::color_scheme_set("brightblue")
plot(model)





posterior_vs_prior(model)
pp_check(model, plotfun = "error_binned")



library(bayesplot)
pp_check(model)

summary(model)

library(rstan) 
monitor(model)
