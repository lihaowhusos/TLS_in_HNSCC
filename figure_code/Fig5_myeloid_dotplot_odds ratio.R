library(dplyr)       # For data manipulation
library(ggplot2)     # For data visualization
library(forcats)     # For dealing with factors
library(RColorBrewer) # For color palette
library(emmeans)     # For estimated marginal means
library(cowplot)     # For combining plots

## Load data --------------------------------------
cell_type = 'Myeloid'
data_file = paste0('G:/TLSdata/TLS_total_data/Processing data/subset/',cell_type,'/',cell_type,'_for_seurat_metadata.csv')
cell_data_tbl <- read.csv(data_file, header = TRUE)
rownames(cell_data_tbl) <- as.vector(as.matrix(cell_data_tbl[, 1]))

# Prepare barcode column and remove first column (ID column)
cell_data_tbl$barcode <- cell_data_tbl[, 1]
cell_data_tbl <- cell_data_tbl[, -1]

## Join and filter data
cell_data_tbl <- cell_data_tbl %>% 
  mutate(
    patient_id = sample_CD45,          # Mapping sample_CD45 to patient_id
    tumor_category = group,            # Renaming group to tumor_category
    cell_identifier = barcode          # Renaming barcode to cell_identifier
  )

# Summarize data by patient, tumor type, and cell type
cell_summary_tbl <- cell_data_tbl %>%
  group_by(patient_id, tumor_category, cell_type) %>%
  summarize(cell_type_count = n_distinct(cell_identifier)) %>% 
  mutate(
    other_cell_count = sum(cell_type_count) - cell_type_count,  # Calculate the count of other cell types
    total_cell_count = sum(cell_type_count),                    # Calculate the total cell count
    cell_type_fraction = cell_type_count / total_cell_count * 100 # Calculate the fraction of each cell type
  ) %>%
  ungroup()

# Define heatmap layers for consistent styling
heatmap_elements <- list(
  scale_fill_manual(values = c(scales::muted("blue"), scales::muted("red"), "grey10")),
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank(),
        plot.margin = margin(0, 0, 0, 0),
        strip.text = element_blank())
)

## Cell type GLM {.tabset}
data_storage <- list()
model_formulas <- list()
model_list <- list()
emmeans_results <- list()

# Prepare data and models for Generalized Linear Model (GLM)
data_storage[["cellData"]] <- cell_summary_tbl
model_formulas[["cellData"]] <- cbind(cell_type_count, other_cell_count) ~ cell_type * tumor_category
model_list[["cellData"]] <- glm(formula = model_formulas[["cellData"]], family = binomial(link = 'logit'), data = data_storage[["cellData"]])

# Calculate estimated marginal means and pairwise comparisons
emmeans_results[["cellData"]] <- emmeans(model_list[["cellData"]], specs = pairwise ~ tumor_category | cell_type, type = "response", tran = "logit", adjust = "bonferroni")
odds_ratios <- as.data.frame(summary(emmeans_results[["cellData"]])$contrasts)

# Redefine data for further analysis
data_storage[["cellData"]] <- cell_summary_tbl
model_formulas[["cellData"]] <- cbind(cell_type_count, other_cell_count) ~ cell_type * tumor_category
model_list[["cellData"]] <- glm(formula = model_formulas[["cellData"]], family = binomial(link = 'logit'), data = data_storage[["cellData"]])

# Estimate marginal means with a different specification
emmeans_results[["cellData"]] <- emmeans(model_list[["cellData"]], specs = eff ~ tumor_category | cell_type * tumor_category, type = "response", tran = "logit", adjust = "bonferroni")
odds_ratios <- emmeans_results[["cellData"]]$contrasts %>%
  rbind() %>%
  as.data.frame()

# Prepare data for plotting
plot_data <- bind_rows(odds_ratios %>% add_column(contrast_variable = "Signature")) %>%
  mutate(p.value = ifelse(p.value < 1E-50, 1E-50, p.value)) %>%
  mutate(contrast = str_remove(contrast, " effect")) %>%
  mutate(log2.odds.ratio = log2(odds.ratio)) %>%
  mutate(log2.odds.ratio = ifelse(log2.odds.ratio > 0.5, 0.5, log2.odds.ratio),
         log2.odds.ratio = ifelse(log2.odds.ratio < -0.5, -0.5, log2.odds.ratio))

## Custom ordering for factors
cell_type_order <- c('cDC_1_15',
                     'cDC_2_06',
                     'DC_LAMP3_11',
                     'pDC_10',
                     'M1_S100A8_07',
                     'M2_CXCL10_12',
                     'M2_MARCO_05',
                     'M2_STAB1_09',
                     'M2_SELENOP_02',
                     'M2_MMP9_08',
                     'M2_COL1A1_04',
                     'Cleaning_macrophage_13',
                     'Cycling_myeloid_cell_14',
                     'Neutrophil_1_01',
                     'Neutrophil_1_03',
                     'Mast_cell_00')
plot_data$cell_type <- factor(plot_data$cell_type, levels = cell_type_order)
plot_data <- plot_data[order(plot_data$cell_type),]

contrast_order <- c("mature", "immature", "none")
plot_data$contrast <- factor(plot_data$contrast, levels = contrast_order)
plot_data <- plot_data[order(plot_data$contrast),]

# Create the main plot for odds ratios
p_odds_ratio <- ggplot(plot_data, aes(x = cell_type, y = forcats::fct_rev(contrast), color = log2.odds.ratio, size = -log10(p.value))) +
  geom_point() +
  scale_color_gradientn(
    colors = rev(brewer.pal(9, "RdBu")),
    limits = c(-0.5, 0.5),
    labels = c("≤-0.5", 0, "≥0.5"),
    breaks = c(-0.5, 0, 0.5)
  ) +
  scale_size_continuous(
    limits = c(0, 50),
    breaks = c(0, 25, 50),
    labels = c(0, 25, 50)
  ) +
  scale_x_discrete(drop = FALSE) +
  theme(
    axis.line = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    legend.title = element_text(size = 14),
    legend.direction = "vertical",
    legend.position = "top",
    legend.box = "vertical",
    legend.box.just = "left",
    panel.background = element_blank(),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "plain", size = 14),
    strip.background = element_blank(),
    strip.text.x = element_blank()
  ) +
  guides(color = guide_colourbar(title.position="top", title.hjust = 0.5, title.vjust = 0),
         size = guide_legend(title.position="top", title.hjust = 0.5, title.vjust = 0)) +
  labs(title = NULL, x = NULL, y = NULL, color = 'log2(odds ratio)', size = "\n-log10(p value)") +
  facet_grid(. ~ cell_type, scales = "free_x", space = "free")

# Create the annotation for signatures
p_signature_site_anno <- ggplot(plot_data %>% distinct(contrast, .keep_all = TRUE)) +
  geom_tile(aes(x = 0, y = forcats::fct_rev(contrast), fill = contrast)) +
  heatmap_elements +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        strip.text = element_blank()) +
  guides(fill = "none")

# Create the color palette for cell types
color_palette <- brewer.pal(11, "Spectral")
color_palette <- rep(color_palette, each = ceiling(22 / 11))[1:22]

# Create the annotation for cell types
p_cluster_anno <- ggplot(plot_data %>% distinct(cell_type, .keep_all = TRUE) %>% mutate(facet_helper = "")) +
  geom_tile(aes(cell_type, facet_helper, fill = cell_type)) +
  heatmap_elements +
  scale_x_discrete(drop = FALSE) +
  facet_grid(~facet_helper, scales = "free", space = "free") +
  theme(axis.text.y = element_blank(), strip.text = element_blank()) +
  guides(fill = "none") + 
  scale_fill_manual(values = color_palette)

# Combine the plots
p_left <- plot_grid(
  p_signature_site_anno + theme(plot.margin = margin(t = 6, r = 0, b = 0, l = 6)),
  ggdraw(),
  nrow = 2,
  align = "hv",
  axis = "tbr",
  rel_heights = c(1, 0.7),
  rel_widths = c(0.75, 1)
)

p_right <- plot_grid(
  p_odds_ratio + theme(legend.position = "none", plot.margin = margin(t = 0, r = 6, b = 0, l = 0)),
  p_cluster_anno + theme(plot.margin = margin(t = 0, r = 6, b = 6, l = 6)),
  nrow = 2,
  align = "v",
  axis = "lr",
  rel_heights = c(1, 0.7),
  rel_widths = c(0.75, 1)
)

p_cell_glm_result <- plot_grid(
  plot_grid(
    p_left,
    p_right,
    nrow = 1,
    align = "hv",
    axis = "tblr",
    rel_widths = c(0.5, 1)
  ),
  get_legend(p_odds_ratio),
  ncol = 1,
  rel_heights = c(1, 0.15)
)

p_cell_glm_result

# Save the final combined plot
ggsave("path/to/figures/cluster_dot_plot.pdf", plot = p_cell_glm_result)

# Convert the data to data.table format for further use if needed
data.table(plot_data)
