library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ggridges)
library(data.table)
library(tidytext)
library(cowplot)
library(ggthemes)
library(ggrepel)
library(corrplot)
library(viridisLite)
library(grid)
library(RColorBrewer)
library(lme4)
library(emmeans)
library(openxlsx)

# Load data
scrna_cell_tbl <- read.csv('XXX/Myeloid_for_seurat_metadata.csv', header = TRUE)
rownames(scrna_cell_tbl) <- as.vector(as.matrix(scrna_cell_tbl[,1]))

# Re-group
scrna_cell_tbl$group <- recode(scrna_cell_tbl$sample, 
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
                               'szj106771' = 'mature')

scrna_cell_tbl$barcode <- scrna_cell_tbl[,1]
scrna_cell_tbl <- scrna_cell_tbl[,-1]

# Join data
scrna_cell_tbl <- scrna_cell_tbl %>% 
  mutate(cell_type = cell_type, 
         patient_id = sample_CD45, 
         tumor_type = group, 
         cell_id = barcode) 

scrna_sample_tbl <- scrna_cell_tbl %>%
  group_by(patient_id, tumor_type, cell_type) %>%
  summarize(cell_type_count = n_distinct(cell_id)) %>% 
  mutate(other_cell_type_count = sum(cell_type_count) - cell_type_count,
         total_cell_type_count = sum(cell_type_count),
         cell_type_proportion = cell_type_count / total_cell_type_count * 100) %>%
  ungroup()

# Define heatmap layers
heatmap_layers <- list(
  scale_fill_manual(values = c(scales::muted("blue"), scales::muted("red"), "grey10")),
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank(),
        plot.margin = margin(0, 0, 0, 0),
        strip.text = element_blank())
)

# Initialize lists
data <- list()
formulas <- list()
models <- list()
TLS_emm <- list()

# Cell type GLM
data[["scRNA"]] <- scrna_sample_tbl
formulas[["scRNA"]] <- cbind(cell_type_count, other_cell_type_count) ~ cell_type * tumor_type
models[["scRNA"]] <- glm(formula = formulas[["scRNA"]], family = binomial(link = 'logit'), data = data[["scRNA"]])

TLS_emm[["scRNA"]] <- emmeans(models[["scRNA"]], specs = eff ~ tumor_type | cell_type * tumor_type, type = "response", tran = "logit", adjust = "bonferroni")
TLS_odds <- TLS_emm[["scRNA"]][["contrasts"]] %>%
  rbind() %>%
  as.data.frame()

# Signature + site contrasts
plot_data <- bind_rows(TLS_odds %>% add_column(contrast_variable = "Signature")) %>%
  mutate(p.value = ifelse(p.value < 1E-50, 1E-50, p.value),
         contrast = str_remove_all(contrast, " effect|\$$|\$$"),
         contrast = ordered(contrast)) %>%
  filter(!str_detect(contrast, "(HRD-Other) effect|Undetermined|Other")) %>%
  mutate(log2.odds.ratio = log2(odds.ratio),
         log2.odds.ratio = ifelse(log2.odds.ratio > 0.5, 0.5, log2.odds.ratio),
         log2.odds.ratio = ifelse(log2.odds.ratio < -0.5, -0.5, log2.odds.ratio))

# Custom order for cell types
custom_order <- c('cDC_1_15', 'cDC_2_06', 'DC_LAMP3_11', 'pDC_10', 'M1_S100A8_07', 
                  'M2_CXCL10_12', 'M2_MARCO_05', 'M2_STAB1_09', 'M2_SELENOP_02', 
                  'M2_MMP9_08', 'M2_COL1A1_04', 'Cleaning_macrophage_13', 
                  'Cycling_myeloid_cell_14', 'Neutrophil_1_01', 'Neutrophil_1_03', 
                  'Mast_cell_00')
plot_data$cell_type <- factor(plot_data$cell_type, levels = custom_order)
plot_data <- plot_data[order(plot_data$cell_type),]

# Plot odds ratio
p_odds_ratio <- plot_data %>%
  ggplot(aes(x = cell_type, y = forcats::fct_rev(contrast), color = log2.odds.ratio, size = -log10(p.value))) +
  geom_point() +
  scale_color_gradientn(colors = rev(brewer.pal(9, "RdBu")), limits = c(-0.5, 0.5), labels = c("≤-0.5", 0, "≥0.5"), breaks = c(-0.5, 0, 0.5)) +
  scale_size_continuous(limits = c(0, 50), breaks = c(0, 25, 50), labels = c(0, 25, 50)) +
  scale_x_discrete(drop = FALSE) +
  theme(axis.line = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), 
        legend.title = element_text(size = 14), legend.direction = "vertical", legend.position = "top", legend.box = "vertical", 
        legend.box.just = "left", panel.background = element_blank(), panel.border = element_blank(), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), plot.title = element_text(face = "plain", size = 14), strip.background = element_blank(), 
        strip.text.x = element_blank()) +
  guides(color = guide_colourbar(title.position = "top", title.hjust = 0.5, title.vjust = 0), size = guide_legend(title.position = "top", title.hjust = 0.5, title.vjust = 0)) +
  labs(title = NULL, x = NULL, y = NULL, color = 'log2(odds ratio)', size = "\n-log10(p value)") +
  facet_grid(contrast_variable ~ ., scales = "free_y", space = "free")

# Plot signature site annotation
p_signature_site_anno <- plot_data %>% 
  distinct(contrast, contrast_variable, .keep_all = TRUE) %>% 
  ggplot() +
  geom_tile(aes(x = 0, y = forcats::fct_rev(contrast), fill = contrast)) +
  heatmap_layers + 
  facet_grid(contrast_variable ~ ., scales = "free", space = "free") +
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(), strip.text = element_blank()) +
  guides(fill = "none")

# Color palette for cluster annotation
color_palette <- brewer.pal(11, "Spectral")
color_palette <- rep(color_palette, each = ceiling(22 / 11))[1:22]

# Plot cluster annotation
p_cluster_anno <- plot_data %>% 
  mutate(facet_helper = "") %>% 
  distinct(cell_type, facet_helper, .keep_all = TRUE) %>% 
  ggplot() +
  geom_tile(aes(cell_type, facet_helper, fill = cell_type)) +
  heatmap_layers + 
  scale_x_discrete(drop = FALSE) +
  facet_grid(~facet_helper, scales = "free", space = "free") +
  theme(axis.text.y = element_blank(), strip.text = element_blank()) +
  guides(fill = "none") + 
  scale_fill_manual(values = color_palette)

# Combine plots
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

p_scrna_odds_ratio <- plot_grid(
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

# Save plot
ggsave("XXX/Myeloid_cluster_dot_plot.pdf", plot = p_scrna_odds_ratio, width = 7)