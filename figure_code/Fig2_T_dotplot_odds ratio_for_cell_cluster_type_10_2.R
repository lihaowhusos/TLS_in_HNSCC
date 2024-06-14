# Load necessary libraries
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ggridges)
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

# Load data
cell_type <- 'T_NK'
data_file <- paste0('XXX/', cell_type, '_for_seurat_metadata.csv')
out_file <- paste0('XXX/', cell_type, '_cluster_2_dot_plot.pdf')

scrna_cell_tbl <- read.csv(data_file, header = TRUE)

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

rownames(scrna_cell_tbl) <- as.vector(as.matrix(scrna_cell_tbl[,1]))
scrna_cell_tbl$barcode <- scrna_cell_tbl[,1]
scrna_cell_tbl <- scrna_cell_tbl[,-1]

# Join data
scrna_cell_tbl <- scrna_cell_tbl %>% 
  mutate(cell_type = cluster_2, 
         patient_id = sample_CD45, 
         tumor_type = group, 
         cell_id = barcode) 

scrna_sample_tbl <- scrna_cell_tbl %>%
  group_by(patient_id, tumor_type, cell_type) %>%
  summarize(cell_type_count = n_distinct(cell_id)) %>% 
  mutate(
    other_cell_type_count = sum(cell_type_count) - cell_type_count,
    total_cell_type_count = sum(cell_type_count),
    cell_type_proportion = cell_type_count / total_cell_type_count * 100
  ) %>%
  ungroup 

heatmap_layers <- list(
  scale_fill_gradient2(low = scales::muted("blue"), high = scales::muted("red"), na.value = "grey10"),
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank(),
        plot.margin = margin(0, 0, 0, 0),
        strip.text = element_blank())
)

data <- list()
formulas <- list()
models <- list()
TLS_emm <- list()

# Cell type GLM
data[["scRNA"]] <- scrna_sample_tbl
formulas[["scRNA"]] <- cbind(cell_type_count, other_cell_type_count) ~ cell_type * tumor_type
models[["scRNA"]] <- glm(formula = formulas[["scRNA"]], family = binomial(link = 'logit'), data = data[["scRNA"]])

TLS_emm[["scRNA"]] <- emmeans(models[["scRNA"]], specs = eff ~ tumor_type | cell_type * tumor_type, type = "response", tran = "logit", adjust = "bonferroni")
TLS_odds <- TLS_emm[["scRNA"]]$contrasts %>%
  rbind() %>%
  as.data.frame()

# Signature + site contrasts
plot_data <- bind_rows(TLS_odds %>% add_column(contrast_variable = "Signature")) %>%
  mutate(p.value = ifelse(p.value < 1E-50, 1E-50, p.value)) %>%
  mutate(contrast = str_remove(contrast, " effect")) %>%
  mutate(contrast = str_remove(contrast, "\\(")) %>%
  mutate(contrast = str_remove(contrast, "\\)")) %>%
  mutate(contrast = ordered(contrast)) %>%
  filter(!str_detect(contrast, "(HRD-Other) effect"),
         !str_detect(contrast, "Undetermined")) %>%
  filter(!str_detect(contrast, "Other")) %>%
  mutate(log2.odds.ratio = log2(odds.ratio)) %>%
  mutate(log2.odds.ratio = ifelse(log2.odds.ratio > 0.5, 0.5, log2.odds.ratio),
         log2.odds.ratio = ifelse(log2.odds.ratio < -0.5, -0.5, log2.odds.ratio))

# Rank
custom_order <- c("CD8_T_effector_memory_01", 
                  "CD8_T_naive_01", 
                  "CD8_T_cytotoxic_02",
                  "CD8_T_late_exhausted_01",
                  "CD8_T_ISG",
                  "Cycling_CD4_T",
                  "Cycling_CD8_T",
                  "Cycling_NK",
                  "CD4_T_fh_01_00",
                  "CD4_T_naive_1_01",
                  "CD4_T_naive_2_03",
                  "CD4_T_early_exhausted_06",
                  "CD4_T_late_exhausted_08",
                  "CD4_T_ISG",
                  "CD4_T_reg_naive_07",
                  "CD4_T_reg_1_04",
                  "CD4_T_reg_2_05",
                  "CD4_T_reg_3_10",
                  "CD4_T_reg_exhausted_01_02",
                  "CD4_T_reg_exhausted_02_09",
                  "CD4_T_reg_ISG",
                  "NKT_1_00",
                  "NKT_2_02",
                  "NK_cytotoxic_04",
                  "NK_reg_CRTAM_03",
                  "NK_reg_KRT86_01",
                  "dg_T_05",
                  "ILC_06")
plot_data$cell_type <- factor(plot_data$cell_type, levels = custom_order)
plot_data <- plot_data[order(plot_data$cell_type),]

custom_order <- c("mature", "immature", "none")
plot_data$contrast <- factor(plot_data$contrast, levels = custom_order)
plot_data <- plot_data[order(plot_data$contrast),]

p_odds_ratio <- plot_data %>%
  ggplot(aes(x = cell_type, y = forcats::fct_rev(contrast), color = log2.odds.ratio, size = -log10(p.value))) +
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
  guides(color = guide_colourbar(title.position = "top", title.hjust = 0.5, title.vjust = 0),
         size = guide_legend(title.position = "top", title.hjust = 0.5, title.vjust = 0)) +
  labs(title = NULL, x = NULL, y = NULL, color = 'log2(odds ratio)', size = "\n-log10(p value)") +
  facet_grid(contrast_variable ~ ., scales = "free_y", space = "free")

p_signature_site_anno <- plot_data %>% 
  distinct(contrast, contrast_variable, .keep_all = TRUE) %>% 
  ggplot() +
  geom_tile(aes(x = 0, y = forcats::fct_rev(contrast), fill = contrast)) +
  heatmap_layers + 
  facet_grid(contrast_variable ~ ., scales = "free", space = "free") +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        strip.text = element_blank()) +
  guides(fill = "none")

# Create color palette
color_palette <- brewer.pal(11, "Spectral")
color_palette <- rep(color_palette, each = 3)

p_cluster_anno <- plot_data %>% 
  mutate(facet_helper = "") %>% 
  distinct(cell_type, facet_helper, .keep_all = TRUE) %>% 
  ggplot() +
  geom_tile(aes(cell_type, facet_helper, fill = cell_type)) +
  heatmap_layers + 
  scale_x_discrete(drop = FALSE) +
  facet_grid(~facet_helper, scales = "free", space = "free") +
  theme(axis.text.y = element_blank(),
        strip.text = element_blank()) +
  guides(fill = "none") + 
  scale_fill_manual(values = color_palette)

p_left <-
  plot_grid(
    p_signature_site_anno + theme(plot.margin = margin(t = 6, r = 0, b = 0, l = 6)),
    ggdraw(),
    nrow = 2,
    align = "hv",
    axis = "tbr",
    rel_heights = c(1, 0.7),
    rel_widths = c(0.75, 1)
  )

p_right <-
  plot_grid(
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

p_scrna_odds_ratio

ggsave(out_file, plot = p_scrna_odds_ratio, width = 10, height = 4)

data.table(plot_data)