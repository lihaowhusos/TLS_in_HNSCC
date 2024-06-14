library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ggridges)
library(emmeans)
library(viridis)
library(RColorBrewer)
library(cowplot)

# Load data
scrna_cell_tbl <- read.csv('XXX/adata.obs.csv', header = TRUE)
rownames(scrna_cell_tbl) <- as.vector(as.matrix(scrna_cell_tbl[,1]))

# Data tidying
scrna_cell_tbl$barcode <- scrna_cell_tbl[,1]
scrna_cell_tbl <- scrna_cell_tbl[,-1]

scrna_cell_tbl <- scrna_cell_tbl %>%
  filter(tissue_type == "tumor") %>%
  mutate(cell_type = cluster_1,
         patient_id = sample,
         tumor_type = group,
         cell_id = barcode) 

scrna_sample_tbl <- scrna_cell_tbl %>%
  group_by(patient_id, tumor_type, cell_type) %>%
  summarize(cell_type_count = n_distinct(cell_id)) %>%
  mutate(other_cell_type_count = sum(cell_type_count) - cell_type_count,
         total_cell_type_count = sum(cell_type_count),
         cell_type_proportion = cell_type_count/total_cell_type_count*100) %>%
  ungroup()

# Heatmap settings
heatmap_layers <- list(
  scale_fill_gradient2(low = scales::muted("blue"), high = scales::muted("red"), 
                       na.value = "grey10"),
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank(),
        plot.margin = margin(0, 0, 0, 0),
        strip.text = element_blank())
)

# GLM model setup and analysis
data <- list()
formulas <- list()
models <- list()
TLS_emm <- list()

data[["scRNA"]] <- scrna_sample_tbl
formulas[["scRNA"]] <- cbind(cell_type_count, other_cell_type_count) ~ cell_type * tumor_type
models[["scRNA"]] <- glm(formula = formulas[["scRNA"]], family = binomial(link = 'logit'), data = data[["scRNA"]])
TLS_emm[["scRNA"]] <- emmeans(models[["scRNA"]], specs = eff ~ tumor_type | cell_type * tumor_type, type = "response", tran = "logit", adjust = "bonferroni")
TLS_odds <- TLS_emm[["scRNA"]]$contrasts %>%
  rbind() %>%
  as.data.frame()

# Plot data preparation
plot_data <- bind_rows(TLS_odds %>% add_column(contrast_variable = "Signature")) %>%
  mutate(p.value = ifelse(p.value < 1E-200, 1E-200, p.value)) %>%
  mutate(contrast = str_remove(contrast, " effect")) %>%
  mutate(contrast = str_remove(contrast, "\\(")) %>%
  mutate(contrast = str_remove(contrast, "\\)")) %>%
  mutate(contrast = ordered(contrast)) %>%
  mutate(log2.odds.ratio = log2(odds.ratio)) %>%
  mutate(log2.odds.ratio = ifelse(log2.odds.ratio > 1.5, 1.5, log2.odds.ratio),
         log2.odds.ratio = ifelse(log2.odds.ratio < -1.5, -1.5, log2.odds.ratio))

# Custom order for factors
custom_order <- c("T_NK", 
                  "B", 
                  "Plasma_cell",
                  "Myeloid_cell",
                  "Neutrophil",
                  "Mast_cell",
                  "Endothelial_cell",
                  "Lymphatic_endothelial_cell",
                  "Fibroblast",
                  "Pericytes",
                  "Cancer_cell")
plot_data$cell_type <- factor(plot_data$cell_type, levels = custom_order)
plot_data <- plot_data[order(plot_data$cell_type),]

custom_order <- c("none", "immature", "mature")
plot_data$contrast <- factor(plot_data$contrast, levels = custom_order)
plot_data <- plot_data[order(plot_data$contrast),]

# Plot generation
p_odds_ratio <- plot_data %>%
  filter(SE < 0.15) %>%
  ggplot(aes(x = cell_type, y = forcats::fct_rev(contrast), color = log2.odds.ratio, size = -log10(p.value))) +
  geom_point() +
  scale_color_gradientn(
    colors = viridis(9, option = "D"),
    limits = c(-1.5, 1.5),
    labels = c("≤-1.5", 0, "≥1.5"),
    breaks = c(-1.5, 0, 1.5)
  ) +
  scale_size_continuous(
    limits = c(0, 200),
    breaks = c(0, 100, 200),
    labels = c(0, 100, 200)
  ) +
  scale_x_discrete(drop=FALSE) +
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
  facet_grid(contrast_variable ~ ., scales = "free_y", space = "free")

p_signature_site_anno <- plot_data %>% 
  distinct(contrast, contrast_variable, .keep_all = T) %>% 
  ggplot() +
  geom_tile(aes(x=0, y=forcats::fct_rev(contrast), fill = contrast)) +
  heatmap_layers + 
  facet_grid(contrast_variable~., scales = "free", space = "free") +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        strip.text = element_blank()) +
  guides(fill = "none")

color_palette <- brewer.pal(11, "Spectral") 
color_palette <- rep(color_palette, each = ceiling(22 / 11))[1:22]

p_cluster_anno <- plot_data %>% 
  mutate(facet_helper = "") %>% 
  distinct(cell_type, facet_helper, .keep_all = T) %>% 
  ggplot() +
  geom_tile(aes(cell_type, facet_helper, fill = cell_type)) +
  heatmap_layers + 
  scale_x_discrete(drop=FALSE) +
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

# Save plot to PDF
ggsave("XXX/cluster_1_dot_plot.pdf", plot = p_scrna_odds_ratio)

# Output data table
data.table(plot_data)