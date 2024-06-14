library(tidyverse)
library(ggplot2)
library(ggpubr)
library(emmeans)
library(viridisLite)
library(RColorBrewer)
library(cowplot)
library(gridExtra)

# Load data
raw_table <- read.csv('XXX/spatial_obs_for_R.csv', header = TRUE)
rownames(raw_table) <- as.vector(as.matrix(raw_table[, 1]))

# Function to create gene dotplot odds ratio
gene_dotplot_odds_ratio <- function(table, gene) {
  scrna_cell_tbl <- raw_table
  scrna_cell_tbl$barcode <- scrna_cell_tbl[, 1]
  scrna_cell_tbl <- scrna_cell_tbl[, -1]
  
  # Add additional columns
  scrna_cell_tbl <- scrna_cell_tbl %>%
    mutate(
      patient_id = sample,
      tumor_type = group,
      cell_id = barcode,
      cell_type = scrna_cell_tbl[, gene]
    )
  
  scrna_sample_tbl <- scrna_cell_tbl %>%
    group_by(patient_id, tumor_type, cell_type) %>%
    summarize(cell_type_count = n_distinct(cell_id)) %>%
    mutate(
      other_cell_type_count = sum(cell_type_count) - cell_type_count,
      total_cell_type_count = sum(cell_type_count),
      cell_type_proportion = cell_type_count / total_cell_type_count * 100
    ) %>%
    ungroup()
  
  # Define model and data
  data <- list(scRNA = scrna_sample_tbl)
  formulas <- list(scRNA = cbind(cell_type_count, other_cell_type_count) ~ cell_type * tumor_type)
  models <- list(scRNA = glm(formula = formulas[["scRNA"]], family = binomial(link = 'logit'), data = data[["scRNA"]]))
  
  TLS_emm <- list(scRNA = emmeans(models[["scRNA"]], specs = eff ~ tumor_type | cell_type * tumor_type, type = "response", tran = "logit", adjust = "bonferroni"))
  TLS_odds <- TLS_emm[["scRNA"]]$contrasts %>%
    rbind() %>%
    as.data.frame()
  
  plot_data <- bind_rows(TLS_odds %>% add_column(contrast_variable = "Signature")) %>%
    mutate(
      p.value = ifelse(p.value < 1E-200, 1E-200, p.value),
      contrast = str_remove_all(contrast, " effect|\$$|\$$"),
      contrast = ordered(contrast),
      log2.odds.ratio = log2(odds.ratio),
      log2.odds.ratio = ifelse(log2.odds.ratio > 1.5, 1.5, log2.odds.ratio),
      log2.odds.ratio = ifelse(log2.odds.ratio < -1.5, -1.5)
    ) %>%
    filter(!str_detect(contrast, "(HRD-Other) effect|Undetermined|Other"))
  
  # Create odds ratio plot
  p_odds_ratio <- plot_data %>%
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
  
  # Create signature and site annotation plot
  p_signature_site_anno <- plot_data %>%
    distinct(contrast, contrast_variable, .keep_all = TRUE) %>%
    ggplot() +
    geom_tile(aes(x = 0, y = forcats::fct_rev(contrast), fill = contrast)) +
    scale_fill_manual(values = c(scales::muted("blue"), scales::muted("red"), "grey10")) +
    theme(
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      strip.text = element_blank()
    ) +
    guides(fill = "none")
  
  # Define color palette
  color_palette <- brewer.pal(11, "Spectral")
  color_palette <- rep(color_palette, each = ceiling(22 / 11))[1:22]
  
  # Create cluster annotation plot
  p_cluster_anno <- plot_data %>%
    mutate(facet_helper = "") %>%
    distinct(cell_type, facet_helper, .keep_all = TRUE) %>%
    ggplot() +
    geom_tile(aes(cell_type, facet_helper, fill = cell_type)) +
    scale_fill_manual(values = color_palette) +
    theme(
      axis.text.y = element_blank(),
      strip.text = element_blank()
    ) +
    guides(fill = "none")
  
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
  
  p <- plot_grid(
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
  
  return(p)
}

# Generate plots for specific genes
p_CD20 <- gene_dotplot_odds_ratio(raw_table, "MS4A1")
p_CD8 <- gene_dotplot_odds_ratio(raw_table, "CD8A")
p_cancer <- gene_dotplot_odds_ratio(raw_table, "KRT14")
p_CD4 <- gene_dotplot_odds_ratio(raw_table, "CD4")
p_TCF1 <- gene_dotplot_odds_ratio(raw_table, "TCF7")
p_KI67 <- gene_dotplot_odds_ratio(raw_table, "MKI67")
p_CD23 <- gene_dotplot_odds_ratio(raw_table, "FCER2")

# Arrange all plots in a grid
p_total <- grid.arrange(p_cancer, p_CD20, p_CD8, p_CD4, p_CD23, p_TCF1, nrow = 1)

# Save the combined plot
ggsave("XXX/gene_spatial.pdf", plot = p_total, width = 15, height = 8)