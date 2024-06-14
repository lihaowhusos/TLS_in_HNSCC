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
library(gridExtra)

## load data --------------------------------------

raw_table <- read.csv('G:/TLSdata/spatial/table/spatial_obs_for_R.csv',header = TRUE)
rownames(raw_table)=as.vector(as.matrix(raw_table[,1]))




gene_dotplot_odds_ratio <-function(table,gene){

  scrna_cell_tbl = raw_table
  scrna_cell_tbl$barcode = scrna_cell_tbl[,1]
  scrna_cell_tbl = scrna_cell_tbl[,-1]

## join data
  scrna_cell_tbl <- scrna_cell_tbl
  scrna_cell_tbl[,"cell_type"] <- scrna_cell_tbl[,gene]
  scrna_cell_tbl <- scrna_cell_tbl %>%
    mutate(patient_id = sample) %>% 
    mutate(tumor_type = group) %>% 
    mutate(cell_id = barcode) 

  scrna_sample_tbl <- scrna_cell_tbl
  scrna_sample_tbl[,"cell_type"] <- scrna_sample_tbl[,gene]
  scrna_sample_tbl <- scrna_sample_tbl %>%
    group_by(patient_id, tumor_type, cell_type) %>%
    summarize(cell_type_count = n_distinct(cell_id)) %>% 
    mutate(
      other_cell_type_count = sum(cell_type_count) - cell_type_count,
      total_cell_type_count = sum(cell_type_count),
      cell_type_proportion = cell_type_count/total_cell_type_count*100
   ) %>%
   ungroup 



#####GPT
heatmap_layers <- list(
  scale_fill_manual(values = c(scales::muted("blue"), scales::muted("red"), "grey10")),
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
consensus_signature_emm <- list()
site_emm <- list()
sort_emm <- list()
TLS_emm <- list()
sort_signature_emm <- list()
sort_site_emm <- list()




## Cell type GLM {.tabset}
data[["scRNA"]] <- scrna_sample_tbl
formulas[["scRNA"]] <-cbind(cell_type_count, other_cell_type_count) ~ cell_type * tumor_type
models[["scRNA"]] <- glm(formula = formulas[["scRNA"]], family = binomial(link = 'logit'), data = data[["scRNA"]])


TLS_emm[["scRNA"]] <- emmeans(models[["scRNA"]], specs = eff ~ tumor_type | cell_type * tumor_type, type = "response", tran = "logit", adjust = "bonferroni")
TLS_odds <- TLS_emm[["scRNA"]]$contrasts %>%
  rbind() %>%
  as.data.frame()


## Signature + site contrasts

plot_data <- bind_rows(TLS_odds %>% add_column(contrast_variable = "Signature")) %>%
  mutate(p.value = ifelse(p.value < 1E-200, 1E-200, p.value)) %>%
  mutate(contrast = str_remove(contrast, " effect")) %>%
  mutate(contrast = str_remove(contrast, "\\(")) %>%
  mutate(contrast = str_remove(contrast, "\\)")) %>%
  mutate(contrast = ordered(contrast)) %>%
  filter(!str_detect(contrast, "(HRD-Other) effect"),
         !str_detect(contrast, "Undetermined")) %>%
  filter(!str_detect(contrast, "Other")) %>%
  mutate(log2.odds.ratio = log2(odds.ratio)) %>%
  mutate(log2.odds.ratio = ifelse(log2.odds.ratio > 1.5, 1.5, log2.odds.ratio),
         log2.odds.ratio = ifelse(log2.odds.ratio < -1.5, -1.5, log2.odds.ratio))

p_odds_ratio <- plot_data %>%
  #filter(p.value < 1E-20) %>%
#  filter(SE < 0.15) %>%
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


library(RColorBrewer)
# 获取一个有22种颜色的调色板
color_palette <- brewer.pal(11, "Spectral") # 获取一个包含11种颜色的调色板
color_palette <- rep(color_palette, each = ceiling(22 / 11))[1:22] # 重复调色板并截断到22种颜色 



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

p_CD20 <- gene_dotplot_odds_ratio(raw_table,"MS4A1")
ggsave("G:/TLSdata/TLS_total_data/figures/dotplot_hotmap/B_spatial.pdf",plot = p_CD20, width = 7, height = 4)


p_CD8 <- gene_dotplot_odds_ratio(raw_table,"CD8A")
#ggsave("G:/TLSdata/TLS_total_data/figures/dotplot_hotmap/CD8_spatial.pdf",plot = p_CD8, width = 7, height = 4)


p_cancer <- gene_dotplot_odds_ratio(raw_table,"KRT14")
#ggsave("G:/TLSdata/TLS_total_data/figures/dotplot_hotmap/B_spatial.pdf",plot = p_CD20, width = 7, height = 4)

p_CD4 <- gene_dotplot_odds_ratio(raw_table,"CD4")

p_TCF1 <- gene_dotplot_odds_ratio(raw_table,"TCF7")

p_KI67 <- gene_dotplot_odds_ratio(raw_table,"MKI67")

p_CD23 <- gene_dotplot_odds_ratio(raw_table,"FCER2")



p_total <- grid.arrange(p_cancer, p_CD20, p_CD8, p_CD4, p_CD23, p_TCF1, nrow = 1)

ggsave("G:/TLSdata/TLS_total_data/figures/dotplot_hotmap/gene_spatial.pdf",plot = p_total, width = 15, height = 8)


