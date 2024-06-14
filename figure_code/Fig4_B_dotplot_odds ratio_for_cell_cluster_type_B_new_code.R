library(dplyr)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(lme4)
library(emmeans)

## 加载并处理数据 --------------------------------------
cell_data <- read.csv('G:/TLSdata/TLS_total_data/Processing data/subset/B_Plasma/B_Plasma_for_seurat_metadata.csv', header = TRUE)
rownames(cell_data) <- as.vector(as.matrix(cell_data[,1]))

cell_data$tumor_group <- recode(cell_data$sample, 
                                'szj105988'= 'immature', 'szj106005'= 'none', 'szj106495'= 'none',
                                'szj106560'= 'immature', 'szj106562'= 'mature', 'szj107010'= 'mature',
                                'szj107145'= 'mature', 'szj107734'= 'mature', 'szj107849'= 'none',
                                'szj108352'= 'immature', 'szj106121'= 'mature', 'szj106138'= 'immature',                  
                                'szj106759'= 'none', 'szj106771'= 'mature')

cell_data$barcode <- cell_data[,1]
cell_data <- cell_data[,-1]

cell_data <- cell_data %>% 
  filter(tissue_type == "tumor") %>% 
  mutate(cell_type = cell_type, patient = sample_CD45, tumor_type = tumor_group, cell_barcode = barcode)

sample_data <- cell_data %>%
  group_by(patient, tumor_type, cell_type) %>%
  summarize(cell_count = n_distinct(cell_barcode),
            other_cell_count = sum(cell_count) - cell_count,
            total_cell_count = sum(cell_count),
            cell_proportion = cell_count/total_cell_count*100) %>%
  ungroup()

## 建模与分析
data_list <- list(scRNA = sample_data)
formula_list <- list(scRNA = cbind(cell_count, other_cell_count) ~ cell_type * tumor_type)
model_list <- list(scRNA = glm(formula = formula_list$scRNA, family = binomial(link = 'logit'), data = data_list$scRNA))

emm_result <- emmeans(model_list$scRNA, specs = eff ~ tumor_type | cell_type * tumor_type, 
                      type = "response", tran = "logit", adjust = "bonferroni")




odds_result <- emm_result$contrasts %>% rbind() %>% as.data.frame() %>% 
  mutate(contrast_var = "Signature",
         p_val = ifelse(p.value < 1E-50, 1E-50, p.value),
         contrast_name = str_remove_all(contrast, " effect|\\(|\\)"),
         contrast_name = ordered(contrast_name),
         log2_odds = log2(odds.ratio),
         log2_odds = ifelse(log2_odds > 1.5, 1.5, log2_odds),
         log2_odds = ifelse(log2_odds < -1.5, -1.5, log2_odds))

## 绘图
color_pal <- rep(brewer.pal(11, "Spectral"), each = ceiling(22 / 11))[1:22]

p_odds <- odds_result %>%
  ggplot(aes(x = cell_type, y = forcats::fct_rev(contrast_name), color = log2_odds, size = -log10(p_val))) +
  geom_point() +
  scale_color_gradientn(colors = rev(brewer.pal(9, "RdBu")), limits = c(-1.5, 1.5), labels = c("≤-1.5", 0, "≥1.5"), breaks = c(-1.5, 0, 1.5)) +
  scale_size_continuous(limits = c(0, 50), breaks = c(0, 25, 50), labels = c(0, 25, 50)) +
  scale_x_discrete(drop=FALSE) +
  theme_minimal() +
  theme(axis.text = element_blank(), axis.ticks = element_blank(), legend.position = "top", 
        panel.grid = element_blank(), strip.text = element_blank()) +
  guides(color = guide_colourbar(title.position="top", title.hjust = 0.5, title.vjust = 0),
         size = guide_legend(title.position="top", title.hjust = 0.5, title.vjust = 0)) +
  labs(title = NULL, x = NULL, y = NULL, color = 'log2(odds ratio)', size = "\n-log10(p value)") +
  facet_grid(contrast_var ~ ., scales = "free_y", space = "free")

p_sig_anno <- odds_result %>% 
  distinct(contrast_name, contrast_var, .keep_all = T) %>% 
  ggplot() +
  geom_tile(aes(x=0, y=forcats::fct_rev(contrast_name), fill = contrast_name)) +
  scale_fill_manual(values = c(scales::muted("blue"), scales::muted("red"), "grey10")) +
  facet_grid(contrast_var~., scales = "free", space = "free") +
  theme_void() + theme(strip.text = element_blank()) +
  guides(fill = "none")

p_cell_anno <- odds_result %>% 
  mutate(facet_var = "") %>% 
  distinct(cell_type, facet_var, .keep_all = T) %>% 
  ggplot() +
  geom_tile(aes(cell_type, facet_var, fill = cell_type)) +
  scale_fill_manual(values = color_pal) +
  scale_x_discrete(drop=FALSE) +
  facet_grid(~facet_var, scales = "free", space = "free") +
  theme_void() + theme(strip.text = element_blank()) +
  guides(fill = "none")

p_left_panel <- plot_grid(p_sig_anno, ggdraw(), nrow = 2, align = "hv", axis = "tbr", rel_heights = c(1, 0.7))
p_right_panel <- plot_grid(p_odds + theme(legend.position = "none"), p_cell_anno, nrow = 2, align = "v", axis = "lr", rel_heights = c(1, 0.7))
p_final <- plot_grid(plot_grid(p_left_panel, p_right_panel, nrow = 1, align = "hv", axis = "tblr", rel_widths = c(0.5, 1)),
                     get_legend(p_odds), ncol = 1, rel_heights = c(1, 0.15))

## 保存图片
ggsave("dot_plot_result.pdf", plot = p_final, width = 10, height = 8)