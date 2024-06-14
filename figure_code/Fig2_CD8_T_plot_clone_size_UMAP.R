####BCR clone size UMAP####
library(ggplot2)
library(ggrepel)  # this line is needed in every new session


CD8_TCR_table <- read.csv('G:/TLSdata/TLS_total_data/table/CD8_T_plot_clone_size_UMAP_for_R.csv',header = TRUE)
CD8_TCR_table$clone_id_size10 = CD8_TCR_table$clone_id_size /10

CD8_TCR_table_immature <- CD8_TCR_table[CD8_TCR_table$group == 'immature', ]
CD8_TCR_table_mature <- CD8_TCR_table[CD8_TCR_table$group == 'mature', ]
CD8_TCR_table_none <- CD8_TCR_table[CD8_TCR_table$group == 'none', ]


color = c('#d4b44eff', '#e9dd7bff', '#f1f1c1ff', '#a87d2bff', '#603b10ff')

p_im = ggplot(CD8_TCR_table_immature, aes(x=X, y=Y, size=clone_id_size, color=cluster_2)) +
  geom_point(alpha=0.25) +
  scale_color_manual(values = color)+
  scale_size_continuous(
    limits = c(0, 250),
    breaks = c(0, 125, 250),
    labels = c(0, 125, 250)
  )+
  theme_bw()


p_m = ggplot(CD8_TCR_table_mature, aes(x=X, y=Y, size=clone_id_size, color=cluster_2)) +
  geom_point(alpha=0.25) +
  scale_color_manual(values = color)+
  scale_size_continuous(
    limits = c(0, 250),
    breaks = c(0, 125, 250),
    labels = c(0, 125, 250)
  )+
  theme_bw()

P_n = ggplot(CD8_TCR_table_none, aes(x=X, y=Y, size=clone_id_size, color=cluster_2)) +
  geom_point(alpha=0.25) +
  scale_color_manual(values = color)+
  scale_size_continuous(
    limits = c(0, 250),
    breaks = c(0, 125, 250),
    labels = c(0, 125, 250)
  )+
  theme_bw()

############# fill color #######
p_im_fill = ggplot(CD8_TCR_table_immature, aes(x=X, y=Y, size=clone_id_size, fill=cluster_2)) +
  geom_point(colour = "black", shape = 21,alpha=1) +
  scale_fill_manual(values = color)+
  scale_size_continuous(
    range =  c(1,9),
    limits = c(0, 250),
    breaks = c(0, 125, 250),
    labels = c(0, 125, 250)
  )+
  theme_bw()


p_m_fill = ggplot(CD8_TCR_table_mature, aes(x=X, y=Y, size=clone_id_size, fill=cluster_2)) +
  geom_point(colour = "black", shape = 21,alpha=1) +
  scale_fill_manual(values = color)+
  scale_size_continuous(
    range =  c(1,9),
    limits = c(0, 250),
    breaks = c(0, 125, 250),
    labels = c(0, 125, 250)
  )+
  theme_bw()


P_n_fill = ggplot(CD8_TCR_table_none, aes(x=X, y=Y, size=clone_id_size, fill=cluster_2)) +
  geom_point(colour = "black", shape = 21,alpha=1) +
  scale_fill_manual(values = color)+
  scale_size_continuous(
    range =  c(1,9),
    limits = c(0, 250),
    breaks = c(0, 125, 250),
    labels = c(0, 125, 250)
  )+
  theme_bw()





ggsave("G:/TLSdata/TLS_total_data/figures/TCR_CD8/CD8_T_clone_type_immature.pdf",plot = p_im_fill, width = 7, height = 5)

ggsave("G:/TLSdata/TLS_total_data/figures/TCR_CD8/CD8_T_clone_type_mature.pdf",plot = p_m_fill, width = 7, height = 5)

ggsave("G:/TLSdata/TLS_total_data/figures/TCR_CD8/CD8_T_clone_type_none.pdf",plot = P_n_fill, width = 7, height = 5)

