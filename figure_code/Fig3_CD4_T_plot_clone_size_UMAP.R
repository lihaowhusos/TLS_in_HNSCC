####BCR clone size UMAP####
library(ggplot2)
library(ggrepel)  # this line is needed in every new session


CD4_TCR_table <- read.csv('G:/TLSdata/TLS_total_data/table/CD4_T_and_Treg_plot_clone_size_UMAP_for_R.csv',header = TRUE)

CD4_TCR_table_immature <- CD4_TCR_table[CD4_TCR_table$group == 'immature', ]
CD4_TCR_table_mature <- CD4_TCR_table[CD4_TCR_table$group == 'mature', ]
CD4_TCR_table_none <- CD4_TCR_table[CD4_TCR_table$group == 'none', ]


color = c('#ea615dff', '#f18988ff', '#6a120c', '#f2cece', '#d8433cff', '#b12921ff', 
          '#f18da7', '#ec6595', '#c72e82', '#f2ced4',
          '#580944', '#e0478c', '#9c196e')

p_im = ggplot(CD4_TCR_table_immature, aes(x=X, y=Y, size=clone_id_size, color=cluster_2)) +
  geom_point(alpha=0.25) +
  scale_color_manual(values = color)+
  scale_size_continuous(
    limits = c(0, 250),
    breaks = c(0, 125, 250),
    labels = c(0, 125, 250)
  )+
  theme_bw()


p_m = ggplot(CD4_TCR_table_mature, aes(x=X, y=Y, size=clone_id_size, color=cluster_2)) +
  geom_point(alpha=0.25) +
  scale_color_manual(values = color)+
  scale_size_continuous(
    limits = c(0, 250),
    breaks = c(0, 125, 250),
    labels = c(0, 125, 250)
  )+
  theme_bw()

P_n = ggplot(CD4_TCR_table_none, aes(x=X, y=Y, size=clone_id_size, color=cluster_2)) +
  geom_point(alpha=0.25) +
  scale_color_manual(values = color)+
  scale_size_continuous(
    limits = c(0, 250),
    breaks = c(0, 125, 250),
    labels = c(0, 125, 250)
  )+
  theme_bw()

############# fill color #######
p_im_fill = ggplot(CD4_TCR_table_immature, aes(x=X, y=Y, size=clone_id_size, fill=cluster_2)) +
  geom_point(colour = "black", shape = 21,alpha=1) +
  scale_fill_manual(values = color)+
  scale_size_continuous(
    range =  c(1,9),
    limits = c(0, 60),
    breaks = c(0, 30, 60),
    labels = c(0, 30, 60)
  )+
  theme_bw()


p_m_fill = ggplot(CD4_TCR_table_mature, aes(x=X, y=Y, size=clone_id_size, fill=cluster_2)) +
  geom_point(colour = "black", shape = 21,alpha=1) +
  scale_fill_manual(values = color)+
  scale_size_continuous(
    range =  c(1,9),
    limits = c(0, 60),
    breaks = c(0, 30, 60),
    labels = c(0, 30, 60)
  )+
  theme_bw()


P_n_fill = ggplot(CD4_TCR_table_none, aes(x=X, y=Y, size=clone_id_size, fill=cluster_2)) +
  geom_point(colour = "black", shape = 21,alpha=1) +
  scale_fill_manual(values = color)+
  scale_size_continuous(
    range =  c(1,9),
    limits = c(0, 60),
    breaks = c(0, 30, 60),
    labels = c(0, 30, 60)
  )+
  theme_bw()





ggsave("G:/TLSdata/TLS_total_data/figures/TCR_CD4/CD4_T_clone_type_immature.pdf",plot = p_im_fill, width = 7, height = 5)

ggsave("G:/TLSdata/TLS_total_data/figures/TCR_CD4/CD4_T_clone_type_mature.pdf",plot = p_m_fill, width = 7, height = 5)

ggsave("G:/TLSdata/TLS_total_data/figures/TCR_CD4/CD4_T_clone_type_none.pdf",plot = P_n_fill, width = 7, height = 5)

