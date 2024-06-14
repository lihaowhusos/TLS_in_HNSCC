####BCR clone size UMAP####
library(ggplot2)
library(ggrepel)  # this line is needed in every new session


BCR_table <- read.csv('G:/TLSdata/TLS_total_data/table/plot_clone_size_UMAP_for_R.csv',header = TRUE)
BCR_table$clone_id_size10 = BCR_table$clone_id_size *10



BCR_table$group <- BCR_table$sample
BCR_table$group <- recode(BCR_table$group, 
                          'szj105988'= 'immature',
                          'szj106005'= 'none',
                          'szj106495'= 'none',
                          'szj106560'= 'immature',                
                          'szj106562'= 'mature',
                          'szj107010'= 'mature',
                          'szj107145'= 'mature',
                          'szj107734'= 'mature',
                          'szj107849'= 'none',
                          'szj108352'= 'immature',
                          'szj106121'= 'mature',                  
                          'szj106138'= 'immature',                  
                          'szj106759'= 'none',
                          'szj106771'= 'mature'
)

BCR_table_immature <- BCR_table[BCR_table$group == 'immature', ]
BCR_table_mature <- BCR_table[BCR_table$group == 'mature', ]
BCR_table_none <- BCR_table[BCR_table$group == 'none', ]



color = c('#6bd0bcff', '#40b7b0ff', '#5e94d0ff', '#adced8ff', '#b3e0d1ff','#228c95ff', '#04263aff', '#263c9aff', '#0f576bff')

color = c('#40B7B0', '#6BD0BC', '#0F576B', '#B3E0D1', '#04263A','#228C95', '#263C9A', '#5E94D0', '#ADCED8')


p_im = ggplot(BCR_table_immature, aes(x=X, y=Y, size=clone_id_size, color=cluster_1)) +
  geom_point(alpha=0.5) +
  scale_color_manual(values = color)+
  theme_bw()


p_m = ggplot(BCR_table_mature, aes(x=X, y=Y, size=clone_id_size, color=cluster_1)) +
  geom_point(alpha=0.5) +
  scale_color_manual(values = color)+
  theme_bw()


P_n = ggplot(BCR_table_none, aes(x=X, y=Y, size=clone_id_size, color=cluster_1)) +
  geom_point(alpha=0.5) +
  scale_color_manual(values = color)+
  theme_bw()

############# fill color #######
p_im_fill = ggplot(BCR_table_immature, aes(x=X, y=Y, size=clone_id_size, fill=cluster_1)) +
  geom_point(colour = "black", shape = 21,alpha=1) +
  scale_fill_manual(values = color)+
  scale_size_continuous(
    range =  c(2,11),
    limits = c(1, 80),
    breaks = c(1, 40, 80),
    labels = c(1, 40, 80)
  )+
  theme_bw()


p_m_fill = ggplot(BCR_table_mature, aes(x=X, y=Y, size=clone_id_size, fill=cluster_1)) +
  geom_point(colour = "black", shape = 21,alpha=1) +
  scale_fill_manual(values = color)+
  scale_size_continuous(
    range =  c(2,11),
    limits = c(1, 80),
    breaks = c(1, 40, 80),
    labels = c(1, 40, 80)
  )+
  theme_bw()


P_n_fill = ggplot(BCR_table_none, aes(x=X, y=Y, size=clone_id_size, fill=cluster_1)) +
  geom_point(colour = "black", shape = 21,alpha=1) +
  scale_fill_manual(values = color)+
  scale_size_continuous(
    range =  c(2,11),
    limits = c(1, 80),
    breaks = c(1, 40, 80),
    labels = c(1, 40, 80)
  )+
  theme_bw()



ggsave("G:/TLSdata/TLS_total_data/figures/B/B_clone_type_immature.pdf",plot = p_im_fill, width = 6, height = 3.9)

ggsave("G:/TLSdata/TLS_total_data/figures/B/B_clone_type_mature.pdf",plot = p_m_fill, width = 6, height = 3.9)

ggsave("G:/TLSdata/TLS_total_data/figures/B/B_clone_type_none.pdf",plot = P_n_fill, width = 6, height = 3.9)

