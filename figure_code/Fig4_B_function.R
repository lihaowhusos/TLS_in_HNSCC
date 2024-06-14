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





ggplot(aes(x = cell_type, y = forcats::fct_rev(contrast), color = log2.odds.ratio, size = -log10(p.value))) +
  geom_point() +
  scale_color_gradientn(
    colors = rev(brewer.pal(9, "RdBu")),
    limits = c(-1, 1),
    labels = c("≤-1", 0, "≥1"),
    breaks = c(-1, 0, 1)
  ) +
  scale_size_continuous(
    limits = c(0, 50),
    breaks = c(0, 25, 50),
    labels = c(0, 25, 50)
  )







