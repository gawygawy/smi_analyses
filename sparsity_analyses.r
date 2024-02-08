library(tidyverse)
library(cowplot)

rm(list = ls())

my_t <- theme_cowplot(font_size=11) + 
    theme(text = element_text(family = "Sans serif"))

S <- readRDS("./data/sparse.rds")

example_ra <- S %>% 
    filter(unit_no == 5, clip_type %in% c("B", "N"), mouse == "2022-09-14")

cc <- c("#264653", "#e76f51", "#2a9d8f")

b_tun <- ggplot(example_ra, aes(x=st, y=fr)) + 
    geom_line(aes(color=clip_type), size=1) +
    geom_point(aes(color=clip_type), size=2) +
    scale_color_manual(values=cc) +
    facet_wrap(~clip_type, nrow=2) + 
    scale_x_continuous(breaks = seq(-6, 4,by=2), expand=c(0.1, 0.1)) + 
    my_t + 
    theme(legend.position="none") + 
    annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf)

ggsave("./thesis_graphs/tuning.pdf", b_tun, width=2, height=2.5, device=cairo_pdf)