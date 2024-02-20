library(tidyverse)
library(ggplot2)
library(ggthemes)
library(RColorBrewer)

rm(list = ls())

dat <- readRDS("./data/sparse.rds")

sparse_calc <- function(x){
    N <- nrow(x)
    s2 <- sum(x$rsquared/N)
    s1 <- (sum(x$firing)/N)**2
    S = (1-(s1/s2))/(1-(1/N))

    x <- x %>% mutate(sparse = round(S, 2), number = N)
    
}

pitch_tuning <- dat %>% 
    ungroup() %>% 
    mutate(clip_name = str_extract(name, "[A-Z]{1}_[0-9]{6}")) %>%
    group_by(mouse, id, session) %>% 
    mutate(norm.fr = (fr - min(fr))/(max(fr)-min(fr)))

clip_sets <- pitch_tuning %>% 
    group_by(mouse, id, session, clip_name) %>%
    select(mouse, id, session, clip_name) %>% 
    distinct() %>% 
    group_by(mouse, id) %>% 
    mutate(clip_no = 1:n())

sparsity <- pitch_tuning %>%
    group_by(mouse, id, clip_type, st) %>%
    mutate(session_no = 1:n(), firing=fr*100, rsquared = firing**2) %>% 
    group_by(mouse, id, clip_type, session) %>%
    group_modify( ~ sparse_calc(.x)) %>% 
    select(mouse, session, unit_no, ch, id, name, clip_name, sparse)
    
sparsity_vals <- sparsity %>% 
    select(mouse, session, id, clip_name, sparse) %>% 
    distinct()

pitch_tuning2 <- pitch_tuning %>% 
    left_join(clip_sets, by=c("mouse", "id", "session", "clip_name")) %>% 
    left_join(sparsity_vals, by=c("mouse", "id", "session", "clip_name", "sparse")) %>% 
    mutate(nnID = interaction(mouse, id)) %>% 
    arrange(nnID, sparse)
    
plt <- ggplot(pitch_tuning2 %>% filter(sparse > 0.2), aes(x=st, y=norm.fr)) + 
    geom_line(show.legend=FALSE, aes(colour=as.factor(nnID))) + 
    facet_wrap(nnID~clip_no) +
    geom_text(aes(label = sparse, y=0.75, x=3), size=4, show.legend=FALSE, colour="grey20", family="Helvetica") + 
    #facet_grid(id+mouse~clip_no) + 
    scale_y_continuous(breaks = c(0, 1), limits = c(0,1.1)) +
    scale_x_continuous(breaks = c(-6, 0, 6), limits = c(-6.1,6.1)) +
    theme(
    strip.background = element_blank(),
    strip.text.x = element_blank(), 
    axis.line = element_line(colour = "grey50"),
    axis.ticks = element_blank(),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    axis.title = element_text(size = 15, family="Helvetica"),
    axis.text = element_text(size = 11, family="Helvetica"))+
    #axis.ticks.y = element_blank()) + 
    annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf)+
    annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf) + 
    ylab("Normalised firing rate") + 
    xlab("Semitone shift")

ggsave(
  filename = "./thesis_graphs/sparse_02.pdf",
  plot = plt,
  width = 8,
  height = 8
)


  theme(
    axis.ticks = element_blank(),
    axis.line = element_line(colour = "grey50"),
    panel.grid = element_line(color = "#b4aea9"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(linetype = "dashed"),
    panel.background = element_rect(fill = "#fbf9f4", color = "#fbf9f4"),
    plot.background = element_rect(fill = "#fbf9f4", color = "#fbf9f4")
  )


