library(tidyverse)
library(ggplot2)
library(ggExtra)
library(GGally)
library(ggpubr)
library(lmerTest)
library(pracma)
library(Convolutioner)

source("selectivity.r")
my_t <- theme_cowplot(font_size=11) + 
    theme(text = element_text(family = "Sans serif"))

rm(list = ls())
pitch_files <- list.files(path = '~/rds/projects/kozlovlab_rds2/live/Grace/data/processed/combined_stims/pitch_shift_data/', pattern="*ch(\\d+).csv", recursive=TRUE, full.names = TRUE)

dat_file <- as_tibble(do.call(rbind, lapply(pitch_files, read.csv)))

# read in ids 

ids <- read.csv("~/rds/projects/kozlovlab_rds2/live/Grace/data/processed/combined_stims/pitch_shift_data/pitch_shift_id.csv") %>% 
    left_join(dat_file, by=c("mouse", "session", "unit_no", "ch"))

#fit_gauss(ff)

pitch_tuning <- ids %>% 
    group_by(mouse, id, clip_type, st) %>% 
    mutate(fr = mean(fr)) %>% 
    group_by(mouse, id, clip_type) %>%
    mutate(avg_fr = mean(fr)) %>% 
    #group_modify(~ fit_gauss(.x)) %>% 
    #group_modify( ~ peak_func(.x)) %>% # doesn't really work - need to think harder 
    mutate(response_selectivity = round(max(fr)/mean(fr), 2), nsig = round(sum(p_val < 0.05)/n(), 2), clip_name = str_extract(name, "[A-Z]{1}_[0-9]{6}")) #%>% 
    #select(mouse, session, unit_no, ch, id, clip_name, response_selectivity, nsig) %>% 
    #distinct()

saveRDS(pitch_tuning, "curves.rds")

sparsity <- ids %>%
    group_by(mouse, id, clip_type, st) %>%
    mutate(session_no = 1:n(), firing=fr*100, rsquared = firing**2) %>% 
    group_by(mouse, id, clip_type, session) %>%
    group_modify( ~ sparse_calc(.x)) %>% 
    select(mouse, session, unit_no, ch, id, name, clip_type, sparse)

saveRDS(sparsity, "sparse.rds")

mice <- unique(sparsity$mouse)

for (m in mice){
    mouse_dat <- sparsity %>% filter(mouse == m)  
    for (unit_id in unique(mouse_dat$id)){    
        p <- ggplot(mouse_dat %>% filter(id == unit_id), aes(x=st, y=firing)) + 
        geom_line() + 
        geom_point()
        geom_text(aes(label = sparse, y=firing[1]), x=0) + 
        facet_wrap(~clip_type+session_no)
        ggsave(file.path(getwd(), "tuning_curve_sparse", paste(m, "_", unit, ".png", sep="")), p)}
}

cc <- c("#264653", "#e76f51", "#2a9d8f")

ggplot(pitch_tuning %>% filter(mouse == "2022-03-25"), aes(x=st, y=fr)) + 
    geom_line() + 
    geom_text(aes(label = nsig, y=avg_fr[1]), x=0) + 
    facet_wrap(~mouse+id+clip_type)


example_ra <- pitch_tuning %>% 
    filter(unit_no == 5, clip_type %in% c("B", "N"), mouse == "2022-09-14")

example_rb <- pitch_tuning %>% 
    filter(unit_no == 5, clip_type %in% "N", mouse == "2022-09-14")

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

n_tun <- ggplot(example_rb, aes(x=st, y=fr, 3)) + 
    geom_line(color=c[3], size=1) +
    geom_point(size=2, color=c[3]) +
    ylim(0, 250) + 
    scale_x_continuous(breaks = seq(-6, 4,by=1), expand=c(0.1, 0.1)) + 
    theme_classic()

pdf("./thesis_graphs/tuning.pdf", width=2, height=2.5)
grid.newpage()
grid.draw(rbind(ggplotGrob(n_tun), ggplotGrob(b_tun), size="last"))
dev.off() 

ggsave("./cosyne_graphs/n_tuning.pdf", n_tun, width=4, height=3)

smooth <- Hamming(dd$fr, 3)

nrow(findpeaks(smooth, nups=0, ndowns=1, minpeakheight=1.5*dd$avg_fr[1], minpeakdistance=1))    




findpeaks(Hamming(dd$fr, 3), nups=0, ndowns=1, minpeakheight=1.5*dd$avg_fr[1], minpeakdistance=2)


tt <- pitch_tuning %>% 
    group_by(mouse, clip_type, id) %>% 
    mutate(step = diff(st)[1])



angles <- ids  %>%    
    group_by(mouse, id, clip_type, st) %>% 
    mutate(fr = mean(fr)) %>% 
    mutate(angle = st*30) %>% 
    mutate(trig = fr*exp(1i*2*angle)) %>% 
    group_by(mouse, id, clip_type) %>% 
    mutate(osi = round(abs(sum(trig))/sum(fr), 3))

for (m in mice){    
    p <- ggplot(angles %>% filter(mouse == m), aes(x=st, y=fr)) + 
    geom_line() + 
    geom_text(aes(label = osi, y=fr[1]), x=0) + 
    facet_wrap(~mouse+id+clip_type)
    ggsave(file.path(getwd(), "tuning_curves_nsig", paste(m, ".png", sep="")), p)
}


x <- example_rb$st
y <- example_rb$fr

mod.smsp <- smooth.spline(x, y, nknots = 6)
plot(x, y)
lines(x, mod.smsp$y, lty = 3, col = 3, lwd = 2)
