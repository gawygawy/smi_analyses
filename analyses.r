library(tidyverse)
library(ggplot2)
library(ggExtra)
library(GGally)
library(ggridges)
library(ggpubr)
library(lmerTest)
library(RColorBrewer)
library(grid)
library(cowplot)

rm(list = ls())

cc <- c("#264653", "#e76f51", "#2a9d8f")

# plot theme 
my_t <- theme_cowplot(font_size=11) + 
    theme(text = element_text(family = "Sans serif"))

dat <- readRDS("./data/filtered_df.rds")
sparse <- readRDS("./sparse.rds")


#dat <- dat %>% 
#    filter(!(nnID == 33 & session == "21-03-30" & base_pair == "AG" ))

## data frame with firing rates and stim frequencies
stim_info <- dat %>% 
    filter(!is_combi) %>%  
    select(mouse, ch, unit_no, stim_pair, base_pair, deviant, standard, meanSpec) %>% 
    distinct() %>% 
    mutate(meanSpec = round(meanSpec/1000, 2))

fr <- dat %>%
    ungroup() %>% 
    select(mouse, session, nnID, ch, unit_no, stim_pair, base_pair, fr, deviant, standard) 
   
merged_fr <- fr %>% 
    left_join(fr, by = c("deviant" = "standard", "stim_pair" = "stim_pair", "base_pair" = "base_pair", "mouse" = "mouse", "nnID" = "nnID", "session" = "session", "ch" = "ch", "unit_no" = "unit_no")) %>% 
    select(!c(deviant.y))

sound_info <- stim_info %>% 
    left_join(stim_info, by = c("deviant" = "standard", "stim_pair" = "stim_pair", "base_pair" = "base_pair", "mouse" = "mouse", "nnID" = "nnID", "ch" = "ch", "unit_no" = "unit_no")) %>% 
    ungroup() %>% 
    select(!c(deviant.y))

adapt_info <- dat %>% select(mouse, session, ch, unit_no, nnID, stim_pair, deviant, standard, base_pair, first_response, second_response, above_threshold, pitch1, pitch2, smi, is_combi, CI)

new_dat <- merged_fr %>% 
    left_join(sound_info, by = c("mouse", "stim_pair", "base_pair", "ch", "unit_no", "deviant", "nnID", "standard")) %>% 
    left_join(adapt_info, by=c("mouse", "stim_pair", "base_pair", "ch", "unit_no", "nnID", "deviant", "standard", "session")) %>% 
    select(mouse, session, nnID, ch, unit_no, base_pair, stim_pair, deviant, standard, is_combi, fr.x, fr.y, smi, above_threshold, first_response, second_response, meanSpec.x, meanSpec.y, pitch1, pitch2, CI)

#  
smis <- new_dat %>% 
    filter(!is_combi) %>% 
    select(mouse, session, nnID, ch, unit_no, base_pair, stim_pair, smi, first_response, second_response, above_threshold) %>% 
    mutate(recovery = (second_response/first_response)*100) %>% 
    group_by(mouse, session, nnID, stim_pair) %>% 
    mutate(max_recovery = max(recovery), stim_no = 1:n()) %>% 
    filter(stim_no == 1) %>% 
    mutate(independent = if_else(above_threshold == 2, if_else(max_recovery > 40, "independent", "overlapping"), "z")) %>% 
    mutate(summation = cut(smi, breaks=c(-Inf, 0.1, 1, Inf), labels=c("MAX", "sublinear", "AND")))

smis_ind_only <- smis %>% 
    filter(independent == "independent")

smis_no_inf <- smis %>% 
    mutate(smi = if_else(smi == Inf, 7, smi)) 

smis_no_inf_multiple <- smis_no_inf %>% 
    group_by(nnID) %>% 
    filter(n() > 2)  %>% 
    mutate(smiM = mean(smi)) %>% 
    distinct() 

# plot to show the variability in SmIs
ridges <- ggplot(smis_no_inf_multiple, aes(x=smi, y=fct_reorder(nnID, smiM))) + 
    geom_density_ridges(
        jittered_points = TRUE, 
        position = position_points_jitter(width = 0.05, height=0), 
        point_shape = '|', point_size = 2.5, point_alpha = 2, alpha=0.5,
    ) + 
    theme_ridges(font_size=15) + 
    scale_x_continuous(breaks=seq(-1, 7), expand=c(0,0)) + 
    scale_y_discrete(expand = c(0.02, 0))

ggsave("./thesis_graphs/ridges.pdf", ridges, dpi=500)

p2 <- ggplot(smis_no_inf, aes(x=smi)) + 
    geom_histogram(boundary=0, fill="grey") + 
    geom_histogram(data = smis_ind_only) + 
    geom_vline(xintercept=c(0, 1), linetype=2, colour="black", size=0.5, alpha=0.8) + 
    scale_x_continuous(breaks = seq(-3, 8, by=1)) +
    scale_y_continuous(expand = c(0, 0)) + 
    theme_cowplot(font_size=11) + 
    ylab("Count") + 
    xlab("SmI") + 
    theme(text = element_text(family = "Sans serif"))
    
ggsave("./thesis_graphs/distribution_adapt.pdf", p2, width=3.5, height=2, device=cairo_pdf)

## show that some units perform more than one computation 
flex_units <- smis %>% 
    ungroup() %>%
    select(nnID, summation) %>% 
    group_by(nnID, summation) %>% 
    summarize(ct=n()) %>%
    mutate(ct = 1) %>% 
    pivot_wider(names_from=summation, values_from=ct) %>% 
    mutate(sublinear = ifelse(!is.na(sublinear), 0.5, sublinear), 
        AND = ifelse(!is.na(AND), 2, AND)) %>% 
    mutate(funcsum= sum(sublinear, MAX, AND, na.rm=TRUE)) 

flex_units$behaviour <- recode(flex_units$funcsum, 
            `0.5` = "sublinear", 
            `3.5` = "MAX/sub/AND", 
            `2` = "AND", 
            `2.5` = "sub/AND",
            `1.5` = "MAX/sub", 
            `1` = "MAX")

flex_units$order <- recode(flex_units$funcsum, 
            `1.5` = 0, 
            `2.5` = 1.0,
            `2` = 5.0)

flex_units <- flex_units %>% 
    mutate(summation = as.factor(behaviour)) %>% 
    mutate(summation = fct_relevel(summation, "AND", "sublinear", "MAX", "sub/AND", "MAX/sub", "MAX/sub/AND"))

p3 <- ggplot(flex_units, aes(x=summation)) + # shows types of summation behaviour in each neuron 
    geom_bar(width=0.7) + 
    coord_flip() + 
    ylab("Number of neurons") + 
    xlab("Summation") + 
    my_t


ggsave("./thesis_graphs/flexible_units.pdf", p3, width=3, height=1.5, device=cairo_pdf)

# merge flex labels and ridge together 

smi_range <- smis_no_inf_multiple %>% 
    left_join(flex_units %>% select(nnID, order, summation), by="nnID")

ridge2 <- ggplot(smi_range, aes(x=smi, y=fct_reorder(nnID, order), fill=summation.y)) + 
    geom_density_ridges(
        jittered_points = TRUE, 
        position = position_points_jitter(width = 0.05, height=0), 
        point_shape = '|', point_size = 1, point_alpha = 2, alpha=0.6,
    ) + 
    scale_fill_manual(values=c("#8c57b5","#69a746","#b54e5d","#669b90","#b5813f")) + 
    scale_x_continuous(breaks=seq(-1, 7), expand=c(0,0)) + 
    scale_y_discrete(expand = c(0.02, 0)) + 
    xlab("SmI") + 
    my_t + 
    theme( 
    axis.ticks.y=element_blank(), 
    axis.text.y = element_blank(), 
    axis.title.y = element_blank())

ggsave("./thesis_graphs/flexible_units_labels.pdf", ridge2, width=3.65, height=4.5, device=cairo_pdf )

# summation labels 
smi_labels <- smis %>% ungroup() %>% select(mouse, nnID, stim_pair, summation)

# merge back onto original data frame 
dat_with_labels <- new_dat %>% 
    ungroup() %>% 
    left_join(smi_labels, by=c("mouse", "nnID", "stim_pair"))

# look at sets with pitch shifts 
pitch_diff <- dat_with_labels %>% 
    mutate(pitch_diff = meanSpec.x - meanSpec.y) %>% 
    group_by(nnID, stim_pair) %>% 
    mutate(set_no = 1:n()) %>% 
    filter(set_no == 1)

pitch_diff_only <- pitch_diff %>% select(nnID, base_pair, stim_pair, pitch_diff)

loong <- dat_with_labels %>% 
    group_by(nnID, base_pair) %>% 
    mutate(no_sets = length(unique(stim_pair))) %>%
    #filter(no_sets > 1) %>% 
    #left_join(pitch_diff_only, by = c("nnID", "base_pair", "stim_pair")) %>% 
    arrange(mouse, nnID, base_pair, is_combi, deviant) %>% 
    select(!c(first_response, second_response)) 

firing_rates_long_fr <- loong %>% 
    select(mouse, session, nnID, ch, unit_no, base_pair, stim_pair, deviant, is_combi, fr.x, smi, above_threshold, summation) %>% 
    #drop_na(pitch_diff) %>% 
    group_by(nnID, stim_pair) %>% 
    mutate(stim_type = 1:n()) %>% 
    select(!c(deviant, is_combi)) %>% 
    pivot_wider(names_from=stim_type, values_from=fr.x) %>% 
    mutate(fr_tot = `1` + `2`) %>% 
    arrange(mouse, nnID, base_pair, fr_tot) %>% 
    group_by(nnID, base_pair) %>% 
    mutate(set = 1:n())

set_change <- firing_rates_long_fr %>% 
    ungroup() %>% 
    filter(nnID == 10, stim_pair %in% c("MQ", "LO")) %>% 
    mutate(set = 1:n(), ra =`1`, rb=`2`, rab=`3`)

rate1 <- ggplot(set_change, aes(x=set)) +
    geom_line(aes(y=ra), color=cc[1], size=1) +
    geom_point(aes(y=ra), size=1.5, color=cc[1]) +
    geom_line(aes(y=rb), color=cc[2], size=1) + 
    geom_point(aes(y=rb), size=1.5, color=cc[2]) +
    geom_line(aes(y=rab), color=cc[3], size=1) +
    geom_point(aes(y=rab), size=1.5, color=cc[3]) +
    #geom_line(aes(y=smi*100), color=c[4], size=1.5) +
    #geom_point(aes(y=smi*100), size=2.5, color=c[4]) +
    #scale_y_continuous(sec.axis = sec_axis(~./100, name="smi")) +
    ylab("Spikes") +
    xlab("Syllable pair") + 
    scale_x_continuous(breaks = seq(1, 2,by=1), expand=c(0.2, 0.2)) + 
    my_t

smi1 <- ggplot(set_change, aes(x=set)) +
    ylab("SmI") +
    geom_line(aes(y=smi), size=1) +
    geom_point(aes(y=smi), size=1.5) + 
    scale_x_continuous(breaks = seq(1, 2,by=1), expand=c(0.2, 0.2)) + 
    my_t + 
    theme(axis.text.x=element_blank(), axis.title.x=element_blank())

plot_grid(smi1, rate1, ncol=1, rel_heights=c(1, 1.5), align="v")

ggsave("./thesis_graphs/set_change.pdf", width=1.3, height=2, device=cairo_pdf)

pitch_change <- firing_rates_long_fr %>% 
    ungroup() %>% 
    filter(base_pair == "BN", nnID == 46) %>% 
    mutate(set = 1:n(), ra =`1`, rb=`2`, rab=`3`)

rate2 <- ggplot(pitch_change, aes(x=set)) +
    geom_line(aes(y=ra), color=cc[1], size=1) +
    geom_point(aes(y=ra), size=1.5, color=cc[1]) + #green
    geom_line(aes(y=rb), color=cc[2], size=1) + 
    geom_point(aes(y=rb), size=1.5, color=cc[2]) + # blue
    geom_line(aes(y=rab), color=cc[3], size=1) +
    geom_point(aes(y=rab), size=1.5, color=cc[3]) + # orange
    #geom_line(aes(y=smi*100), color=c[4], size=1.5) +
    #geom_point(aes(y=smi*100), size=2.5, color=c[4]) + # pink
    #scale_y_continuous(sec.axis = sec_axis(~./100, name="smi")) +
    ylab("Spikes") +
    xlab("Syllable pair") +
    scale_x_continuous(breaks = seq(1, 4,by=1), expand=c(0.1, 0.2)) + 
    my_t 

smi2 <- ggplot(pitch_change, aes(x=set)) + 
    geom_line(aes(y=smi), size=1) +
    geom_point(aes(y=smi), size=1.5) + 
    scale_x_continuous(breaks = seq(1, 4,by=1), expand=c(0.2, 0.2)) + 
    my_t + 
    ylab("SmI") + 
    theme(axis.text.x=element_blank(), axis.title.x=element_blank())

plot_grid(smi2, rate2, ncol=1, rel_heights=c(1, 1.5), align="v")

ggsave("./thesis_graphs/pitch_change.pdf", width=2, height=2, device=cairo_pdf)
