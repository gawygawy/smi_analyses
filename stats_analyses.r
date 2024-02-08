library(tidyverse)
library(lmerTest)
library(cowplot)
library(ggpubr)

rm(list=ls())

smi_dat <- readRDS("./data/filtered_df.rds")

# theme for plots 
my_t <- theme_cowplot(font_size=11) + 
    theme(text = element_text(family = "Sans serif"))

# calculate simple summary statistics 
fr <- smi_dat %>%
    ungroup() %>% 
    select(mouse, session, nnID, ch, unit_no, stim_pair, base_pair, fr, deviant, standard) 

merged_fr <- fr %>% 
    left_join(fr, by = c("deviant" = "standard", "stim_pair" = "stim_pair", "base_pair" = "base_pair", "mouse" = "mouse", "nnID" = "nnID", "session" = "session", "ch" = "ch", "unit_no" = "unit_no")) %>% 
    select(!c(deviant.y))

adapt_info <- filtered_dat %>% select(mouse, session, ch, unit_no, nnID, stim_pair, deviant, standard, base_pair, first_response, second_response, above_threshold, smi, is_combi, CI)

new_dat <- merged_fr %>% 
    left_join(adapt_info, by=c("mouse", "stim_pair", "base_pair", "ch", "unit_no", "nnID", "deviant", "standard", "session")) %>% 
    select(mouse, session, nnID, ch, unit_no, base_pair, stim_pair, deviant, standard, is_combi, fr.x, fr.y, smi, above_threshold, first_response, second_response, CI)

firing_metrics <- new_dat %>% 
    filter(!is_combi) %>% 
    select(mouse, session, nnID, ch, unit_no, base_pair, stim_pair, smi, first_response, second_response, above_threshold, fr.x, fr.y) %>% 
    mutate(recovery = (second_response/first_response), recovery = if_else(is.infinite(recovery), 0, recovery), recovery = replace_na(recovery, 0)) %>% 
    group_by(mouse, session, nnID, stim_pair) %>% 
    mutate(max_recovery = max(recovery), stim_no = 1:n(), tot_fr = fr.x + fr.y) %>% 
    filter(stim_no == 1) %>% 
    select(mouse, session, nnID, stim_pair, base_pair, smi, above_threshold, tot_fr, max_recovery)

# correlate smi with various stimulus features 
euclidean_dist <- read_csv("./data/euclidean.csv") %>% # this is the euclidean distance measured between 
    select(mouse, nnID, stim_pair, dists) %>% 
    mutate(nnID = as.factor(nnID), mouse = as.character(mouse)) %>% 
    rename(euclidean_d = dists)

windowed_dist <- read_csv("./data/dist_short_specs.csv") %>% # distance metrics for stimuli within 55ms of spiking 
    select(mouse, nnID, stim_pair, d, eu, temp, freq) %>% 
    mutate(nnID = as.factor(nnID), mouse = as.character(mouse)) %>% 
    rename(eu_win = eu, d_win = d)

firingNfeatures <- firing_metrics %>% 
    left_join(euclidean_dist, by=c("nnID", "stim_pair", "mouse")) %>%
    left_join(windowed_dist, by=c("nnID", "stim_pair", "mouse")) %>% 
    distinct() %>% 
    ungroup() %>% 
    filter(!(nnID == 49 & stim_pair == "(('K', '-4'), ('M', '-1'))")) %>% 
    mutate(smi = if_else(smi == Inf, 7, smi)) %>% 
    mutate(summation = cut(smi, breaks=c(-Inf, 0.1, 1, Inf), labels=c("MAX", "sublinear", "AND"))) %>% 
    mutate(smi.t = smi + 0.8800001) %>% # https://sscc.wisc.edu/sscc/pubs/RegDiag-R/normality.html#fn2 scale smi to positive values for boxcox transformation
    filter(complete.cases(.)) %>% 
    filter(smi < 7) 

cor.test(formula = ~ tot_fr + smi, data=firingNfeatures, method="spearman", exact=FALSE) # R_a + R_b vs SmI
cor.test(formula = ~ max_recovery + smi, data=firingNfeatures, method="spearman", exact=FALSE) # recovery from adaptation vs SmI
cor.test(formula = ~ temp + smi, data=firingNfeatures, method="spearman", exact=FALSE) # temporal correlation vs SmI
cor.test(formula = ~ euclidean_d + smi, data=firingNfeatures, method="spearman", exact=FALSE) # euclidean distance vs SmI
cor.test(formula = ~ d_win + smi, data=firingNfeatures, method="spearman", exact=FALSE) # spectral feature distance 
cor.test(formula = ~ d_win + smi, data=firingNfeatures %>% filter(temp < 0), method="spearman", exact=FALSE) # for uncorrelated inputs 

# plots for Supplementary figure 3
boxdat <- firingNfeatures %>% 
    mutate(summation=as.factor(summation)) %>% 
    dplyr::select(euclidean_d, max_recovery, summation, nnID, d_win, nnID, tot_fr, smi) %>% 
    rename(AB_input=tot_fr, spectral_feature_dist = d_win, adaptation=max_recovery) %>% 
    pivot_longer(!c(summation, nnID, smi), names_to = "variable", values_to = "val")

boxdat$summation <- recode(boxdat$summation, sublinear = "sub")

boxplots <- ggplot(boxdat %>% filter(variable %in% c("AB_input", "adaptation", "euclidean_d", "spectral_feature_dist")), aes(summation, val)) + 
    geom_boxplot(outlier.shape = NA) + 
    #geom_jitter(fill="white", size=0.6, alpha=0.8, width=0.2, shape=21) + 
    facet_wrap(~variable, scales="free", nrow=4) + 
    xlab("") +
    coord_flip()+ 
    my_t + 
    theme(panel.spacing=unit(5, "lines"))

corrplots <- ggplot(boxdat %>% filter(variable %in% c("AB_input", "adaptation", "euclidean_d", "spectral_feature_dist")), aes(val, smi)) + 
    geom_point(shape=21, size=0.6) + 
    facet_wrap(~variable, scales="free",nrow=4) + 
    stat_cor(method="spearman") + 
    geom_hline(yintercept = c(0, 1), linetype='dashed',col="#b4aea9") +
    ylab("SmI") + 
    xlab("") + 
    my_t

analyses <- plot_grid(corrplots, boxplots)
#ggsave("box_plots.pdf", analyses, width=5, height=7.5, device=cairo_pdf)

# smi distribution is non-normal 
#https://www.css.cornell.edu/faculty/dgr2/_static/files/R_html/Transformations.html
BCTransform <- function(y, lambda=0) {
    if (lambda == 0L) { log(y) }
    else { (y^lambda - 1) / lambda }
    }

BCTransformInverse <- function(yt, lambda=0) {
    if (lambda == 0L) { exp(yt) }
    else { exp(log(1 + lambda * yt)/lambda) }
    }

transform.mod <- lm(smi.t  ~ tot_fr + max_recovery + temp + d_win, data=firingNfeatures %>% filter(smi < 7))
bc <- MASS::boxcox(transform.mod) 
bc.power <- (bc.power <- bc$x[which.max(bc$y)])
firingNfeatures$smi.bc <- BCTransform(firingNfeatures$smi.t, bc.power)
par(mfrow=c(2,2))
hist(firingNfeatures$smi.bc, breaks=30);

bc.mixed.mod <- lmer(smi.bc ~ d_win + tot_fr + max_recovery + temp + (1|nnID), data=firingNfeatures %>% drop_na()) # slightly better qq norm
qqnorm(resid(bc.mixed.mod))
qqline(resid(bc.mixed.mod))

library(partR2)

R2_BMa <- partR2(bc.mixed.mod, partvars = c("tot_fr", "d_win", "temp"), 
                  R2_type = "marginal", nboot = 100)

library("performance")
r2_nakagawa(bc.mixed.mod)

sparsity <- readRDS("./data/sparse.rds")

sparsity_dat <- filtered_dat %>%
    ungroup() %>% 
    dplyr::select(mouse, id, nnID) %>% 
    distinct() %>%
    right_join(sparsity, by=c("mouse","id")) %>% 
    dplyr::select(!id) %>% 
    mutate(clip_name = str_extract(name, "[A-Z]{1}_[0-9]{6}")) %>% 
    dplyr::select(mouse, nnID, clip_name, sparse) %>%
    distinct() %>% 
    drop_na()

sparse_smi <- filtered_dat %>% 
    mutate(clip_name = str_extract(deviant, "[A-Z]{1}_[0-9]{6}")) %>% 
    left_join(sparsity_dat, by=c("mouse", "nnID", "clip_name")) %>% 
    dplyr::select(!"clip_name") %>% 
    group_by(nnID, stim_pair) %>% 
    filter(!is_combi) %>%
    dplyr::select(mouse, stim_pair, base_pair, deviant, smi, sparse) %>% 
    drop_na(sparse, smi) %>%  
    arrange(mouse, nnID, stim_pair, base_pair, sparse) %>% 
    mutate(mean_sparse = mean(sparse), order=1:n()) %>% 
    filter(order < 3) %>% 
    dplyr::select(!deviant) %>%  
    pivot_wider(names_from=order, values_from=sparse) %>% 
    rename(tun1 = `1`, tun2 = `2`) %>% 
    group_by(nnID, base_pair) %>% 
    mutate(mean_smi = mean(smi), tun_both = tun1+tun2) %>% 
    mutate(summation = as.factor(cut(smi, breaks=c(-Inf, 0.1, 1, Inf), labels=c("MAX", "sublinear", "AND"))))

art.mod2 <- art(tun1 ~ summation + (1|nnID), sparse_smi %>% drop_na()) # significant
contrast(emmeans(artlm(art.mod2, "summation"), ~summation), method="pairwise")
art.con(art.mod2, "summation", adjust="holm") %>%  summary()



ggplot(sparse_smi, aes(x=as.factor(summation), y=tun_both)) + 
    geom_boxplot(shape=21, size=0.8, outlier.shape=NA) + 
    geom_jitter(fill="white", size=2, alpha=0.8, width=0.2, shape=21)


cor.test(sparse_smi$smi, sparse_smi$tun1, method="spearman", exact=FALSE)

ccc <- sparse_smi %>% select(mean_smi, tun_both) %>% distinct()
cor.test(ccc$mean_smi, ccc$tun_both, method="spearman", exact=FALSE)


ggplot(sparse_smi, aes(x=tun2, y=smi))+ 
    geom_point()

sp <- ggplot(sparse_smi, aes(x=tun1, y=tun2))+ 
    geom_point(shape=21) + 
    theme_cowplot(font_size=11) + 
    theme(text = element_text(family = "Sans serif")) + 
    ylab("Higher selectivity index") +
    xlab("Smaller selectivity index")

ggsave("selectivity.pdf", sp, device=cairo_pdf, width=2.5, height=2.5)

smiVsparse <- sparse_smi  %>% 
    dplyr::select(smi, tun1, tun2) %>% 
    pivot_longer(cols=c("tun1", "tun2"))


smivsparse <- ggplot(smiVsparse, aes(x=value, y=smi))+ 
    geom_point(shape=21, size=0.8) + 
    facet_wrap(~ name, nrow=2, scales="free") + 
    scale_x_continuous(limits=c(0, 0.85))+
    theme_cowplot(font_size=11) + 
    theme(text = element_text(family = "Sans serif"), axis.line=element_line(color="black")) + 
    ylab("SmI") +
    xlab("Selectivity index")

ggsave("smiVtun.pdf", smivsparse, device=cairo_pdf, width=2.5, height=5.3)