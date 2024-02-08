library(tidyverse)
library(ggplot2)
library(ggExtra)
library(GGally)
library(ggpubr)
library(lmerTest)
library(ggeffects)
library(ARTool)
library(ggthemes)
library(GGally)
library(cowplot)
library(patchwork)
#library(MASS) #if doing boxcox transformations 

rm(list = ls())

filtered_dat <- readRDS("./data/filtered_df.rds")

my_t <- theme_cowplot(font_size=11) + 
    theme(text = element_text(family = "Sans serif"))

# smi groupings 

smis <- filtered_dat%>% 
    filter(!is_combi) %>% 
    mutate( stim_no = 1:n()) %>% 
    filter(stim_no == 1) %>% 
     %>% 
    select(mouse, nnID, session, smi, summation)

# adaptation and firing rate metrics 
fr <- filtered_dat %>%
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

euclidean_dist <- read_csv("euclidean.csv") %>% # this is the euclidean distance measured between 
    select(mouse, nnID, stim_pair, dists) %>% 
    mutate(nnID = as.factor(nnID), mouse = as.character(mouse)) %>% 
    rename(euclidean_d = dists)

windowed_dist <- read_csv("dist_short_specs.csv") %>% # distance metrics for stimuli within 55ms of spiking 
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
    #mutate(smi.t = 1/(smi + 1.88)) # https://sscc.wisc.edu/sscc/pubs/RegDiag-R/normality.html#fn2 
    mutate(smi.t = smi + 0.8800001) %>% 
    filter(complete.cases(.)) %>% 
    filter(smi < 7) 

cor.test(formula = ~ tot_fr + smi, data=firingNfeatures, method="spearman", exact=FALSE)
cor.test(formula = ~ max_recovery + smi, data=firingNfeatures, method="spearman", exact=FALSE)
cor.test(formula = ~ temp + smi, data=firingNfeatures, method="spearman", exact=FALSE)
cor.test(formula = ~ euclidean_d + smi, data=firingNfeatures, method="spearman", exact=FALSE)
cor.test(formula = ~ d_win + smi, data=firingNfeatures, method="spearman", exact=FALSE)
cor.test(formula = ~ d_win + smi, data=firingNfeatures %>% filter(temp < 0), method="spearman", exact=FALSE)

boxdat <- firingNfeatures %>% 
    mutate(summation=as.factor(summation)) %>% 
    dplyr::select(euclidean_d, max_recovery, summation, nnID, d_win, nnID, tot_fr, smi) %>% 
    rename(AB_input=tot_fr, spectral_feature_dist = d_win, adaptation=max_recovery) %>% 
    pivot_longer(!c(summation, nnID, smi), names_to = "variable", values_to = "val")

boxdat$summation <- recode(boxdat$summation, sublinear = "sub")

b <- ggplot(boxdat %>% filter(variable %in% c("AB_input", "adaptation", "euclidean_d", "spectral_feature_dist")), aes(summation, val)) + 
    geom_boxplot(outlier.shape = NA) + 
    #geom_jitter(fill="white", size=0.6, alpha=0.8, width=0.2, shape=21) + 
    facet_wrap(~variable, scales="free", nrow=4) + 
    xlab("") +
    coord_flip()+ 
    my_t + 
    theme(panel.spacing=unit(5, "lines"))

d <- ggplot(boxdat %>% filter(variable %in% c("AB_input", "adaptation", "euclidean_d", "spectral_feature_dist")), aes(val, smi)) + 
    geom_point(shape=21, size=0.6) + 
    facet_wrap(~variable, scales="free",nrow=4) + 
    stat_cor(method="spearman") + 
    geom_hline(yintercept = c(0, 1), linetype='dashed',col="#b4aea9") +
    ylab("SmI") + 
    xlab("") + 
    my_t

analyses <- plot_grid(d, b)

ggsave("box_plots.pdf", analyses, width=5, height=7.5, device=cairo_pdf)

# compare euclidean distances 
library(emmeans)
art.mod <- art(val ~ summation +(1|nnID), boxdat %>% filter(variable == "euclidean_d")) # significant
contrast(emmeans(artlm(art.mod, "summation"), ~summation), method="pairwise")
%art.con(art.mod, "summation", adjust="holm") %>%  summary()

stat.test <- tibble::tribble(
  ~group1, ~group2,   ~p.adj,
    "MAX", "sub", 0.013,
  )

anova(art(val ~ summation + (1|nnID), boxdat %>% filter(variable == "adaptation"))) # n.s. 
anova(art(val ~ summation + (1|nnID), boxdat %>% filter(variable == "AB_input") %>% drop_na(val)))

p1 <- ggplot(boxdat %>% filter(variable == "euclidean_d"), aes(summation, val)) + 
    geom_boxplot(outlier.shape = NA) + 
    #geom_jitter(fill="white", size=0.6, alpha=0.8, width=0.2, shape=21) + 
    facet_wrap(~variable, scales="free") + 
    ylab("euclidean dist between spectrograms") + 
    theme(
    text=element_text(family="Helvetica", size=12), 
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    axis.text.y=element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_line(colour = "grey50"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(linetype = "dashed"),
    panel.background = element_rect(fill='transparent'),
    plot.background = element_rect(fill='transparent', color=NA)
    ) + 
    ylim(0, 90) + 
    ylab("") + 
    xlab("") + 
    stat_pvalue_manual(
    stat.test, 
    y.position = 80, step.increase = 0.3,
    label = "p.adj"
    )

dtplot <- ggplot(boxdat %>% filter(variable %in% c("AB_input", "adaptation", "euclidean_d", "spectral_feature_dist")), aes(val, smi)) + 
    geom_point(shape=21, size=0.6) + 
    facet_wrap(~variable, scales="free",ncol=2) + 
    stat_cor(method="spearman", label.x.npc=0.7, aes(label = ..r.label..)) + 
    geom_hline(yintercept = c(0, 1), linetype='dashed',col="#b4aea9") +
    theme(
    aspect.ratio = 1,
    text=element_text(family="Helvetica", size=18), 
    axis.ticks = element_blank(),
    axis.line = element_line(colour = "grey50"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.background = element_rect(fill='transparent'),
    plot.background = element_rect(fill='transparent', color=NA),
    ) + 
    ylab("SmI") + 
    xlab("") 

p2 <- dtplot + inset_element(p1, 0.1, 0.16, 0.25, 0.46, clip=FALSE)

ggsave("scatter_plots.pdf", p2)

# transform max_recovery and euclidean_d
#https://www.css.cornell.edu/faculty/dgr2/_static/files/R_html/Transformations.html
BCTransform <- function(y, lambda=0) {
    if (lambda == 0L) { log(y) }
    else { (y^lambda - 1) / lambda }
    }

BCTransformInverse <- function(yt, lambda=0) {
    if (lambda == 0L) { exp(yt) }
    else { exp(log(1 + lambda * yt)/lambda) }
    }

# transform smi instead 
 
transform.mod <- lm(smi.t  ~ tot_fr + max_recovery + temp + d_win, data=firingNfeatures %>% filter(smi < 7))
bc <- boxcox(transform.mod) # add 0.880000001 
bc.power <- (bc.power <- bc$x[which.max(bc$y)])
firingNfeatures$smi.bc <- BCTransform(firingNfeatures$smi.t, bc.power)
par(mfrow=c(2,2))
hist(firingNfeatures$smi.bc, breaks=30);

hist(firingNfeatures$tot_fr)

bc.mod2 <- lm(smi.bc ~ tot_fr + max_recovery + d.au, data=firingNfeatures) # euclidean distance wasn't significant 
bc.mixed.mod <- lmer(smi.bc ~ d_win + tot_fr + max_recovery + temp + (1|nnID), data=firingNfeatures %>% drop_na()) # slightly better qq norm 
par(mfrow=c(2,2))
plot(bc.mod2)

qqnorm(resid(bc.mixed.mod))
qqline(resid(bc.mixed.mod))

library(partR2)

R2_BMa <- partR2(bc.mixed.mod, partvars = c("tot_fr", "d_win", "temp"), 
                  R2_type = "marginal", nboot = 100) # distance measures only account for 2% of the variance 

#kruskal.test(val ~ summation, boxdat %>% filter(variable == "spectral_feature_dist") %>% drop_na())
# cannot do kruskal test because observations are correlated within each neuron 

library("performance")
r2_nakagawa(bc.mixed.mod)

# try doing bootstrapping 

linear.mixed.mod <- lmer(smi.bc ~ tot_fr + d.au + (1+tot_fr|nnID), data=firingNfeatures %>% filter(smi < 7) %>% drop_na())

mbootstraps <- bootMer(linear.mixed.mod, FUN=fixef, nsim=1000, verbose=T)

check_model(bc.mixed.mod)
    model_performance(linear.mixed.mod)


###### data frame with wavnames to compare with the model 

spec_features <- read_csv("euclidean.csv") %>%
    select(mouse, nnID, stim_pair, A, B, delay) %>% 
    mutate(nnID = as.factor(nnID), mouse = as.character(mouse)) %>% 
    left_join(firingNfeatures, by=c("nnID", "stim_pair", "mouse")) %>% 
    rename(smi_exp = smi)

write.csv(spec_features, "specfeatures.csv", row.names=FALSE)

### sparsity 

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
