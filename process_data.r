library(tidyverse)

rm(list = ls())

summary_files <- list.files(path = '~/rds/projects/kozlovlab_rds2/live/Grace/data/processed/combined_stims/units_merged/csv_files', pattern="*_metrics.csv", recursive=TRUE, full.names = TRUE)

raw_dat <- as_tibble(do.call(rbind, lapply(summary_files, read.csv))) %>% 
    mutate(is_combi = as.logical(is_combi))

dat <- raw_dat %>% 
    mutate(sig_fr = if_else(p_val <= 0.05, TRUE, FALSE)) %>% 
    mutate(sig_fr = if_else(is.na(sig_fr), ifelse(fr > 0.15, TRUE, FALSE), sig_fr)) %>% # for instances when standalone clips weren't included 
    group_by(mouse, id, stim_pair) %>% 
    mutate(valid = if_else(all(sig_fr == FALSE), FALSE, TRUE)) %>% 
    filter(valid, !is.na(smi), smi > -Inf) %>% 
    group_by(mouse, id) %>% 
    mutate(nnID = factor(group_indices())) %>% 
    group_by(nnID, stim_pair) %>% 
    mutate(above_threshold = sum(sig_fr)-1)

# examine responses where there was one or two stims below threshold 
# keep set if Rab exceeds Ra by 25%, and Ra is more than 0.2
dominant <- dat %>% 
    filter(above_threshold < 2, !is_combi) %>%
    mutate(Ra = max(fr)) %>% 
    select(mouse, session, nnID, ch, unit_no, stim_pair, base_pair, Ra, above_threshold)

combined <- dat %>% 
    filter(above_threshold < 2, is_combi) %>% 
    mutate(Rab = fr) %>% 
    select(mouse, session, nnID, ch, unit_no, stim_pair, base_pair, Rab, above_threshold)

exclude_stims <- combined %>% 
    left_join(dominant, by=c("mouse", "session", "nnID", "ch", "unit_no", "stim_pair", "base_pair", "above_threshold")) %>% 
    distinct() %>% 
    mutate(frIncrease = (Rab - Ra)/Ra) %>% 
    filter(frIncrease < 0.25 | Ra < 0.2) %>% 
    select(mouse, session, nnID, ch, unit_no, stim_pair, base_pair)

write.csv(exclude_stims, "./exclusions.csv", row.names=FALSE)

# but include these ones that have passed manual inspection 
keep <- read.csv("./retain.csv") %>% 
    mutate(nnID = as.factor(nnID))

exclude_stims2 <- exclude_stims %>% 
    anti_join(keep, by=c("mouse", "session", "nnID", "ch", "unit_no", "stim_pair", "base_pair"))

filtered_dat <- dat %>%
    anti_join(exclude_stims2, by=c("mouse", "session", "nnID", "ch", "unit_no", "stim_pair", "base_pair")) 

write_rds(filtered_dat, "./data/filtered_df.rds")

