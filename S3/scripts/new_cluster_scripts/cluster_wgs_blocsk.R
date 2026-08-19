library(dplyr)
library(tidyr)
library(magrittr)
library(purrr)
library(tibble)
library(tidyverse)

# Read set of wgs, cluster codes, and summed pairwise JSD
allflu_wgs_ref <- read.csv("S3\\data\\segmentwise\\Raw_Sequences\\allflu_wgs_df.csv") 
clust_codes <- read.csv("S3\\data\\segmentwise\\cds_wgs_clust_codes.csv")
JSD_wgs_matrix <- readRDS("S3\\data\\segmentwise\\clust_dist\\pairwise_JSD_wgs.rds")

# Read avian/zoon mixed status for clusters of all segments of wgs
mix_files <- list.files(cluster_dir, pattern = "_members\\.csv$", full.names = TRUE)
mix_list <- lapply(seq_along(mix_files), function(x) {
  read.csv(mix_files[x]) %>%
    
    #####
  # temporary fix to make sure each gid only appears in one cluster. preferentially take mixes to calc worst case scenario
  #####
  arrange(gid, -mix) %>%
    group_by(gid) %>% 
    dplyr::slice(1) %>%
    ungroup() %>%
    #####
  select(gid, mix) %>%
    rename_with(~ paste0("mix_", x), .cols = mix)
})

reduced_mix_list <- reduce(mix_list, left_join, by = "gid") %>%
  distinct()

# Combine mixed status with wgs gids
allflu_wgs_ref <- allflu_wgs_ref %>%
  select(gid, subtype, label) %>%
  distinct() %>%
  mutate(label = factor(case_when(label == "zoon" ~ "hzoon", label == "nz" ~ "nz"))) %>% # Rearrange factor levels for better compatibility with model functions
  left_join(reduced_mix_list, by = "gid") %>%
  mutate(any_mix = case_when(
    if_all(starts_with("mix_"), ~ .x == 1)  ~ 0,
    if_any(starts_with("mix_"), ~ .x != 1)  ~ 1
  ))

# Check table of mixed status
allflu_wgs_ref %>% with(table(label, any_mix))

# Cut to those with unmixed clusters only
allflu_wgs_ref_nomix <- allflu_wgs_ref %>% filter(!(label == "nz" & any_mix == 1))
JSD_wgs_matrix_nomix <- JSD_wgs_matrix[allflu_wgs_ref_nomix$gid, allflu_wgs_ref_nomix$gid]




## Hierarchical clustering and fold definition

# Hierarchically from hclust package
hclust_wgs <- JSD_wgs_matrix_nomix %>%
  as.dist(diag = TRUE) %>%
  hclust(method="ward.D2") # use Ward's minimum variance method to force more compact clustering


# Cut into different values of j folds and check distribution (remembering no-preclustering has been done here, so lots of very similar)
cut_k <- hclust_wgs %>% cutree(k = c(10)) 

cut_df <- data.frame(gid = hclust_wgs$labels) %>% 
  left_join(allflu_wgs_ref_nomix %>% select(gid, label) %>% distinct)

# At this stage review manually

table(fold = cut_k, label = cut_df$label) %>% 
  addmargins(1) 

table(fold = cut_k, label = cut_df$label) %>% 
  prop.table(margin = 2) %>%
  round(2) %>%
  addmargins(1) %>% 
  addmargins(2)

data.frame(cut_df, fold = cut_k) %>% left_join(allflu_wgs_ref_nomix %>% select(gid, subtype) %>% distinct) %>% with(., table(subtype, fold))

data.frame(cut_df, fold = cut_k) %>% left_join(allflu_wgs_ref_nomix %>% select(gid, subtype, label) %>% distinct) %>% filter(fold == 1) %>% with(., table(subtype, label))

# Manually designate test fold

folds_wgs <- hclust_wgs %>%
  cutree(k = 10) %>%
  as.data.frame() %>%
  bind_cols(gid = hclust_wgs$labels) %>%
  set_colnames(c("fold", "gid")) %>%
  select(gid, fold)

folds_wgs %<>%
  mutate(fold = as.character(fold)) %>%
  mutate(fold = ifelse(fold == "1", "test", fold)) %>%
  mutate(fold = ifelse(!is.na(as.numeric(fold)), as.character(dense_rank(fold)), fold)) # Renumber non-test folds from 1

folds_wgs %>% saveRDS("S3\\data\\segmentwise\\clust_dist\\folds_wgs.rds")
