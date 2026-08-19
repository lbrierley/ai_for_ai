library(tidyverse)

or <- read.csv("S3\\data\\full\\mapping\\kmer_or_results\\prot\\HA_position_weighted_mean_abs_logOR.csv")

atom <- read.csv("S3\\data\\full\\mapping\\alphafold_struct\\cif_csv\\HA_consensus_base.csv", header = F)

atom %>% 
  left_join(or %>% select(position, w_mean_abs_logOR), by = c("V29" = "position")) %>% 
  select(-V27) %>%
  dplyr::rename(V27 = w_mean_abs_logOR) %>%
  select(paste0("V", seq(1:33))) %>% # arrange columns as in original %>%
  select_if(~ !any(is.na(.))) %>%
  write.table("S3\\data\\full\\mapping\\alphafold_struct\\cif_csv\\HA_consensus_new.csv", sep=",", col.names=FALSE, row.names=FALSE)



atom <- read.csv("S3\\data\\full\\mapping\\alphafold_struct\\cif_csv\\HA_consensus_trimer_base.csv", header = F)

atom %>% 
  left_join(or %>% select(position, w_mean_abs_logOR), by = c("V29" = "position")) %>% 
  select(-V27) %>%
  dplyr::rename(V27 = w_mean_abs_logOR) %>%
  select(paste0("V", seq(1:33))) %>% # arrange columns as in original %>%
  select_if(~ !any(is.na(.))) %>%
  write.table("S3\\data\\full\\mapping\\alphafold_struct\\cif_csv\\HA_consensus_trimer_new.csv", sep=",", col.names=FALSE, row.names=FALSE)
