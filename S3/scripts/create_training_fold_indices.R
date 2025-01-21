#######################
# Load packages, data #
#######################

library(caret)
library(magrittr)
library(janitor)
library(dplyr)
library(tidyr)

#########
# Setup #
#########

holdout_cluster_grid <- list.files(path = "S3\\data\\full\\holdout_clusters\\", pattern = "labels.csv") %>%
  gsub("ex_|_labels.csv", "", .) %>%
  stringr::str_split(., "_") %>% 
  do.call(rbind.data.frame, .) %>%
  magrittr::set_colnames(c("subtype", "minseqid", "C"))

cluster_set <- "70_7"

fold_indices_list <- list()

#############################################################
# Extract training folds from previously run model and save #
#############################################################

model_list <- list.files(path = "S3/analysis", pattern = ".rds", recursive = TRUE, full.names = TRUE) %>%
  .[grepl(cluster_set, .)] %>%
  .[grepl("results_14_02_24", .)] %>%
  .[1] %>%
  readRDS

for(i in 1:length(unique(holdout_cluster_grid$subtype))){
  
  fold_indices_list[[i]] <- model_list[[i]]$control$index
  
}
  
# #############################################
# # Generate training folds a priori and save #
# #############################################
# 
# fold_fun <- function(subtype){
#   
#   labels <- read.csv(paste0("S3/data/full/holdout_clusters/ex_", subtype, "_", cluster_set, "_labels.csv")) %>% 
#     select(cluster_rep, label) %>%
#     mutate(label = factor(case_when(label == "zoon" ~ "hzoon", label == "nz" ~ "nz")) # Rearrange factor levels for better compatibility with model functions
#     )
#   
#   # Create folds for 5-fold cross-validation
#   set.seed(1657)
#   createMultiFolds(labels$label, k = 5, times = 1) %>% 
#     return()
# }
# 
# fold_indices_list <- Map(f = fold_fun, subtype = unique(holdout_cluster_grid$subtype)) 
names(fold_indices_list) <- unique(holdout_cluster_grid$subtype)
fold_indices_list %>% saveRDS("S3/data/fold_indices_list.rds")