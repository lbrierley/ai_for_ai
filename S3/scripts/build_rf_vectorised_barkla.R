#######################
# Load packages, data #
#######################

library(caret)
library(magrittr)
library(janitor)
library(dplyr)
library(tidyr)
library(forcats)
library(stringr)
library(parallel)
library(doParallel)
library(foreach)
library(ranger)

####################################################################################
# Options and global definitions used in all runs to keep training sets consistent #
####################################################################################

# Set parallelisation
cl <- makePSOCKcluster(detectCores() - 1)
registerDoParallel(cl)
clusterSetRNGStream(cl, 1429)

holdout_cluster_grid <- list.files(path = "S3/data/full/holdout_clusters/", pattern = "labels.csv") %>%
  gsub("ex_|_labels.csv", "", .) %>%
  str_split(., "_") %>% 
  do.call(rbind.data.frame, .) %>%
  set_colnames(c("subtype", "minseqid", "C"))

run_date <- "14_02_24"

dir.create(paste0("results_", run_date))

cluster_sets <- holdout_cluster_grid %>% 
  select(minseqid, C) %>% 
  distinct %>% 
  mutate(cluster_set = paste0(minseqid, "_", C)) %>% 
  pull(cluster_set)

cluster_sets <- c("70_7")

#######################
# Define ML procedure #
#######################

rf_fun <- function(subtype){
  
  ####################################################
  # Train models on each feature set on each protein #
  ####################################################
  
  labels <- read.csv(paste0("S3/data/full/holdout_clusters/ex_", subtype, "_", cluster_set, "_labels.csv")) %>% 
    select(cluster_rep, label) %>%
    mutate(label = factor(case_when(label == "zoon" ~ "hzoon", label == "nz" ~ "nz")) # Rearrange factor levels for better compatibility with model functions
    )
  
  # Load feature sets for training data clusters, rename features to indicate gene/protein being modelled
  train <- left_join(labels,
                     readRDS(paste0("S3/data/full/mlready/allflu_", featset, "_pt_", focgene, ".rds")) %>% 
                       select(-any_of(c("segment", "cds_id", "enc", "GC_content"))),
                     by = c("cluster_rep" = "gid")) %>%
    rename_with(~paste(., focgene, sep = "_"), -c(cluster_rep, label))
  
  # Specify variables used
  preds <- train %>% select(-label, -cluster_rep) %>% remove_constant %>% names
  
  # Create folds for 5-fold cross-validation
  set.seed(1657)
  fold_indices <- createMultiFolds(train$label, k = 5, times = 1)
  
  #################
  # Run ML models #
  #################
  
  # Train and validate RF (tuning mtry, min.node.size parameters) through 5-fold cross-validation using ranger
  # Store result as list of n ensemble models
  train(x = train %>% select(all_of(preds)),
        y = train %>% pull(label),
        method = "ranger",
        preProc = c("center", "scale"),
        metric = "Kappa",
        num.trees = 1000,
        weights = ifelse(train$label == "nz",
                         (1/table(train$label)[2]) * 0.5,
                         (1/table(train$label)[1]) * 0.5),
        importance = "impurity",
        trControl = trainControl(method = "repeatedcv", 
                                 index = fold_indices,
                                 number = 5,
                                 repeats = 1,
                                 #verboseIter = TRUE,
                                 classProbs = TRUE,
								 savePredictions = TRUE),
        tuneGrid = expand.grid(
          .splitrule = "gini",
          #.min.node.size = 5,
          .min.node.size = seq(from = 5, to = 45, length = 3),
          .mtry = round(sqrt(length(preds))))) %>%
    return()
}

######################################################
# Repeat for each cluster set, gene, and feature set #
######################################################

foreach (cluster_set = cluster_sets) %:% 
  foreach (focgene = c("HA", "M1", "NA", "NP", "NS1", "PA", "PB1", "PB2")) %:% 
  foreach (featset = list.files(path = "S3/data/full/mlready", pattern = focgene) %>% gsub("allflu_|_pt.*.rds", "", .),
           .packages = c("caret",
                         "dplyr",
                         "janitor",
                         "magrittr", 
                         "ranger")) %dopar% {
    
    ##############################################################
    # Save list of ML models, each holding out a holdout subtype #
    ##############################################################
    
    if (!dir.exists(paste0("results_", run_date, "/", cluster_set))){
      dir.create(paste0("results_", run_date, "/", cluster_set))
    }
    
    Map(f = rf_fun, subtype = unique(holdout_cluster_grid$subtype)) %>% 
      suppressWarnings() %>% 
      saveRDS(file=paste0("results_", run_date, "/", cluster_set, 
                          "/rf_list_", gsub("allflu_|_pt.*.rds", "", featset), "_pt_", focgene, ".rds"))
    
  }

stopCluster(cl)