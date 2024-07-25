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
# library(parallel)
# library(doParallel)
library(foreach)
library(xgboost)

####################################################################################
# Options and global definitions used in all runs to keep training sets consistent #
####################################################################################

## FOR XGBOOST, USE NATIVE PARALLELISATION INSTEAD
# # Set parallelisation
# cl <- makePSOCKcluster(detectCores() - 1)
# registerDoParallel(cl)
# clusterSetRNGStream(cl, 1429)

holdout_cluster_grid <- list.files(path = "holdout_clusters/", pattern = "labels.csv") %>%
  gsub("ex_|_labels.csv", "", .) %>%
  str_split(., "_") %>% 
  do.call(rbind.data.frame, .) %>%
  set_colnames(c("subtype", "minseqid", "C"))

run_date <- "16_02_24"

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

xgb_fun <- function(subtype){
  
  ####################################################
  # Train models on each feature set on each protein #
  ####################################################
  
  labels <- read.csv(paste0("holdout_clusters/ex_", subtype, "_", cluster_set, "_labels.csv")) %>% 
    select(cluster_rep, label) %>%
    mutate(label = factor(case_when(label == "zoon" ~ "hzoon", label == "nz" ~ "nz")) # Rearrange factor levels for better compatibility with model functions
    )
  
  # Load feature sets for training data clusters, rename features to indicate gene/protein being modelled
  train <- left_join(labels,
                     readRDS(paste0("/users/lbrier/mlready/allflu_", featset, "_pt_", focgene, ".rds")) %>% 
                       select(-any_of(c("segment", "cds_id", "enc", "GC_content"))),
                     by = c("cluster_rep" = "gid")) %>%
    rename_with(~paste(., focgene, sep = "_"), -c(cluster_rep, label))
  
  # Specify variables used
  preds <- train %>% select(-label, -cluster_rep) %>% remove_constant %>% names
  
  # Create folds for 5-fold cross-validation
  fold_indices <- readRDS("S3/data/fold_indices_list.rds") %>% .[[subtype]]
  
  ##############
  # Run models #
  ##############
  
  # Train and validate xgboost (tuning max_depth, eta, min_child_weight parameters) through 5-fold cross-validation
  # Store result as list of n ensemble models
  train(form = label ~ .,
        data = train %>% select(-cluster_rep),
        method = "xgbTree",
        preProc = c("center", "scale"),
        metric = "Kappa",
        colsample_bylevel = c(2/5),
        scale_pos_weight = 10,
        # weights = ifelse(train$label == "nz",                 # CARET CLASS WEIGHTS DOES NOT PLAY WELL WITH XGBOOST
        #                  (1/table(train$label)[2]) * 0.5,
        #                  (1/table(train$label)[1]) * 0.5),
        trControl = trainControl(method = "repeatedcv", 
                                 index = fold_indices,
                                 number = 5,
                                 repeats = 1,
                                 #verboseIter = TRUE,
                                 classProbs = TRUE,
								 savePredictions = TRUE),
        tuneGrid = expand.grid(
          nrounds = c(2000),
          max_depth = c(3,4,5), 
          eta = c(0.1,0.3,1), 
          gamma = c(0),          
          subsample = c(1), 
          colsample_bytree = c(4/5),  
          min_child_weight = c(1,3,5))    
  ) %>%
    return()
}

######################################################
# Repeat for each cluster set, gene, and feature set #
######################################################

foreach (cluster_set = cluster_sets) %:% 
  foreach (focgene = c("HA", "M1")) %:% 
  foreach (featset = list.files(path = "mlready", pattern = focgene) %>% gsub("allflu_|_pt.*.rds", "", .),
           .packages = c("caret",
                         "dplyr",
                         "janitor",
                         "magrittr", 
                         "xgboost")) %do% {
    
    ##############################################################
    # Save list of ML models, each holding out a holdout subtype #
    ##############################################################
    
    if (!dir.exists(paste0("results_", run_date, "/", cluster_set))){
      dir.create(paste0("results_", run_date, "/", cluster_set))
    }
    
    Map(f = xgb_fun,  subtype = unique(holdout_cluster_grid$subtype)) %>% 
      suppressWarnings() %>% 
      saveRDS(file=paste0("results_", run_date, "/", cluster_set, 
                          "/xgb_list_", gsub("allflu_|_pt.*.rds", "", featset), "_pt_", focgene, ".rds"))
    
  }

# stopCluster(cl)