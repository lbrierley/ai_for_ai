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
library(kernlab)

####################################################################################
# Options and global definitions used in all runs to keep training sets consistent #
####################################################################################

# Set parallelisation
cl <- makePSOCKcluster(detectCores() - 1)
registerDoParallel(cl)
clusterSetRNGStream(cl, 1429)

holdout_cluster_grid <- list.files(path = "S3/data/full/holdout_clusters", pattern = "labels.csv") %>%
  gsub("ex_|_labels.csv", "", .) %>%
  str_split(., "_") %>% 
  do.call(rbind.data.frame, .) %>%
  set_colnames(c("subtype", "minseqid", "C"))

run_date <- format(Sys.time(), "%Y_%m_%d")

dir.create(paste0("results_", run_date), showWarnings = FALSE, recursive = TRUE)

cluster_sets <- holdout_cluster_grid %>% 
  select(minseqid, C) %>% 
  distinct %>% 
  mutate(cluster_set = paste0(minseqid, "_", C)) %>% 
  pull(cluster_set)
  
cluster_sets <- "70_7"

#######################
# Define ML procedure #
#######################

svm_fun <- function(subtype){
  
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
  
  # Train and validate SVM (tuning mtry, min.node.size parameters) through 5-fold cross-validation
  # Store result as list of n ensemble models
  train(form = label ~ .,
        data = train %>% select(-cluster_rep),
        method = "svmRadialWeights",
        preProc = c("center", "scale"),
        metric = "Kappa",
        trControl = trainControl(method = "repeatedcv", 
                                 index = fold_indices,
                                 number = 5,
                                 repeats = 1,
                                 #verboseIter = TRUE,
                                 classProbs = TRUE,
								 savePredictions = TRUE),
        tuneGrid = expand.grid(
          sigma = c(0.001, 0.01, 0.1),
          C = c(1, 3, 9),
          Weight = c(3, 9))
        ) %>%
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
                         "kernlab")) %dopar% {
    
    ##############################################################
    # Save list of ML models, each holding out a holdout subtype #
    ##############################################################
    
    if (!dir.exists(paste0("results_", run_date, "/", cluster_set))){
      dir.create(paste0("results_", run_date, "/", cluster_set))
    }
    
    Map(f = svm_fun,  subtype = unique(holdout_cluster_grid$subtype)) %>% 
      suppressWarnings() %>% 
      saveRDS(file=paste0("results_", run_date, "/", cluster_set, 
                          "/svm_list_", gsub("allflu_|_pt.*.rds", "", featset), "_pt_", focgene, ".rds"))
    
  }

stopCluster(cl)