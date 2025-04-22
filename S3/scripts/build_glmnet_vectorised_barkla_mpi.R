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
library(doParallel)
library(foreach)
library(glmnet)
library(Rmpi)

####################################################################################
# Options and global definitions used in all runs to keep training sets consistent #
####################################################################################

# Set parallelisation
workers <- mpi.universe.size() - 1
cat("Number of workers = ", workers, "\n")
cl <- parallelly::makeClusterMPI(workers, autoStop = TRUE)
registerDoParallel(cl)

holdout_cluster_grid <- list.files(path = "S3/data/full/holdout_clusters", pattern = "labels.csv") %>%
  gsub("ex_|_labels.csv", "", .) %>%
  str_split(., "_") %>% 
  do.call(rbind.data.frame, .) %>%
  set_colnames(c("subtype", "minseqid", "C"))

# run_date <- format(Sys.time(), "%Y_%m_%d")

# dir.create(paste0("results_", run_date), showWarnings = FALSE, recursive = TRUE)

cluster_sets <- holdout_cluster_grid %>% 
  select(minseqid, C) %>% 
  distinct %>% 
  mutate(cluster_set = paste0(minseqid, "_", C)) %>% 
  pull(cluster_set)
  
cluster_sets <- c("70_7")

#######################
# Define ML procedure #
#######################

glmnet_fun <- function(subtype){
  
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
  
  # Read-in predefined folds for 5-fold cross-validation
  fold_indices <- readRDS("S3/data/fold_indices_list.rds") %>% .[[subtype]]
  
  #################
  # Run ML models #
  #################
  
  # Train and validate glmnet (tuning mtry, min.node.size parameters) through 5-fold cross-validation
  # Store result as list of n ensemble models
  train(x = train %>% select(all_of(preds)),
        y = train %>% pull(label),
        method = "glmnet",
        preProc = c("center", "scale"),
        metric = "Kappa",
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
      alpha = c(0,0.5,1),         # mixing parameter (0 = ridge, 0.5 = elastic, 1 = lasso)
      lambda = c(0.001, 0.01, 0.1, 1, 10))     # regularisation parameter
  )%>%
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
                         "glmnet")) %dopar% {
    
    ##############################################################
    # Save list of ML models, each holding out a holdout subtype #
    ##############################################################
    
    if (!dir.exists(paste0("results_", run_date, "/", cluster_set))){
      dir.create(paste0("results_", run_date, "/", cluster_set))
    }
    
    Map(f = glmnet_fun,  subtype = unique(holdout_cluster_grid$subtype)) %>% 
      suppressWarnings() %>% 
      saveRDS(file=paste0("results_", run_date, "/", cluster_set, 
                          "/glmnet_list_", gsub("allflu_|_pt.*.rds", "", featset), "_pt_", focgene, ".rds"))
    
  }

stopCluster(cl)
mpi.exit()