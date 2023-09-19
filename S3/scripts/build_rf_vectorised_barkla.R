#######################
# Load packages, data #
#######################

library(caret)
library(e1071)
library(matrixStats)
library(magrittr)
library(pROC)
library(randomForest)
library(janitor)
library(dplyr)
library(tidyr)
library(readr)
library(purrr)
library(forcats)
library(stringr)
library(tibble)
library(ranger)
library(parallel)
library(doParallel)

for (focgene in c("HA", "M1", "NA", "NP", "NS1", "PA", "PB1", "PB2")){
  
  featsets <- list.files(path = "mlready", pattern = focgene)
  
  for (i in 1:length(featsets)) {
    
    # Define feature name
    featname <- featsets[i] %>% gsub("allflu_|_pt.*.rds", "", .)
    
    # Load in data
    train <- inner_join(read.csv("cluster_rep_labels.csv") %>% 
                          select(cluster_rep, label),
                        readRDS(paste0("/users/lbrier/mlready/",featsets[i])) %>% 
                          select(-any_of(c("segment", "cds_id", "enc", "GC_content"))),
                        by = c("cluster_rep" = "gid")) %>%
      mutate(label = factor(label))
    
    ##################################################################################
    # Prepare data frame for modelling and define set of variables used in the model #
    ##################################################################################
    
    # Set options for local analyses
    set.seed(1657)
    
    # Specify variables used
    preds <- train %>% select(-label, -cluster_rep) %>% remove_constant %>% names
    
    # Create outer folds for 5-fold cross validation
    outer_fold_dfs <- createMultiFolds(train$label, k = 5, times = 1)
    outer_fold_dfs %<>% lapply(function(x)
      train %>% slice(x)
    )
    
    # Create inner folds for 10-fold cross-validation
    inner_fold_indices <- lapply(outer_fold_dfs, function(x)
      createMultiFolds(x$label, k = 10, times = 1)
    )
    
    ######################
    # Run random forests #
    ######################
    
    # Set parallelisation
    
    cl <- makePSOCKcluster(5)
    registerDoParallel(cl)
    clusterSetRNGStream(cl, 1429)
    
    # Train and validate RF (tuning mtry, min.node.size parameters) through 10-fold cross-validation using ranger
    
    timer_start <- Sys.time()
    
    # Store result as list of n ensemble models
    rf_list <- Map(function(outer, inner) 
      
      train(x = outer %>% select(all_of(preds)),
            y = outer %>% pull(label),
            method = "ranger",
            preProc = c("center", "scale"),
            metric = "Accuracy",
            num.trees = 1000,
            importance = "impurity",
            trControl = trainControl(method = "repeatedcv", 
                                     index = inner,
                                     number = 10,
                                     repeats = 1,
                                     #verboseIter = TRUE,
                                     classProbs = TRUE),
            tuneGrid = expand.grid(
              .splitrule = "gini",
              #.min.node.size = 5,
              .min.node.size = seq(from = 5, to = 45, length = 3),
              .mtry = round(sqrt(length(preds))))
      ),
      outer = outer_fold_dfs,
      inner = inner_fold_indices
    )
    
    timer_end <- Sys.time()
    
    ##################################
    # Save models and computing time #
    ##################################
    
    write(paste0(focgene, ",", featname, ",", as.numeric(timer_end - timer_start, units = "hours")), "ml_times_hours.txt", append = TRUE)
    
    rf_list %>% saveRDS(file=paste0("rf_list_", featname, "_", focgene, "_", format(Sys.time(), "%d_%m_%y"), ".RData"))
    
  }
}