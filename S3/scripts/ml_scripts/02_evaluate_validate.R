#################
# Load packages #
#################

library(caret)
library(e1071)
library(matrixStats)
library(magrittr)
library(pROC)
library(janitor)
library(dplyr)
library(tidyr)
library(purrr)
library(forcats)
library(stringr)
library(tibble)
library(parallel)
library(doParallel)
library(foreach)
library(kernlab)
library(xgboost)
library(ranger)
library(glmnet)

######################
# Set options, paths #
######################

results_date <- "17_07_26"
results_path <- paste0("S3\\analysis\\results_", results_date)

cluster_dir <- "S3/data/segmentwise/clust_processed/"
mlready_dir <- "S3/data/segmentwise/mlready/"

####################
# Define functions #
####################

# helper function
canon_id <- function(x) {
  as.character(x) %>%
    trimws() %>%
    str_replace_all("\\s+", " ")
}

#######################
# Load model pointers #
#######################

model_list <- list.files(path = paste0("S3\\analysis\\results_", results_date))
methods <- model_list %>% gsub("_.*", "", .) %>% unique
featsets <- model_list %>% gsub("^[^_]*_(.*)_pt.*$", "\\1", .) %>% unique
genes <- c("PB2","PB1","PA","HA","NP","NA","M1","NS1")

# Specific feature set files
feature_set_files <- list.files(mlready_dir, pattern = "\\.rds$", full.names = TRUE)

# Specify fold ids for all segments up front
folds_list <- readRDS("S3/data/segmentwise/clust_dist/folds_list.rds")

folds_list <- lapply(seq_along(folds_list), function(x) {
  fold_df <- folds_list[[x]]
  fold_df$gid <- canon_id(fold_df$gid)
  return(fold_df)
})

# Specify labels for all segments up front
label_files <- list.files("S3\\data\\segmentwise\\clust_processed", pattern = "_labels\\.csv$", full.names = TRUE)

labels_list <- lapply(seq_along(label_files), function(x) {
  read.csv(label_files[x]) %>%
    select(cluster_rep, label, mix) %>%
    mutate(
      cluster_rep = canon_id(cluster_rep), # Rearrange factor levels for better compatibility with model functions
      label = factor(case_when(
        label == "zoon" ~ "hzoon",
        label == "nz" ~ "nz",
        TRUE ~ as.character(label)
      ))
    ) %>%
    distinct(cluster_rep, .keep_all = TRUE) %>%
    rename(gid = cluster_rep)
})

names(folds_list) <- genes
names(labels_list) <- genes

############################################################
# Load in models, extract and store optimal parameter sets #
############################################################

foreach (method = methods) %do% {
  gridsearch <- foreach (focgene = genes, .combine = bind_rows) %:%
    foreach (featset = featsets, .combine = bind_rows,
             .packages = c("caret","magrittr","pROC","dplyr")) %do% {
               
               # Load in ML model
               model_file <- paste0(results_path, "\\", method, "_", featset, "_pt_", focgene, ".rds")
               if (!file.exists(model_file)) {
                 message("Note missing model: ", model_file)
                 return(NULL)
               }
               model <- readRDS(model_file)
               # Grid search parameter optimisation during cross-validation
               gridsearch <- model$results %>%
                 mutate(method = method,
                        focgene = focgene,
                        featset = featset) %>%
                 relocate(method, focgene, featset)
               return(gridsearch)
             }
  gridsearch %>%
    write.table(file=paste0("S3\\analysis\\gridsearch_", results_date, "_", method, ".csv"),
                sep=',', row.names=F, col.names=T)
}

####################################################################################################
# Load in models, save training set cross-validation predictions and performance #
####################################################################################################

foreach (focgene = genes, .combine = bind_rows) %do% {
  result_all <- foreach (method = methods, .combine = bind_rows) %:%
    foreach (featset = featsets, .combine = bind_rows,
             .packages = c("caret","e1071","matrixStats","magrittr","pROC","janitor","dplyr","tidyr","purrr","forcats","stringr","tibble","kernlab","xgboost","ranger","glmnet")) %do% {
               
               # Load in ML model
               model_file <- paste0(results_path, "\\", method, "_", featset, "_pt_", focgene, ".rds")
               if (!file.exists(model_file)) {
                 message("Note missing model: ", model_file)
                 return(NULL)
               }
               model <- readRDS(model_file)

               # Calculate receiver operating curve
               ROC = roc(response = model$pred$obs, 
                         predictor = model$pred$hzoon,
                         direction = ">") # note that because pROC::roc assumes the FIRST factor level is the control, we want control > cases
               
               pred_raw <- data.frame(label = model$pred$obs,
                                      hzoon = model$pred$hzoon,
                                      nz = model$pred$nz) %>%
                 mutate(pred = factor(ifelse(hzoon > coords(ROC, "best", best.method="closest.topleft")$threshold, "hzoon", "nz")))
               
               # # Save raw predictions per sequence
               # write.csv(pred_raw, file=paste0("S3\\analysis\\pred_raw\\train_pred_raw_", method, "_", featset, "_", focgene, "_", results_date, ".csv"))
               
               # Calculate confusion matrix
               matrix_train <- confusionMatrix(data = pred_raw$pred, 
                                              reference = pred_raw$label, 
                                              positive = "hzoon")
               
               # Save selected performance metrics
               line <- data.frame(method = method,
                                  featset = featset,
                                  focgene = focgene, 
                                  threshold = coords(ROC, "best", best.method="closest.topleft")$threshold,
                                  matrix_train$overall %>% t(),
                                  AUC = ROC$auc %>% as.numeric(),
                                  matrix_train$byClass %>% t()) %>%
                 mutate(across(where(is.numeric), round, 3))
               
               return(line)
               
             }
  result_all %>%
    write.table(file=paste0("S3\\analysis\\train_results_", results_date, "_", focgene, ".csv"),
                sep=',', row.names=F, col.names=T)
}


####################################################################################################
# Load in models, predict on held out test subtype sequences, and save predictions and performance #
####################################################################################################

foreach (focgene = genes, .combine = bind_rows) %do% {
  result_all <- foreach (method = methods, .combine = bind_rows) %:%
    foreach (featset = featsets, .combine = bind_rows,
             .packages = c("caret","e1071","matrixStats","magrittr","pROC","janitor","dplyr","tidyr","purrr","forcats","stringr","tibble","kernlab","xgboost","ranger","glmnet")) %do% {
               
               # Load in ML model
               model_file <- paste0(results_path, "\\", method, "_", featset, "_pt_", focgene, ".rds")
               if (!file.exists(model_file)) {
                 message("Note missing model: ", model_file)
                 return(NULL)
               }
               model <- readRDS(model_file)
               
               # Load in feature files
               pattern <- paste0("allflu_.*", featset, ".*_pt_", focgene, "\\.rds$")
               feat_file <- feature_set_files[grep(pattern, feature_set_files)]
               
               feat_df <- readRDS(feat_file) %>%
                 mutate(gid = canon_id(gid)) %>%
                 distinct(gid, .keep_all = TRUE)
               
               # Select held out test fold
               test_set_df <- folds_list[[focgene]] %>% 
                 filter(fold == "test") %>%
                 left_join(labels_list[[focgene]], by = "gid") %>% 
                 left_join(feat_df, by = "gid") %>% 
                 filter(!is.na(label)) %>%
                 mutate(across(where(is.numeric), ~ tidyr::replace_na(.x, 0))) # error fixing
               
               feats <- model$trainingData %>% select(-any_of(c(".outcome", ".weights"))) %>% names
               test_set_reduced <- test_set_df[, feats, drop = FALSE]
               
               predict_prob_test <- predict(model, newdata=test_set_df, type="prob")
               
               # Calculate receiver operating curve
               ROC = roc(response = test_set_df$label, 
                         predictor = predict_prob_test$hzoon,
                         direction = ">") # note that because pROC::roc assumes the FIRST factor level is the control, we want control > cases
               
               pred_raw <- data.frame(gid = test_set_df$gid,
                                      label = test_set_df$label,
                                      predict_prob_test) %>%
                 
                 # Could instead use the threshold determined from the training here
                 
                 mutate(pred = factor(ifelse(hzoon > coords(ROC, "best", best.method="closest.topleft")$threshold, "hzoon", "nz"))) 
               
               # Save raw predictions per sequence
               write.csv(pred_raw, file=paste0("S3\\analysis\\pred_raw\\test_pred_raw_", method, "_", featset, "_", focgene, "_", results_date, ".csv"))
               
               # Calculate confusion matrix
               matrix_test <- confusionMatrix(data = pred_raw$pred, 
                                              reference = pred_raw$label, 
                                              positive = "hzoon")
               
               # Save selected performance metrics
               line <- data.frame(method = method,
                                 featset = featset,
                                 focgene = focgene, 
                                 threshold = coords(ROC, "best", best.method="closest.topleft")$threshold,
                                 matrix_test$overall %>% t(),
                                 AUC = ROC$auc %>% as.numeric(),
                                 matrix_test$byClass %>% t()) %>%
                 mutate(across(where(is.numeric), round, 3))

               return(line)
         
             }
  result_all %>%
    write.table(file=paste0("S3\\analysis\\test_results_", results_date, "_", focgene, ".csv"),
                sep=',', row.names=F, col.names=T)
}
