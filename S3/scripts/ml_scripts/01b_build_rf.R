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
# Global definitions
####################################################################################

##inputs
cluster_dir <- "S3/data/segmentwise/clust_processed/"
mlready_dir <- "S3/data/segmentwise/mlready/"
folds_list <- readRDS("S3/data/segmentwise/clust_dist/folds_list.rds")

#run date for file name
run_date <- format(Sys.Date(), "%d_%m_%y")

##define output directory
out_dir <- paste0("S3/analysis/results_",run_date)

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

# Set order of genes relative to segments
focgene <- c("PB2", "PB1", "PA", "HA", "NP", "NA", "M1", "NS1")


# helper function
canon_id <- function(x) {
  as.character(x) %>%
    trimws() %>%
    str_replace_all("\\s+", " ")
}

# Specify labels for all segments up front
label_files <- list.files(cluster_dir, pattern = "_labels\\.csv$", full.names = TRUE)
labels_list <- lapply(seq_along(label_files), function(x) {
  read.csv(label_files[x]) %>%
    select(cluster_rep, label, mix) %>%
    mutate(
      cluster_rep = canon_id(cluster_rep),
      label = factor(case_when(
        label == "zoon" ~ "hzoon",
        label == "nz" ~ "nz",
        TRUE ~ as.character(label)
      ))
    ) %>%
    distinct(cluster_rep, .keep_all = TRUE) %>%
    rename(gid = cluster_rep)
})

# Specify features
feature_set_files <- list.files(mlready_dir, pattern = "\\.rds$", full.names = TRUE)
featsets <- unique(gsub(".*allflu_(.*)_pt.*", "\\1", feature_set_files))

# Specify fold ids for all segments up front
folds_list <- lapply(seq_along(folds_list), function(x) {
  fold_df <- folds_list[[x]]
  fold_df$gid <- canon_id(fold_df$gid)
  return(fold_df)
})

# Set parallelisation - to disable, comment out, and change %dopar% in line 107 to %do%
cl <- makePSOCKcluster(detectCores() - 1)
registerDoParallel(cl)
clusterSetRNGStream(cl, 1429)

#######################
# Define ML procedure #
#######################

rf_fun <- function(seg, featset){
  
  ####################################################
  # Features
  ####################################################
  
  # Extract feature files
  pattern <- paste0("allflu_.*", featset, ".*_pt_", focgene[seg], "\\.rds$")
  feat_file <- feature_set_files[grep(pattern, feature_set_files)]
  
  if (is.na(feat_file)) {
    cat("No feature file found\n")
    return(NULL)
  }
  
  cat("Using feature file:\n", feat_file, "\n")
  
  feat_df <- readRDS(feat_file) %>%
    mutate(gid = canon_id(gid)) %>%
    distinct(gid, .keep_all = TRUE)
  
  # Build training set
  train <- labels_list[[seg]] %>%
    left_join(feat_df, by = "gid") %>%
    filter(!is.na(label)) %>%
    mutate(across(where(is.numeric), ~ tidyr::replace_na(.x, 0))) # error fixing
  
  preds_df <- train %>% select(!any_of(c("label", "segment", "string", "gid", "mix")))
  
  # Remove constant predictors
  nzv <- nearZeroVar(preds_df)
  if (length(nzv) > 0) preds_df <- preds_df[, -nzv, drop = FALSE]
  
  preds <- names(preds_df)
  
  if (length(preds) < 2) {
    cat("too few predictors\n")
    return(NULL)
  }
  
  ####################################################
  # Folds
  ####################################################
  
  fold_df <- folds_list[[seg]]
  fold_vec <- fold_df$fold[match(train$gid, fold_df$gid)]
  
  # # Checking progress of run
  # cat("Fold coverage:\n")
  # print(table(is.na(fold_vec)))
  
  keep <- !is.na(fold_vec) & fold_vec != "test"
  train <- train[keep, ]
  preds_df <- preds_df[keep, , drop = FALSE]
  fold_vec <- fold_vec[keep]
  
  cat("Class balance (after folds):\n")
  print(table(train$label))
  
  fold_levels <- 1:9 # 9-fold, excluding one fold for validation
  fold_indices <- lapply(fold_levels, function(f) {
    which(fold_vec != f)
  })
  
  names(fold_indices) <- paste0("Fold", seq_along(fold_indices))
  
  ####################################################
  # Weights
  ####################################################
  
  label_table <- table(train$label)
  w_nz <- 1 / as.numeric(label_table["nz"]) * 0.5
  w_hzoon <- 1 / as.numeric(label_table["hzoon"]) * 0.5
  weights_vect <- ifelse(train$label == "nz", w_nz, w_hzoon)
  # sum(weights_vect)   # should sum to 1
  
  ####################################################
  # Train model
  ####################################################
  
  model <- caret::train(
    x = preds_df,
    y = train$label,
    method = "ranger",
    preProc = c("center", "scale"),
    metric = "Kappa",
    weights = weights_vect,
    importance = "impurity",
    trControl = trainControl(
      method = "repeatedcv",
      index = fold_indices,
      number = length(fold_indices),
      repeats = 1,
      classProbs = TRUE,
      savePredictions = "final" # only save predictions for optimal tuning
    ),
    tuneGrid = expand.grid(
      .splitrule = "gini",
      .min.node.size = seq(from = 5, to = 45, length = 3),
      .mtry = round(sqrt(length(preds)))
    )
  )
  
  return(model)
}

######################################################
# Run models
######################################################

#runs models through each gene, each cluster set and each feature 
foreach (seg = 1:8) %:% 
  foreach (featset = featsets,
           .packages = c("caret","dplyr","ranger","tidyr","stringr")) %dopar% {
             
             model <- rf_fun(seg, featset)
             
             if (!is.null(model)) {
               # save outputs
               saveRDS(
                 model,
                 file = file.path(out_dir, paste0("rf_", featset, "_pt_", focgene[seg], ".rds"))
               )
             }
             
             NULL
           }

stopCluster(cl)