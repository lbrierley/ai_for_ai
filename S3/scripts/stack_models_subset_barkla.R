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
library(ranger)
library(xgboost)
library(glmnet)
library(caretEnsemble)
library(pROC)

#########################################
# Load required data and model pointers #
#########################################

# Set parallelisation
cl <- makePSOCKcluster(3)
registerDoParallel(cl)
clusterSetRNGStream(cl, 1429)

allflu_wgs_ref <- read.csv("allflu_wgs_ref.csv") %>%
  mutate(label = factor(case_when(label == "zoon" ~ "hzoon", label == "nz" ~ "nz"))) # Rearrange factor levels for better compatibility with model functions

holdout_cluster_grid <- list.files(path = "holdout_clusters/", pattern = "labels.csv") %>%
  gsub("ex_|_labels.csv", "", .) %>%
  str_split(., "_") %>% 
  do.call(rbind.data.frame, .) %>%
  magrittr::set_colnames(c("subtype", "minseqid", "C"))

# cluster_sets <- holdout_cluster_grid %>% 
#   select(minseqid, C) %>% 
#   distinct %>% 
#   mutate(cluster_set = paste0(minseqid, "_", C)) %>% 
#   pull(cluster_set)

cluster_sets <- "70_7"

holdouts <- unique(holdout_cluster_grid$subtype) %>% as.character
holdout_zoon <- c("H7N9", "H5N1", "H9N2", "H5N6", "H10N8", "H7N3", "H3N8", "H7N7", "H7N4")
holdout_nz <- c("H4N6", "H16N3", "H4N8", "H8N4")

model_files <- list.files(pattern = ".rds", recursive = TRUE, full.names = TRUE) %>%
  .[grepl(cluster_sets, .)] %>%
  .[grepl("results_1.*_02_24", .)]

# Set up result list
result <- list()

#################################################
# Apply a stacked model to each holdout subtype #
#################################################

foreach (alphaparam = c(0,0.5,1),
           .packages = c("caret","caretEnsemble","e1071","matrixStats","magrittr","pROC","janitor","dplyr","tidyr","purrr","forcats","stringr","tibble","kernlab","xgboost","ranger","glmnet")) %dopar% {
  
  # Read ALL models, each holding out the given subtype (requires large workspace)
  model_list <- purrr::map(model_files %>% as.list(), 
                           function (x) readRDS(x) %>% .[[which(holdouts == "H7N9")]])
  
  names(model_list) <- model_files %>% gsub(".*/|.rds|_list|_pt", "", .)
  
  
  
  
  
  
  # Retain only models with positive discriminatory power (AUC > 0.5)
  
  all_res <- read.csv("results_all_methods.csv", na.strings = "NaN") %>%
    mutate(modname = paste0(method, "_", featset, "_", focgene))
  
  model_list <- model_list[names(model_list) %in% (all_res %>% filter(AUC > 0.5) %>% pull(modname))]
  
  # Retain only best models for each protein * feature set combination (i.e., eliminating dimension of method)
  
  model_list <- model_list[names(model_list) %in% (all_res %>% group_by(featset, focgene) %>% slice_max(AUC) %>% pull(modname))]

  
  
  
  # model_list %>% as.caretList %>% resamples %>% bwplot
  
  set.seed(1146)
  plr_stack <- caretStack(model_list %>% as.caretList(),
                          method = "glmnet",      # Can change as needed
                          preProc = c("center", "scale"),
                          intercept = FALSE,
                          metric = "ROC",      # Can change to sensitivity if needed
                          trControl = trainControl(
                            method = "repeatedcv", 
                            number = 5,
                            repeats = 5,
                            savePredictions = "final",
                            classProbs = TRUE,
                            summaryFunction = twoClassSummary
                          ),
                          tuneGrid = expand.grid(
                            alpha = alphaparam,         # mixing parameter (0 = ridge, 0.5 = elastic, 1 = lasso)
                            lambda = c(0.0001, 0.001, 0.01, 0.1, 0.5))     # regularisation parameter
  )
  
  
 png(file=paste0("stacks/stack_test_", alphaparam, ".png"), width=800, height=480)
 plot(plr_stack$ens_model$finalModel, xvar="lambda", label=TRUE)
 dev.off()

  saveRDS(plr_stack, file=paste0("stacks/stack_test_", alphaparam, ".rds"))
  
}

stopCluster(cl)