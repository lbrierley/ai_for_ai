#################
# Load packages #
#################

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

###############
# Set options #
###############

# Only attempt to build stack models with positive discriminatory power (AUC > 0.5)?
# This slightly cuts down computing time by pre-filtering some candidate models out.
auc_filter <- TRUE

# Set parallelisation - to disable, comment out, and change %dopar% in line 67 to %do%
cl <- makePSOCKcluster(detectCores() - 1)
registerDoParallel(cl)
clusterSetRNGStream(cl, 1429)

#######################
# Load model pointers #
#######################

allflu_wgs_ref <- read.csv("S3\\data\\full\\allflu_wgs_ref.csv") %>%
  mutate(label = factor(case_when(label == "zoon" ~ "hzoon", label == "nz" ~ "nz"))) # Rearrange factor levels for better compatibility with model functions

holdout_cluster_grid <- list.files(path = "S3\\data\\full\\holdout_clusters", pattern = "labels.csv") %>%
  gsub("ex_|_labels.csv", "", .) %>%
  str_split(., "_") %>% 
  do.call(rbind.data.frame, .) %>%
  magrittr::set_colnames(c("subtype", "minseqid", "C"))

cluster_chosen <- "70_7"

holdouts <- unique(holdout_cluster_grid$subtype) %>% as.character
holdout_zoon <- c("H7N9", "H5N1", "H9N2", "H5N6", "H10N8", "H7N3", "H3N8", "H7N7", "H7N4")
holdout_nz <- c("H4N6", "H16N3", "H4N8", "H8N4")

#########################################
# Load required data and model pointers #
#########################################

# Define list of individual model files
model_files <- list.files(path = "S3\\analysis\\", pattern = ".rds", recursive = TRUE, full.names = TRUE) %>%
  .[grepl(cluster_chosen, .)] %>%
  .[grepl("results_1.*_02_24", .)]

# Set up result list
result <- list()

#################################################
# Apply a stacked model to each holdout subtype #
#################################################

result <- foreach (subtypepicked = holdouts,
                   .packages = c("caret","caretEnsemble","e1071","matrixStats","magrittr","pROC","janitor","dplyr","tidyr","purrr","forcats","stringr","tibble","kernlab","xgboost","ranger","glmnet")) %dopar% {
                     
                     # Determine weights
                     labels <- read.csv(paste0("S3\\data\\full\\holdout_clusters\\ex_", subtypepicked, "_", cluster_chosen, "_labels.csv")) %>% 
                       select(cluster_rep, label) %>%
                       mutate(label = factor(case_when(label == "zoon" ~ "hzoon", label == "nz" ~ "nz")) # Rearrange factor levels for better compatibility with model functions
                       )
                     
                     # Read ALL models, each holding out the given subtype (requires large workspace)
                     model_list <- purrr::map(model_files %>% as.list(), 
                                              function (x) readRDS(x) %>% .[[which(holdouts == subtypepicked)]])
                     
                     names(model_list) <- model_files %>% gsub(".*/|.rds|_list|_pt", "", .)
                     
                     if(auc_filter == TRUE){
                       
                       # Load in results of all individual models
                       results_rf <- read.csv(paste0("S3\\analysis\\results_", "14_02_24", ".csv"), na.strings = "NaN") %>% filter(cluster_set == "70_7") %>% mutate(method = "rf") 
                       results_plr <- read.csv(paste0("S3\\analysis\\results_", "15_02_24", ".csv"), na.strings = "NaN") %>% filter(cluster_set == "70_7") %>% mutate(method = "glmnet")
                       results_xgb <- read.csv(paste0("S3\\analysis\\results_", "16_02_24", ".csv"), na.strings = "NaN") %>% filter(cluster_set == "70_7") %>% mutate(method = "xgb")
                       results_svmlin <- read.csv(paste0("S3\\analysis\\results_", "17_02_24", ".csv"), na.strings = "NaN") %>% filter(cluster_set == "70_7") %>% mutate(method = "svmlin")
                       results_svmrad <- read.csv(paste0("S3\\analysis\\results_", "18_02_24", ".csv"), na.strings = "NaN") %>% filter(cluster_set == "70_7") %>% mutate(method = "svm")
                       
                       all_res <- bind_rows(results_rf,
                                            results_svmlin,
                                            results_svmrad,
                                            results_xgb,
                                            results_plr) %>%
                         mutate(modname = paste0(method, "_", featset, "_", focgene))
                       
                       
                       # Filter model list to positive discriminatory models
                       model_list <- model_list[names(model_list) %in% (all_res %>% filter(AUC > 0.5) %>% pull(modname))]
                       
                     }
                     
                     # Filter to only best models for each protein * feature set combination (i.e., eliminating dimension of method)
                     model_list <- model_list[names(model_list) %in% (all_res %>% group_by(featset, focgene) %>% slice_max(AUC) %>% pull(modname))]
                     
                     # Fit initial model and extract coefficients
                     set.seed(1146)
                     temp_stack <- caretStack(model_list %>% as.caretList(),
                                              method = "glmnet",      
                                              preProc = c("center", "scale"),
                                              metric = "ROC",      
                                              weights = ifelse(labels$label == "nz",
                                                               (1/table(labels$label)[2]) * 0.5,
                                                               (1/table(labels$label)[1]) * 0.5),
                                              trControl = trainControl(
                                                method = "repeatedcv", 
                                                number = 5,
                                                repeats = 5,
                                                savePredictions = "final",
                                                classProbs = TRUE,
                                                summaryFunction = twoClassSummary
                                              ),
                                              tuneGrid = expand.grid(
                                                alpha = c(1),                                  # mixing parameter (0 = ridge, 0.5 = elastic, 1 = lasso)
                                                lambda = c(0.0001, 0.001, 0.01, 0.1, 0.5))     # regularisation parameter
                     )
                     
                     
                     keep_coefs <- coef(temp_stack$ens_model$finalModel, 
                                        s = temp_stack$ens_model$finalModel$lambdaOpt) %>% 
                       as.matrix %>% 
                       as.data.frame %>%
                       tibble::rownames_to_column(var = "param") %>% 
                       filter(param != "(Intercept)" & s1 != 0) %>%
                       pull(param)
                     
                     rm(temp_stack)
                     
                     # Refit stack models filtering to only individual models that have non-zero parameters (i.e., remove models not selected for the stack). This makes the overall stack object smaller and more manageable for permutation variable importance.
                     model_list <- model_list[names(model_list) %in% keep_coefs]
                     
                     set.seed(1146)
                     plr_stack <- caretStack(model_list %>% as.caretList(),
                                             method = "glmnet",   
                                             preProc = c("center", "scale"),
                                             metric = "ROC",   
                                             weights = ifelse(labels$label == "nz",
                                                              (1/table(labels$label)[2]) * 0.5,
                                                              (1/table(labels$label)[1]) * 0.5),
                                             trControl = trainControl(
                                               method = "repeatedcv", 
                                               number = 5,
                                               repeats = 5,
                                               savePredictions = "final",
                                               classProbs = TRUE,
                                               summaryFunction = twoClassSummary
                                             ),
                                             tuneGrid = expand.grid(
                                               alpha = c(1),                                  # mixing parameter (0 = ridge, 0.5 = elastic, 1 = lasso)
                                               lambda = c(0.0001, 0.001, 0.01, 0.1, 0.5))     # regularisation parameter
                     )
                     
                     # Read in EVERY feature set of EVERY protein for test set subtype (53184 features)
                     allfeats <- purrr::map(list.files(path = "S3\\data\\full\\mlready", full.names = TRUE), 
                                            function (x) 
                                              readRDS(x) %>%
                                              select(-any_of(c("segment", "cds_id", "enc", "GC_content"))) %>%
                                              rename_with(~paste(., gsub(".*_pt_|.rds", "", x), sep = "_"), -c(gid)) %>%
                                              right_join(allflu_wgs_ref %>% 
                                                           filter(subtype == subtypepicked) %>%
                                                           filter(subtype %in% holdout_zoon & label == "hzoon"|subtype %in% holdout_nz) %>%  # Only consider zoonotic sequences for zoonotic holdouts
                                                           select(gid),
                                                         by = c("gid")) %>%
                                              arrange(gid)
                     ) %>% 
                       purrr::list_cbind(name_repair = "unique_quiet") %>%
                       select(-contains("gid..")) %>%
                       bind_cols(allflu_wgs_ref %>% 
                                   filter(subtype == subtypepicked) %>%
                                   filter(subtype %in% holdout_zoon & label == "hzoon"|subtype %in% holdout_nz) %>%  # Only consider zoonotic sequences for zoonotic holdouts
                                   select(gid, subtype, label, src) %>% 
                                   arrange(gid),
                                 .)
                     
                     # Generate predictions on test set subtype
                     element <- list(main = data.frame(hzoon = predict(plr_stack, newdata=allfeats, type = "prob"), 
                                                       label = allfeats$label,
                                                       subtype = subtypepicked),
                                     coef = coef(plr_stack$ens_model$finalModel, 
                                                 s = plr_stack$ens_model$finalModel$lambdaOpt) %>% 
                                       as.matrix %>% 
                                       as.data.frame %>%
                                       tibble::rownames_to_column(var = "param") %>%
                                       bind_rows(data.frame(s1 = plr_stack$ens_model$finalModel$lambdaOpt, param = "lambda")) %>%
                                       mutate(subtype = subtypepicked))
                     
                     # Save the stack for this test set subtype
                     saveRDS(plr_stack, file=paste0("S3\\analysis\\stacks_weight\\stack_", subtypepicked, ".rds"))
                     rm(model_list, plr_stack, allfeats)
                     gc()
                     
                     return(element)
                   }

stopCluster(cl)

##########################################################################
# Summarise aggregated performance across held out subtypes as test sets #
##########################################################################

result_all <- result %>% purrr::transpose() %>% .[["main"]] %>% bind_rows

# Calculate receiver operating curve
ROC = roc(response = result_all$label,
          predictor = result_all$hzoon,
          direction = ">")

result_all %<>% mutate(pred = factor(ifelse(hzoon > coords(ROC, "best", best.method="closest.topleft")$threshold, "hzoon", "nz")))

# Save raw predictions per sequence
write.csv(result_all, "S3\\analysis\\stack_weight_subtypeacc_raw.csv")

# Calculate confusion matrix
matrix_test <- confusionMatrix(data = result_all$pred, 
                               reference = result_all$label, 
                               positive = "hzoon")

# Save selected performance metrics
line <- bind_cols(threshold = coords(ROC, "best", best.method="closest.topleft")$threshold,
                  matrix_test$overall %>% t(),
                  AUC = ROC$auc %>% as.numeric(),
                  matrix_test$byClass %>% t()) %>%
  mutate(across(where(is.numeric), round, 3))

write.csv(line, "S3\\analysis\\stack_weight_results.csv")

# Save coefficients on each individual model from the glmnet stack model
result_coefs <- result %>% purrr::transpose() %>% .[["coef"]] %>% bind_rows
write.csv(result_coefs, "S3\\analysis\\stack_weight_coef.csv")