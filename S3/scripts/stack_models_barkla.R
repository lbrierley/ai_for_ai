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
cl <- makePSOCKcluster(detectCores() - 1)
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

irat_df <- read.csv("cdc_irat.csv") %>% filter(incomplete != "Y")

# Set up result list
result <- list()

#################################################
# Apply a stacked model to each holdout subtype #
#################################################

result <- foreach (subtypepicked = holdouts,
                   .packages = c("caret","caretEnsemble","e1071","matrixStats","magrittr","pROC","janitor","dplyr","tidyr","purrr","forcats","stringr","tibble","kernlab","xgboost","ranger","glmnet")) %dopar% {
                     
                     # Read ALL models, each holding out the given subtype (requires large workspace)
                     model_list <- purrr::map(model_files %>% as.list(), 
                                              function (x) readRDS(x) %>% .[[which(holdouts == subtypepicked)]])
                     
                     names(model_list) <- model_files %>% gsub(".*/|.rds|_list|_pt", "", .)
                     
                     # Retain only models with positive discriminatory power (AUC > 0.5)
                     
                     all_res <- read.csv("results_all_methods.csv", na.strings = "NaN") %>%
                       mutate(modname = paste0(method, "_", featset, "_", focgene))
                     
                     model_list <- model_list[names(model_list) %in% (all_res %>% filter(AUC > 0.5) %>% pull(modname))]
                     
                     # Retain only best models for each protein * feature set combination (i.e., eliminating dimension of method)
                     
                     model_list <- model_list[names(model_list) %in% (all_res %>% group_by(featset, focgene) %>% slice_max(AUC) %>% pull(modname))]
                     
                     # Fit initial model and extract coefficients
                     set.seed(1146)
                     temp_stack <- caretStack(model_list %>% as.caretList(),
                                              method = "glmnet",      # Can change as needed
                                              preProc = c("center", "scale"),
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
                                                alpha = c(1),         # mixing parameter (0 = ridge, 0.5 = elastic, 1 = lasso)
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
                     
                     # Refit models keeping only non-zero parameters
                     model_list <- model_list[names(model_list) %in% keep_coefs]
                     
                     set.seed(1146)
                     plr_stack <- caretStack(model_list %>% as.caretList(),
                                             method = "glmnet",      # Can change as needed
                                             preProc = c("center", "scale"),
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
                                               alpha = c(1),         # mixing parameter (0 = ridge, 0.5 = elastic, 1 = lasso)
                                               lambda = c(0.0001, 0.001, 0.01, 0.1, 0.5))     # regularisation parameter
                     )
                     
                     # coef(plr_stack$ens_model$finalModel)
                     # plot(plr_stack$ens_model$finalModel, xvar="lambda", label=TRUE)
                     
                     # Read in EVERY feature set of EVERY protein for test set subtype (53184 features)
                     allfeats <- purrr::map(list.files(path = "mlready", full.names = TRUE), 
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
                     
                     # Read in EVERY feature set of EVERY protein for CDC IRAT (53184 features)
                     
                     iratfeats <- purrr::map(list.files(path = "mlready_irat", full.names = TRUE),
                                             function (x) 
                                               readRDS(x) %>%
                                               select(-any_of(c("segment", "cds_id", "enc", "GC_content"))) %>%
                                               rename_with(~paste(., gsub(".*_pt_|.rds", "", x), sep = "_"), -c(gid)) %>%
                                               arrange(gid) %>%
                                               filter(gid %in% irat_df$gid)      # Only use IRAT where full genome sequences (possibly allow partial sequences later, but would need to impute or similar before modelling)
                     ) %>% 
                       purrr::list_cbind(name_repair = "unique_quiet") %>%
                       select(-contains("gid..")) %>%
                       bind_cols(irat_df %>% 
                                   arrange(gid),
                                 .)
                     
                     element <- list(main = data.frame(hzoon = predict(plr_stack, newdata=allfeats, type = "prob"), 
                                                       label = allfeats$label,
                                                       subtype = subtypepicked),
                                     irat = data.frame(hzoon = predict(plr_stack, newdata=iratfeats, type = "prob"),
                                                       rawpred = predict(plr_stack, newdata=iratfeats),
                                                       gid = irat_df$gid),
                                     coef = coef(plr_stack$ens_model$finalModel, 
                                                 s = plr_stack$ens_model$finalModel$lambdaOpt) %>% 
                                       as.matrix %>% 
                                       as.data.frame %>%
                                       tibble::rownames_to_column(var = "param") %>%
                                       bind_rows(data.frame(s1 = plr_stack$ens_model$finalModel$lambdaOpt, param = "lambda")) %>%
                                       mutate(subtype = subtypepicked))
                     
                     saveRDS(plr_stack, file=paste0("stacks/stack_", subtypepicked, ".rds"))
                     rm(model_list, plr_stack, allfeats, iratfeats)
                     gc()
                     
                     return(element)
                   }

stopCluster(cl)

result_all <- result %>% purrr::transpose() %>% .[["main"]] %>% bind_rows
result_irat <- result %>% purrr::transpose() %>% .[["irat"]] %>% bind_rows %>% group_by(gid) %>%    
  summarise_at(vars("hzoon"), list(med = median, upper = ~quantile(., probs = 0.25), lower = ~quantile(., probs = 0.75)))

ROC = roc(response = result_all$label,
          predictor = result_all$hzoon,
          direction = ">")

result_all %<>% mutate(pred = factor(ifelse(hzoon > coords(ROC, "best", best.method="closest.topleft")$threshold, "hzoon", "nz")))
result_irat %<>% mutate(pred = factor(ifelse(med > coords(ROC, "best", best.method="closest.topleft")$threshold, "hzoon", "nz")))

write.csv(result_all, "stack_subtypeacc_raw.csv")

matrix_test <- confusionMatrix(data = result_all$pred, 
                               reference = result_all$label, 
                               positive = "hzoon")

line <- bind_cols(threshold = coords(ROC, "best", best.method="closest.topleft")$threshold,
                  matrix_test$overall %>% t(),
                  AUC = ROC$auc %>% as.numeric(),
                  matrix_test$byClass %>% t()) %>%
  mutate(across(where(is.numeric), round, 3))

result_coefs <- result %>% purrr::transpose() %>% .[["coef"]] %>% bind_rows

write.csv(line, "stack_results.csv")
write.csv(result_irat, "stack_irat.csv")
write.csv(result_coefs, "stack_coef.csv")

saveRDS(result %>% purrr::transpose() %>% .[["irat"]], "result_irat_raw.rds")