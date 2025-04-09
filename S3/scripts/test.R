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

###################################################
# Options and global definitions used in all runs #
###################################################

dir.create("subtyperaw", showWarnings = FALSE, recursive = TRUE)

# Set parallelisation
cores <- 41
cat("cores = ", cores, "\n")

cl <- makePSOCKcluster(cores)
registerDoParallel(cl)
clusterSetRNGStream(cl, 1429)

results_date <- "2025_03_10"
method <- "svm"

allflu_wgs_ref <- read.csv("S3/data/full/allflu_wgs_ref.csv") %>%
  mutate(label = factor(case_when(label == "zoon" ~ "hzoon", label == "nz" ~ "nz"))) # Rearrange factor levels for better compatibility with model functions

holdout_cluster_grid <- list.files(path = "S3/data/full/holdout_clusters/", pattern = "labels.csv") %>%
  gsub("ex_|_labels.csv", "", .) %>%
  str_split(., "_") %>% 
  do.call(rbind.data.frame, .) %>%
  magrittr::set_colnames(c("subtype", "minseqid", "C"))

# cluster_sets <- holdout_cluster_grid %>% 
# select(minseqid, C) %>% 
# distinct %>% 
# mutate(cluster_set = paste0(minseqid, "_", C)) %>% 
# pull(cluster_set)

cluster_sets <- "70_7"

holdouts <- unique(holdout_cluster_grid$subtype) %>% as.character
holdout_zoon <- c("H7N9", "H5N1", "H9N2", "H5N6", "H10N8", "H7N3", "H3N8", "H7N7", "H7N4")
holdout_nz <- c("H4N6", "H16N3", "H4N8", "H8N4")

#########################f
# Load ML model objects #
#########################

gridsearch <- foreach (cluster_set = cluster_sets) %:%
  foreach (focgene = c("HA", "M1", "NA", "NP", "NS1", "PA", "PB1", "PB2")) %:%
  foreach (featset = list.files(path = "S3/data/full/mlready", pattern = focgene) %>% gsub("allflu_|_pt.*.rds", "", .),
           .packages = c("caret","magrittr","pROC","dplyr")) %dopar% {

             # Load in ML model
             model_list <- readRDS(paste0("results_", results_date, "/", cluster_set, "/", method, "_list_", featset, "_pt_", focgene, ".rds"))


             # Grid search parameter optimisation on validation sets
             gridsearch <- lapply(model_list, function(x)  x$results) %>%
               bind_rows() %>%
               mutate(cluster_set = cluster_set,
                      featset = featset,
                      focgene = focgene)

             gridsearch %<>%
               mutate(subtype = rep(holdouts, each=nrow(gridsearch)/length(holdouts))) %>%
               relocate(cluster_set, featset, focgene, subtype)

             return(gridsearch)
           }

gridsearch %>%
  unlist(recursive=FALSE) %>%
  unlist(recursive=FALSE) %>%
  bind_rows() %>%
  write.table(file=paste0("gridsearch_", results_date, ".csv"),
              sep=',', row.names=F, col.names=T)

result_all <- foreach (cluster_set = cluster_sets) %:% 
  foreach (focgene = c("HA", "M1", "NA", "NP", "NS1", "PA", "PB1", "PB2")) %:% 
  foreach (featset = list.files(path = "S3/data/full/mlready", pattern = focgene) %>% gsub("allflu_|_pt.*.rds", "", .),
           .packages = c("caret","e1071","matrixStats","magrittr","pROC","janitor","dplyr","tidyr","purrr","forcats","stringr","tibble","kernlab","xgboost","ranger","glmnet")) %dopar% {
             
             # Load in ML model
             model_list <- readRDS(paste0("results_", results_date, "/", cluster_set, "/", method, "_list_", featset, "_pt_", focgene, ".rds"))
             
             # Set up result list
             result_all <- list()
             
             # Set up corresponding holdout test set list
             test_set_list <- lapply(holdouts, function(x)   
               left_join(allflu_wgs_ref %>% 
                           filter(subtype == x) %>%
                           filter(subtype %in% holdout_zoon & label == "hzoon"|subtype %in% holdout_nz) %>%  # Only consider zoonotic sequences for zoonotic holdouts
                           select(gid, subtype, label, src),
                         readRDS(paste0("S3/data/full/mlready/allflu_", featset, "_pt_", focgene, ".rds")) %>% 
                           select(-any_of(c("segment", "cds_id", "enc", "GC_content"))) %>%
                           rename_with(~paste(., focgene, sep = "_"), -c(gid)),
                         by = c("gid")))
             
             predict_prob_test <- Map(function(model, newdata)
               
               if (nrow(newdata)>0){
                 predict(model, newdata=newdata, type="prob")
               },
               model = model_list,
               newdata = test_set_list
             ) %>% bind_rows
             
             ROC = roc(response = test_set_list %>% bind_rows() %>% pull(label), 
                       predictor = predict_prob_test %>% pull(hzoon),
                       direction = ">")
             
             subtypeacc <- data.frame(cluster_set = cluster_set,
                                      featset = featset,
                                      focgene = focgene,
                                      test_set_list %>%
                                        bind_rows %>%
                                        select(subtype, label),
                                      predict_prob_test) %>%
               mutate(pred = factor(ifelse(hzoon > coords(ROC, "best", best.method="closest.topleft")$threshold, "hzoon", "nz")))
             
             write.csv(subtypeacc, file=paste0("subtyperaw/subtypeacc_raw_", cluster_set, "_", featset, "_", focgene, "_", results_date, ".csv"))
             
             matrix_test <- confusionMatrix(data = subtypeacc$pred, 
                                            reference = subtypeacc$label, 
                                            positive = "hzoon")
             
             line <- bind_cols(cluster_set = cluster_set,
                               featset = featset,
                               focgene = focgene, 
                               threshold = coords(ROC, "best", best.method="closest.topleft")$threshold,
                               matrix_test$overall %>% t(),
                               AUC = ROC$auc %>% as.numeric(),
                               matrix_test$byClass %>% t()) %>%
               mutate(across(where(is.numeric), round, 3))
             
             subtypeacc_reduced <- subtypeacc %>%
               group_by(cluster_set, featset, focgene, subtype) %>%
               summarise(accuracy = sum(label == pred) / n())
             
             result_all[["line"]] <- line
             result_all[["subtypeacc_reduced"]] <- subtypeacc_reduced
             return(result_all)
             
           }

result_all %>% 
  unlist(recursive=FALSE) %>% 
  unlist(recursive=FALSE) %>% 
  purrr::transpose(.) %>% 
  extract2("subtypeacc_reduced") %>% 
  bind_rows() %>%
  write.table(file=paste0("resultsbysubtype_", results_date, ".csv"), 
              sep=',', row.names=F, col.names=T)

result_all %>% 
  unlist(recursive=FALSE) %>% 
  unlist(recursive=FALSE) %>% 
  purrr::transpose(.) %>% 
  extract2("line") %>% 
  bind_rows() %>%
  write.table(file=paste0("results_", results_date, ".csv"), 
              sep=',', row.names=F, col.names=T)

stopCluster(cl)
