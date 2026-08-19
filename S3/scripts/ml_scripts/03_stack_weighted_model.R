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
library(purrr)

###############
# Set options #
###############

# Only attempt to build stack models with positive discriminatory power (AUC > 0.5)?
# This slightly cuts down computing time by pre-filtering some candidate models out.
auc_filter <- TRUE

######################
# Global definitions
######################

#run date for file name
run_date <- "17_07_26"

##inputs
cluster_dir <- "S3/data/segmentwise/clust_processed/"
mlready_dir <- "S3/data/segmentwise/mlready/"
model_dir <- paste0("S3/analysis/results_",run_date)
results_dir <- "S3/analysis/"
folds_wgs <- readRDS("S3/data/segmentwise/clust_dist/folds_wgs.rds")

##define output directory
out_dir <- paste0("S3/analysis/stack_weight_",run_date)

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

###############
# Load models 
###############

# Read all constituent model results
results_list <- list.files(results_dir, pattern = paste0("test_results_", run_date), full.names = TRUE)
all_res <- map_df(results_list, read.csv, na.strings = "") %>% 
  mutate(modpath = paste0(model_dir, "/", method, "_", featset, "_pt_", focgene, ".rds"))


# Filter to only best models for each protein * feature set combination (i.e., choose only one ML method)
all_res <- all_res %>% 
  group_by(featset, focgene) %>% 
  slice_max(AUC)

if(auc_filter == TRUE){
  # Filter model list to positive discriminatory models
  all_res <- all_res %>% filter(AUC > 0.5)
}

# Read models to be placed into the stack
model_list <- purrr::map(all_res$modpath %>% as.list(), 
                         function (x) readRDS(x))

names(model_list) <- all_res$modpath %>% gsub(".*/|.rds|_list|_pt", "", .)

############
# Load data 
############

allflu_wgs_ref <- read.csv("S3\\data\\segmentwise\\Raw_Sequences\\allflu_wgs_df.csv", na.strings = "") %>%
  select(gid, subtype, label) %>%
  distinct() %>%
  filter(gid %in% folds_wgs$gid) %>%
  mutate(label = factor(case_when(label == "zoon" ~ "hzoon", label == "nz" ~ "nz"))) # Rearrange factor levels for better compatibility with model functions

# Specify features
feature_set_files <- list.files(mlready_dir, pattern = "\\.rds$", full.names = TRUE)
featsets <- unique(gsub(".*allflu_(.*)_pt.*", "\\1", feature_set_files))

# Read in every feature set of every segment for these wgs (53184 features)
allfeats <- purrr::map(feature_set_files, 
                       function (x) 
                         readRDS(x) %>%
                         filter(gid %in% allflu_wgs_ref$gid) %>%
                         select(-any_of(c("segment", "cds_id", "enc", "GC_content"))) %>%
                         rename_with(~paste(., gsub(".*_pt_|.rds", "", x), sep = "_"), -c(gid))
) %>% 
  purrr::list_cbind(name_repair = "unique_quiet") %>%
  select(-contains("gid..")) %>%
  bind_cols(allflu_wgs_ref %>% 
              filter(subtype == subtypepicked) %>%
              filter(subtype %in% holdout_zoon & label == "hzoon"|subtype %in% holdout_nz) %>%  # Only consider zoonotic sequences for zoonotic holdouts
              select(gid, subtype, label, src) %>% 
              arrange(gid),
            .)


###################################################################
# Apply a manual stacked model for the chosen external holdout set 
###################################################################

#Generate predictions from each independent base model on that holdout set.
#Bind those predictions into a new dataframe (where each column is a model's prediction).
#Train a meta-model (e.g., via standard caret::train()) using the predictions to predict the final target.


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