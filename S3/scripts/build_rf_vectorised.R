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

# Load in data
cluster_cds_wgs_df <- get(load(file = "cluster_cds_wgs_df_16_5_23.RData"))

# Set options for global analyses
outcome_name <- "label"
use_stop_codons <- TRUE

##################################################################################
# Prepare data frame for modelling and define set of variables used in the model #
##################################################################################

# Set options for local analyses
set.seed(1657)

# Prepare data frame for modelling
# If outcome not in the sequence-level dataset, merge it in from the species-level dataset
 model_df <- cluster_cds_wgs_df %>%
    mutate(outcome = factor(!!sym(outcome_name))) %>%
    filter(!is.na(outcome))

# Include nucleotide, dincucleotide and codon usage composition bias measurements, with option to exclude stop codon RSCU

## BUILD IN OPTION FOR CODON PAIRS TOO

if (use_stop_codons == TRUE){
  model_df %<>% select(matches("^[A|C|G|T]_Bias$|^[A|C|G|T][A|C|G|T]_p[1|2|3]_Bias$|^[A|C|G|T][A|C|G|T][A|C|G|T]_Bias$"), outcome, gid, subtype, date, label, src)
} else {
  model_df %<>% select(matches("^[A|C|G|T]_Bias$|^[A|C|G|T][A|C|G|T]_p[1|2|3]_Bias$|^[A|C|G|T][A|C|G|T][A|C|G|T]_Bias$"), outcome, gid, subtype, date, label, src) %>% select(-c(TAG_Bias, TAA_Bias, TGA_Bias))
}

# Specify variables used
preds <- model_df %>% select(-outcome, -gid, -subtype, -date, -label, -src) %>% remove_constant %>% names

# Create outer folds for hold-one-out validation
outer_fold_dfs <- lapply(unique(model_df$gid), function(x)
  model_df %>% filter(gid != x)
)

# Create inner folds for 1 x 10-fold cross-validation
inner_fold_indices <- lapply(outer_fold_dfs, function(x)
  createMultiFolds(x$outcome, k = 10, times = 1)
)

######################
# Run random forests #
######################

# Train and validate RF (tuning mtry, min.node.size parameters) through 10-fold cross-validation using ranger

timer_start <- Sys.time()

# Store result as list of n ensemble models where n = number of clusters, holding one out each time
rf_list <- Map(function(outer, inner) 
  
  train(x = outer %>% select(preds),
        y = outer %>% pull(outcome),
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
          .min.node.size = seq(from = 5, to = 20, length = 3),
          .mtry = seq(from = 5, to = 20, length = 3))
  ),
  outer = outer_fold_dfs,
  inner = inner_fold_indices
)

timer_end <- Sys.time()

#####################
# Store all outputs #
#####################

save(model_df, use_stop_codons,
     timer_start, timer_end, rf_list,
     file=paste0("listresults_ml_vector_cds_", format(Sys.time(), "%d_%m_%y"), ".RData"))

varimp_order_cds <- lapply(rf_list, function(x)
  varImp(x)$importance %>% rownames_to_column("name") %>% mutate(Overall = Overall/100) %>% rename(relGini = Overall)
) %>% 
  bind_rows() %>%
  group_by(name) %>% 
  summarise(mean = mean(relGini)) %>% 
  arrange(-mean) %>%
  pull(name)

saveRDS(varimp_order_cds, "varimp_order_cds.rds")