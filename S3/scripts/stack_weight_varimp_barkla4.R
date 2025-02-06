#######################
# Load packages, data #
#######################

library(caret)
library(magrittr)
library(dplyr)
library(tidyr)
library(stringr)
library(parallel)
library(doParallel)
library(foreach)
library(glmnet)
library(caretEnsemble)
library(pROC)

###############
# Set options #
###############

# Number of permutations
n_perms <- 1

## Set parallelisation
cl <- makePSOCKcluster(6)
registerDoParallel(cl)
clusterSetRNGStream(cl, 1429)

#######################
# Load model pointers #
#######################

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

########################################
# Load in stacked models and test data #
########################################

stack_list <- holdouts %>% 
  purrr::map(function(x)
    readRDS(paste0("stacks_weight/stack_", x, ".rds")))

feats_list <- holdouts %>% 
  purrr::map(function(y)
    purrr::map(list.files(path = "mlready", full.names = TRUE), 
               function(x) 
                 readRDS(x) %>%
                 select(-any_of(c("segment", "cds_id", "enc", "GC_content"))) %>%
                 rename_with(~paste(., gsub(".*_pt_|.rds", "", x), sep = "_"), -c(gid)) %>%
                 right_join(allflu_wgs_ref %>% 
                              filter(subtype == y) %>%
                              filter(subtype %in% holdout_zoon & label == "hzoon"|subtype %in% holdout_nz) %>%  # Only consider zoonotic sequences for zoonotic holdouts
                              select(gid),
                            by = c("gid")) %>%
                 arrange(gid)
    ) %>% 
      purrr::list_cbind(name_repair = "unique_quiet") %>%
      select(-contains("gid..")) %>%
      bind_cols(allflu_wgs_ref %>% 
                  filter(subtype == y) %>%
                  filter(subtype %in% holdout_zoon & label == "hzoon"|subtype %in% holdout_nz) %>%  # Only consider zoonotic sequences for zoonotic holdouts
                  select(gid, subtype, label, src) %>% 
                  arrange(gid),
                .)
  )

###################################
# Permutation variable importance #
###################################

set.seed(1717)

# Restrict variables to only those used in models
varnames <- feats_list %>% 
  bind_rows %>% 
  select(-gid, -label, -subtype, -src) %>% 
  names

stacked_coef <- read.csv("stack_weight_coef.csv") %>%
  select(-X) %>%
  filter(param != "(Intercept)" & param != "lambda")  %>%
  separate_wider_delim(param, delim = "_", names = c("method", "feat_a", "feat_b", "gene")) %>%
  mutate(featset = paste(feat_a, feat_b, sep = "_")) %>%
  select(featset, gene) %>%
  distinct() %>% 
  arrange(featset, gene)

# Filter out feature names sequentially based on retained models

for (focgene in  c("HA", "M1", "NA", "NP", "NS1", "PA", "PB1", "PB2")){
  
  feats_in <- stacked_coef %>% filter(gene == focgene) %>% pull(featset)
  
  if (!("nuc_2mer" %in% feats_in)){
    varnames <- varnames[!grepl(paste0("^", paste0(rep("[A|C|G|T]", 2), collapse = ""), "_", focgene), varnames)]
  }
  
  if (!("nuc_3mer" %in% feats_in)){
    varnames <- varnames[!grepl(paste0("^", paste0(rep("[A|C|G|T]", 3), collapse = ""), "_", focgene), varnames)]
  }
  
  if (!("nuc_4mer" %in% feats_in)){
    varnames <- varnames[!grepl(paste0("^", paste0(rep("[A|C|G|T]", 4), collapse = ""), "_", focgene), varnames)]
  }
  
  if (!("nuc_5mer" %in% feats_in)){
    varnames <- varnames[!grepl(paste0("^", paste0(rep("[A|C|G|T]", 5), collapse = ""), "_", focgene), varnames)]
  }
  
  if (!("nuc_6mer" %in% feats_in)){
    varnames <- varnames[!grepl(paste0("^", paste0(rep("[A|C|G|T]", 6), collapse = ""), "_", focgene), varnames)]
  }
  
  if (!("cds_compbias" %in% feats_in)){
    varnames <- varnames[!grepl(paste0("Bias_", focgene), varnames)]
  }
  
  if (!("prot_2mer" %in% feats_in)){
    varnames <- varnames[!grepl(paste0("^DPC_.*_", focgene), varnames)]
  }
  
  if (!("prot_pseeac" %in% feats_in)){
    varnames <- varnames[!grepl(paste0("^PAAC_.*_", focgene), varnames)]
  }
  
  if (!("prot_ctriad" %in% feats_in)){
    varnames <- varnames[!grepl(paste0("^CTriad_.*_", focgene), varnames)]
  }
  
  if (!("prot_ctdc" %in% feats_in)){
    varnames <- varnames[!grepl(paste0("^CTDC_.*_", focgene), varnames)]
  }
  
  if (!("prot_ctdt" %in% feats_in)){
    varnames <- varnames[!grepl(paste0("^CTDT_.*_", focgene), varnames)]
  }
  
  if (!("prot_ctdd" %in% feats_in)){
    varnames <- varnames[!grepl(paste0("^CTDD_.*_", focgene), varnames)]
  }
}

varnames <- rep(varnames, 
                each = n_perms)

# feats_list %<>% lapply(function(x)
#   x %>% select(gid, subtype, label, src, any_of(varnames)) # Drop unused features
# )

# Apply stack to make predictions

predict_prob_test <- Map(function(model, newdata)
  data.frame(hzoon = predict(model, newdata=newdata, type = "prob"), 
             label = newdata$label),
  model = stack_list,
  newdata = feats_list
) %>% bind_rows

AUC_base = roc(response = feats_list %>% bind_rows() %>% pull(label),
               predictor = predict_prob_test %>% pull(hzoon),
               direction = ">",
               quiet = TRUE)$auc %>% 
  as.numeric

varimp_perm <- foreach (varname = varnames[30001:length(varnames)],
.packages = c("caret","caretEnsemble","matrixStats","magrittr","pROC","dplyr","tidyr","purrr","stringr","tibble","kernlab","xgboost","ranger","glmnet"),
.inorder = FALSE) %dopar% {
  
  start <- Sys.time()
  
  feats_list_perm <- feats_list %>%
    bind_rows %>% # bind all subtype data frames together
    mutate_at(vars(varname), sample) %>% # permute given column
    group_split(subtype) # reform into list of subtype data frames
  
  perm_test <- Map(function(model, newdata)
    data.frame(hzoon = predict(model, newdata=newdata, type = "prob"), 
               label = newdata$label),
    model = stack_list,
    newdata = feats_list_perm
  ) %>% bind_rows
  
  AUC_perm = roc(response = feats_list_perm %>% bind_rows() %>% pull(label),
                 predictor = perm_test %>% pull(hzoon),
                 direction = ">",
                 quiet = TRUE)$auc %>% 
    as.numeric
    
  end <- Sys.time()
 
  print(end - start)
 
  return(data.frame(var = varname, AUC_loss = AUC_base - AUC_perm))
  
}


stopCluster(cl)

varimp_perm %>% bind_rows %>% write.csv("varimp_perm_weight_pt4.csv")