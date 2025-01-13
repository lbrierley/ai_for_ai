#######################
# Load packages, data #
#######################

holdout_cluster_grid <- list.files(path = "E:\\Working\\ai_for_ai\\S3\\data\\full\\holdout_clusters\\", pattern = "labels.csv") %>%
  gsub("ex_|_labels.csv", "", .) %>%
  str_split(., "_") %>% 
  do.call(rbind.data.frame, .) %>%
  set_colnames(c("subtype", "minseqid", "C"))

cluster_set <- "70_7"

fold_indices_list <- list()

model_list <- list.files(path = "S3/analysis", pattern = ".rds", recursive = TRUE, full.names = TRUE) %>%
  .[grepl(cluster_sets, .)] %>%
  .[grepl("results_1.*_02_24", .)] %>%
  .[1] %>%
  readRDS

for(i in 1:length(unique(holdout_cluster_grid$subtype))){
  
  ####################################################
  # Train models on each feature set on each protein #
  ####################################################
  
  # labels <- read.csv(paste0("S3\\data\\full\\holdout_clusters\\ex_", "H10N8", "_", cluster_set, "_labels.csv")) %>% 
  #   select(cluster_rep, label) %>%
  #   mutate(label = factor(label, levels = c("zoon", "nz"))
  #   )
  # 
  # # Load feature sets for training data clusters
  # train <- left_join(labels,
  #                    readRDS(paste0("S3\\data\\full\\mlready\\", "allflu_prot_2mer_pt_HA.rds")) %>% 
  #                      select(-any_of(c("segment", "cds_id", "enc", "GC_content"))),
  #                    by = c("cluster_rep" = "gid"))
  # 
  # # Specify variables used
  # preds <- train %>% select(-label, -cluster_rep) %>% remove_constant %>% names
  # 
  # # Create folds for 5-fold cross-validation
  # set.seed(1657)
  # fold_indices_list[[subtype]] <- createMultiFolds(train$label, k = 5, times = 1)
  
  fold_indices_list[[i]] <- model_list[[i]]$control$index
  
}

fold_indices_list %>% saveRDS("S3/data/fold_indices_list.rds")
