library(dplyr)
library(purrr)
library(magrittr)
library(stringr)

source("S3/scripts/data_challenge_dataset_functions.R")


# define parameters -----------------------------------------------------------


holdout_zoon <- c(
  "H7N9", 
  "H5N1", 
  "H9N2", 
  "H5N6", 
  "H10N8", 
  "H7N3", 
  "H3N8", 
  "H7N7", 
  "H7N4"
)

holdout_nz <- c("H4N6", "H16N3", "H4N8", "H8N4")

root <- "S3/data/full/holdout_clusters/ex_"

cluster_set <- "70_7"

genes <- c("HA", "M1", "NA", "NP", "NS1", "PA", "PB1", "PB2")

feature_sets <- list.files(
  path = "S3/data/full/mlready", 
  pattern = genes[1]) %>% 
  gsub("allflu_|_pt.*.rds", "", .)


set.seed(16354)


# load data -------------------------------------------------------------------


holdout_cluster_grid <- list.files(
  path = "S3/data/full/holdout_clusters", 
  pattern = "labels.csv") %>%
  gsub("ex_|_labels.csv", "", .) %>%
  str_split(., "_") %>% 
  do.call(rbind.data.frame, .) %>%
  set_colnames(c("subtype", "minseqid", "C"))

subtypes <- unique(holdout_cluster_grid$subtype)

labels <- read.csv(paste0(root, subtypes[8], "_", cluster_set, "_labels.csv")) %>% 
  select(cluster_rep, label) %>%
  # Rearrange factor levels for better compatibility with model functions
  mutate(label = factor(case_when(label == "zoon" ~ "hzoon", label == "nz" ~ "nz")) 
  )

allflu_wgs_ref <- read.csv("S3/data/full/allflu_wgs_ref.csv") %>%
  # Rearrange factor levels for better compatibility with model functions
  mutate(label = factor(case_when(label == "zoon" ~ "hzoon", label == "nz" ~ "nz"))) %>%
  select(gid, subtype, label, src)


# processing ------------------------------------------------------------------


# read and concatenate n_col random features for all genes 
all_featsets_ls <- map(genes, read_featset, n_col = 10)
all_featsets <- reduce(all_featsets_ls, left_join, by = "gid")

# join features to labels (of the clustered dataset where one subtype was held out for test) 
featsets_labels_cluster <- left_join(labels, all_featsets, by = c("cluster_rep" = "gid"))

# join features to labels (of the full dataset) 
featsets_labels_full <- left_join(allflu_wgs_ref, all_featsets, by = "gid")

# the clustered data set is to be tested on a separate set
test_set <- left_join(
  allflu_wgs_ref %>%
    filter(subtype == subtypes[8]),
    # Only consider zoonotic sequences for zoonotic holdouts
    # filter(subtype %in% holdout_zoon & label == "hzoon"|subtype %in% holdout_nz)
  all_featsets,
  by = c("gid"))


# save ------------------------------------------------------------------------


write.csv(
  featsets_labels_cluster,
  file = "zoonosis_dataset_train_eval.csv",
  row.names = FALSE)

write.csv(
  test_set,
  file = "zoonosis_dataset_test.csv",
  row.names = FALSE)

write.csv(
  featsets_labels_full,
  file = "zoonosis_dataset_full.csv", 
  row.names = FALSE)
