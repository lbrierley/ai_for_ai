library(dplyr)
library(purrr)
library(magrittr)
library(stringr)


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

read_RDS <- function(x, focgene) {
  
  readRDS(x) %>% 
    select(-any_of(c("segment", "cds_id", "enc", "GC_content"))) %>%
    rename_with(~paste(., focgene, sep = "_"), -c(gid))
  
}

read_featset <- function(focgene) {
  
  feature_sets_files <- list.files(
    path = "S3/data/full/mlready", 
    pattern = focgene, 
    full.names = TRUE)
  
  feature_sets <- map(feature_sets_files, read_RDS, focgene = focgene) %>%
    reduce(left_join, by = "gid")

}


# load data -------------------------------------------------------------------


holdout_cluster_grid <- list.files(
  path = "S3/data/full/holdout_clusters", 
  pattern = "labels.csv") %>%
  gsub("ex_|_labels.csv", "", .) %>%
  str_split(., "_") %>% 
  do.call(rbind.data.frame, .) %>%
  set_colnames(c("subtype", "minseqid", "C"))

subtypes <- unique(holdout_cluster_grid$subtype)

labels <- read.csv(paste0(root, subtypes[1], "_", cluster_set, "_labels.csv")) %>% 
  select(cluster_rep, label) %>%
  # Rearrange factor levels for better compatibility with model functions
  mutate(label = factor(case_when(label == "zoon" ~ "hzoon", label == "nz" ~ "nz")) 
  )

allflu_wgs_ref <- read.csv("S3/data/full/allflu_wgs_ref.csv") %>%
  # Rearrange factor levels for better compatibility with model functions
  mutate(label = factor(case_when(label == "zoon" ~ "hzoon", label == "nz" ~ "nz"))) 


# processing ------------------------------------------------------------------


# read and concatenate features for one gene only 
# all_featsets <- read_featset(focgene = genes[1])

# read and concatenate features for more genes 
all_featsets_ls <- map(genes[c(1,3)], read_featset)
all_featsets <- reduce(all_featsets_ls, left_join, by = "gid")

# join to labels 
all_featsets_labels <- left_join(labels, all_featsets, by = c("cluster_rep" = "gid"))

# test set
test_set <- left_join(
  allflu_wgs_ref %>% 
    filter(subtype == subtypes[1]) %>%
    # Only consider zoonotic sequences for zoonotic holdouts
    # filter(subtype %in% holdout_zoon & label == "hzoon"|subtype %in% holdout_nz) %>%  
    select(gid, subtype, label, src),
  all_featsets,
  by = c("gid"))


# save ------------------------------------------------------------------------


write.csv(
  all_featsets_labels, 
  file = sprintf("train_eval_set_%s_%s.csv", genes[1], genes[3]), 
  row.names = FALSE) 

write.csv(
  test_set,
  file = sprintf("test_set_%s_%s.csv", genes[1], genes[3]), 
  row.names = FALSE)


# option 1 provide all sequences and let them do the split (and potentially clustering?)
# 19531 records x 13303 features 
# 2.5G csv file 

# option 2 provide only the cluster representative sequences and associated test set 
# 2123 records x 13303 features
# 250K file

# either ways I am only providing features for two genes only (too many features otherwise)
