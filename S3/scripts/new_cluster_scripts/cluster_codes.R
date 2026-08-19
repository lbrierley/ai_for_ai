library(dplyr)
library(tidyr)
library(magrittr)
library(purrr)
library(tidyverse)


# Read reference set of clusters
cluster_mems_list <- NULL
cluster_reps_list <- NULL

for(i in 1:8){
  cluster_mems_list[[i]] <- read.csv(list.files("S3\\data\\segmentwise\\clust_processed", full.names = TRUE, pattern = "members")[i])
  cluster_reps_list[[i]] <- read.csv(list.files("S3\\data\\segmentwise\\clust_processed", full.names = TRUE, pattern = "labels")[i])
}

folds_list <- readRDS("S3\\data\\segmentwise\\clust_dist\\folds_list.rds")
pairwise_JSD <- readRDS("S3\\data\\segmentwise\\clust_dist\\pairwise_JSD_.rds")

####################################################################
# Look at representation of subtypes among cluster representatives #
####################################################################

# Take specific clusters and match them to GID to find subtype breakdown (regardless of stripped-out test set for n)
## Filter to only cluster representatives and calculate unique representation

cluster_meta_list <- map2(cluster_reps_list, folds_list, .f = left_join, by = c("cluster_rep" = "gid"))

# Check n subtypes for each label in each segment cluster representatives
lapply(cluster_meta_list, function(x) x %>%
         group_by(label) %>% 
         summarise(n_distinct(subtype)))

# Check n subtypes for each fold in each segment cluster repres entatives (bearing in mind folds are very different sizes)
lapply(cluster_meta_list, function(x) x %>%
         group_by(fold) %>% 
         summarise(n_distinct(subtype)))

# How many wgs would be available for the stack if you took any wgs with at least one segment selected as a cluster-representative as the data scope?
bind_rows(cluster_meta_list) %>% group_by(cluster_rep) %>% tally %>% pull(n) %>% table
bind_rows(cluster_meta_list) %>% group_by(cluster_rep) %>% tally %>% pull(n) %>% table %>% barplot
bind_rows(cluster_meta_list) %>% select(cluster_rep) %>% distinct %>% nrow
bind_rows(cluster_meta_list) %>% select(cluster_rep, label) %>% distinct %>% with(., table(label))


#############################################################################
# Identify available whole genomes (all 8 segments) and their cluster codes #
#############################################################################

# Select whole genome sequences suitable for stack model - all 8 segments and not a duplicate

allflu_cds_ref <- read.csv("S3\\data\\segmentwise\\Raw_Sequences\\allflu_cds_df.csv", na.strings = "") %>%
  select(-fastahead, -string) %>%
  filter(n >= 8 & wgs_dup == 0) 

# allflu_cds_ref %>% 
#   group_by(gid) %>%
#   select(gid, n) %>%
#   distinct %>%
#   with(., table(n))

allflu_wgs_ref <- allflu_cds_ref %>%
  group_by(gid) %>%
  filter(all(c("PB2","PB1","PA","HA","NP","NA","M1","NS1") %in% gene)) %>% # remove gids that are missing any of the genes
  ungroup %>%
  group_by(gid, gene) %>%
  filter(!(seg_dup == 1 & any(seg_dup == 0))) %>% # remove cases with duplicate seg where another seg exists for the same gene
  slice(1)   # otherwise just take the first entry per gene 

allflu_wgs_ref %>% write.csv("S3\\data\\segmentwise\\Raw_Sequences\\allflu_wgs_df.csv", row.names = FALSE) 

# Assign numeric codes to each cluster based on frequency, per segment
cluster_mems_list <- lapply(seq_along(cluster_mems_list), function(x) 
  cluster_mems_list[[x]] %>%
    group_by(cluster_rep) %>%
    add_count %>%
    ungroup %>%
    arrange(-n) %>%
    mutate(segment = x,
           cluster_id = as.numeric(
             factor(cluster_rep, 
                    levels = unique(cluster_rep)))) %>%
    select(-n, -mix)
)

cluster_mems_list %>% saveRDS("S3\\data\\segmentwise\\clust_dist\\cluster_codes.rds")

# Match each wgs's segment sequences to the respective clusters
clust_codes <- allflu_wgs_ref %>%
  left_join(bind_rows(cluster_mems_list), by = c("gid", "segment")) %>%
  group_by(gid, gene) %>%
  slice(1) %>%  # to fix issues of few gids being represented multiple times in clustering
  ungroup %>%
  select(segment, cluster_id, gid) %>%
  pivot_wider(names_from = segment, values_from = cluster_id) %>%
  mutate(code = paste(`1`, `2`, `3`, `4`, `5`, `6`, `7`, `8`, sep = "-")) %>%
  select(gid, `1`, `2`, `3`, `4`, `5`, `6`, `7`, `8`, code)

clust_codes %>% write.csv("S3\\data\\segmentwise\\cds_wgs_clust_codes.csv", row.names = FALSE)

##########################################################################
# Calculate weighted sum of JSD to capture dissimilarity between all WGS #
##########################################################################

# Set weighting per segment
JSD_wgs_seg_weights <- rep(1, 8)

# For easier processing give each segment gid an index
gid_index <- lapply(cluster_mems_list, function(df) {
  distinct_df <- df %>% select(cluster_rep, cluster_id) %>% distinct()
  setNames(distinct_df$cluster_id, distinct_df$cluster_rep)
})

# Setup empty matrix of all pairwise comparisons of wgs
JSD_wgs_matrix <- matrix(0, nrow = nrow(clust_codes), ncol = nrow(clust_codes),
                         dimnames = list(clust_codes$gid, clust_codes$gid))

# Calculate weighted sum of all segmentwise JSD distances per wgs pair
for (seg in 1:8) {
  
  col_name <- paste0("X", seg)
  seg_gids <- gid_index[[seg]][(clust_codes[[col_name]])]
  seg_matrix <- pairwise_JSD[[seg]][seg_gids, seg_gids]

  JSD_wgs_matrix <- JSD_wgs_matrix + (seg_matrix * JSD_wgs_seg_weights[seg]) # Iteratively sum contributions from each segment
}

JSD_wgs_matrix %>% saveRDS("S3\\data\\segmentwise\\clust_dist\\pairwise_JSD_wgs.rds")

