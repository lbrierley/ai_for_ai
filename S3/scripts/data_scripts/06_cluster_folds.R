################################################################
# Designate training and test fold membership using clustering #
################################################################
###########
# Options #
###########

# Set k for k-mers and j for j blocks
k <- 3
j <- 10

###########################
# Read in and set up data #
###########################

cluster_reps_list <- NULL
fasta_list <- NULL

for(i in 1:8){
  # Read reference set of clusters
  cluster_reps_list[[i]] <- read.csv(list.files("H:\\Working\\ai_for_ai\\S3\\data\\segmentwise\\clust_processed", full.names = TRUE, pattern = "labels")[i])
  
  fasta_list[[i]] <- readSet(file = paste0("S3\\data\\segmentwise\\allflu_cds_seg_", i, "_full.fasta"))
}

# Read in sequences and generate k-mer counts at desired k
reps_kmer_list <- Map(function(x,y) 
  x[y$cluster_rep] %>%
    oligonucleotideFrequency(width = 3, step = 1) %>% 
    as.data.frame %>% 
    mutate(across(where(is.numeric))/rowSums(across(where(is.numeric)))) %>%
    bind_cols(gid = y$cluster_rep, .),
  fasta_list, cluster_reps_list)

reps_kmer_list %>% saveRDS("S3\\data\\segmentwise\\clust_dist\\reps_kmer_list.rds")
#reps_kmer_list <- readRDS("S3\\data\\segmentwise\\clust_dist\\reps_kmer_list.rds") # For loading at intermediate stage

#######################################
# Calculate pairwise distance metrics #
#######################################

# Initialise empty list
pairwise_JSD <- vector(mode = "list", length = 8)

# Calc Jensen-Shannon distance in parallel
cl <- makePSOCKcluster(8)
registerDoParallel(cl)
clusterSetRNGStream(cl, 1545)

pairwise_JSD <- foreach(seg = 1:8,
                        .packages = c("dplyr", "philentropy")) %dopar% {
                          
                          kmer_df <- reps_kmer_list[[seg]] %>% distinct()
                          kmer_matrix <- as.matrix(select(kmer_df, -gid)) %>% 
                            set_rownames(kmer_df$gid)
                          
                          JSD_matrix <- kmer_matrix %>%
                            distance(method = "jensen-shannon") %>% 
                            set_rownames(kmer_df$gid) %>% 
                            set_colnames(kmer_df$gid)
                          
                          return(JSD_matrix)
                        }

pairwise_JSD %>% saveRDS("S3\\data\\segmentwise\\clust_dist\\pairwise_JSD_.rds")
# pairwise_JSD <- readRDS("S3\\data\\segmentwise\\clust_dist\\pairwise_JSD_.rds") # For loading at intermediate stage

########################################################
# Hierarchical clustering and testing various fold j's #
########################################################


# Hierarchically from hclust package
hclust_list <- lapply(pairwise_JSD, function(x) x %>%
                        as.dist(diag = TRUE) %>%
                        hclust(method="ward.D2")) # use Ward's minimum variance method to force more compact clustering

hclust_list %>% saveRDS("S3\\data\\segmentwise\\clust_dist\\hclust_list.rds")
#hclust_list <- readRDS("S3\\data\\segmentwise\\clust_dist\\hclust_list.rds") # For loading at intermediate stage

# Cut into different values of j folds and check distribution
folds_list <- foreach(i = 1:8,
                      .packages = c("dplyr",
                                    "tibble",
                                    "magrittr")) %dopar% {
                                      
                                      tree <- hclust_list[[i]]
                                      cut_k <- tree %>% cutree(k = c(5, 7, 10)) 
                                      
                                      cut_df <- data.frame(gid = tree$labels) %>% 
                                        left_join(cluster_reps_list[[i]] %>% select(cluster_rep, label) %>% distinct,
                                                  by = c("gid" = "cluster_rep"))
                                      
                                      list(
                                        k5 = table(fold = cut_k[, "5"],  label = cut_df$label) %>% 
                                          prop.table(margin = 2) %>%
                                          round(2) %>%
                                          addmargins(1) %>% 
                                          addmargins(2),
                                        k7 = table(fold = cut_k[, "7"],  label = cut_df$label) %>% 
                                          prop.table(margin = 2) %>%
                                          round(2) %>%
                                          addmargins(1) %>% 
                                          addmargins(2),
                                        k10 = table(fold = cut_k[, "10"], label = cut_df$label) %>% 
                                          prop.table(margin = 2) %>%
                                          round(2) %>%
                                          addmargins(1) %>% 
                                          addmargins(2))
                                    }

stopCluster(cl)

# At this stage review manually and check gene-by-gene
# Manually designate test fold per segment
folds_list <- vector(mode = "list", length = 8)
test_folds <- c("3", "4", "3", "2", "5", "1", "2", "7")

for (i in 1:8) {
  folds_list[[i]] <- hclust_list[[i]] %>%
    cutree(k = 10) %>%
    as.data.frame() %>%
    rownames_to_column() %>%
    set_colnames(c("gid", "fold"))
  
  folds_list[[i]] %<>%
    mutate(fold = as.character(fold)) %>%
    mutate(fold = ifelse(fold == test_folds[i], "test", fold)) %>%
    mutate(fold = ifelse(!is.na(as.numeric(fold)), as.character(dense_rank(fold)), fold)) # Renumber non-test folds from 1
}

folds_list %>% saveRDS("S3\\data\\segmentwise\\clust_dist\\folds_list.rds")
