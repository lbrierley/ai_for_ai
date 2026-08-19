library(foreach)
library(parallel)
library(doParallel)

# Read in reference of all sequences used
allflu_cds_ref <- read.csv("S3\\data\\segmentwise\\Raw_Sequences\\allflu_cds_df.csv", na.strings = "")

# Set extended MMseqs2 Linclust parameters to try
cluster_minseqid <- c(0.85, 0.90, 0.925, 0.95, 0.975, 0.99)
cluster_C <- c(0.7, 0.8, 0.90, 0.95)
combos <- crossing(seg = 1:8, minseqid = cluster_minseqid, C = cluster_C)

# # Save each segment sequences for clustering from old wgs version
# for (seg in 1:8){
#   seq_df <- allflu_cds_df %>% filter(gid %in% allflu_wgs_ref$gid & segment == seg) %>% as.data.frame
#   seq <- seq_df$string %>% DNAStringSet()
#   names(seq) <- seq_df$gid
#   writeXStringSet(seq, filepath = paste0("S3\\data\\segmentwise\\allflu_cds_seg_", seg, "_full.FASTA"))
# }

# # Set parallelisation
# cl <- makePSOCKcluster(4)
# registerDoParallel(cl)
# clusterSetRNGStream(cl, 1429)

# Call MMseqs2 14-7e284 to cluster on similarity score - replace with your own MMseqs2 path if needed
foreach(seg = 1:8) %:%
  foreach(minseqid_param = cluster_minseqid) %:%
  foreach(C_param = cluster_C) %do% {
    
    if(!file.exists(paste0("S3\\data\\segmentwise\\Clustering\\clust_", seg, "_", minseqid_param*100, "_", C_param*10, "_cluster.tsv"))){
      system(paste0("S3\\scripts\\mmseqs-win64\\mmseqs.bat easy-linclust S3\\data\\segmentwise\\allflu_cds_seg_", seg, "_full.FASTA S3\\data\\segmentwise\\Clustering\\clust_", seg, "_", minseqid_param*100, "_", C_param*10, " tmp --min-seq-id ", minseqid_param ," -c ", C_param, " --cov-mode 0"),
             intern = TRUE, show.output.on.console = TRUE, ignore.stdout=FALSE, ignore.stderr=FALSE, wait=FALSE)
      
      file.remove(list.files(path = "S3\\data\\segmentwise\\Clustering\\", pattern = "*.fasta", full.names = TRUE)) # Keep only the .tsv output from MMseqs2
    }
  }
stopCluster(cl)





# Take clusters and match them to GID to find dataset sizes of cluster-centroids (regardless of stripped-out test set for  n)
tabulate_n_cluster_combos <- function(file){
  cluster_ref <- read.table(file, sep = '\t', quote = "\"",  encoding="UTF-8", comment.char = '@', header = FALSE)
  
  # Identify clusters having mixed host labels (i.e., both avian and zoonotic sequences within cluster)
  cluster_ref %<>% left_join((allflu_cds_ref %>% select(gid, label) %>% distinct), by = c("V2" = "gid")) %>% filter(complete.cases(label))
  cluster_ref %<>% left_join(cluster_ref %>% group_by(V1) %>% summarise(mix = n_distinct(label)))
  
  # cluster_ref %>% distinct(V1, mix) %>% with(., table(mix)) %>% prop.table() %>% print()            # Check % of mixed-label clusters
  # cluster_ref %>% group_by(V1) %>% summarise(n = n_distinct(V2)) %>% with(., table(n)) %>% print()  # Check distribution of cluster sizes
  # cluster_ref %>% filter(mix == 2)                                                                  # Check identity of mixed-label clusters
  
  # For mixed-label clusters, if cluster representative is not zoonotic, then select a random zoonotic representative instead
  cluster_ref %<>% mutate(manual_cluster_rep = case_when(
    mix == 2 & !(V1 %in% (cluster_ref %>% filter(mix == 2 & V1 == V2 & label == "zoon") %>% pull(V1))) ~ 1,
    TRUE ~ 0)
  )
  
  set.seed(1516)
  manual_cluster_df <- cluster_ref %>% filter(manual_cluster_rep == 1 & label == "zoon") %>% group_by(V1) %>% slice_sample(n = 1) %>% rename(cluster_rep = V2) %>% ungroup
  
  cluster_ref %<>% left_join(manual_cluster_df %>% select(V1, cluster_rep), by = "V1") %>% mutate(cluster_rep = coalesce(cluster_rep, V1)) %>% select(cluster_rep, V2, mix)

  # Pull reference set of cluster-representative sequences
  cluster_ref %>% select(cluster_rep, mix) %>% distinct %>% 
    left_join((allflu_cds_ref %>% select(gid, src, label, subtype) %>% distinct),  by = c("cluster_rep" = "gid")) %>%
    with(., table(label)) %>%
    addmargins() %>%
    as.vector

}

results <- combos %>%
  mutate(
    file = paste0("S3\\data\\segmentwise\\Clustering\\clust_", seg, "_", minseqid*100, "_", C*10, "_cluster.tsv"),
    output = map(file, possibly(tabulate_n_cluster_combos, otherwise = NULL))
  ) %>%
  unnest_wider(output, names_sep = "temp") %>%
  select(-file) %>%
  set_colnames(c("seg", "minseqid", "C", "nz", "zoon", "total"))

results %>% write.csv("S3\\data\\segmentwise\\cds_clust_sizes.csv")

for(i in 1:8){
  results %>%
    filter(seg == i) %>%
    mutate(cell = paste0(zoon, "  /  ", total, "  (", round((total-zoon)/zoon, 1), ")" )) %>%
    select(minseqid, C, cell) %>%
    pivot_wider(names_from = C, values_from = cell) %>% write.csv(paste0("temp_",i,".csv"))
}



# Save chosen clustering set as processed labels

# Directories
cluster_dir <- "S3\\data\\segmentwise\\clust_final/"
outdir <- "S3\\data\\segmentwise\\clust_processed/"


cluster_files <- list.files(path = cluster_dir,pattern = "_cluster.tsv", full.names = TRUE)

cluster_chosen_list <- NULL

for (i in 1:length(cluster_files)) {
  cluster_ref <- read.table(cluster_files[i], sep = "\t",quote = "\"",encoding = "UTF-8",comment.char = "@",header = FALSE)
  
  # Identify clusters having mixed host labels (i.e., both avian and zoonotic sequences within cluster)
  cluster_ref %<>% left_join((allflu_cds_ref %>% select(gid, label) %>% distinct), by = c("V2" = "gid")) %>% filter(complete.cases(label))
  cluster_ref %<>% left_join(cluster_ref %>% group_by(V1) %>% summarise(mix = n_distinct(label)))
  
  # cluster_ref %>% distinct(V1, mix) %>% with(., table(mix)) %>% prop.table() %>% print()            # Check % of mixed-label clusters
  # cluster_ref %>% group_by(V1) %>% summarise(n = n_distinct(V2)) %>% with(., table(n)) %>% print()  # Check distribution of cluster sizes
  # cluster_ref %>% filter(mix == 2)                                                                  # Check identity of mixed-label clusters
  
  # For mixed-label clusters, if cluster representative is not zoonotic, then select a random zoonotic representative instead
  cluster_ref %<>% mutate(manual_cluster_rep = case_when(
    mix == 2 & !(V1 %in% (cluster_ref %>% filter(mix == 2 & V1 == V2 & label == "zoon") %>% pull(V1))) ~ 1,
    TRUE ~ 0)
  )
  
  set.seed(1516)
  manual_cluster_df <- cluster_ref %>% filter(manual_cluster_rep == 1 & label == "zoon") %>% group_by(V1) %>% slice_sample(n = 1) %>% rename(cluster_rep = V2) %>% ungroup
  
  cluster_ref %<>% left_join(manual_cluster_df %>% select(V1, cluster_rep), by = "V1") %>% mutate(cluster_rep = coalesce(cluster_rep, V1)) %>% select(cluster_rep, V2, mix)
  
  # Save reference set of clusters
  cluster_chosen_list[[i]] <- cluster_ref
  
  # Save cluster members
  write.csv(
    cluster_ref %>% rename(gid = V2) %>% distinct(),
    file = paste0(outdir, basename(cluster_files[i]), "_members.csv"),
    row.names = FALSE
  )
  
  write.csv(cluster_ref %>% select(cluster_rep, mix) %>% distinct() %>%
              left_join(allflu_cds_ref %>% select(gid, src, label, subtype) %>% distinct(),
        by = c("cluster_rep" = "gid")
      ),
    file = paste0(outdir, basename(cluster_files[i]), "_labels.csv"),
    row.names = FALSE
  )

}

















