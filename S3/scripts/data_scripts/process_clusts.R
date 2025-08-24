####################################################################
# Examine and process sequence clustering with selected parameters #
####################################################################
#########################################################
# Select clustering and associated labels for ML models #
#########################################################

# Read in clustering from MMseqs2 and select representative sequences

for(clusterset in list.files(path = "S3\\data\\full\\nuc\\", pattern = paste0("*", cluster_chosen, "_cluster.tsv"))){
  
  cluster_ref <- read.table(paste0("S3\\data\\full\\nuc\\", clusterset), sep = '\t', quote = "\"",  encoding="UTF-8", comment.char = '@', header = FALSE)
  
  # Identify clusters having mixed host labels (i.e., both avian and zoonotic sequences within cluster)
  cluster_ref %<>% left_join((allflu_wgs_df %>% select(gid, label)), by = c("V2" = "gid")) 
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
  
  # Save reference set of which sequences belong to which cluster
  cluster_ref %>% rename(gid = V2) %>% write.csv(paste0("S3\\data\\full\\holdout_clusters\\", gsub("wgs_|cluster.tsv", "", clusterset), "members.csv"))
  
  # Save reference set of cluster-representative sequences
  cluster_ref %>% select(cluster_rep, mix) %>% distinct %>% 
    left_join(allflu_wgs_df %>% select(gid, src, label, subtype, date),  by = c("cluster_rep" = "gid")) %>% 
    write.csv(paste0("S3\\data\\full\\holdout_clusters\\", gsub("wgs_|cluster.tsv", "", clusterset), "labels.csv"))
  
}

####################################################################
# Save chosen clusters protein-by-protein for genome mapping later #
####################################################################

for(foclabel in c("nz", "zoon")){
  
  seqs_to_save <- read.csv(paste0("S3\\data\\full\\holdout_clusters\\full\\ex_full_", cluster_chosen, "_labels.csv")) %>% 
    filter(label == foclabel)
  
  for(focgene in c("HA", "M1", "NA", "NP", "NS1", "PA", "PB1", "PB2")){
    
    # Save nuc sequences for alignment and mapping
    nuc_temp <- seqs_to_save %>% 
      select(cluster_rep) %>% 
      left_join(allflu_nuc_df, by = c("cluster_rep" = "gid")) %>%
      mutate(gene = case_when(gene == "MP" ~ "M1",
                              gene == "NS" ~ "NS1",
                              TRUE ~ gene)) %>%
      filter(gene == focgene)
    
    nuc <- nuc_temp$string %>% DNAStringSet()
    names(nuc) <- nuc_temp$cluster_rep
    writeXStringSet(nuc, filepath = paste0("S3\\data\\full\\mapping\\nuc\\",foclabel,"_clusterreps_",focgene,".FASTA"))
    rm(nuc_temp, nuc)
    
    # Save cds sequences for alignment and mapping
    cds_temp <- seqs_to_save %>%
      select(cluster_rep) %>%
      left_join(allflu_cds_df, by = c("cluster_rep" = "gid")) %>%
      filter(gene == focgene)
    
    cds <- cds_temp$string %>% DNAStringSet()
    names(cds) <- cds_temp$cluster_rep
    writeXStringSet(cds, filepath = paste0("S3\\data\\full\\mapping\\cds\\",foclabel,"_clusterreps_",focgene,".FASTA"))
    rm(cds_temp, cds)
    
    # Save prot sequences for alignment and mapping 
    prot_temp <- seqs_to_save %>%
      select(cluster_rep) %>%
      left_join(bind_rows(GISAID_avian_prot_df, GISAID_human_prot_df, NCBI_avian_prot_df, NCBI_human_prot_df), by = c("cluster_rep" = "gid")) %>%
      filter(gene == focgene)
    
    prot_strings <- c(readAAStringSet("S3\\data\\full\\prot\\GISAID_avian_prot.FASTA"),
                      readAAStringSet("S3\\data\\full\\prot\\GISAID_human_prot.FASTA"),
                      readAAStringSet("S3\\data\\full\\prot\\NCBI_avian_prot.FASTA"),
                      readAAStringSet("S3\\data\\full\\prot\\NCBI_human_prot.FASTA"))
    
    prot_temp <- left_join(prot_temp,
                           data.frame(string = prot_strings %>% as.character(use.names=FALSE),
                                      fastahead = names(prot_strings)))
    
    prot <- prot_temp$string %>% AAStringSet()
    names(prot) <- prot_temp$cluster_rep
    writeXStringSet(prot, filepath = paste0("S3\\data\\full\\mapping\\prot\\",foclabel,"_clusterreps_",focgene,".FASTA"))
    rm(prot_temp, prot, prot_strings)
    
  }
}