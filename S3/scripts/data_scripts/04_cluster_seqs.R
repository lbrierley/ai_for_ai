###########################################
# Cluster to collapse redundant sequences #
###########################################
##########################################################################
# Calculate clustering based on k-mer overlaps, considering holdout sets #
##########################################################################

wgs_all <- allflu_wgs_df %>% 
  filter(wgs_dup == 0) %>% 
  inner_join(allflu_cds_df %>% filter(n >= 8) %>% select(gid) %>% distinct())   # Select only wgs having all 8 cds

wgs_all %>%
  select(gid, title, subtype, label, src, date) %>%
  write.csv("S3\\data\\full\\allflu_wgs_ref.csv")

## Cluster complete dataset 

# Save whole genome sequences (complete coding sequences for all 8 segments) for clustering
wgs <- wgs_all$wgs_string %>% DNAStringSet()
names(wgs) <- wgs_all$gid
writeXStringSet(wgs, filepath = paste0("S3\\data\\full\\allflu_nuc_wgs_full.FASTA"))
rm(wgs_all, wgs)

# Call MMseqs2 14-7e284 to cluster on similarity score - replace with your own MMseqs2 path if needed
for (minseqid_param in cluster_minseqid){
  
  for(C_param in cluster_C){
    
    system(paste0("mmseqs.bat easy-linclust S3\\data\\full\\allflu_nuc_wgs_full.FASTA S3\\data\\full\\nuc\\wgs_full_", minseqid_param*100, "_", C_param*10, " tmp --min-seq-id ", minseqid_param ," -c ", C_param, " --cov-mode 0"),
           intern = TRUE, show.output.on.console = TRUE, ignore.stdout=FALSE, ignore.stderr=FALSE, wait=FALSE)
    
    file.remove(list.files(path = "S3\\data\\full\\nuc\\", pattern = "*.fasta", full.names = TRUE)) # Keep only the .tsv output from MMseqs2
    
  }
}

## Cluster datasets for holdout training, holding out specific test subtype each time

for(holdout in c(holdout_nz, holdout_zoon)){
  
  # Save whole genome sequences with 8 coding sequences for clustering
  wgs_temp <- allflu_wgs_df %>% 
    filter(wgs_dup == 0 & subtype != holdout) %>% 
    inner_join(allflu_cds_df %>% filter(n >= 8) %>% select(gid) %>% distinct()) %>%   # Select only wgs having all 8 cds
    select(gid, wgs_string)
  wgs <- wgs_temp$wgs_string %>% DNAStringSet()
  names(wgs) <- wgs_temp$gid
  writeXStringSet(wgs, filepath = paste0("S3\\data\\full\\allflu_nuc_wgs_ex_",holdout,".FASTA"))
  rm(wgs_temp, wgs)
  
}

# Call MMseqs2 14-7e284 to cluster on similarity score - replace with your own MMseqs2 path if needed
for(holdout in c(holdout_nz, holdout_zoon)){
  
  for (minseqid_param in cluster_minseqid){
    
    for(C_param in cluster_C){
      
      system(paste0("mmseqs.bat easy-linclust S3\\data\\full\\allflu_nuc_wgs_ex_", holdout, ".FASTA S3\\data\\full\\nuc\\wgs_ex_", holdout, "_", minseqid_param*100, "_", C_param*10, " tmp --min-seq-id ", minseqid_param ," -c ", C_param, " --cov-mode 0"),
             intern = TRUE, show.output.on.console = TRUE, ignore.stdout=FALSE, ignore.stderr=FALSE, wait=FALSE)
      
      file.remove(list.files(path = "S3\\data\\full\\nuc\\", pattern = "*.fasta", full.names = TRUE)) # Keep only the .tsv output from MMseqs2
      
    }
  }
  
}
