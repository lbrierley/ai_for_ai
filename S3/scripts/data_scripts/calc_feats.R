#######################################################
# Calculate and process machine learning feature sets #
#######################################################
#######################################
# Calculate and save features: k-mers #
#######################################

for(i in 2:6){
  for(j in 1:8) {
    allflu_nuc_df %>%
      filter(segment == j) %>%
      calc_kmer_counts(k = i, overlap = TRUE) %>% 
      saveRDS(paste0("S3\\data\\full\\mlready\\allflu_nuc_",i,"mer_pt_",c("PB2", "PB1", "PA", "HA", "NP", "NA", "M1", "NS1")[j],".rds"))
  }
}

###################################################
# Calculate and save features: genome composition #
###################################################

# Calculate genomic composition for individual genes - restrict to whole genomes

# Calculate composition counts and biases for individual protein cds (can take a long time)
for(i in (1:length(unique(allflu_cds_df$gene)))){
  
  cds <- allflu_cds_df %>% filter(gene == unique(allflu_cds_df$gene)[i]) %>% pull(string) %>% DNAStringSet()
  names(cds) <- allflu_cds_df %>% filter(gene == unique(allflu_cds_df$gene)[i]) %>% pull(cds_id)
  
  temp_df <- cds %>%
    calc_composition_counts(codonpairs = TRUE)
  
  rm(cds)
  
  temp_df %>% 
    saveRDS(paste0("S3\\data\\full\\cds\\allflu_cds_compcounts_pt_",unique(allflu_cds_df$gene)[i],".rds"))
  
  temp_df %>%
    calc_composition_bias(codonpairs = FALSE) %>%                          # Do not calculate codon pair biases for individual protein cds
    select(cds_id, matches("_Bias$")) %>%
    left_join(allflu_cds_df %>% select(cds_id, gid)) %>%
    relocate(gid) %>%
    saveRDS(paste0("S3\\data\\full\\mlready\\\allflu_cds_compbias_pt_",unique(allflu_cds_df$gene)[i],".rds"))
  
  rm(temp_df)
}

gc()

# Calculate composition biases for whole genome sequences based off of all proteins - limiting to complete whole genomes only
list.files(path = "S3\\data\\full\\cds\\", pattern = "allflu_cds_compcounts_pt_.*\\.rds", full.names = TRUE) %>%
  map_dfr(readRDS) %>%   
  left_join(allflu_cds_df %>% select(cds_id, gid, n, wgs_dup)) %>%
  filter(n >= 8 & wgs_dup == 0) %>% 
  select(-n, -wgs_dup) %>%
  group_by(gid) %>%
  summarise_if(is.numeric, sum) %>%
  ungroup %>%
  distinct %>%
  calc_composition_bias(codonpairs = TRUE) %>%                          # Calculate codon pair biases for wgs
  saveRDS("S3\\data\\full\\cds\\allflu_wgs_compbias.rds")


###############################################################################
# Convert .csv outputs from iFeatureOmega to model-ready protein feature sets #
# See: protein_feat_extract.py                                                #
###############################################################################

for (feat in c("2mer", "ctriad", "ctdc", "ctdt", "ctdd", "pseaac")){
  x <- list.files("S3\\data\\full\\prot\\", pattern = feat, full.names = TRUE) %>%
    map_dfr(read.csv) %>%
    bind_cols(bind_rows(GISAID_avian_prot_df, GISAID_human_prot_df, NCBI_avian_prot_df, NCBI_human_prot_df)) %>% 
    filter(!(grepl("\\|N40\\||\\|M42\\|", fastahead))) %>%
    filter(src == "NCBI" | src == "GISAID" & protaccession %in% meta_ref$value)  # Use canonical GISAID metadata to define which prot sequences belong to which wgs
  
  for (j in c("HA", "M1", "NA", "NP", "NS1", "PA", "PB1", "PB2")){
    x %>%
      filter(gene == j) %>%
      select(-X, -title, -UID, -subtype, -date, -protINSDC, -protaccession, -gene, -length, -label, -src, -fastahead, -prot_id, -segment, -accession) %>%
      relocate(gid) %>%
      saveRDS(paste0("S3\\data\\full\\mlready\\allflu_prot_",feat,"_pt_",j,".rds"))
  }
  
  rm(x)
}
