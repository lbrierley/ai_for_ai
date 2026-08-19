#######################################################
# Calculate and process machine learning feature sets #
#######################################################

##################################################################
# Parallelise: to revert, comment out and change %dopar% to %do% #
##################################################################

# Set parallelisation
cl <- makePSOCKcluster(12)
registerDoParallel(cl)
clusterSetRNGStream(cl, 1429)

#######################################
# Calculate and save features: k-mers #
#######################################

# Calculate k-mers of entire segment nucleotide sequences

foreach (k_val = 2:6) %:%
  foreach (j = 1:8,
           .packages = c("Biostrings",
                         "coRdon",
                         "dplyr",
                         "magrittr")
  ) %dopar% {
    allflu_nuc_df %>%
      filter(segment == j) %>%
      calc_kmer_counts(k = k_val, overlap = TRUE) %>% 
      saveRDS(paste0("S3\\data\\segmentwise\\mlready/allflu_nuc_",k_val,"mer_pt_",c("PB2", "PB1", "PA", "HA", "NP", "NA", "M1", "NS1")[j],".rds"))
  }

###################################################
# Calculate and save features: genome composition #
###################################################

# Calculate genome composition for individual gene coding sequences

genes <- c("PB2", "PB1", "PA", "HA", "NP", "NA", "M1", "NS1")
foreach(i = seq_along(genes),
        .packages = c("Biostrings",
                      "coRdon",
                      "dplyr",
                      "stringr",
                      "magrittr")
) %dopar% {
  
  source("S3\\scripts\\data_scripts\\01_functions.R")
  
  gene_i <- genes[i]             # for each gene in turn
  if (is.na(gene_i)) {
    df_sub <- allflu_cds_df %>%
      filter(is.na(gene)) %>%
      distinct(gid, .keep_all = TRUE) #If multiple rows contain the same gid, this keeps the first occurance based on row order
    gene_label <- "NA_gene"    #Changing NA gene to NA_gene to avoid error when NA was processed as missing
  } else {
    df_sub <- allflu_cds_df %>%
      filter(gene == gene_i) %>%
      distinct(gid, .keep_all = TRUE)  #Keeping only distinct genes 
    gene_label <- gene_i
  }
  if (nrow(df_sub) == 0) return(NULL)
  
  gids <- df_sub$gid
  cds  <- DNAStringSet(df_sub$string)
  names(cds) <- gids
  temp_df <- cds %>% calc_composition_counts(codonpairs = FALSE) # Calling function
  rm(cds)
  
  temp_df %>%
    saveRDS(paste0( "S3\\data\\segmentwise\\cds/allflu_cds_compcounts_pt_",gene_label,".rds"))
  
  bias_raw <- temp_df %>%
    calc_composition_bias(codonpairs = FALSE)   #Calling function
  
  bias_cols <- colnames(bias_raw)[grepl("_Bias$", colnames(bias_raw))]
  bias_df   <- bind_cols(gid = gids, bias_raw[, bias_cols, drop = FALSE])
  
  bias_df %>%
    saveRDS(paste0("S3\\data\\segmentwise\\cds/allflu_cds_compbias_pt_", gene_label,".rds"))
  
  rm(temp_df, bias_raw, bias_df, df_sub)
  gc()
}

###############################################################################
# Convert .csv outputs from iFeatureOmega to model-ready protein feature sets #
# See: protein_feat_extract.py                                                #
###############################################################################

out_dir <- "S3\\data\\segmentwise\\mlready"
in_dir  <- "S3\\data\\segmentwise\\prot"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
gisaid_keep <- unique(trimws(as.character(meta_ref$Isolate_Id)))

meta_all <- bind_rows(GISAID_avian_prot_df,GISAID_human_prot_df,NCBI_avian_prot_df,NCBI_human_prot_df) %>%
  mutate(gid = trimws(as.character(gid))) %>%
  distinct(gid, .keep_all = TRUE)

foreach (feat = c("2mer", "ctriad", "ctdc", "ctdt", "ctdd", "pseaac")) %:%
  foreach (j = c("PB2", "PB1", "PA", "HA", "NP", "NA", "M1", "NS1"),
           .packages = "tidyverse") %dopar% {
             
             x <- list.files(
               in_dir,
               pattern = paste0(j, "_", feat, "\\.csv$"),
               full.names = TRUE
             ) %>%
               map_dfr(~ read.csv(.x, stringsAsFactors = FALSE))
             
             names(x)[1] <- "gid"
             
             x %>%
               mutate(gid = trimws(as.character(gid))) %>%
               inner_join(meta_all, by = "gid") %>%
               filter(!(grepl("\\|N40\\||\\|M42\\|", fastahead))) %>%
               filter(!grepl("^GISAID", src) | gid %in% gisaid_keep) %>%
               select(-any_of(c(
                 "X", "title", "UID", "subtype", "date", "protINSDC",
                 "protaccession", "gene", "length", "label", "src",
                 "fastahead", "prot_id", "segment", "accession"
               ))) %>%
               relocate(gid) %>%
               saveRDS(file.path(out_dir, paste0("allflu_prot_", feat, "_pt_", j, ".rds")))
           }


stopCluster(cl)
