#######################################
# Process other sequences of interest #
#######################################

library(caretEnsemble)

# List available external seqs
list.files("S3\\data\\irat\\", pattern = ".fasta", ignore.case=TRUE) %>% gsub("GISAID_|NCBI_|_nuc|_cds|_prot|.fasta|.FASTA", "", .) %>% unique

# Define seqs to be used
set <- "irat"

################################
# Read-in and process metadata #
################################

# Read-in and construct summary df, tidy up where known to be mislabelled
# Nucleotide sequences of segments
g_nuc <- process_GISAID_seq(x = readSet(file = paste0("S3\\data\\irat\\GISAID_", set, "_nuc.fasta")), label = set, type = "nuc")
n_nuc <- process_NCBI_seq(x = readSet(file = paste0("S3\\data\\irat\\NCBI_", set, "_nuc.fasta")), label = set, type = "nuc")

# Coding sequences for proteins
g_orfs <- g_nuc %>% pull(string) %>% lapply(bind_ORF) %>% bind_rows() # Create GISAID coding sequences from transcripts
g_cds <- g_nuc %>% bind_cols(g_orfs) %>% rowwise %>% mutate(string = substr(string, start, end)) %>% as.data.frame

n_cds <- process_NCBI_seq(x = readSet(file =paste0("S3\\data\\irat\\NCBI_", set, "_cds.fasta")), label = "irat", type = "cds")

# Clean names of protein sequence FASTA files and resave for iFeatureOmega
files <- list.files(path = "S3\\data\\irat\\", pattern = "prot.fasta", full.names = TRUE, ignore.case = TRUE)
for(i in 1:length(files)){
  prot_fasta_name_clean(files[i])
}

g_prot <- process_GISAID_seq(x = readAAStringSet(file = paste0("S3\\data\\irat\\GISAID_", set, "_prot.fasta")), label = set, type = "prot")
n_prot <- process_NCBI_seq(x = readAAStringSet(file = paste0("S3\\data\\irat\\NCBI_", set, "_prot.fasta")), label = set, type = "prot") %>%
  select(-string, -gid) %>%
  left_join(n_nuc %>% select(gid, title) %>% distinct, by = "title")

##################################################
# Process nucleotide, coding sequences, proteins #
##################################################

nuc_df <- bind_rows(g_nuc, 
                    n_nuc)

cds_df <- bind_rows(g_cds, 
                    n_cds) %>%  # Relabel GISAID genes to reflect likely ORF captured (M1 for MP, NS1 for NS), filter out NCBI genes not present in GISAID (M2, NS2, PA-X, PB1-F2)
  mutate(gene = case_when(
    gene == "MP" ~ "M1",
    gene == "NS" ~ "NS1",
    TRUE ~ gene
  )) %>%
  filter(!(gene %in% c("M2", "NS2", "PA-X", "PB1-F2"))) %>%
  select(-start, -end) %>%
  mutate(cds_id = accession) %>%
  select(-INSDC, -accession)

prot_df <- bind_rows(g_prot,
                     n_prot) %>%
  mutate(gene = case_when(
    gene == "MP" ~ "M1",
    gene == "NS" ~ "NS1",
    TRUE ~ gene
  ))
 
#######################################################
# Calculate and process machine learning feature sets #
#######################################################
#######################################
# Calculate and save features: k-mers #
#######################################

for(i in 2:6){
  for(j in 1:8) {
    nuc_df %>%
      filter(segment == j) %>%
      calc_kmer_counts(k = i, overlap = TRUE) %>% 
      saveRDS(paste0("S3\\data\\irat\\mlready\\", set, "\\allflu_nuc_",i,"mer_pt_",c("PB2", "PB1", "PA", "HA", "NP", "NA", "M1", "NS1")[j],".rds"))
  }
}

###################################################
# Calculate and save features: genome composition #
###################################################

# Calculate genomic composition for individual genes - restrict to whole genomes

# Calculate composition counts and biases for individual protein cds
for(i in (1:length(unique(cds_df$gene)))){
  
  cds <- cds_df %>% filter(gene == unique(cds_df$gene)[i]) %>% pull(string) %>% DNAStringSet()
  names(cds) <- cds_df %>% filter(gene == unique(cds_df$gene)[i]) %>% pull(gid)
  
  temp_df <- cds %>%
    calc_composition_counts(codonpairs = TRUE)
  
  rm(cds)
  
  temp_df %>% 
    saveRDS(paste0("S3\\data\\irat\\cds\\allflu_cds_compcounts_pt_", set, "_", unique(cds_df$gene)[i],".rds"))
  
  temp_df %>%
    calc_composition_bias(codonpairs = FALSE) %>%                     # Do not calculate codon pair biases for individual protein cds
    select(cds_id, matches("_Bias$")) %>%
    rename(gid = cds_id) %>%
    saveRDS(paste0("S3\\data\\irat\\mlready\\", set, "\\allflu_cds_compbias_pt_",unique(cds_df$gene)[i],".rds"))
  
  rm(temp_df)
  
}

gc()

# # Calculate composition biases for whole genome sequences based off of all proteins - limiting to complete whole genomes only
# list.files(path = "S3\\data\\irat\\cds\\", pattern = paste0("allflu_cds_compcounts_pt_", set, "_.*\\.rds"), full.names = TRUE) %>%
#   map_dfr(readRDS) %>%
#   left_join(cds_df %>% select(cds_id, gid)) %>%
#   group_by(gid) %>%
#   summarise_if(is.numeric, sum) %>%
#   ungroup %>%
#   distinct %>%
#   calc_composition_bias(codonpairs = TRUE) %>%                          # Calculate codon pair biases for wgs
#   saveRDS(paste0("S3\\data\\irat\\cds\\allflu_wgs_", set, "_compbias.rds"))

######################################################################
# Check pre-calculated and saved protein features from iFeatureOmega # Note these need to run on WSL installation of Python, and after cleaning FASTA with custom function
######################################################################

for (feat in c("2mer", "ctriad", "ctdc", "ctdt", "ctdd", "pseaac")){
  x <- list.files("S3\\data\\irat\\prot\\", full.names = TRUE) %>% .[grepl(feat, .)] %>%
    map_dfr(read.csv) %>%
    bind_cols(prot_df) %>% 
    filter(!(gene %in% c("M2", "NS2", "PA-X", "PB1-F2"))) %>%
    filter(!(grepl("\\|N40\\||\\|M42\\|", fastahead)))
  
  for (j in c("HA", "M1", "NA", "NP", "NS1", "PA", "PB1", "PB2")){
    x %>%
      filter(gene == j) %>%
      select(-X, -title, -UID, -subtype, -date, -protINSDC, -protaccession, -gene, -length, -label, -src, -fastahead) %>%
      relocate(gid) %>%
      saveRDS(paste0("S3\\data\\irat\\mlready\\", set, "\\allflu_prot_",feat,"_pt_",j,".rds"))
  }
  
  rm(x)
}

##############
# Stack test #
##############

# Restrict to only complete sequences
nuc_df %<>% filter(gid %in% (nuc_df %>% group_by(gid) %>% tally %>% filter(n == 8) %>% pull(gid)))
cds_df %<>% filter(gid %in% (cds_df %>% group_by(gid) %>% tally %>% filter(n == 8) %>% pull(gid)))
prot_df %<>% filter(!(gene %in% c("M2", "NS2", "PA-X", "PB1-F2"))) %>%
  filter(gid %in% (prot_df %>% filter(!(gene %in% c("M2", "NS2", "PA-X", "PB1-F2"))) %>% group_by(gid) %>% tally %>% filter(n == 8) %>% pull(gid)))

# Formalise list of complete sequences
if (length(unique(nuc_df$gid)) <= length(unique(prot_df$gid))){
  canon <- nuc_df %>% select(gid) %>% distinct
} else {
  canon <- prot_df %>% select(gid) %>% distinct
}

irat_test <- purrr::map(list.files(path = paste0("H:\\Working\\ai_for_ai\\S3\\data\\irat\\mlready\\", set), full.names = TRUE),
                        function (x) 
                          readRDS(x) %>%
                          select(-any_of(c("segment", "cds_id", "enc", "GC_content", "prot_id", "segment", "accession"))) %>%
                          rename_with(~paste(., gsub(".*_pt_|.rds", "", x), sep = "_"), -c(gid)) %>%
                          right_join(canon) %>%                                                            
                          arrange(gid) %>%
                          select(-gid)
) %>% 
  purrr::list_cbind(name_repair = "unique_quiet")

# Check failing n rows, and which columns with NAs
which(is.na(irat_test), arr.ind=TRUE) %>% as.data.frame() %>% pull(col) %>% unique %>% names(irat_test)[.]
irat_test %>% filter(!(complete.cases(.))) %>% nrow

# Bind in gid and remove rows with NAs
irat_test <- bind_cols(gid = canon, irat_test)
irat_test %<>% filter(complete.cases(.))

# Generate predictions over all stacks
pred <- lapply(list.files("S3\\analysis\\stacks\\", pattern = ".rds", full.names = TRUE),
               function(x)
                 
                 data.frame(set = set,
                            stack = gsub("S3\\\\analysis\\\\stacks\\\\stack_|.rds", "", x),
                            hzoon = predict(readRDS(x), newdata=irat_test, type = "prob"),
                            gid = irat_test$gid) %>%
                 left_join(g_nuc %>% select(UID, date) %>% distinct, by = c("gid" = "UID")) 
               
               
) %>%
  bind_rows

pred %>% write.csv(paste0("S3\\data\\irat\\preds_", set, "_raw.csv"), row.names=FALSE)

pred_weight <- lapply(list.files("S3\\analysis\\stacks_weight\\", pattern = ".rds", full.names = TRUE),
               function(x)
                 
                 data.frame(set = set,
                            stack = gsub("S3\\\\analysis\\\\stacks_weight\\\\stack_|.rds", "", x),
                            hzoon = predict(readRDS(x), newdata=irat_test, type = "prob"),
                            gid = irat_test$gid) %>%
                 left_join(g_nuc %>% select(UID, date) %>% distinct, by = c("gid" = "UID")) 
               
               
) %>%
  bind_rows

pred_weight %>% write.csv(paste0("S3\\data\\irat\\preds_weight_", set, "_raw.csv"), row.names=FALSE)

## Simple plot for first visualisation
# Point estimate for a given stack

# # Determine a threshold in classification task: zoonotic H5N1 vs all non-zoonotic subtypes held out
# stacked_raw <- read.csv("S3/analysis/stack_subtypeacc_raw.csv") %>%
#   filter(subtype %in% c("H5N1", "H4N6", "H4N8", "H8N4", "H16N3"))
# 
# new_threshold <- pROC::roc(response = stacked_raw$label,
#                  predictor = stacked_raw$hzoon,
#                  direction = ">") %>%
#   pROC::coords("best", best.method="closest.topleft") %>%
#   .$threshold
