###############################################################
# Process sequences manually downloaded from NCBI Flu, GISAID #
###############################################################
################################
# Read-in and process metadata #
################################

# Read-in
# Nucleotide sequences of segments

GISAID_avian_nuc <- readSet(file = "S3\\data\\full\\GISAID_avian_nuc.fasta")
GISAID_human_nuc <- readSet(file = "S3\\data\\full\\GISAID_human_nuc.fasta")

NCBI_avian_nuc <- readSet(file = "S3\\data\\full\\NCBI_avian_nuc.fasta")
NCBI_human_nuc <- readSet(file = "S3\\data\\full\\NCBI_human_nuc.fasta")

# Coding sequences for proteins

NCBI_avian_cds <- readSet(file = "S3\\data\\full\\NCBI_avian_cds.fasta")
NCBI_human_cds <- readSet(file = "S3\\data\\full\\NCBI_human_cds.fasta")

#   # Proteins
#
# GISAID_avian_prot <- read.fasta(file = "S3\\data\\full\\GISAID_avian_prot.fasta")
# GISAID_human_prot <- read.fasta(file = "S3\\data\\full\\GISAID_human_prot.fasta")
#
# NCBI_avian_prot <- read.fasta(file = "S3\\data\\full\\NCBI_avian_prot.fasta")
# NCBI_human_prot <- read.fasta(file = "S3\\data\\full\\NCBI_human_prot.fasta")
#
# # Metadata
#
# GISAID reads in from file
GISAID_avian_meta <- readxl::read_excel("S3/data/full/GISAID_avian_meta.xls", guess_max = 1048576)
GISAID_human_meta <- readxl::read_excel("S3/data/full/GISAID_human_meta.xls", guess_max = 1048576)

meta_ref <- bind_rows(
  GISAID_avian_meta %>%
    select(1:9) %>% reshape2::melt(id.vars = "Isolate_Id") %>% mutate(
      value = gsub("EPI", "", gsub("\\|.*", "", value)),
      label = "nz"
    ),
  GISAID_human_meta %>% select(1:9) %>% reshape2::melt(id.vars = "Isolate_Id") %>% mutate(
    value = gsub("EPI", "", gsub("\\|.*", "", value)),
    label = "zoo"
  )
)

# Clean metadata

# Construct summary df for individual genes from particular FASTA headers and tidy up where known to be mislabelled

NCBI_human_nuc_df <- process_NCBI_DNA(x = NCBI_human_nuc, label = "zoon", type = "nuc")
NCBI_human_cds_df <- process_NCBI_DNA(x = NCBI_human_cds, label = "zoon", type = "cds")

NCBI_avian_nuc_df <- process_NCBI_DNA(x = NCBI_avian_nuc, label = "nz", type = "nuc")
NCBI_avian_cds_df <- process_NCBI_DNA(x = NCBI_avian_cds, label = "nz", type = "cds")

GISAID_human_nuc_df <- process_GISAID_DNA(x = GISAID_human_nuc, label = "zoon", type = "nuc")
GISAID_avian_nuc_df <- process_GISAID_DNA(x = GISAID_avian_nuc, label = "nz", type = "nuc")

# Adjust where FASTA headers have pulled wrong information directly
GISAID_avian_nuc_fix_ref <- GISAID_avian_nuc[!(GISAID_avian_nuc_df$gene %in% c("PB1", "PB2", "PA", "HA", "NP", "NA", "MP", "NS")), ] %>%
  names() %>%
  as.data.frame() %>%
  setNames("title") %>%
  tidyr::separate(title, sep = "\\|", into = c("nn1", "nn9", "nn2", "nn3", "nn4", "INSDC", "nn6", "nn7", "nn8", "gene", "segment"), extra = "drop") %>%
  select(-contains("nn"))

GISAID_avian_nuc_df %<>%
  rows_update(., GISAID_avian_nuc_fix_ref, by = "INSDC")


#############################################################
# Process whole genome sequences, i.e. full-length segments #
#############################################################

# Combine all data - nucleotide sequences of segments
allflu_nuc_df <- bind_rows(GISAID_human_nuc_df, 
                           GISAID_avian_nuc_df, 
                           NCBI_human_nuc_df, 
                           NCBI_avian_nuc_df) %>%
  
  # Remove mixed infections (usually signified by 'MIXED' in NCBI data but also slashes, periods, commas and by H0N0 in GISAID data but some labelled as such in title), (removes 14,733 individual segment sequences)
  filter(!grepl("MIXED|\\,|\\.|\\/", subtype)) %>%
  filter(subtype != "H0N0" & !(grepl("Mixed", title))) %>%
  
  # Remove sequences known a priori to be mislabelled (removes 24 individual segment sequences, 3 wgs)
  filter(title != "A/goose/Chiayi/18020014-1/2018" & title != "H5N6")

# Check to see if multiple of the same segment per genome id:
allflu_nuc_df %>%
  group_by(gid) %>%
  summarise(seqs = n(), segs = n_distinct(segment), genes = n_distinct(gene)) %>%
  filter(seqs != segs | seqs != genes)

# How many sequences are there for each?
allflu_nuc_df %>%
  group_by(gid, label, src) %>%
  summarise(seqs = n(), segs = n_distinct(segment), genes = n_distinct(gene)) %>%
  with(., table(seqs, paste(label, src)))

# Retain canonical sequences belonging to each whole genome according to GISAID metadata (and keep all NCBI) (removes 575 individual segment sequences)
allflu_nuc_df %<>% 
  filter(is.na(UID) | INSDC %in% meta_ref$value)

# Assign whole genomes (all segments)
# Where there is a multiple of 8 segments, assign each run as an individual wgs
ids <- allflu_nuc_df %>%
  group_by(gid) %>%
  summarise(seqs = n()) %>%
  filter(seqs%%8 == 0 & seqs > 8) %>% pull(gid)

temp_8s <- allflu_nuc_df %>% filter(gid %in% ids) %>% group_by(gid,segment) 
temp_8s <- bind_rows(temp_8s %>% slice(1), # Arrange such that we get segment 1-8 in each instance, regardless of original row order. 
                     temp_8s %>% slice(2),
                     temp_8s %>% slice(3),
                     temp_8s %>% slice(4)) # Max seqs per gid = 32, so 4 sets
temp_8s %<>% 
  group_by(gid) %>% 
  mutate(tempalpha = rep(LETTERS, each=8, length=n())) %>%
  mutate(gid = paste0(gid, "_", tempalpha)) %>%
  select(-tempalpha) %>%
  ungroup()

allflu_nuc_df %<>% 
  filter(!(gid %in% ids)) %>%
  bind_rows(temp_8s)

# For NCBI data, where there is > 8 segments and not a multiple of 8, take the first instance of each segment as the wgs (removes 339 individual segment sequences)
allflu_nuc_df %<>% 
  group_by(gid, segment) %>%
  slice(1) %<>%
  ungroup()

# Remove based on acceptable lower bounds of segment lengths (removes 46 individual segment sequences and 38 wgs)
allflu_nuc_df %<>% filter(segment == 1 & length >= 2100 |
                            segment == 2 & length >= 2200 |
                            segment == 3 & length >= 2000 |
                            segment == 4 & length >= 1650 |
                            segment == 5 & length >= 1450 |
                            segment == 6 & length >= 1300 |
                            segment == 7 & length >= 900 |
                            segment == 8 & length >= 800) 

# Work out duplicated whole genome sequences and duplicated individual segments and label
allflu_wgs_df <- allflu_nuc_df %>% 
  arrange(gid, segment) %>% # arrange segments in order
  group_by(src, label, subtype, title, date, gid) %>%
  summarise(n = n(), wgs_string = paste(string, collapse="")) %>%
  ungroup() %>%
  arrange(desc(label), desc(src), title, gid) %>%  # arrange such that when calculating duplicates, zoonotic seqs and NCBI seqs are preferentially retained
  mutate(wgs_dup = if_else(duplicated(wgs_string) == TRUE, 1, 0))

allflu_nuc_df %<>% 
  left_join(allflu_wgs_df %>% select(-wgs_string)) %>%
  arrange(desc(label), desc(src), title, gid, segment) %>%  # arrange such that when calculating duplicates, zoonotic seqs and NCBI seqs are preferentially retained
  mutate(seg_dup = case_when(duplicated(string) == TRUE ~ 1, 
                             duplicated(accession) == TRUE ~ 1,
                             TRUE ~ 0))

# # Plot lengths
# # All retained segment sequences
# allflu_nuc_df %>%
#   filter(seg_dup == 0) %>%
#   mutate(group = paste(src, label)) %>%
#   ggplot(aes(x = length, fill = group)) +
#   geom_density(alpha = 0.25) +
#   facet_wrap(~segment, ncol = 2, scales = "free")
# 
# # Plot segment sequences belonging to whole genomes only
# allflu_nuc_df %>%
#   filter(gid %in% (allflu_nuc_df %>% filter(wgs_dup == 0) %>% group_by(gid) %>% count() %>% arrange(-n) %>% filter(n == 8) %>% pull(gid))) %>%
#   mutate(group = paste(src, label)) %>%
#   ggplot(aes(x = length, fill = group)) +
#   geom_density(alpha = 0.25) +
#   facet_wrap(~segment, ncol = 2, scales = "free")

# Create combined DNAStringSet and remove unfiltered sets from working memory
rm(GISAID_avian_nuc, GISAID_human_nuc, NCBI_avian_nuc, NCBI_human_nuc)


## CALCULATE K-MERS OF WGS (I.E. STITCHING SEGMENTS)



######################################################
# Process gene sequences, i.e. only coding sequences #
######################################################

# Create GISAID coding sequences from transcripts
# GISAID_avian_orfs <- allflu_nuc_df %>% filter(src == "GISAID" & label == "nz") %>% pull(string) %>% lapply(bind_ORF) %>% bind_rows()
# GISAID_human_orfs <- allflu_nuc_df %>% filter(src == "GISAID" & label == "zoon") %>% pull(string) %>% lapply(bind_ORF) %>% bind_rows()

# write.csv(GISAID_avian_orfs, "S3\\data\\full\\GISAID_avian_orfs_temp.csv")
# write.csv(GISAID_human_orfs, "S3\\data\\full\\GISAID_human_orfs_temp.csv")

GISAID_avian_orfs <- read.csv("S3\\data\\full\\GISAID_avian_orfs_temp.csv")
GISAID_human_orfs <- read.csv("S3\\data\\full\\GISAID_human_orfs_temp.csv")

GISAID_avian_cds_df <- allflu_nuc_df %>% filter(src == "GISAID" & label == "nz") %>% bind_cols(GISAID_avian_orfs) %>% rowwise %>% mutate(string = substr(string, start, end)) %>% as.data.frame
GISAID_human_cds_df <- allflu_nuc_df %>% filter(src == "GISAID" & label == "zoon") %>% bind_cols(GISAID_human_orfs) %>% rowwise %>% mutate(string = substr(string, start, end)) %>% as.data.frame

# Combine all data and remove duplicates - coding sequences
allflu_cds_df <- bind_rows(GISAID_human_cds_df, 
                           GISAID_avian_cds_df, 
                           NCBI_human_cds_df %>% 
                             select(-gid) %>% inner_join(allflu_nuc_df %>%                  # Overwrite gids with updated ones from allflu_nuc_df, which accounts for multiple sets of 8 within same title
                                                           filter(src == "NCBI") %>%
                                                           select(gid, accession, wgs_dup, seg_dup), by = "accession"),
                           NCBI_avian_cds_df %>%
                             select(-gid) %>% inner_join(allflu_nuc_df %>% 
                                                           filter(src == "NCBI") %>%
                                                           select(gid, accession, wgs_dup, seg_dup), by = "accession")
) %>%
  
  # Coalesce unique ID columns (INSDC for GISAID, cds_id for NCBI)
  mutate(cds_id = coalesce(cds_id, INSDC)) %>%
  select(-INSDC) %>%
  
  # Remove coding sequences that do not start ATG (removes 32 individual coding sequences)
  filter(grepl("^ATG", string)) %>%
  
  # Remove partial sequences that should not be present (removes 74 individual coding sequences)
  filter(!(grepl(">|<", fastahead))) %>%
  
  # Remove mixed infections (usually signified by 'MIXED' in NCBI data but also slashes, periods, commas and by H0N0 in GISAID data but some labelled as such in title), (removes 19,187 individual coding sequences)
  filter(!grepl("MIXED|\\,|\\.|\\/", subtype)) %>%
  filter(subtype != "H0N0" & !(grepl("Mixed", title))) %>%
  
  # Remove sequences known a priori to be mislabelled (removes 34 individual coding sequences)
  filter(title != "A/goose/Chiayi/18020014-1/2018" & title != "H5N6")

# Relabel GISAID genes to reflect likely ORF captured (M1 for MP, NS1 for NS), filter out NCBI genes not present in GISAID (M2, NS2, PA-X, PB1-F2)
allflu_cds_df %<>% 
  mutate(gene = case_when(
    gene == "MP" ~ "M1",
    gene == "NS" ~ "NS1",
    TRUE ~ gene
  )) %>%
  filter(!(gene %in% c("M2", "NS2", "PA-X", "PB1-F2"))) %>%
  select(-X, -start, -end)

# Check to see if multiple of the same gene per genome id:
allflu_cds_df %>%
  group_by(gid) %>%
  summarise(seqs = n(), segs = n_distinct(segment), genes = n_distinct(gene)) %>%
  filter(seqs != genes)

# How many sequences are there for each?
allflu_cds_df %>%
  group_by(gid, label, src) %>%
  summarise(seqs = n(), segs = n_distinct(segment), genes = n_distinct(gene)) %>%
  with(., table(seqs, paste(label, src)))

# Remove based on acceptable lower bounds of gene lengths (removes 7 individual coding sequences)
allflu_cds_df %<>% filter(gene == "PB2" & length >= 2200 |
                            gene == "PB1-F2" & length >= 150 |
                            gene == "PB1" & length >= 2250 |
                            gene == "PA-X" & length >= 700 |
                            gene == "PA" & length >= 2100 |
                            gene == "HA" & length >= 1600 |
                            gene == "NP" & length >= 1400 |
                            gene == "NA" & length >= 1300 |
                            gene == "M2" & length >= 250 |
                            gene == "M1" & length >= 750 |
                            gene == "NS2" & length >= 350 |
                            gene == "NS1" & length >= 600)

# Work out duplicated coding sequences and label
allflu_cds_df %<>%
  arrange(desc(label), desc(src), title, gid, segment) %>%  # arrange such that when calculating duplicates, zoonotic seqs are preferentially retained
  mutate(cds_dup = case_when(duplicated(string) == TRUE ~ 1,
                             seg_dup == 1 ~ 1,             # if the segment is duplicated, the cds must also be (except in cases of alternative splicing/missing data, but consider those close enough to be worth dropping also)
                             TRUE ~ 0))   %>%
  select(-n) %>%
  add_count(gid)

# # Plot lengths 
# # All retained segment sequences
# allflu_cds_df %>%
#   mutate(group = paste(src, label)) %>%
#   ggplot(aes(x = length, fill = group)) +
#   geom_density(alpha = 0.25) +
#   facet_wrap(~gene, ncol = 2, scales = "free")
# 
# # Coding sequences belonging to complete sets only
# allflu_cds_df %>%
#   filter(gid %in% (allflu_cds_df %>% group_by(gid) %>% count() %>% arrange(-n) %>% filter(n == 12) %>% as.data.frame() %>% pull(gid))) %>%
#   mutate(group = paste(src, label)) %>%
#   ggplot(aes(x = length, fill = group)) +
#   geom_density(alpha = 0.25) +
#   facet_wrap(~gene, ncol = 2, scales = "free")

# Create combined DNAStringSet and remove unfiltered sets from working memory
rm(NCBI_avian_cds, NCBI_human_cds)

#######################################
# Calculate and save features: k-mers #
#######################################

cds %>% calc_kmer_counts(k = 9)



###################################################
# Calculate and save features: genome composition #
###################################################

# Calculate genomic composition for individual genes - restrict to whole genomes

if(load_prev_calcs == TRUE) {
  
  
  allflu_cds_compbias_feats <- list.files(path = "S3\\data\\full\\", pattern = "allflu_cds_compbias_pt_.*\\.rds", full.names = TRUE) %>%
    map_dfr(readRDS)
  
  allflu_wgs_compbias_feats <- readRDS("S3\\data\\full\\allflu_wgs_compbias.rds")
  
  
} else {
  
  cds_temp <- allflu_cds_df %>%
    filter(n >= 8 & wgs_dup == 0) %>%
    arrange(gid, cds_id)
  cds <- cds_temp$string %>% DNAStringSet()
  names(cds) <- cds_temp$cds_id
  rm(cds_temp)
  
  cds_list <- batch(cds, 40000) # Split into batch of 40,000 sequences
  rm(cds)
  
  # Calculate composition counts and biases for individual protein cds (~4 hrs work desktop)
  for(i in (1:length(cds_list))){
    temp_df <- cds_list[[i]] %>%
      calc_composition_counts(codonpairs = TRUE)
    
    temp_df %>% 
      saveRDS(paste0("S3\\data\\full\\allflu_cds_compcounts_pt_",i,".rds"))
    
    temp_df %>%
      calc_composition_bias(codonpairs = FALSE) %>%                          # Do not calculate codon pair biases for individual protein cds
      saveRDS(paste0("S3\\data\\full\\allflu_cds_compbias_pt_",i,".rds"))
    
    rm(temp_df)
  }
  
  rm(cds_list)
  gc()
  
  # Calculate composition biases for whole genome sequences based off of all proteins (~1hrs work desktop)
  list.files(path = "S3\\data\\full\\", pattern = "allflu_cds_compcounts_pt_.*\\.rds", full.names = TRUE) %>%
    map_dfr(readRDS) %>%   
    left_join(allflu_cds_df %>% select(cds_id, gid)) %>%
    group_by(gid) %>%
    summarise_if(is.numeric, sum) %>%
    ungroup %>%
    distinct %>%
    calc_composition_bias(codonpairs = TRUE) %>%                          # Calculate codon pair biases for wgs
    saveRDS("S3\\data\\full\\allflu_wgs_compbias.rds")
  
  rm(allflu_wgs_compcounts)
  
}


################################################################################
# Calculate clustering based on k-mer overlaps to collapse redundant sequences #
################################################################################

# Save whole genome sequences with 8 coding sequences for clustering
wgs_temp <- allflu_wgs_df %>% 
  filter(wgs_dup == 0) %>% 
  inner_join(allflu_cds_df %>% filter(n >= 8) %>% select(gid) %>% distinct()) %>%   # Select only wgs having all 8 cds
  select(gid, wgs_string)
wgs <- wgs_temp$wgs_string %>% DNAStringSet()
names(wgs) <- wgs_temp$gid
writeXStringSet(wgs, filepath = "S3\\data\\full\\nuc\\allflu_nuc_wgs.FASTA")
rm(wgs_temp, wgs)

# Call mmSeq2 14-7e284 to cluster on similarity score
system("E:\\Working\\ai_for_ai\\S3\\scripts\\mmseqs-win64\\mmseqs.bat easy-linclust E:\\Working\\ai_for_ai\\S3\\data\\full\\nuc\\allflu_nuc_wgs.FASTA E:\\Working\\ai_for_ai\\S3\\data\\full\\nuc\\wgs_90_8 tmp --min-seq-id 0.90 -c 0.8 --cov-mode 0",
       intern = TRUE, show.output.on.console = TRUE, ignore.stdout=FALSE, ignore.stderr=FALSE, wait=FALSE)

system("E:\\Working\\ai_for_ai\\S3\\scripts\\mmseqs-win64\\mmseqs.bat easy-linclust E:\\Working\\ai_for_ai\\S3\\data\\full\\nuc\\allflu_nuc_wgs.FASTA E:\\Working\\ai_for_ai\\S3\\data\\full\\nuc\\wgs_80_8 tmp --min-seq-id 0.80 -c 0.8 --cov-mode 0",
       intern = TRUE, show.output.on.console = TRUE, ignore.stdout=FALSE, ignore.stderr=FALSE, wait=FALSE)

# system("H:\\Working\\ai_for_ai\\S3\\scripts\\mmseqs-win64\\mmseqs.bat easy-linclust H:\\Working\\ai_for_ai\\S3\\data\\full\\nuc\\allflu_nuc_wgs.FASTA H:\\Working\\ai_for_ai\\S3\\data\\full\\nuc\\wgs_90_8 tmp --min-seq-id 0.90 -c 0.8 --cov-mode 0",
#        intern = TRUE, show.output.on.console = TRUE, ignore.stdout=FALSE, ignore.stderr=FALSE, wait=FALSE)
# 
# system("H:\\Working\\ai_for_ai\\S3\\scripts\\mmseqs-win64\\mmseqs.bat easy-linclust H:\\Working\\ai_for_ai\\S3\\data\\full\\nuc\\allflu_nuc_wgs.FASTA H:\\Working\\ai_for_ai\\S3\\data\\full\\nuc\\wgs_80_8 tmp --min-seq-id 0.80 -c 0.8 --cov-mode 0",
#        intern = TRUE, show.output.on.console = TRUE, ignore.stdout=FALSE, ignore.stderr=FALSE, wait=FALSE)

# # Test pairwise alignment identity of selected sequences to check clustering
# allflu_wgs_df %>% filter(gid == "A/China/ZMD-22-2/2022O") %>% pull(wgs_string) %>%
#   pairwiseAlignment(allflu_wgs_df %>% filter(gid == "A/Anhui/2/2005H") %>% pull(wgs_string)) %>% pid


##################################################################
# Select clustering run and features to prepare df for ML models #
##################################################################

# Read in clustering from Mmseq2 and select indicative sequences: select random zoonotic if zoonotic in cluster, else select centroid
cluster_ref <- read.table("S3\\data\\full\\nuc\\wgs_90_8_cluster.tsv", sep = '\t', quote = "\"",  encoding="UTF-8", comment.char = '@', header = FALSE)

# Check mixed labels within clusters
cluster_ref %<>% left_join((allflu_wgs_df %>% select(gid, label)), by = c("V2" = "gid")) 
cluster_ref %<>% left_join(cluster_ref %>% group_by(V1) %>% summarise(mix = n_distinct(label)))

cluster_ref %>% group_by(V1) %>% summarise(mix = n_distinct(label)) %>% with(., table(mix)) %>% prop.table()  # Most (98%) of clusters not mixed
cluster_ref %>% group_by(V1) %>% summarise(n = n_distinct(V2)) %>% with(., table(n)) # Clusters range in size from 1 to ~300, though over half are singletons
cluster_ref %>% filter(mix == 2)

# For mixed clusters where cluster representative is not zoonotic, select a random zoonotic representative instead
cluster_ref %<>% mutate(manual_cluster_rep = case_when(
  mix == 2 & !(V1 %in% (cluster_ref %>% filter(mix == 2 & V1 == V2 & label == "zoon") %>% pull(V1))) ~ 1,
  TRUE ~ 0)
)

set.seed(1516)
manual_cluster_df <- cluster_ref %>% filter(manual_cluster_rep == 1 & label == "zoon") %>% group_by(V1) %>% slice_sample(n = 1) %>% rename(cluster_rep = V2) %>% ungroup

cluster_ref %<>% left_join(manual_cluster_df %>% select(V1, cluster_rep), by = "V1") %>% mutate(cluster_rep = coalesce(cluster_rep, V1)) %>% select(cluster_rep, V2, mix)

cluster_ref %>% select(cluster_rep) %>% distinct %>% left_join((allflu_wgs_df %>% select(gid, src, label)), by = c("cluster_rep" = "gid")) %>%  with(., table(src, label)) %>% addmargins # Most (92%; 4165) of the 4536 clusters NCBI (though NCBI was read first?)

#
#
# SELECT FEATURES FOR GIVEN CLUSTER REPRESENTATIVES N' SAVE
#
#
#

##########################################
# Matrix plot of cluster representatives #
##########################################

png("S3\\figures_tables\\heatmap.png", width = 15, height = 9, units = "in", res = 320)
heatmap_cds <- cluster_cds_wgs_df %>%
  as.data.frame() %>%
  set_rownames(.$gid) %>%
  select(gid, label, matches("_Bias")) %>%
  rename_with(., ~ str_replace_all(., c("_Bias" = "", "T" = "U"))) %>%
  select(matches("^[A|C|G|U][A|C|G|U][A|C|G|U]$")) %>%
  janitor::remove_constant() %>%
  as.matrix() %>%
  t() %>%
  heatmap.2(
    density.info = "none",
    trace = "none",
    margins = c(1, 5),
    dendrogram = "col",
    Rowv = "NA",
    labCol = NA,
    ColSideColors = cluster_cds_wgs_df %>% mutate(colsidecol = case_when( # Set side colours using same genus colours as ggplots elsewhere
      label == "nz" ~ "grey80",
      label == "zoon" ~ "gray10"
    )) %>% pull(colsidecol),
    cexRow = 1.1,
    col = colorRampPalette(c("dodgerblue", "gray10", "firebrick2"))(n = 65),
    breaks = c(seq(0, 0.95, length = 28), seq(0.951, 1.05, length = 3), seq(1.051, 3.2, length = 35)),
    keysize = 0.75, key.par = list(cex = 0.5), key.title = NA, key.xlab = "RSCU",
    lhei = c(2, 10), lwid = c(2, 10)
  )
par(lend = 1) # square line ends for the color legend
legend(
  cex = 0.8, x = -0.05, y = 0.9, xpd = TRUE,
  legend = cluster_cds_wgs_df$label %>% unique() %>% sort(), fill = c("grey80", "gray10"), ncol = 1
)
dev.off()




####################
# Checks and scrap #
####################

###### Quick check for possible duplicate sequences and check what unique data GISAID contributes and feasibility of extracting cds myself
test <- bind_rows(
  calc_composition_counts(NCBI_avian_nuc_selected),
  calc_composition_counts(NCBI_human_nuc_selected),
  calc_composition_counts(GISAID_avian_nuc_selected),
  calc_composition_counts(GISAID_human_nuc_selected))

test %<>%
  mutate(dup =ifelse(duplicated(test[,-1]) | duplicated(test[,-1], fromLast = TRUE), 1,0))

test %>% select(dup, fastahead) %>% right_join(allflu_nuc_df) %>% with(., table(dup, paste(label, src)))
test %>% select(dup, fastahead) %>% right_join(allflu_nuc_df) %>% filter(is.na(accession)) %>% with(., table(dup, paste(label, src)))
test %>% select(dup, fastahead) %>% right_join(allflu_nuc_df) %>% filter(!is.na(accession)) %>% with(., table(dup, paste(label, src)))

# Ok so these are the unique ones.. check if they begin/end with a start/stop codon?
zoonlab <- test %>% select(dup, fastahead) %>% right_join(allflu_nuc_df) %>% filter(is.na(accession) & dup == 0 & label == "zoon" & src == "GISAID") %>% pull(fastahead)

nzlab <- test %>% select(dup, fastahead) %>% right_join(allflu_nuc_df) %>% filter(is.na(accession) & dup == 0 & label == "nz" & src == "GISAID") %>% pull(fastahead)

GISAID_avian_nuc[GISAID_avian_nuc_df$fastahead %in% nzlab] %>% subseq(start = 1, end = 3) %>% as.data.frame %>% with(., table(x))
GISAID_avian_nuc[GISAID_avian_nuc_df$fastahead %in% nzlab] %>% subseq(start = -3, end = -1) %>% as.data.frame %>% with(., table(x))

GISAID_human_nuc[GISAID_human_nuc_df$fastahead %in% zoonlab] %>% subseq(start = 1, end = 3) %>% as.data.frame %>% with(., table(x))
GISAID_human_nuc[GISAID_human_nuc_df$fastahead %in% zoonlab] %>% subseq(start = -3, end = -1)  %>% as.data.frame %>% with(., table(x))

# GISAID_avian_nuc %>% subseq(start = 1, end = 3) %>% as.data.frame %>% with(., table(x))
# GISAID_avian_nuc %>% subseq(start = -3, end = -1) %>% as.data.frame %>% with(., table(x))
# 
# GISAID_human_nuc %>% subseq(start = 1, end = 3) %>% as.data.frame %>% with(., table(x))
# GISAID_human_nuc %>% subseq(start = -3, end = -1)  %>% as.data.frame %>% with(., table(x))

#######




# # Filter DNAStringsSet objects to chosen segment sequences (currently ALL segments, regardless of whether complete whole genomes)
# GISAID_avian_nuc_index <- GISAID_avian_nuc_df %>%
#   mutate(index = row_number()) %>%
#   inner_join(allflu_nuc_df %>% select(-any_of(c("INSDC", "UID"))), by = "fastahead") %>%
#   pull(index)
# GISAID_human_nuc_index <- GISAID_human_nuc_df %>%
#   mutate(index = row_number()) %>%
#   inner_join(allflu_nuc_df %>% select(-any_of(c("INSDC", "UID"))), by = "fastahead") %>%
#   pull(index)
# NCBI_avian_nuc_index <- NCBI_avian_nuc_df %>%
#   mutate(index = row_number()) %>%
#   inner_join(allflu_nuc_df %>% select(-any_of(c("INSDC", "UID"))), by = "fastahead") %>%
#   pull(index)
# NCBI_human_nuc_index <- NCBI_human_nuc_df %>%
#   mutate(index = row_number()) %>%
#   inner_join(allflu_nuc_df %>% select(-any_of(c("INSDC", "UID"))), by = "fastahead") %>%
#   pull(index)
# 
# GISAID_avian_nuc_selected <- GISAID_avian_nuc[GISAID_avian_nuc_index]
# GISAID_human_nuc_selected <- GISAID_human_nuc[GISAID_human_nuc_index]
# NCBI_avian_nuc_selected <- NCBI_avian_nuc[NCBI_avian_nuc_index]
# NCBI_human_nuc_selected <- NCBI_human_nuc[NCBI_human_nuc_index]




# list.files("1") %>% paste0("ai_for_ai/data/full/nuc/seg/1/", .) %>% writeLines("reference_list_1_b.txt") # Above code has errors for some reason? Use list.files instead

# 
# for (i in c(1:8)){
#   nuc <- allflu_nuc_df %>% filter(segment == as.character(i)) %>% pull(string) %>% DNAStringSet
#   names(nuc) <- allflu_nuc_df %>% filter(segment == as.character(i)) %>% pull(fastahead)
#   assign(paste0("allflu_nuc_", i), nuc)
#   writeXStringSet(nuc, filepath = paste0("S3\\data\\full\\nuc\\allflu_nuc_fastani_",as.character(i),".FASTA"))
# }
# 
# 
# for (i in c(1:8)){
#   nuc <- allflu_nuc_df %>% filter(segment == as.character(i)) %>% pull(string) %>% DNAStringSet
#   names(nuc) <- allflu_nuc_df %>% filter(segment == as.character(i)) %>% pull(fastahead)
#   assign(paste0("allflu_nuc_", i), nuc)
#   writeXStringSet(nuc, filepath = paste0("S3\\data\\full\\nuc\\allflu_nuc_fastani_",as.character(i),".FASTA"))
# }

## Save each segment of complete whole genomes to one collective FASTA file
# 
# for (i in c(1:8)){
#   nuc <- allflu_nuc_df %>% filter(segment == as.character(i)) %>% pull(string) %>% DNAStringSet
#   names(nuc) <- allflu_nuc_df %>% filter(segment == as.character(i)) %>% pull(fastahead)
#   assign(paste0("allflu_nuc_", i), nuc)
#   writeXStringSet(nuc, filepath = paste0("S3\\data\\full\\nuc\\allflu_nuc_fastani_",as.character(i),".FASTA"))
# }
## Remove duplicates at whole genome level (removes 59,091 individual segments), save each segment of complete whole genomes to individual FASTA files
#allflu_nuc_df %<>% filter(wgs_dup == 0)
#
# for (i in c(1:8)){
#   nuc_df <- allflu_nuc_df %>%
#     filter(gid %in% (allflu_nuc_df %>% group_by(gid) %>% count() %>% arrange(-n) %>% filter(n == 8) %>% pull(gid)))
#   nuc <- nuc_df %>% 
#     filter(segment == as.character(i)) %>% pull(string) %>% DNAStringSet
#   names(nuc) <- nuc_df %>% filter(segment == as.character(i)) %>% pull(fastahead)
#   assign(paste0("allflu_nuc_", i), nuc)
#   dir.create(file.path("S3\\data\\full\\nuc\\wgs\\", as.character(i)))
#   for(j in 1:nrow(nuc_df %>% filter(segment == as.character(i)))){
#     gid <- nuc_df %>% filter(segment == as.character(i)) %>% pull(gid)
#     writeXStringSet(nuc[j], paste0("S3\\data\\full\\nuc\\wgs\\",i,"\\",fs::path_sanitize(gid[j]),".gz"), compress=TRUE)
#     cat(paste0("ai_for_ai/data/full/nuc/wgs",i,"/",fs::path_sanitize(gid[j]),".gz"), file = paste0("S3\\data\\full\\nuc\\wgs\\reference_list_",i,".txt"), sep = "\n", append = TRUE)
#   }
# }
#
## Remove duplicates in segments(removes 45,964 individual segments), save each segment (regardless of whether part of complete wgs) to individual FASTA files
# allflu_nuc_df %<>% filter(seg_dup == 0)
# 
# for (i in c(1:8)){
#   nuc <- allflu_nuc_df %>% 
#     filter(segment == as.character(i)) %>% pull(string) %>% DNAStringSet
#   names(nuc) <- allflu_nuc_df %>% filter(segment == as.character(i)) %>% pull(fastahead)
#   assign(paste0("allflu_nuc_", i), nuc)
#   dir.create(file.path("S3\\data\\full\\nuc\\seg\\", as.character(i)))
#   for(j in 1:nrow(allflu_nuc_df %>% filter(segment == as.character(i)))){
#     gid <- allflu_nuc_df %>% filter(segment == as.character(i)) %>% pull(gid)
#     writeXStringSet(nuc[j], paste0("S3\\data\\full\\nuc\\seg\\",i,"\\",fs::path_sanitize(gid[j]),".gz"), compress=TRUE)
#     cat(paste0("ai_for_ai/data/full/nuc/seg/",i,"/",fs::path_sanitize(gid[j]),".gz"), file = paste0("S3\\data\\full\\nuc\\seg\\reference_list_",i,".txt"), sep = "\n", append = TRUE)
#   }
# }
#
## FastANI read-in, no longer using
#
# for (i in 1:length(unique(allflu_cds_df$gene))){
#   cds <- allflu_cds_df %>% filter(gene == unique(allflu_cds_df$gene)[i]) %>% pull(string) %>% DNAStringSet
#   names(cds) <- allflu_cds_df %>% filter(gene == unique(allflu_cds_df$gene)[i]) %>% pull(fastahead)
#   assign(paste0("allflu_cds_", i), cds)
#   writeXStringSet(cds, filepath = paste0("S3\\data\\full\\cds\\allflu_cds_fastani_",unique(allflu_cds_df$gene)[i],".FASTA"))
# }


