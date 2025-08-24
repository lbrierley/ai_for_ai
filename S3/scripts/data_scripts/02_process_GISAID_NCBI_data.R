######################################################
# Process sequences downloaded from NCBI Flu, GISAID #
######################################################
################################
# Read-in and process metadata #
################################

# Read-in and construct summary data frame, tidy up where known to be mislabelled
# Nucleotide sequences of whole segments
GISAID_avian_nuc_df <- process_GISAID_seq(x = readSet(file = "S3\\data\\full\\GISAID_avian_nuc.fasta"), label = "nz", type = "nuc")
GISAID_human_nuc_df <- process_GISAID_seq(x = readSet(file = "S3\\data\\full\\GISAID_human_nuc.fasta"), label = "zoon", type = "nuc")
NCBI_avian_nuc_df <- process_NCBI_seq(x = readSet(file = "S3\\data\\full\\NCBI_avian_nuc.fasta"), label = "nz", type = "nuc")
NCBI_human_nuc_df <- process_NCBI_seq(x = readSet(file = "S3\\data\\full\\NCBI_human_nuc.fasta"), label = "zoon", type = "nuc")

# Coding sequences only
NCBI_avian_cds_df <- process_NCBI_seq(x = readSet(file = "S3\\data\\full\\NCBI_avian_cds.fasta"), label = "nz", type = "cds")
NCBI_human_cds_df <- process_NCBI_seq(x = readSet(file = "S3\\data\\full\\NCBI_human_cds.fasta"), label = "zoon", type = "cds")

# Adjust where FASTA headers have pulled wrong information directly
GISAID_avian_nuc_fix_ref <- readSet(file = "S3\\data\\full\\GISAID_avian_nuc.fasta")[!(GISAID_avian_nuc_df$gene %in% c("PB1", "PB2", "PA", "HA", "NP", "NA", "MP", "NS")), ] %>%
  names() %>%
  as.data.frame() %>%
  setNames("title") %>%
  tidyr::separate(title, sep = "\\|", into = c("nn1", "nn9", "nn2", "nn3", "nn4", "INSDC", "nn6", "nn7", "nn8", "gene", "segment"), extra = "drop") %>%
  select(-contains("nn"))

GISAID_avian_nuc_df %<>%
  rows_update(., GISAID_avian_nuc_fix_ref, by = "INSDC")

# GISAID metadata - reshape, as this comes in long format data specifying IDs for each respective protein of each sequence
meta_ref <- bind_rows(
  readxl::read_excel("S3/data/full/GISAID_avian_meta.xls", guess_max = 1048576) %>%
    select(1:9) %>%
    as.data.frame %>%
    mutate_at(vars(matches("Id$")), ~ gsub("\\|.*$", "", .)) %>%
    reshape2::melt(id.vars = "Isolate_Id") %>%
    mutate(variable = gsub(" Segment_Id", "", variable)),
  readxl::read_excel("S3/data/full/GISAID_human_meta.xls", guess_max = 1048576) %>%
    select(1:9) %>%
    as.data.frame %>%
    mutate_at(vars(matches("Id$")), ~ gsub("\\|.*$", "", .)) %>%
    reshape2::melt(id.vars = "Isolate_Id") %>%
    mutate(variable = gsub(" Segment_Id", "", variable))
)

# Clean names of protein sequence FASTA files and resave for compatibility with iFeatureOmega package
files <- list.files(path = "S3\\data\\full\\prot\\", pattern = "prot.FASTA", full.names = TRUE)
for(i in 1:length(files)){
 prot_fasta_name_clean(files[i])
}

########################################################################
# Filter and process whole genome sequences, i.e. full-length segments #
########################################################################

# Combine all data in nucleotide sequences of segments
allflu_nuc_df <- bind_rows(GISAID_human_nuc_df, 
                           GISAID_avian_nuc_df, 
                           NCBI_human_nuc_df, 
                           NCBI_avian_nuc_df) %>%
  
# Remove mixed infections (usually signified by 'MIXED' in NCBI data but also slashes, periods, commas and by H0N0 in GISAID data but some labelled as such in title), (removes 14,733 individual segment sequences)
  filter(!grepl("MIXED|\\,|\\.|\\/", subtype)) %>%
  filter(subtype != "H0N0" & !(grepl("Mixed", title))) %>%
  
  # Remove sequences known a priori to be mislabelled (removes 24 individual segment sequences, 3 wgs)
  filter(title != "A/goose/Chiayi/18020014-1/2018" & title != "H5N6")

# Filter canonical sequences belonging to each whole genome according to GISAID metadata (removes 575 individual segment sequences)
allflu_nuc_df %<>% 
  filter(is.na(INSDC) | INSDC %in% gsub("EPI", "", meta_ref$value))

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

rm(temp_8s)

# For NCBI data, where there is > 8 segments and not a multiple of 8, consider the first instance of each segment as the canonical one of the wgs (removes 339 individual segment sequences)
allflu_nuc_df %<>% 
  group_by(gid, segment) %>%
  slice(1) %<>%
  ungroup()

# Filter based on predefined lower bounds of segment lengths (removes 46 individual segment sequences and 38 wgs)
allflu_nuc_df %<>% filter(segment == 1 & length >= 2100 |
                            segment == 2 & length >= 2200 |
                            segment == 3 & length >= 2000 |
                            segment == 4 & length >= 1650 |
                            segment == 5 & length >= 1450 |
                            segment == 6 & length >= 1300 |
                            segment == 7 & length >= 900 |
                            segment == 8 & length >= 800) 

# Identify duplicated whole genome sequences and duplicated individual segment sequences
allflu_wgs_df <- allflu_nuc_df %>% 
  arrange(gid, segment) %>% # arrange segments in order
  mutate(date = date %>% gsub("--", "-06-", .) %>% # assume midpoint for missing months
           gsub("-$", "-15", .) %>% # assume midpoint for missing days
           as.Date(format = "%Y-%m-%d")) %>%
  group_by(src, label, subtype, title, gid) %>%
  summarise(n = n(), wgs_string = paste(string, collapse=""), date=as.Date(min(date))) %>%
  ungroup() %>%
  arrange(desc(label), desc(src), title, gid) %>%  # arrange such that when calculating duplicates, zoonotic seqs and NCBI seqs are preferentially retained
  mutate(wgs_dup = if_else(duplicated(wgs_string) == TRUE, 1, 0))

allflu_nuc_df %<>% 
  select(-date) %>%
  left_join(allflu_wgs_df %>% select(-wgs_string)) %>%
  arrange(desc(label), desc(src), title, gid, segment) %>%  # arrange such that when calculating duplicates, zoonotic seqs and NCBI seqs are preferentially retained
  mutate(seg_dup = case_when(duplicated(string) == TRUE ~ 1, 
                             duplicated(accession) == TRUE ~ 1,
                             TRUE ~ 0))

######################################################
# Process gene sequences, i.e. only coding sequences #
######################################################

# Identify GISAID coding sequences, defined as longest ORFs within segment sequences
GISAID_avian_orfs <- allflu_nuc_df %>% filter(src == "GISAID" & label == "nz") %>% pull(string) %>% lapply(bind_ORF) %>% bind_rows()
GISAID_human_orfs <- allflu_nuc_df %>% filter(src == "GISAID" & label == "zoon") %>% pull(string) %>% lapply(bind_ORF) %>% bind_rows()

write.csv(GISAID_avian_orfs, "S3\\data\\full\\GISAID_avian_orfs.csv")
write.csv(GISAID_human_orfs, "S3\\data\\full\\GISAID_human_orfs.csv")

## (as this step takes a while, optionally read in previously saved ORFs)
# GISAID_avian_orfs <- read.csv("S3\\data\\full\\GISAID_avian_orfs.csv")
# GISAID_human_orfs <- read.csv("S3\\data\\full\\GISAID_human_orfs.csv")

GISAID_avian_cds_df <- allflu_nuc_df %>% filter(src == "GISAID" & label == "nz") %>% bind_cols(GISAID_avian_orfs) %>% rowwise %>% mutate(string = substr(string, start, end)) %>% as.data.frame
GISAID_human_cds_df <- allflu_nuc_df %>% filter(src == "GISAID" & label == "zoon") %>% bind_cols(GISAID_human_orfs) %>% rowwise %>% mutate(string = substr(string, start, end)) %>% as.data.frame


# Combine all data and filter out duplicates in coding sequences
allflu_cds_df <- bind_rows(GISAID_human_cds_df, 
                           GISAID_avian_cds_df, 
                           NCBI_human_cds_df %>% 
                             select(-gid) %>% 
                             mutate(date = date %>% gsub("--", "-06-", .) %>% gsub("-$", "-15", .) %>% as.Date(format = "%Y-%m-%d")) %>%  # fix date errors
                             inner_join(allflu_nuc_df %>%           # Overwrite gids with updated ones from allflu_nuc_df, which accounts for multiple sets of 8 within same title
                                          filter(src == "NCBI") %>%
                                          select(gid, accession, wgs_dup, seg_dup), by = "accession"),
                           NCBI_avian_cds_df %>%
                             select(-gid) %>% 
                             mutate(date = date %>% gsub("--", "-06-", .) %>% gsub("-$", "-15", .) %>% as.Date(format = "%Y-%m-%d")) %>% # fix date errors
                             inner_join(allflu_nuc_df %>% 
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
  
  # Remove mixed infections (usually signified by 'MIXED' in NCBI data but also slashes, periods, commas and by "H0N0" in GISAID data but some labelled as such in title), (removes 19,187 individual coding sequences)
  filter(!grepl("MIXED|\\,|\\.|\\/", subtype)) %>%
  filter(subtype != "H0N0" & !(grepl("Mixed", title))) %>%
  
  # Remove sequences known a priori to be mislabelled (removes 34 individual coding sequences)
  filter(title != "A/goose/Chiayi/18020014-1/2018" & title != "H5N6")

# Relabel GISAID genes to reflect the ORF captured (M1 for MP, NS1 for NS), filter out NCBI genes not present in GISAID (M2, NS2, PA-X, PB1-F2)
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

# Filter based on predefined lower bounds of segment lengths (removes 7 individual coding sequences)
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

# Identify duplicated coding sequences and label
allflu_cds_df %<>%
  arrange(desc(label), desc(src), title, gid, segment) %>%  # arrange such that when calculating duplicates, zoonotic seqs are preferentially retained
  mutate(cds_dup = case_when(duplicated(string) == TRUE ~ 1,
                             seg_dup == 1 ~ 1,             # if the segment is duplicated, the cds must also be (except in cases of alternative splicing/missing data, but consider those close enough to be worth dropping also)
                             TRUE ~ 0))   %>%
  select(-n) %>%
  add_count(gid)

#############################
# Process protein sequences #
#############################

# Combine all data - protein sequences
# For NCBI, use pre-made gids (from nucleotide sequences) to assign which prot sequences belong to which wgs
NCBI_avian_prot_df <- process_NCBI_seq(x = readAAStringSet(file = "S3\\data\\full\\prot\\NCBI_avian_prot.fasta"), label = "nz", type = "prot") %>% 
  select(-string, -gid) %>%
  left_join(allflu_nuc_df %>% filter(src == "NCBI") %>% select(gid, accession), by = "accession") 
NCBI_human_prot_df <- process_NCBI_seq(x = readAAStringSet(file = "S3\\data\\full\\prot\\NCBI_human_prot.fasta"), label = "zoon", type = "prot") %>% 
  select(-string, -gid) %>%
  left_join(allflu_nuc_df %>% filter(src == "NCBI") %>% select(gid, accession), by = "accession") 

# For GISAID, use metadata to assign which prot sequences belong to which wgs
GISAID_avian_prot_df <- process_GISAID_seq(x = readAAStringSet(file = "S3\\data\\full\\prot\\GISAID_avian_prot.fasta"), label = "nz", type = "prot")
GISAID_human_prot_df <- process_GISAID_seq(x = readAAStringSet(file = "S3\\data\\full\\prot\\GISAID_human_prot.fasta"), label = "zoon", type = "prot")