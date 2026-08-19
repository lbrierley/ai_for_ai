###############################################
# Process sequences from NCBI Genbank, GISAID #
###############################################
################################
# Read-in and process metadata #
################################

# Read-in and construct summary df, tidy up where known to be mislabelled
# Nucleotide sequences of segments
nuc_results <- process_segmented_sequences(
  path       = "//filepor10/DOP$/X037_AVI_GeneticMarkers/S3/data/Segmented/Raw_Sequences/nuc",
  ending     = "_Nuc",
  type_label = "nuc"
)

GISAID_avian_nuc_df <- nuc_results$GISAID_avian_df
GISAID_human_nuc_df <- nuc_results$GISAID_human_df
NCBI_avian_nuc_df   <- nuc_results$NCBI_avian_df
NCBI_human_nuc_df   <- nuc_results$NCBI_human_df

# Coding sequences of genes
cds_results <- process_segmented_sequences(
  path       = "//filepor10/DOP$/X037_AVI_GeneticMarkers/S3/data/Segmented/Raw_Sequences/cds",
  ending     = "_Cds",
  type_label = "cds"
)

NCBI_avian_cds_df <- cds_results$NCBI_avian_df
NCBI_human_cds_df <- cds_results$NCBI_human_df

# Adjust where FASTA headers have pulled wrong information directly
GISAID_avian_nuc_fix_ref <- names(
  nuc_results$raw_GISAID_avian[
    !(GISAID_avian_nuc_df$gene %in% c("PB1","PB2","PA","HA","NP","NA","MP","NS"))]) %>%
  as.data.frame() %>%
  setNames("fastahead") %>%
  mutate(fastahead_keep=fastahead) %>%
  tidyr::separate(
    fastahead, sep = "\\|",
    into = c("nn1","nn9","nn2","nn3","nn4","INSDC","nn6","nn7","nn8","gene","segment"),
    extra = "drop") %>%
  mutate(fastahead=fastahead_keep) %>%
  select(-contains("nn"), -fastahead_keep)

GISAID_avian_nuc_df %<>%
  rows_update(., GISAID_avian_nuc_fix_ref, by = "fastahead")

# Create combined segment metadata files 
metadata_dir <- "//filepor10/DOP$/X037_AVI_GeneticMarkers/S3/data/Segmented/Raw_Sequences/metadata"
meta_results <- build_metadata_excels(metadata_dir)

# GISAID metadata - in long format specifying IDs for each respective protein of each sequence
meta_ref <- bind_rows(
  readxl::read_excel("//filepor10/DOP$/X037_AVI_GeneticMarkers/S3/data/Segmented/Raw_Sequences/metadata/GISAID_Avian_Metadata.xlsx", guess_max = 1048576) %>%
    select(1:9) %>%
    as.data.frame %>%
    mutate_at(vars(matches("Id$")), ~ gsub("\\|.*$", "", .)) %>%
    reshape2::melt(id.vars = "Isolate_Id") %>%
    mutate(variable = gsub(" Segment_Id", "", variable)),
  readxl::read_excel("//filepor10/DOP$/X037_AVI_GeneticMarkers/S3/data/Segmented/Raw_Sequences/metadata/GISAID_Human_Metadata.xlsx", guess_max = 1048576) %>%
    select(1:9) %>%
    as.data.frame %>%
    mutate_at(vars(matches("Id$")), ~ gsub("\\|.*$", "", .)) %>%
    reshape2::melt(id.vars = "Isolate_Id") %>%
    mutate(variable = gsub(" Segment_Id", "", variable))
)


# # Clean names of protein sequence FASTA files and resave for iFeatureOmega
files <- list.files(path = "//filepor10/DOP$/X037_AVI_GeneticMarkers/S3/data/Segmented/Raw_Sequences/prot",
                    pattern = "_Pro\\.(fa|fasta)$",
                    full.names = TRUE)
for(i in 1:length(files)){
  prot_fasta_name_clean(files[i])
}

############################################################
# Filter and process complete segment nucleotide sequences #
############################################################

# Combine all data - nucleotide sequences of segments
allflu_nuc_df <- bind_rows(GISAID_human_nuc_df, 
                           GISAID_avian_nuc_df, 
                           NCBI_human_nuc_df, 
                           NCBI_avian_nuc_df) %>%
  
  # Remove mixed infections (usually signified by 'MIXED' in NCBI data but also slashes, periods, commas and by H0N0 in GISAID data but some labelled as such in title), (removes 14,733 individual segment sequences)
  filter(!grepl("MIXED|\\,|\\.|\\/", subtype)) %>%
  filter(subtype != "H0N0" & !(grepl("Mixed", title))) %>%
  
  # Remove sequences known a priori to be mislabelled (removes 24 individual segment sequences, 3 wgs)
  filter(title != "A/goose/Chiayi/18020014-1/2018" & title != "H5N6") %>%
  
  # Ensure subtype is not blank or missing
  filter(!is.na(subtype), str_trim(subtype) != "") %>%
  
  # Keep only subtypes that contain BOTH H and N
  filter(!grepl("^H0N[1-9][0-9]*$|^H[1-9][0-9]*N0$|^H0N0$", subtype, ignore.case = TRUE)) %>%
  filter(grepl("H[0-9]+N[0-9]+", subtype, ignore.case = TRUE))

# Fixing where gene labels are missing or switched with segment
allflu_nuc_df <- allflu_nuc_df %>%
  mutate(
    gene = if_else(is.na(gene), NA_character_, as.character(gene)),
    segment = if_else(is.na(segment), NA_character_, as.character(segment)),
    gene = str_trim(gene),
    segment = str_trim(segment),
    
    # Treat missing indicators as missing
    gene = case_when(
      gene %in% c("", ".", "unknown", "Unknown") ~ NA_character_,
      TRUE ~ gene
    ),
    segment = case_when(
      segment %in% c("", ".", "unknown", "Unknown") ~ NA_character_,
      TRUE ~ segment
    ),
    
    # Detect where the segment number is in gene column
    gene_is_segment = !is.na(gene) & str_detect(gene, "^[1-8]$")
  ) %>%
  mutate(
    segment = if_else(gene_is_segment & !segment %in% as.character(1:8),
                      gene, segment),
    gene = if_else(gene_is_segment, NA_character_, gene)
  ) %>%
  select(-gene_is_segment) %>%
  
  # Infer gene from segment when gene is missing
  mutate(
    gene = case_when(
      is.na(gene) & segment == "1" ~ "PB2",
      is.na(gene) & segment == "2" ~ "PB1",
      is.na(gene) & segment == "3" ~ "PA",
      is.na(gene) & segment == "4" ~ "HA",
      is.na(gene) & segment == "5" ~ "NP",
      is.na(gene) & segment == "6" ~ "NA",
      is.na(gene) & segment == "7" ~ "MP",
      is.na(gene) & segment == "8" ~ "NS",
      TRUE ~ gene
    )
  ) %>%
  
  # Remove broken entries
  filter(!is.na(gene) & !is.na(segment))

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

#  # Retain canonical sequences belonging to each whole genome according to GISAID metadata (and keep all NCBI) (removes 575 individual segment sequences)
#  allflu_nuc_df %<>% 
#    filter(is.na(INSDC) | INSDC %in% gsub("EPI", "", meta_ref$value))

#  # Assign whole genomes (all segments)
#  # Where there is a multiple of 8 segments, assign each run as an individual wgs
#  ids <- allflu_nuc_df %>%
#    group_by(gid) %>%
#    summarise(seqs = n()) %>%
#    filter(seqs%%8 == 0 & seqs > 8) %>% pull(gid)

#  temp_8s <- allflu_nuc_df %>% filter(gid %in% ids) %>% group_by(gid,segment) 
#  temp_8s <- bind_rows(temp_8s %>% slice(1), # Arrange such that we get segment 1-8 in each instance, regardless of original row order. 
#                       temp_8s %>% slice(2),
#                       temp_8s %>% slice(3),
#                       temp_8s %>% slice(4)) # Max seqs per gid = 32, so 4 sets
#  temp_8s %<>% 
#    group_by(gid) %>% 
#    mutate(tempalpha = rep(LETTERS, each=8, length=n())) %>%
#    mutate(gid = paste0(gid, "_", tempalpha)) %>%
#    select(-tempalpha) %>%
#    ungroup()

#  allflu_nuc_df %<>% 
#    filter(!(gid %in% ids)) %>%
#    bind_rows(temp_8s)

#  rm(temp_8s)

#  # For NCBI data, where there is > 8 segments and not a multiple of 8, take the first instance of each segment as the wgs (removes 339 individual segment sequences)
#  allflu_nuc_df %<>% 
#    group_by(gid, segment) %>%
#    slice(1) %<>%
#    ungroup()

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

# allflu_wgs_df %>%
#   mutate(date = date %>% gsub("--", "-06-", .) %>% # assume midpoint for missing months
#            gsub("-$", "-15", .) %>% # assume midpoint for missing days
#            as.Date(format = "%Y-%m-%d")) %>%
#   filter(date > as.Date("2016-12-31")) %>%
#   nrow()

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

######################################################
# Process gene sequences, i.e. only coding sequences #
######################################################

# #Extract Open Reading Frames (ORFs) for GISAID avian
# GISAID_avian_df <- allflu_nuc_df %>%
# filter(src == "GISAID", label == "nz")

# GISAID_avian_orfs <- GISAID_avian_df$string %>%
# lapply(bind_ORF) %>%
# bind_rows()

# #Extract ORFs for GISAID human
# GISAID_human_df <- allflu_nuc_df %>%
#  filter(src == "GISAID", label == "zoon")

# GISAID_human_orfs <- GISAID_human_df$string %>%
#  lapply(bind_ORF) %>%
#  bind_rows()

# #Build CDS sequences by slicing the ORF out of the nucleotide string
# GISAID_avian_cds_df <- GISAID_avian_df %>%
#  bind_cols(GISAID_avian_orfs) %>%
#  rowwise() %>%
#  mutate(string = substr(string, start, end)) %>%
#  ungroup()

#  GISAID_human_cds_df <- GISAID_human_df %>%
#  bind_cols(GISAID_human_orfs) %>%
#  rowwise() %>%
#  mutate(string = substr(string, start, end)) %>%
#  ungroup()

#Load in saved files
GISAID_avian_orfs <- read.csv("//filepor10/DOP$/X037_AVI_GeneticMarkers/S3/data/Segmented/Raw_Sequences/GISAID_avian_orfs.csv")
GISAID_human_orfs <- read.csv("//filepor10/DOP$/X037_AVI_GeneticMarkers/S3/data/Segmented/Raw_Sequences/GISAID_human_orfs.csv")
GISAID_avian_cds_df <- read.csv("//filepor10/DOP$/X037_AVI_GeneticMarkers/S3/data/Segmented/Raw_Sequences/GISAID_avian_cds.csv", na.strings = "")
GISAID_human_cds_df <- read.csv("//filepor10/DOP$/X037_AVI_GeneticMarkers/S3/data/Segmented/Raw_Sequences/GISAID_human_cds.csv", na.strings = "")

# Making segment characters
GISAID_human_cds_df$segment <- as.character(GISAID_human_cds_df$segment)
GISAID_avian_cds_df$segment <- as.character(GISAID_avian_cds_df$segment)
NCBI_human_cds_df$segment  <- as.character(NCBI_human_cds_df$segment)
NCBI_avian_cds_df$segment  <- as.character(NCBI_avian_cds_df$segment)

# Making dates characters
GISAID_human_cds_df$date <- as.character(GISAID_human_cds_df$date)
GISAID_avian_cds_df$date <- as.character(GISAID_avian_cds_df$date)
NCBI_human_cds_df$date   <- as.character(NCBI_human_cds_df$date)
NCBI_avian_cds_df$date   <- as.character(NCBI_avian_cds_df$date)

# Ensure NCBI CDS frames have wgs_dup and seg_dup columns
if (!"wgs_dup" %in% names(NCBI_human_cds_df)) NCBI_human_cds_df$wgs_dup <- NA_integer_
if (!"seg_dup" %in% names(NCBI_human_cds_df)) NCBI_human_cds_df$seg_dup <- NA_integer_

if (!"wgs_dup" %in% names(NCBI_avian_cds_df)) NCBI_avian_cds_df$wgs_dup <- NA_integer_
if (!"seg_dup" %in% names(NCBI_avian_cds_df)) NCBI_avian_cds_df$seg_dup <- NA_integer_


# Build mapping from nucleotide dataframe (NCBI only)
nuc_map <- allflu_nuc_df %>%
  filter(src == "NCBI") %>%
  mutate(accession = as.character(accession) %>% str_trim()) %>%
  group_by(accession) %>%
  summarise(
    gid_nuc = dplyr::first(gid[!is.na(gid)]),
    wgs_dup_nuc = dplyr::first(wgs_dup),
    seg_dup_nuc = dplyr::first(seg_dup),
    .groups = "drop"
  )

# Combine all data and remove duplicates - coding sequences
allflu_cds_df <- bind_rows(GISAID_human_cds_df, 
                           GISAID_avian_cds_df, 
                           
                           NCBI_human_cds_df %>% 
                             
                             # NEW: recover accession from fastahead if missing
                             mutate(
                               accession = as.character(accession) %>% str_trim(),
                               accession = if_else(
                                 is.na(accession),
                                 str_extract(fastahead, "(?<=cds:)[A-Za-z0-9]+"),
                                 accession
                               ),
                               date = date %>% gsub("--", "-06-", .) %>% gsub("-$", "-15", .)
                             ) %>%
                             
                             # NEW: join to recover gid (no select(-gid))
                             left_join(nuc_map, by = "accession") %>%
                             
                             mutate(
                               gid = coalesce(as.character(gid), gid_nuc),
                               wgs_dup = coalesce(wgs_dup, wgs_dup_nuc),
                               seg_dup = coalesce(seg_dup, seg_dup_nuc)
                             ) %>%
                             mutate(gid = if_else(src == "NCBI" & str_detect(gid, "NA$"), paste0(title, str_sub(accession, 1, 1)),gid)) %>%
                             select(-gid_nuc, -wgs_dup_nuc, -seg_dup_nuc)
                           
                           ,
                           
                           NCBI_avian_cds_df %>%
                             
                             # NEW: recover accession from fastahead if missing
                             mutate(
                               accession = as.character(accession) %>% str_trim(),
                               accession = if_else(
                                 is.na(accession),
                                 str_extract(fastahead, "(?<=cds:)[A-Za-z0-9]+"),
                                 accession
                               ),
                               date = date %>% gsub("--", "-06-", .) %>% gsub("-$", "-15", .)
                             ) %>%
                             
                             # NEW: join to recover gid
                             left_join(nuc_map, by = "accession") %>%
                             
                             mutate(
                               gid = coalesce(as.character(gid), gid_nuc),
                               wgs_dup = coalesce(wgs_dup, wgs_dup_nuc),
                               seg_dup = coalesce(seg_dup, seg_dup_nuc)
                             ) %>%
                             mutate(gid = if_else(src == "NCBI" & str_detect(gid, "NA$"),paste0(title, str_sub(accession, 1, 1)),gid)) %>%
                             select(-gid_nuc, -wgs_dup_nuc, -seg_dup_nuc)
                           
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
  filter(title != "A/goose/Chiayi/18020014-1/2018" & title != "H5N6") %>%
  
  #ensure subtype is not blank or missing
  filter(!is.na(subtype), str_trim(subtype) != "") %>%
  
  #keep only subtypes that contain BOTH H and N, and remove H0 or N0's
  filter(!grepl("^H0N[1-9][0-9]*$|^H[1-9][0-9]*N0$|^H0N0$", subtype, ignore.case = TRUE)) %>%
  filter(grepl("H[0-9]+N[0-9]+", subtype, ignore.case = TRUE))

# Relabel GISAID genes to reflect likely ORF captured (M1 for MP, NS1 for NS), filter out NCBI genes not present in GISAID (M2, NS2, PA-X, PB1-F2)
allflu_cds_df %<>% 
  mutate(gene = case_when(
    gene == "MP" ~ "M1",
    gene == "NS" ~ "NS1",
    TRUE ~ gene
  )) %>%
  filter(!(gene %in% c("M2", "NS2", "PA-X", "PB1-F2"))) %>%
  select(-start, -end)

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
  arrange(desc(label), desc(src), title, gid, segment) %>%
  mutate(cds_dup = case_when(
    duplicated(string) == TRUE ~ 1,
    seg_dup == 1 ~ 1,
    TRUE ~ 0
  ))   %>%
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

#############################
# Process protein sequences #
#############################

#Read in FASTA files and combined segment results 
prot_results <- process_segmented_sequences_pro(
  path       = "//filepor10/DOP$/X037_AVI_GeneticMarkers/S3/data/Segmented/Raw_Sequences/prot",
  ending     = "_Pro",
  type_label = "prot"
)

GISAID_avian_prot_df <- prot_results$GISAID_avian_df %>%
  left_join(raw_to_lookup(prot_results$raw_GISAID_avian), by = "fastahead") %>%
  mutate(
    string = raw_string,
    accession = na_if(protaccession, "")
  ) %>%
  select(-raw_string)

GISAID_human_prot_df <- prot_results$GISAID_human_df %>%
  left_join(raw_to_lookup(prot_results$raw_GISAID_human), by = "fastahead") %>%
  mutate(
    string = raw_string,
    accession = na_if(protaccession, "")
  ) %>%
  select(-raw_string)

NCBI_avian_prot_df <- prot_results$NCBI_avian_df %>%
  mutate(fastahead = gsub(" ", "_", fastahead)) %>%
  left_join(extract_ncbi_accession(prot_results$raw_NCBI_avian), by = "fastahead") %>%
  mutate(
    accession = coalesce(na_if(prot_id, ""), raw_accession, accession),
    gid = paste0(title, str_sub(accession, 1, 1))  
  ) %>%
  select(-raw_accession)

NCBI_human_prot_df <- prot_results$NCBI_human_df %>%
  mutate(fastahead = gsub(" ", "_", fastahead)) %>%
  left_join(extract_ncbi_accession(prot_results$raw_NCBI_human), by = "fastahead") %>%
  mutate(
    accession = coalesce(na_if(prot_id, ""), raw_accession, accession),
    gid = paste0(title, str_sub(accession, 1, 1))  
  ) %>%
  select(-raw_accession)

#############################################
# Save outputs as .csv and processed fastas #
#############################################
#######################################################
# Save segment cds data as fastas for clustering #
#######################################################
for (seg in 1:8){
  seq_df <- allflu_cds_df %>% filter(segment == seg) %>% as.data.frame
  seq <- seq_df$string %>% DNAStringSet()
  names(seq) <- seq_df$gid #changed
  writeXStringSet(seq, filepath = paste0("//filepor10/DOP$/X037_AVI_GeneticMarkers/S3/data/Segmented/Clustering/allflu_cds_seg_", seg, "_full.FASTA"))
}

###########################
# Creating protein fastas #
###########################

# Output directory
output_dir <- "//filepor10/DOP$/X037_AVI_GeneticMarkers/S3/data/Segmented/Raw_Sequences/protein_fastas"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Listing protein dataframes 
prot_dfs <- list(
  GISAID_avian_prot_df = GISAID_avian_prot_df %>% filter(gene %in% c("PB2", "PB1", "PA", "HA", "NP", "NA", "M1", "NS1")),
  GISAID_human_prot_df  = GISAID_human_prot_df %>% filter(gene %in% c("PB2", "PB1", "PA", "HA", "NP", "NA", "M1", "NS1")),
  NCBI_avian_prot_df    = NCBI_avian_prot_df %>% filter(gene %in% c("PB2", "PB1", "PA", "HA", "NP", "NA", "M1", "NS1")),
  NCBI_human_prot_df    = NCBI_human_prot_df %>% filter(gene %in% c("PB2", "PB1", "PA", "HA", "NP", "NA", "M1", "NS1"))
)

# Write all FASTA files (should create 32 FASTA files)
purrr::iwalk(prot_dfs, ~ write_gene_fastas(.x, .y, output_dir))

# Check for illegal characters
files <- list.files(path = output_dir,
                    pattern = ".(fa|fasta)$",
                    full.names = TRUE)
lapply(files, function(x) x %>% 
         readAAStringSet() %>% 
         alphabetFrequency %>% 
         colSums() %>% 
         t() %>% 
         as.data.frame) %>% 
  bind_rows %>% 
  set_rownames(list.files(path = output_dir, pattern = ".(fa|fasta)$", full.names = FALSE))

######################
# Write main outputs #
######################

# Write ORFs
write.csv(GISAID_avian_orfs, "//filepor10/DOP$/X037_AVI_GeneticMarkers/S3/data/Segmented/Raw_Sequences/GISAID_avian_orfs.csv")
write.csv(GISAID_human_orfs, "//filepor10/DOP$/X037_AVI_GeneticMarkers/S3/data/Segmented/Raw_Sequences/GISAID_human_orfs.csv")

# Write csv files
write.csv(allflu_nuc_df,"//filepor10/DOP$/X037_AVI_GeneticMarkers/S3/data/Segmented/Raw_Sequences/allflu_nuc_df.csv",row.names = FALSE)

write.csv(GISAID_avian_cds_df,"//filepor10/DOP$/X037_AVI_GeneticMarkers/S3/data/Segmented/Raw_Sequences/GISAID_avian_cds.csv",row.names = FALSE)
write.csv(GISAID_human_cds_df,"//filepor10/DOP$/X037_AVI_GeneticMarkers/S3/data/Segmented/Raw_Sequences/GISAID_human_cds.csv",row.names = FALSE)
write.csv(allflu_cds_df, "//filepor10/DOP$/X037_AVI_GeneticMarkers/S3/data/Segmented/Raw_Sequences/allflu_cds_df.csv",row.names = FALSE)

write.csv(GISAID_avian_prot_df,"//filepor10/DOP$/X037_AVI_GeneticMarkers/S3/data/Segmented/Raw_Sequences/GISAID_avian_prot_df.csv",row.names = FALSE)
write.csv(NCBI_avian_prot_df,"//filepor10/DOP$/X037_AVI_GeneticMarkers/S3/data/Segmented/Raw_Sequences/NCBI_avian_prot_df.csv",row.names = FALSE)
write.csv(GISAID_human_prot_df,"//filepor10/DOP$/X037_AVI_GeneticMarkers/S3/data/Segmented/Raw_Sequences/GISAID_human_prot_df.csv",row.names = FALSE)
write.csv(NCBI_human_prot_df,"//filepor10/DOP$/X037_AVI_GeneticMarkers/S3/data/Segmented/Raw_Sequences/NCBI_human_prot_df.csv",row.names = FALSE)


#######################################################
# Calculate and process machine learning feature sets #
#######################################################
#######################################
# Calculate and save features: k-mers #
#######################################

if(load_prev_calcs == FALSE) {
  
  # Calculate k-mers of entire segment nucleotide sequences
  
  # Set parallelisation
  cl <- makePSOCKcluster(12)
  registerDoParallel(cl)
  clusterSetRNGStream(cl, 1429)
  
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
        saveRDS(paste0("//filepor10/DOP$/X037_AVI_GeneticMarkers/S3/data/Segmented/mlready/allflu_nuc_",k_val,"mer_pt_",c("PB2", "PB1", "PA", "HA", "NP", "NA", "M1", "NS1")[j],".rds"))
    }
  
  stopCluster(cl)
}


###################################################
# Calculate and save features: genome composition #
###################################################

# Calculate genome composition for individual gene coding sequences

if(load_prev_calcs == TRUE) {
  
  allflu_cds_compbias_feats <- list.files(path = "S3\\data\\Segmented\\cds", pattern = "allflu_cds_compbias_pt_.*\\.rds", full.names = TRUE) %>%
    map_dfr(readRDS)
 
  
  
} else {
  
  # Calculate composition counts and biases for individual protein cds
  
  genes <- c("PB2", "PB1", "PA", "HA", "NP", "NA", "M1", "NS1")
  
  # Set parallelisation
  cl <- makePSOCKcluster(8)
  registerDoParallel(cl)
  clusterSetRNGStream(cl, 1429)
  
  foreach(i = seq_along(genes),
          .packages = c("Biostrings",
                        "coRdon",
                        "dplyr",
                        "stringr",
                        "magrittr")
  ) %dopar% {
    
    source("S3\\scripts\\functions.R")
    
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
      saveRDS(paste0( "//filepor10/DOP$/X037_AVI_GeneticMarkers/S3/data/Segmented/cds/allflu_cds_compcounts_pt_",gene_label,".rds"))
    
    bias_raw <- temp_df %>%
      calc_composition_bias(codonpairs = FALSE)   #Calling function
    
    bias_cols <- colnames(bias_raw)[grepl("_Bias$", colnames(bias_raw))]
    bias_df   <- bind_cols(gid = gids, bias_raw[, bias_cols, drop = FALSE])
    
    bias_df %>%
      saveRDS(paste0("//filepor10/DOP$/X037_AVI_GeneticMarkers/S3/data/Segmented/cds/allflu_cds_compbias_pt_", gene_label,".rds"))
    
    rm(temp_df, bias_raw, bias_df, df_sub)
    gc()
  }
  
  stopCluster(cl)
  
  # # Calculate composition biases for whole genome sequences based off of all proteins (~1hrs work desktop) - limiting to complete whole genomes only
  # list.files(path = "//filepor10/DOP$/X037_AVI_GeneticMarkers/S3/data/Segmented/cds/", pattern = "allflu_cds_compcounts_pt_.*\\.rds",full.names = TRUE) %>%
  #   map_dfr(readRDS) %>%
  #   left_join(allflu_cds_df %>% select(gid, n, wgs_dup) %>% 
  #               distinct(gid, .keep_all = TRUE),by = c("cds_id" = "gid")) %>%
  #   filter(n >= 8, wgs_dup == 0) %>%
  #   select(-n, -wgs_dup) %>%
  #   group_by(cds_id) %>%
  #   summarise(across(where(is.numeric), sum, na.rm = TRUE), .groups = "drop") %>%
  #   distinct() %>%
  #   calc_composition_bias(codonpairs = TRUE) %>%
  #   saveRDS("//filepor10/DOP$/X037_AVI_GeneticMarkers/S3/data/Segmented/cds/allflu_wgs_compbias.rds")
  
  # # Arrange protein-wise features in separate columns of one training data set - UNUSED
  # 
  # list.files(path = "S3\\data\\full\\", pattern = "allflu_cds_compbias_pt_.*\\.rds", full.names = TRUE) %>%
  #   map_dfr(readRDS) %>%   
  #   left_join(allflu_cds_df %>% select(cds_id, gid, gene, n, wgs_dup)) %>% 
  #   filter(n >= 8 & wgs_dup == 0) %>% 
  #   select(-n, -wgs_dup) %>%
  #   pivot_wider(names_from = gene, values_from = matches("_Bias"), names_glue = "{gene}_{.value}") %>%
  #   saveRDS(paste0("S3\\data\\full\\allflu_cds_compbias_proteinwise.rds"))
  
  
###############################################################################
# Convert .csv outputs from iFeatureOmega to model-ready protein feature sets #
# See: protein_feat_extract.py                                                #
###############################################################################
  
  out_dir <- "//filepor10/DOP$/X037_AVI_GeneticMarkers/S3/data/Segmented/mlready"
  in_dir  <- "//filepor10/DOP$/X037_AVI_GeneticMarkers/S3/data/Segmented/prot"
  
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  gisaid_keep <- unique(trimws(as.character(meta_ref$Isolate_Id)))
  
  meta_all <- bind_rows(GISAID_avian_prot_df,GISAID_human_prot_df,NCBI_avian_prot_df,NCBI_human_prot_df) %>%
    mutate(gid = trimws(as.character(gid))) %>%
    distinct(gid, .keep_all = TRUE)
  
  cl <- makePSOCKcluster(8)
  registerDoParallel(cl)
  
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
  
  rm(x)
  
  stopCluster(cl)
}

###########################################
# Cluster to collapse redundant sequences #
###########################################
##########################################################################
# Calculate clustering based on k-mer overlaps, considering holdout sets #
##########################################################################

if(load_prev_calcs == FALSE) {
  
  # Read in reference of all sequences used
allflu_cds_ref <- read.csv("//filepor10/DOP$/X037_AVI_GeneticMarkers/S3/data/Segmented/Raw_Sequences\\allflu_cds_df.csv", na.strings = "")

# Parameter sets to test
combos <- crossing(seg = 1:8, minseqid = cluster_minseqid, C = cluster_C)

# Call MMseqs2 14-7e284 to cluster on similarity score - replace with your own MMseqs2 path if needed
foreach(seg = 1:8) %:%
  foreach(minseqid_param = cluster_minseqid) %:%
  foreach(C_param = cluster_C) %do% {
    
    if(!file.exists(paste0("//filepor10/DOP$/X037_AVI_GeneticMarkers/S3/data/Segmented/Clustering\\clust_", seg, "_", minseqid_param*100, "_", C_param*10, "_cluster.tsv"))){
      system(paste0("S3\\scripts\\mmseqs-win64\\mmseqs.bat easy-linclust //filepor10/DOP$/X037_AVI_GeneticMarkers/S3/data/Segmented/allflu_cds_seg_", seg, "_full.FASTA //filepor10/DOP$/X037_AVI_GeneticMarkers/S3/data/Segmented/Clustering\\clust_", seg, "_", minseqid_param*100, "_", C_param*10, " tmp --min-seq-id ", minseqid_param ," -c ", C_param, " --cov-mode 0"),
             intern = TRUE, show.output.on.console = TRUE, ignore.stdout=FALSE, ignore.stderr=FALSE, wait=FALSE)
      
      file.remove(list.files(path = "//filepor10/DOP$/X037_AVI_GeneticMarkers/S3/data/Segmented/Clustering\\", pattern = "*.fasta", full.names = TRUE)) # Keep only the .tsv output from MMseqs2
    }
  }
    
}


###########################################
# Examine and process sequence clustering #
###########################################
##########################################################
# Examine clusters under different clustering parameters #
##########################################################

# Take clusters and match them to GID to find dataset sizes of cluster-centroids
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

# Check clustering sizes for different parameters

results <- combos %>%
  mutate(
    file = paste0("//filepor10/DOP$/X037_AVI_GeneticMarkers/S3/data/Segmented/Clustering\\clust_", seg, "_", minseqid*100, "_", C*10, "_cluster.tsv"),
    output = map(file, possibly(tabulate_n_cluster_combos, otherwise = NULL))
  ) %>%
  unnest_wider(output, names_sep = "temp") %>%
  select(-file) %>%
  set_colnames(c("seg", "minseqid", "C", "nz", "zoon", "total"))

results %>% write.csv("//filepor10/DOP$/X037_AVI_GeneticMarkers/S3/data/Segmented/cds_clust_sizes.csv")

for(i in 1:8){
  results %>%
    filter(seg == i) %>%
    mutate(cell = paste0(zoon, "  /  ", total, "  (", round((total-zoon)/zoon, 1), ")" )) %>%
    select(minseqid, C, cell) %>%
    pivot_wider(names_from = C, values_from = cell) %>% write.csv(paste0("temp_",i,".csv"))
}

############################################################
# Save chosen clustering parameter set as processed labels #
############################################################

# Directories
cluster_dir <- "//filepor10/DOP$/X037_AVI_GeneticMarkers/S3/data/Segmented/clust_final/"
outdir <- "//filepor10/DOP$/X037_AVI_GeneticMarkers/S3/data/Segmented/clust_processed/"

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
    cluster_ref %>% rename(gid = V2),
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
  cluster_reps_list[[i]] <- read.csv(list.files("//filepor10/DOP$/X037_AVI_GeneticMarkers/S3/data/Segmented/clust_processed", full.names = TRUE, pattern = "labels")[i])
  
  fasta_list[[i]] <- readSet(file = paste0("//filepor10/DOP$/X037_AVI_GeneticMarkers/S3/data/Segmented/allflu_cds_seg_", i, "_full.fasta"))
}

# Read in sequences and generate k-mer counts at desired k
reps_kmer_list <- Map(function(x,y) 
  x[y$cluster_rep] %>%
    oligonucleotideFrequency(width = 3, step = 1) %>% 
    as.data.frame %>% 
    mutate(across(where(is.numeric))/rowSums(across(where(is.numeric)))) %>%
    bind_cols(gid = y$cluster_rep, .),
  fasta_list, cluster_reps_list)

reps_kmer_list %>% saveRDS("//filepor10/DOP$/X037_AVI_GeneticMarkers/S3/data/Segmented/clust_dist\\reps_kmer_list.rds")
#reps_kmer_list <- readRDS("//filepor10/DOP$/X037_AVI_GeneticMarkers/S3/data/Segmented/clust_dist\\reps_kmer_list.rds") # For loading at intermediate stage

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

pairwise_JSD %>% saveRDS("//filepor10/DOP$/X037_AVI_GeneticMarkers/S3/data/Segmented/clust_dist\\pairwise_JSD_.rds")
# pairwise_JSD <- readRDS("//filepor10/DOP$/X037_AVI_GeneticMarkers/S3/data/Segmented/clust_dist\\pairwise_JSD_.rds") # For loading at intermediate stage

########################################################
# Hierarchical clustering and testing various fold j's #
########################################################


# Hierarchically from hclust package
hclust_list <- lapply(pairwise_JSD, function(x) x %>%
                        as.dist(diag = TRUE) %>%
                        hclust(method="ward.D2")) # use Ward's minimum variance method to force more compact clustering

hclust_list %>% saveRDS("//filepor10/DOP$/X037_AVI_GeneticMarkers/S3/data/Segmented/clust_dist\\hclust_list.rds")
#hclust_list <- readRDS("//filepor10/DOP$/X037_AVI_GeneticMarkers/S3/data/Segmented/clust_dist\\hclust_list.rds") # For loading at intermediate stage

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

folds_list %>% saveRDS("//filepor10/DOP$/X037_AVI_GeneticMarkers/S3/data/Segmented/clust_dist\\folds_list.rds")


####################
# Checks and scrap #
####################

# Clustering plots

# # Plot cluster dates
# g2 <- cluster_ref %>% select(cluster_rep) %>% distinct %>% 
#   left_join(allflu_wgs_df %>% select(gid, src, label, subtype, date), by = c("cluster_rep" = "gid")) %>%
#   filter(!is.na(date) & date > as.Date("1990-01-01")) %>%
#   add_count(subtype, name = "sub_n") %>%
#   filter(subtype %in% c("H1N1", "H3N8", "H4N6", "H5N1", "H5N2", "H5N6", "H6N2", "H7N3", "H7N4", "H7N9", "H9N2", "H10N7", "H10N8")) %>%
#   ggplot(aes(x = as.Date(date), fill = subtype)) +
#   geom_histogram(position = "stack", binwidth=365) +
#   scale_fill_manual(values = rev(c(RColorBrewer::brewer.pal(12, "Paired"), "black"))) +
#   #scale_x_date(limits = c(as.Date("1990-01-01"), as.Date("2022-12-31"))) +
#   scale_x_date(date_breaks = "1 year", date_labels =  "%y") +
#   facet_grid(rows = vars(label), cols = vars(src), scales = "free_y")
# 
# ggsave(paste0("S3\\figures_tables\\time_dist_clusters_", minseqid_param*100, "_", C_param*10, ".png"), plot = g2, width = 18, height = 6)
}

# # Read cluster ref back in
# cluster_ref <-  read.csv("S3\\data\\full\\cluster_rep_labels.csv")

# png("S3\\figures_tables\\heatmap.png", width = 15, height = 9, units = "in", res = 320)
# heatmap_cds <- cluster_cds_wgs_df %>%
  # as.data.frame() %>%
  # set_rownames(.$gid) %>%
  # select(gid, label, matches("_Bias")) %>%
  # rename_with(., ~ str_replace_all(., c("_Bias" = "", "T" = "U"))) %>%
  # select(matches("^[A|C|G|U][A|C|G|U][A|C|G|U]$")) %>%
  # janitor::remove_constant() %>%
  # as.matrix() %>%
  # t() %>%
  # heatmap.2(
    # density.info = "none",
    # trace = "none",
    # margins = c(1, 5),
    # dendrogram = "col",
    # Rowv = "NA",
    # labCol = NA,
    # ColSideColors = cluster_cds_wgs_df %>% mutate(colsidecol = case_when( # Set side colours using same genus colours as ggplots elsewhere
      # label == "nz" ~ "grey80",
      # label == "zoon" ~ "gray10"
    # )) %>% pull(colsidecol),
    # cexRow = 1.1,
    # col = colorRampPalette(c("dodgerblue", "gray10", "firebrick2"))(n = 65),
    # breaks = c(seq(0, 0.95, length = 28), seq(0.951, 1.05, length = 3), seq(1.051, 3.2, length = 35)),
    # keysize = 0.75, key.par = list(cex = 0.5), key.title = NA, key.xlab = "RSCU",
    # lhei = c(2, 10), lwid = c(2, 10)
  # )
# par(lend = 1) # square line ends for the color legend
# legend(
  # cex = 0.8, x = -0.05, y = 0.9, xpd = TRUE,
  # legend = cluster_cds_wgs_df$label %>% unique() %>% sort(), fill = c("grey80", "gray10"), ncol = 1
# )
# dev.off()

# 
# Save chosen clusters protein-by-protein for genome mapping later #
# 
# 
# for (foclabel in c("nz", "zoon")) {
#   label_files <- list.files(outdir, pattern = "_labels.csv", full.names = TRUE)
#   
#   for (label_file in label_files) {
#     seqs_to_save <- read.csv(label_file, stringsAsFactors = FALSE) %>%
#       mutate(cluster_rep = as.character(cluster_rep),label = as.character(label)) %>%
#       filter(label == foclabel)
#     if (nrow(seqs_to_save) == 0) next
#     for (focgene in gene_order) {
#       
#       # Nucleotide
#       nuc_temp <- seqs_to_save %>%
#         select(cluster_rep) %>%
#         left_join(allflu_nuc_df %>%
#                     mutate(gid = as.character(gid),gene = as.character(gene),gene = if_else(is.na(gene), "NA", gene),gene = case_when(gene == "MP" ~ "M1",gene == "NS" ~ "NS1",TRUE ~ gene)),by = c("cluster_rep" = "gid")) %>%
#         filter(gene == focgene)
#       
#       if (nrow(nuc_temp) > 0) {
#         nuc <- DNAStringSet(nuc_temp$string)
#         names(nuc) <- nuc_temp$cluster_rep
#         writeXStringSet(nuc,filepath = paste0(nuc_outdir, foclabel, "_clusterreps_", focgene, ".FASTA"))
#       }
#       
#       # Coding Region
#       cds_temp <- seqs_to_save %>%
#         select(cluster_rep) %>%
#         left_join(allflu_cds_df %>%
#                     mutate(gid = as.character(gid),gene = as.character(gene),gene = if_else(is.na(gene), "NA", gene)),by = c("cluster_rep" = "gid")) %>%
#         filter(gene == focgene)
#       
#       if (nrow(cds_temp) > 0) {
#         cds <- DNAStringSet(cds_temp$string)
#         names(cds) <- cds_temp$cluster_rep
#         
#         writeXStringSet(cds,filepath = paste0(cds_outdir, foclabel, "_clusterreps_", focgene, ".FASTA"))
#       }
#       
#       # Protein
#       prot_fasta_dir <- "//filepor10/DOP$/X037_AVI_GeneticMarkers/S3/data/Segmented/Raw_Sequences/protein_fastas"
#       prot_files <- list.files(prot_fasta_dir,pattern = paste0("_", focgene, "\\.fasta$"),full.names = TRUE)
#       
#       if (length(prot_files) == 0) next
#       
#       prot_strings <- do.call(c, lapply(prot_files, readAAStringSet))
#       
#       prot_lookup <- data.frame(gid = sub(" .*", "", names(prot_strings)),string = as.character(prot_strings),stringsAsFactors = FALSE) %>%
#         distinct(gid, .keep_all = TRUE)
#       
#       prot_temp <- seqs_to_save %>%
#         select(cluster_rep) %>%
#         left_join(prot_lookup, by = c("cluster_rep" = "gid")) %>%
#         filter(!is.na(string))
#       
#       if (nrow(prot_temp) > 0) {
#         prot <- AAStringSet(prot_temp$string)
#         names(prot) <- prot_temp$cluster_rep
#         
#         writeXStringSet(prot,filepath = paste0(prot_outdir, foclabel, "_clusterreps_", focgene, ".FASTA"))
#       }
#     }
#   }
# }


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


