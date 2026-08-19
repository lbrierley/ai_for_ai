###############################################
# Process sequences from NCBI Genbank, GISAID #
###############################################
################################
# Read-in and process metadata #
################################

# Read-in and construct summary df, tidy up where known to be mislabelled
# Nucleotide sequences of segments
nuc_results <- process_segmented_sequences(
  path       = "S3\\data\\segmentwise\\Raw_Sequences/nuc",
  ending     = "_Nuc",
  type_label = "nuc"
)

GISAID_avian_nuc_df <- nuc_results$GISAID_avian_df
GISAID_human_nuc_df <- nuc_results$GISAID_human_df
NCBI_avian_nuc_df   <- nuc_results$NCBI_avian_df
NCBI_human_nuc_df   <- nuc_results$NCBI_human_df

# Coding sequences of genes
cds_results <- process_segmented_sequences(
  path       = "S3\\data\\segmentwise\\Raw_Sequences/cds",
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
metadata_dir <- "S3\\data\\segmentwise\\Raw_Sequences/metadata"
meta_results <- build_metadata_excels(metadata_dir)

# GISAID metadata - in long format specifying IDs for each respective protein of each sequence
meta_ref <- bind_rows(
  readxl::read_excel("S3\\data\\segmentwise\\Raw_Sequences/metadata/GISAID_Avian_Metadata.xlsx", guess_max = 1048576) %>%
    select(1:9) %>%
    as.data.frame %>%
    mutate_at(vars(matches("Id$")), ~ gsub("\\|.*$", "", .)) %>%
    reshape2::melt(id.vars = "Isolate_Id") %>%
    mutate(variable = gsub(" Segment_Id", "", variable)),
  readxl::read_excel("S3\\data\\segmentwise\\Raw_Sequences/metadata/GISAID_Human_Metadata.xlsx", guess_max = 1048576) %>%
    select(1:9) %>%
    as.data.frame %>%
    mutate_at(vars(matches("Id$")), ~ gsub("\\|.*$", "", .)) %>%
    reshape2::melt(id.vars = "Isolate_Id") %>%
    mutate(variable = gsub(" Segment_Id", "", variable))
)


# # Clean names of protein sequence FASTA files and resave for iFeatureOmega
files <- list.files(path = "S3\\data\\segmentwise\\Raw_Sequences/prot",
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

######################################################
# Process gene sequences, i.e. only coding sequences #
######################################################

# #Extract Open Reading Frames (ORFs) for GISAID avian - this step takes a while, so un-comment to run for first time.
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

# Load in previously saved files
GISAID_avian_orfs <- read.csv("S3\\data\\segmentwise\\Raw_Sequences/GISAID_avian_orfs.csv")
GISAID_human_orfs <- read.csv("S3\\data\\segmentwise\\Raw_Sequences/GISAID_human_orfs.csv")
GISAID_avian_cds_df <- read.csv("S3\\data\\segmentwise\\Raw_Sequences/GISAID_avian_cds.csv", na.strings = "")
GISAID_human_cds_df <- read.csv("S3\\data\\segmentwise\\Raw_Sequences/GISAID_human_cds.csv", na.strings = "")

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

#############################
# Process protein sequences #
#############################

#Read in FASTA files and combined segment results 
prot_results <- process_segmented_sequences_pro(
  path       = "S3\\data\\segmentwise\\Raw_Sequences/prot",
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
    accession = coalesce(na_if(prot_id, ""), raw_accession, accession)
  ) %>%
  select(-raw_accession)

NCBI_human_prot_df <- prot_results$NCBI_human_df %>%
  mutate(fastahead = gsub(" ", "_", fastahead)) %>%
  left_join(extract_ncbi_accession(prot_results$raw_NCBI_human), by = "fastahead") %>%
  mutate(
    accession = coalesce(na_if(prot_id, ""), raw_accession, accession)
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
  writeXStringSet(seq, filepath = paste0("S3\\data\\segmentwise\\Clustering/allflu_cds_seg_", seg, "_full.FASTA"))
}

###########################
# Creating protein fastas #
###########################

# Output directory
output_dir <- "S3\\data\\segmentwise\\Raw_Sequences/protein_fastas"
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
write.csv(GISAID_avian_orfs, "S3\\data\\segmentwise\\Raw_Sequences/GISAID_avian_orfs.csv")
write.csv(GISAID_human_orfs, "S3\\data\\segmentwise\\Raw_Sequences/GISAID_human_orfs.csv")

# Write csv files
write.csv(allflu_nuc_df,"S3\\data\\segmentwise\\Raw_Sequences/allflu_nuc_df.csv",row.names = FALSE)

write.csv(GISAID_avian_cds_df,"S3\\data\\segmentwise\\Raw_Sequences/GISAID_avian_cds.csv",row.names = FALSE)
write.csv(GISAID_human_cds_df,"S3\\data\\segmentwise\\Raw_Sequences/GISAID_human_cds.csv",row.names = FALSE)
write.csv(allflu_cds_df, "S3\\data\\segmentwise\\Raw_Sequences/allflu_cds_df.csv",row.names = FALSE)

write.csv(GISAID_avian_prot_df,"S3\\data\\segmentwise\\Raw_Sequences/GISAID_avian_prot_df.csv",row.names = FALSE)
write.csv(NCBI_avian_prot_df,"S3\\data\\segmentwise\\Raw_Sequences/NCBI_avian_prot_df.csv",row.names = FALSE)
write.csv(GISAID_human_prot_df,"S3\\data\\segmentwise\\Raw_Sequences/GISAID_human_prot_df.csv",row.names = FALSE)
write.csv(NCBI_human_prot_df,"S3\\data\\segmentwise\\Raw_Sequences/NCBI_human_prot_df.csv",row.names = FALSE)