###############################################################################
# Process sequences manually downloaded from NCBI Flu, GISAID : whole genomes #
###############################################################################

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

NCBI_human_nuc_df <- data.frame(title = NCBI_human_nuc %>% names(), length = NCBI_human_nuc %>% width())
NCBI_human_nuc_df %<>%
  tidyr::separate(title, sep = "\\|", into = c("title", "subtype", "date", "accession", "gene", "segment"), extra = "drop") %>%
  mutate(
    gene = str_sub(str_sub(gene, 4, -1), 1, -2),
    label = "zoon",
    src = "NCBI",
    title = gsub("_", " ", title),
    gid = paste0(title, str_sub(accession, 0, -5)), # Create genome id that can be used for both NCBI and GISAID data
    fastahead = NCBI_human_nuc %>% names()
  )

NCBI_human_cds_df <- data.frame(title = NCBI_human_cds %>% names(), length = NCBI_human_cds %>% width())
NCBI_human_cds_df %<>%
  tidyr::separate(title, sep = "\\|", into = c("title", "subtype", "date", "cds_id", "gene", "segment", "gb", "accession"), extra = "drop") %>%
  mutate(
    fastahead = NCBI_human_cds %>% names(),
    accession = gsub("\\:.*","",accession),
    gene = case_when(
      segment == 1 ~ "PB2",
      segment == 2 & length < 500  ~ "PB1-F2",                               # Distinguish PB1-F2 based on size
      segment == 2 ~ "PB1",                                                  # Else set to PB1
      segment == 3 & grepl("(gb|.*\\:.*\\,).*)",fastahead) ~ "PA-X",         # Distinguish PA-X based on joining ORFs
      segment == 3 ~ "PA",                                                   # Else set to PA
      segment == 4 ~ "HA",
      segment == 5 ~ "NP",
      segment == 6 ~ "NA",
      segment == 7 & grepl("(gb|.*\\:.*\\,).*)",fastahead) ~ "M2",           # Distinguish M2 based on joining ORFs
      segment == 7 ~ "M1",                                                   # Else set to M1
      segment == 8 & grepl("(gb|.*\\:.*\\,).*)",fastahead) ~ "NS2",          # Distinguish NS2 based on joining ORFs
      segment == 8 ~ "NS1"                                                   # Else set to NS1
    ),
    label = "zoon",
    src = "NCBI",
    title = gsub("_", " ", title),
    gid = paste0(title, str_sub(accession, 0, -5)), # Create genome id that can be used for both NCBI and GISAID data
  ) %>%
  select(-gb) %>%
  filter(gene != "N40") %>%
  relocate(fastahead, .after = last_col())

NCBI_avian_nuc_df <- data.frame(title = NCBI_avian_nuc %>% names(), length = NCBI_avian_nuc %>% width())
NCBI_avian_nuc_df %<>%
  tidyr::separate(title, sep = "\\|", into = c("title", "subtype", "date", "accession", "gene", "segment"), extra = "drop") %>%
  mutate(
    gene = str_sub(str_sub(gene, 4, -1), 1, -2),
    label = "nz",
    src = "NCBI",
    subtype = case_when(
      subtype == "H3N6,H3" ~ "H3N6",
      subtype == "H6N1,H6" ~ "H6N1",
      subtype == "H01N2" ~ "H1N2",
      TRUE ~ toupper(subtype)
    ),
    title = case_when(
      title == "Influenza A virus" ~ paste0(title, subtype),
      title %in% c(
        "A/mallard/Italy/3401/2005",
        "A/Avian/Viet Nam/Egg/2014",
        "A/Avian/Viet Nam/Egg/2017",
        "A/chicken/China/embryonated chicken eggs/2014",
        "A/chicken/MUWRP-Uganda/853/2018",
        "A/wild duck/South Korea/KNU2018-26/2020"
      ) ~ paste0(title, date),
      TRUE ~ gsub("_", " ", title)
    ),
    gid = paste0(title, str_sub(accession, 0, -5)), # Create genome id that can be used for both NCBI and GISAID data
    fastahead = NCBI_avian_nuc %>% names()
  )

NCBI_avian_cds_df <- data.frame(title = NCBI_avian_cds %>% names(), length = NCBI_avian_cds %>% width())
NCBI_avian_cds_df %<>%
  tidyr::separate(title, sep = "\\|", into = c("title", "subtype", "date", "cds_id", "gene", "segment", "gb", "accession"), extra = "drop") %>%
  mutate(
    fastahead = NCBI_avian_cds %>% names(),
    accession = gsub("\\:.*","",accession),
    gene = case_when(
      segment == 1 ~ "PB2",
      segment == 2 & length < 500  ~ "PB1-F2",                               # Distinguish PB1-F2 based on size
      segment == 2 ~ "PB1",                                                  # Else set to PB1
      segment == 3 & grepl("(gb|.*\\:.*\\,).*)",fastahead) ~ "PA-X",         # Distinguish PA-X based on joining ORFs
      segment == 3 ~ "PA",                                                   # Else set to PA
      segment == 4 ~ "HA",
      segment == 5 ~ "NP",
      segment == 6 ~ "NA",
      segment == 7 & grepl("(gb|.*\\:.*\\,).*)",fastahead) ~ "M2",           # Distinguish M2 based on joining ORFs
      segment == 7 ~ "M1",                                                   # Else set to M1
      segment == 8 & grepl("(gb|.*\\:.*\\,).*)",fastahead) ~ "NS2",          # Distinguish NS2 based on joining ORFs
      segment == 8 ~ "NS1"                                                   # Else set to NS1
    ),
    label = "nz",
    src = "NCBI",
    subtype = case_when(
      subtype == "H3N6,H3" ~ "H3N6",
      subtype == "H6N1,H6" ~ "H6N1",
      subtype == "H01N2" ~ "H1N2",
      TRUE ~ toupper(subtype)
    ),
    title = case_when(
      title == "Influenza A virus" ~ paste0(title, subtype),
      title %in% c(
        "A/mallard/Italy/3401/2005",
        "A/Avian/Viet Nam/Egg/2014",
        "A/Avian/Viet Nam/Egg/2017",
        "A/chicken/China/embryonated chicken eggs/2014",
        "A/chicken/MUWRP-Uganda/853/2018",
        "A/wild duck/South Korea/KNU2018-26/2020"
      ) ~ paste0(title, date),
      TRUE ~ gsub("_", " ", title)
    ),
    gid = paste0(title, str_sub(accession, 0, -5)) # Create genome id that can be used for both NCBI and GISAID data
  ) %>% 
  select(-gb) %>%
  filter(gene != "N40") %>%
  relocate(fastahead, .after = last_col())

GISAID_avian_nuc_df <- data.frame(title = GISAID_avian_nuc %>% names(), length = GISAID_avian_nuc %>% width())
GISAID_avian_nuc_df %<>%
  tidyr::separate(title, sep = "\\|", into = c("title", "UID", "subtype", "null", "date", "INSDC", "accession", "title2", "gene", "segment"), extra = "drop") %>%
  mutate(
    accession = ifelse(accession %in% "", NA, accession),
    segment = case_when(
      gene == "PB2" ~ "1",
      gene == "PB1" ~ "2",
      gene == "PA" ~ "3",
      gene == "HA" ~ "4",
      gene == "NP" ~ "5",
      gene == "NA" ~ "6",
      gene == "MP" ~ "7",
      gene == "NS" ~ "8"
    ),
    subtype = str_sub(subtype, 5, -1),
    label = "nz",
    src = "GISAID",
    title = gsub("_", " ", title),
    gid = UID,
    fastahead = GISAID_avian_nuc %>% names()
  ) %>%
  select(-null, -title2)

# Adjust where FASTA headers have pulled wrong information directly
GISAID_avian_nuc_fix_ref <- GISAID_avian_nuc[!(GISAID_avian_nuc_df$gene %in% c("PB1", "PB2", "PA", "HA", "NP", "NA", "MP", "NS")), ] %>%
  names() %>%
  as.data.frame() %>%
  setNames("title") %>%
  tidyr::separate(title, sep = "\\|", into = c("nn1", "nn9", "nn2", "nn3", "nn4", "INSDC", "nn6", "nn7", "nn8", "gene", "segment"), extra = "drop") %>%
  select(-contains("nn"))

GISAID_avian_nuc_df %<>%
  rows_update(., GISAID_avian_nuc_fix_ref, by = "INSDC")

GISAID_human_nuc_df <- data.frame(title = GISAID_human_nuc %>% names(), length = GISAID_human_nuc %>% width())
GISAID_human_nuc_df %<>%
  tidyr::separate(title, sep = "\\|", into = c("title", "UID", "subtype", "null", "date", "INSDC", "accession", "title2", "gene", "segment"), extra = "drop") %>%
  mutate(
    accession = ifelse(accession %in% "", NA, accession),
    segment = case_when(
      gene == "PB2" ~ "1",
      gene == "PB1" ~ "2",
      gene == "PA" ~ "3",
      gene == "HA" ~ "4",
      gene == "NP" ~ "5",
      gene == "NA" ~ "6",
      gene == "MP" ~ "7",
      gene == "NS" ~ "8"
    ),
    subtype = str_sub(subtype, 5, -1),
    label = "zoon",
    src = "GISAID",
    title = gsub("_", " ", title),
    gid = UID,
    fastahead = GISAID_human_nuc %>% names()
  ) %>%
  select(-null, -title2)

# Combine all data and remove duplicates - nucleotide sequences of segments
allflu_nuc_df <- bind_rows(GISAID_avian_nuc_df, NCBI_avian_nuc_df, GISAID_human_nuc_df, NCBI_human_nuc_df) %>%
  
  # Remove mixed infections (usually signified by 'MIXED' in NCBI data but also slashes, periods, commas and by H0N0 in GISAID data but some labelled as such in title), (removes 11,788 individual segment sequences)
  filter(!grepl("MIXED|\\,|\\.|\\/", subtype)) %>%
  filter(subtype != "H0N0" & !(grepl("Mixed", title))) %>%
 
   # Remove sequences known a priori to be mislabelled (removes 17 individual segment sequenecs)
  filter(title != "A/goose/Chiayi/18020014-1/2018" & title != "H5N6") %>%
  
  # Remove duplicate records by GenBank accession number but keep records with NA (removes 40,659 individual segment sequences)
  filter(!duplicated(accession, incomparables = NA)) %>%
 
   # Retain canonical sequences belonging to each whole genome according to GISAID metadata (and keep all NCBI) (removes 575 individual segment sequences)
  filter(is.na(UID) | INSDC %in% meta_ref$value) %>%
 
   # Remove duplicate records inferred by having same title, segment and length (removes 55 individual segment sequences)
  arrange(accession) %>%
  distinct(gid, segment, length, .keep_all = TRUE) %>%
  
  # Remove records known to be clones or part of incomplete genomes (removes 8 individual segment sequences)
  filter(!(accession %in% c("CY089569", "KJ907476", "KJ907477", "KJ907481", "KX638412", "LC101678", "MF682884", "MH135674"))) %>%
  
  # Remove based on acceptable lower bounds of segment lengths (removes 78 individual segment sequences)
  filter(segment == 1 & length >= 2280 |
           segment == 2 & length >= 2274 |
           segment == 3 & length >= 2145 |
           segment == 4 & length >= 1659 |
           segment == 5 & length >= 1491 |
           segment == 6 & length >= 1308 |
           segment == 7 & length >= 959 |
           segment == 8 & length >= 820)

# Check to see if multiple of the same gene per genome id:
allflu_nuc_df %>%
  group_by(gid) %>%
  summarise(seqs = n(), segs = n_distinct(segment), genes = n_distinct(gene)) %>%
  filter(seqs != segs | seqs != genes)

allflu_nuc_df %>%
  group_by(gid, label, src) %>%
  summarise(seqs = n(), segs = n_distinct(segment), genes = n_distinct(gene)) %>%
  with(., table(seqs, paste(label, src)))


# Plot lengths
# All retained segment sequences
allflu_nuc_df %>%
  mutate(group = paste(src, label)) %>%
  ggplot(aes(x = length, fill = group)) +
  geom_density(alpha = 0.25) +
  facet_wrap(~segment, ncol = 2, scales = "free")

# Segment sequences belonging to whole genomes only
allflu_nuc_df %>%
  filter(gid %in% (allflu_nuc_df %>% group_by(gid) %>% count() %>% arrange(-n) %>% filter(n == 8) %>% as.data.frame() %>% pull(gid))) %>%
  mutate(group = paste(src, label)) %>%
  ggplot(aes(x = length, fill = group)) +
  geom_density(alpha = 0.25) +
  facet_wrap(~segment, ncol = 2, scales = "free")





# Combine all data and remove duplicates - coding sequences
allflu_cds_df <- bind_rows(NCBI_avian_cds_df, NCBI_human_cds_df) %>%
  
  # Remove mixed infections (usually signified by 'MIXED' in NCBI data but also slashes, periods, commas and by H0N0 in GISAID data but some labelled as such in title), (removes 10,205 individual coding sequences)
  filter(!grepl("MIXED|\\,|\\.|\\/", subtype)) %>%
  filter(subtype != "H0N0" & !(grepl("Mixed", title))) %>%
  
  # Remove sequences known a priori to be mislabelled (removes 25 individual coding sequenecs)
  filter(title != "A/goose/Chiayi/18020014-1/2018" & title != "H5N6") %>%

   # Remove duplicate records by GenBank accession number but keep records with NA (removes XXXX individual coding sequences)
  filter(!duplicated(paste(accession,gene), incomparables = NA)) %>%
  
  # Retain canonical sequences belonging to each whole genome according to GISAID metadata (and keep all NCBI) (removes XXXXX individual coding sequences)
 # filter(is.na(UID) | INSDC %in% meta_ref$value) %>%
  
  # Remove duplicate records inferred by having same title, gene and length (removes XXXX individual coding sequences)
  arrange(accession) %>%
  distinct(gid, gene, length, .keep_all = TRUE) %>%
  
  # Remove records known to be clones or part of incomplete genomes (removes XXXX individual coding sequences)
  filter(!(accession %in% c("CY089569", "KJ907476", "KJ907477", "KJ907481", "KX638412", "LC101678", "MF682884", "MH135674", "FJ384753", "CY039745", "CY039747", "CY039749", "CY038354", "CY038356"))) %>%
  
  # Remove based on acceptable lower bounds of gene lengths (removes XXXX individual coding sequences)
  # REDO IF INCLUDING GISAID SEQUENCES
  filter(gene == "PB2" & length >= 2180 |
           gene == "PB1-F2" & length >= 180 |
           gene == "PB1" & length >= 2250 |
           gene == "PA-X" & length >= 690 |
           gene == "PA" & length >= 1900 |
           gene == "HA" & length >= 1650 |
           gene == "NP" & length >= 1400 |
           gene == "NA" & length >= 1300 |
           gene == "M2" & length >= 250 |
           gene == "M1" & length >= 750 |
           gene == "NS2" & length >= 350 |
           gene == "NS1" & length >= 600)

# Check to see if multiple of the same gene per genome id:
allflu_cds_df %>%
  group_by(gid) %>%
  summarise(seqs = n(), segs = n_distinct(segment), genes = n_distinct(gene)) %>%
  filter(seqs != genes)

allflu_cds_df %>%
  group_by(gid, label, src) %>%
  summarise(seqs = n(), segs = n_distinct(segment), genes = n_distinct(gene)) %>%
  with(., table(seqs, paste(label, src)))


##### FILTER ONES THAT DON'T START ATG?




# Plot lengths
# All retained segment sequences
allflu_cds_df %>%
  mutate(group = paste(src, label)) %>%
  ggplot(aes(x = length, fill = group)) +
  geom_density(alpha = 0.25) +
  facet_wrap(~gene, ncol = 2, scales = "free")

# Coding sequences belonging to complete sets only
allflu_cds_df %>%
  filter(gid %in% (allflu_cds_df %>% group_by(gid) %>% count() %>% arrange(-n) %>% filter(n == 12) %>% as.data.frame() %>% pull(gid))) %>%
  mutate(group = paste(src, label)) %>%
  ggplot(aes(x = length, fill = group)) +
  geom_density(alpha = 0.25) +
  facet_wrap(~gene, ncol = 2, scales = "free")



# Filter DNAStringsSet objects to chosen segment sequences

GISAID_avian_nuc_index <- GISAID_avian_nuc_df %>%
  mutate(index = row_number()) %>%
  inner_join(allflu_nuc_df %>% select(-INSDC, -UID), by = "fastahead") %>%
  pull(index)
GISAID_human_nuc_index <- GISAID_human_nuc_df %>%
  mutate(index = row_number()) %>%
  inner_join(allflu_nuc_df %>% select(-INSDC, -UID), by = "fastahead") %>%
  pull(index)
NCBI_avian_nuc_index <- NCBI_avian_nuc_df %>%
  mutate(index = row_number()) %>%
  inner_join(allflu_nuc_df %>% select(-INSDC, -UID), by = "fastahead") %>%
  pull(index)
NCBI_human_nuc_index <- NCBI_human_nuc_df %>%
  mutate(index = row_number()) %>%
  inner_join(allflu_nuc_df %>% select(-INSDC, -UID), by = "fastahead") %>%
  pull(index)

GISAID_avian_nuc_selected <- GISAID_avian_nuc[GISAID_avian_nuc_index]
GISAID_human_nuc_selected <- GISAID_human_nuc[GISAID_human_nuc_index]
NCBI_avian_nuc_selected <- NCBI_avian_nuc[NCBI_avian_nuc_index]
NCBI_human_nuc_selected <- NCBI_human_nuc[NCBI_human_nuc_index]




# Filter DNAStringsSet objects to chosen coding sequences















# Calculate genomic composition for individual genes


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


# Select a training set with OAT, and use lineages??
