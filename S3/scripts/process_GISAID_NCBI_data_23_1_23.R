###############################################################################
# Process sequences manually downloaded from NCBI Flu, GISAID : whole genomes #
###############################################################################

# Read-in
# Nucleotides (of the genes, not particular coding sequences as some genes encode multiple proteins, e.g. matrix protein 1&2)

GISAID_avian_nuc <- readSet(file = "S3\\data\\full\\GISAID_avian_nuc.fasta")
GISAID_human_nuc <- readSet(file = "S3\\data\\full\\GISAID_human_nuc.fasta")

NCBI_avian_nuc <- readSet(file = "S3\\data\\full\\NCBI_avian_nuc.fasta")
NCBI_human_nuc <- readSet(file = "S3\\data\\full\\NCBI_human_nuc.fasta")

#   # Proteins
# 
# GISAID_avian_prot <- read.fasta(file = "S3\\data\\full\\GISAID_avian_prot.fasta")
# GISAID_human_prot <- read.fasta(file = "S3\\data\\full\\GISAID_human_prot.fasta")
# 
# NCBI_avian_prot <- read.fasta(file = "S3\\data\\full\\NCBI_avian_prot.fasta")
# NCBI_human_prot <- read.fasta(file = "S3\\data\\full\\NCBI_human_prot.fasta")

# Metadata

# GISAID reads in from file
GISAID_avian_meta <- readxl::read_excel("S3/data/full/GISAID_avian_meta.xls", guess_max = 1048576)
GISAID_human_meta <- readxl::read_excel("S3/data/full/GISAID_human_meta.xls", guess_max = 1048576)



### CHECK PROTEIN LENGTHS

# Construct summary df for individual genes from particular FASTA headers

NCBI_human_nuc_df <- data.frame(title = NCBI_human_nuc %>% names(), length = NCBI_human_nuc %>% width())
NCBI_human_nuc_df %<>%
  tidyr::separate(title, sep = "\\|", into =  c("title", "subtype", "date", "accession","gene", "segment"))  %>%
  mutate(gene = str_sub(str_sub(gene,4,-1), 1, -2),
         label = "zoon",
         src = "NCBI")


NCBI_avian_nuc_df <- data.frame(title = NCBI_avian_nuc %>% names(), length = NCBI_avian_nuc %>% width())
NCBI_avian_nuc_df %<>%
  tidyr::separate(title, sep = "\\|", into =  c("title", "subtype", "date", "accession","gene", "segment"))  %>%
  mutate(gene = str_sub(str_sub(gene,4,-1), 1, -2),
         label = "nz",
         src = "NCBI")


GISAID_avian_nuc_df <- data.frame(title = GISAID_avian_nuc %>% names(), length = GISAID_avian_nuc %>% width())
GISAID_avian_nuc_df %<>%
  tidyr::separate(title, sep = "\\|", into =  c("title", "UID", "subtype", "null", "date", "INSDC", "accession","title2", "gene", "segment")) %>%
  mutate(accession = ifelse(accession %in% "", NA, accession),
         segment = case_when(gene == "PB2" ~ "1",
                             gene == "PB1" ~ "2",
                             gene == "PA" ~ "3",
                             gene == "HA" ~ "4",
                             gene == "NP" ~ "5",
                             gene == "NA" ~ "6",
                             gene == "MP" ~ "7",
                             gene == "NS" ~ "8"
         ),
         subtype = str_sub(subtype,5,-1),
         label = "nz",
         src = "GISAID") %>%
  select(-null, -title2)

#### Tidy up mislabelled sequence FASTA headers
GISAID_avian_nuc_fix_ref <- GISAID_avian_nuc[!(GISAID_avian_nuc_df$gene %in% c("PB1", "PB2", "PA", "HA", "NP", "NA", "MP", "NS")),] %>% 
  names %>%
  as.data.frame %>%
  setNames("title") %>%
  tidyr::separate(title, sep = "\\|", into =  c("nn1", "nn9", "nn2", "nn3", "nn4", "INSDC", "nn6","nn7", "nn8", "gene", "segment")) %>%
  select(-contains("nn"))

GISAID_avian_nuc_df %<>% 
  rows_update(., GISAID_avian_nuc_fix_ref, by = "INSDC")
####

GISAID_human_nuc_df <- data.frame(title = GISAID_human_nuc %>% names(), length = GISAID_human_nuc %>% width())
GISAID_human_nuc_df %<>%
  tidyr::separate(title, sep = "\\|", into =  c("title", "UID", "subtype", "null", "date", "INSDC", "accession","title2", "gene", "segment")) %>%
  mutate(accession = ifelse(accession %in% "", NA, accession),
         segment = case_when(gene == "PB2" ~ "1",
                             gene == "PB1" ~ "2",
                             gene == "PA" ~ "3",
                             gene == "HA" ~ "4",
                             gene == "NP" ~ "5",
                             gene == "NA" ~ "6",
                             gene == "MP" ~ "7",
                             gene == "NS" ~ "8"
         ),
         subtype = str_sub(subtype,5,-1),
         label = "zoon",
         src = "GISAID") %>%
  select(-null, -title2)


## little surprising there aren't equal numbers of each gene... should all be complete so check from metadata?
NCBI_avian_nuc_df %>% group_by(title) %>% count %>% pull(n) %>% table
GISAID_avian_nuc_df %>% group_by(title) %>% count %>% pull(n) %>% table
NCBI_human_nuc_df %>% group_by(title) %>% count %>% pull(n) %>% table
GISAID_human_nuc_df %>% group_by(title) %>% count %>% pull(n) %>% table
# May just clean when compiling whole genomes
# Seen at least one instance of an NCBI record (!) mislabelled in the FASTA.. MW334541

NCBI_avian_nuc_df %>% 
  group_by(title) %>%
  summarise(seqs = n(), segs = n_distinct(segment), genes = n_distinct(gene)) %>% filter(seqs != segs|seqs != genes)

GISAID_avian_nuc_df %>% 
  group_by(title) %>%
  summarise(seqs = n(), segs = n_distinct(segment), genes = n_distinct(gene)) %>% filter(seqs != segs|seqs != genes)

NCBI_human_nuc_df %>% 
  group_by(title) %>%
  summarise(seqs = n(), segs = n_distinct(segment), genes = n_distinct(gene)) %>% filter(seqs != segs|seqs != genes)

GISAID_human_nuc_df %>% 
  group_by(title) %>%
  summarise(seqs = n(), segs = n_distinct(segment), genes = n_distinct(gene)) %>% filter(seqs != segs|seqs != genes)


# Check lengths

bind_rows(NCBI_avian_nuc_df, GISAID_avian_nuc_df, NCBI_human_nuc_df, GISAID_human_nuc_df) %>%
  filter(length > 750) %>%
  mutate(group = paste(src, label)) %>%
  ggplot(aes(x=length, fill=group)) + geom_density(alpha=0.25) + facet_wrap(~ gene, ncol=4, scales="free")




# rm everything under 750 in the first instance
# Maybe remove all GISAIDs < min GenBank in each gene? or < min reported as complete? or < min consensus?
# Too demanding to get "complete/partial" from metadata, 88hr search..
# bind_rows(NCBI_avian_nuc_df, GISAID_avian_nuc_df, NCBI_human_nuc_df, GISAID_human_nuc_df) %>% filter(segment==1) %>% arrange(length) %>% .[1:9]
# PB2
# PB1
# PA
# HA
# NP
# NA
# MP
# NS:  consensus min seems to be 823




# Find duplicates in GISAID nuc sequences

GISAID_avian_nuc %>% duplicated %>% grep(TRUE, .)
GISAID_human_nuc %>% duplicated %>% grep(TRUE, .)

## Remove duplicates - stitch into whole genomes first??
# GISAID_avian_nuc %<>% .[-(GISAID_avian_nuc %>% duplicated %>% grep(TRUE, .))]
# GISAID_human_nuc %<>% .[-(GISAID_human_nuc %>% duplicated %>% grep(TRUE, .))]

# Should be no duplicates within NCBI nuc sequences

NCBI_avian_nuc %>% duplicated %>% grep(TRUE, .)
NCBI_human_nuc %>% duplicated %>% grep(TRUE, .)

# Find duplicates between GISAID and NCBI nuc sequences (could do by accession, but seems safest to manually check)
# Quite a lot of overlap

c(GISAID_avian_nuc, NCBI_avian_nuc) %>% duplicated %>% grep(TRUE, .)
c(GISAID_human_nuc, NCBI_human_nuc) %>% duplicated %>% grep(TRUE, .)

## Remove duplicates - stitch into whole genomes first??
# all_avian_nuc <- c(GISAID_avian_nuc, NCBI_avian_nuc)
# all_human_nuc <- c(GISAID_human_nuc, NCBI_human_nuc)
# all_avian_nuc %<>% .[-(all_avian_nuc %>% duplicated %>% grep(TRUE, .))]
# all_human_nuc %<>% .[-(all_human_nuc %>% duplicated %>% grep(TRUE, .))]




### FILTER H0N0 - THIS IS A SIGNIFIER FOR A MIXED INFECTION
# Fix H01N2 replace with H1N2


bind_rows(NCBI_avian_nuc_df, GISAID_avian_nuc_df, NCBI_human_nuc_df, GISAID_human_nuc_df) %>% 
  filter(grepl("H0N0", subtype)) %>% 
  select(title, subtype) %>% 
  unique


### MERGE WHOLE GENOMES



























# Calculate genomic composition counts for individual CDS
cov_enc_df %<>% cbind(data.frame(
  enc = cov_cord %>% codonTable() %>% ENC(stop.rm = FALSE), # Calculate Effective Number of Codons (including STOP codons)
  cov_cord %>% letterFrequency("GC", as.prob = TRUE) * 100 %>% as.vector(), # Calculate % GC content
  cov_cord %>% letterFrequency(c("A", "C", "G", "T")), # Nucleotide counts
  cov_cord %>% dinucleotideFrequency(), # Dinucleotide counts
  cov_cord %>%
    DNAStringSet(start = 1) %>%
    dinucleotideFrequency(step = 3) %>%
    as.data.frame() %>%
    rename_all(., ~ paste0(., "_p1")), # Dinucleotide counts between positions 1-2 only
  cov_cord %>%
    DNAStringSet(start = 2) %>%
    dinucleotideFrequency(step = 3) %>%
    as.data.frame() %>%
    rename_all(., ~ paste0(., "_p2")), # Dinucleotide counts between positions 2-3 only
  cov_cord %>%
    DNAStringSet(start = 3) %>%
    dinucleotideFrequency(step = 3) %>%
    as.data.frame() %>%
    rename_all(., ~ paste0(., "_p3")), # Dinucleotide counts between positions 3-1 only
  cov_cord %>%
    codonTable() %>%
    codonCounts() # Codon counts
  #                        cov_cord %>% oligonucleotideFrequency(6, step=3)                   # Codon pair counts - NOT USING FOR NOW
)) %>% rename_at(vars(G.C), ~"GC_content")


merge(cov_wg_df, allcov_meta_df, all.x = TRUE, by = "accessionversion") %>%
  # filter(length > 10000) %>%
  ggplot(aes(x = meta_length, y = CDS_length)) +
  # scale_x_continuous(limits = c(10000,50000)) +
  # scale_y_continuous(limits = c(10000, 50000)) +
  geom_hline(yintercept = min_wgs_length, color = "dodgerblue", alpha = 0.4, size = 1.2) +
  geom_hline(yintercept = max_wgs_length, color = "dodgerblue", alpha = 0.4, size = 1.2) +
  geom_point()

cov_enc_df %>%
  filter(is.na(protein)) %>%
  nrow() # n missing protein annotation
cov_enc_df %>%
  filter(is.na(gene)) %>%
  nrow() # n missing gene annotation

# Clean up metadata and assign each cds as whole S, S1, S2 or other, and highlight problematic cds
cov_enc_df %<>% mutate(seqtype = case_when(
  grepl("^env", protein, ignore.case = TRUE) ~ "E",
  grepl("^s1| s1|^s-1| s-1|subunit 1$|spike 1|spike protein 1|spike glycoprotein 1", protein, ignore.case = TRUE) ~ "S1",
  grepl("^s2| s2|^s-2| s-2|subunit 2$|spike 2|spike protein 2|spike glycoprotein 2", protein, ignore.case = TRUE) ~ "S2",
  grepl("spike|^surface|surface gly|s gly|s prot|peplom|protein S$|^S$", protein, ignore.case = TRUE) ~ "S"
)) %>% replace_na(list(seqtype = "other"))

# Override with gene information where available
cov_enc_df %<>% mutate(seqtype = case_when(
  grepl("^E|env", protein, ignore.case = TRUE) ~ "E",
  grepl("^s1| s1|^s-1| s-1|subunit 1$|spike 1|spike protein 1|spike glycoprotein 1", gene, ignore.case = TRUE) ~ "S1",
  grepl("^s2| s2|^s-2| s-2|subunit 2$|spike 2|spike protein 2|spike glycoprotein 2", gene, ignore.case = TRUE) ~ "S2",
  TRUE ~ seqtype
))
















### Need to combine each nuc/protein into one wgs record to work with, may need coRdon?
### May need to convert GISAID into long format?




# # NCBI may have to have programmatic access, interface consistently failing to fetch metadata, but v slow so do we really need it??
# Seq_summary <- function(x){
#   
#   query_index <- split(seq(1,length(x)), ceiling(seq_along(seq(1,length(x)))/300))
#   Seq_result <- vector("list", length(query_index))
#   
#   for (i in 1:length(query_index)) {
#     Seq_result[[i]] <- entrez_summary(db = "nuccore",id = x[unlist(query_index[i])])
#     Sys.sleep(5)
#   }
#   
#   if(length(x) == 1){
#     return(Seq_result)
#   } else {
#     return(Seq_result %>% flatten %>% unname)
#   }
# }

# # Fetch metadata information for NUCLEOTIDE sequences only (title, accession). NB entrez_summary doesn't accept empty objects, so filtering out those IDs with no sequences
# if (load_prev_seqs == TRUE) {
#   load(file = "S3\\data\\full\\NCBI_avian_meta.RData")
#   load(file = "S3\\data\\full\\NCBI_human_meta.RData")
# } else {
#   NCBI_human_meta <- pblapply(names(NCBI_human_nuc), function(x) Seq_summary(x)) # ETA 2H20 FOR HUMAN
#   save(NCBI_human_meta, file = "data\\NCBI_human_meta.RData")
#   NCBI_avian_meta <- pblapply(names(NCBI_avian_nuc), function(x) Seq_summary(x)) # ETA 88H BASIC CALC FOR AVIAN
#   save(NCBI_avian_meta, file = "data\\NCBI_avian_meta.RData")
# }

# Check lengths?

# Select a training set with OAT, and use lineages??

