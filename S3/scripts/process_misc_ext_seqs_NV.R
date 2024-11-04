#######################################
# Process other sequences of interest #
#######################################

library(caretEnsemble)

# List available external seqs
list.files("S3\\data\\ext\\", pattern = ".fasta", ignore.case=TRUE) %>% gsub("GISAID_|NCBI_|_nuc|_cds|_prot|.fasta|.FASTA", "", .) %>% unique

# Define seqs to be used
set <- "bat"

################################
# Read-in and process metadata #
################################

# Read-in and construct summary df, tidy up where known to be mislabelled
# Nucleotide sequences of segments
v_nuc <- process_NV_seq(x = readSet(file = paste0("S3\\data\\ext\\NCBI_", set, "_nuc.fasta")), label = set, type = "nuc")

# Coding sequences for proteins
v_cds <- process_NV_seq(x = readSet(file =paste0("S3\\data\\ext\\NCBI_", set, "_cds.fasta")), label = set, type = "cds")

# Clean names of protein sequence FASTA files and resave for iFeatureOmega
files <- list.files(path = "S3\\data\\ext\\", pattern = "prot.fasta", full.names = TRUE)
for(i in 1:length(files)){
  prot_fasta_name_clean(files[i])
}

v_prot <- process_GISAID_seq(x = readAAStringSet(file = paste0("S3\\data\\ext\\NCBI_", set, "_prot.fasta")), label = set, type = "prot")

##################################################
# Process nucleotide, coding sequences, proteins #
##################################################

nuc_df <- v_nuc

cds_df <- v_cds %>%  # Relabel GISAID genes to reflect likely ORF captured (M1 for MP, NS1 for NS), filter out NCBI genes not present in GISAID (M2, NS2, PA-X, PB1-F2)
  mutate(gene = case_when(
    gene == "MP" ~ "M1",
    gene == "NS" ~ "NS1",
    TRUE ~ gene
  )) %>%
  filter(!(gene %in% c("M2", "NS2", "PA-X", "PB1-F2"))) %>%
  select(-start, -end) %>%
  mutate(cds_id = accession) %>%
  select(-INSDC, -accession)

prot_df <- v_prot %>%
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
      saveRDS(paste0("S3\\data\\ext\\mlready\\", set, "\\allflu_nuc_",i,"mer_pt_",c("PB2", "PB1", "PA", "HA", "NP", "NA", "M1", "NS1")[j],".rds"))
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
    saveRDS(paste0("S3\\data\\ext\\cds\\allflu_cds_compcounts_pt_", set, "_", unique(cds_df$gene)[i],".rds"))
  
  temp_df %>%
    calc_composition_bias(codonpairs = FALSE) %>%                     # Do not calculate codon pair biases for individual protein cds
    select(cds_id, matches("_Bias$")) %>%
    rename(gid = cds_id) %>%
    saveRDS(paste0("S3\\data\\ext\\mlready\\", set, "\\allflu_cds_compbias_pt_",unique(cds_df$gene)[i],".rds"))
  
  rm(temp_df)
  
}

gc()

# # Calculate composition biases for whole genome sequences based off of all proteins - limiting to complete whole genomes only
# list.files(path = "S3\\data\\ext\\cds\\", pattern = paste0("allflu_cds_compcounts_pt_", set, "_.*\\.rds"), full.names = TRUE) %>%
#   map_dfr(readRDS) %>%   
#   left_join(cds_df %>% select(cds_id, gid)) %>%
#   group_by(gid) %>%
#   summarise_if(is.numeric, sum) %>%
#   ungroup %>%
#   distinct %>%
#   calc_composition_bias(codonpairs = TRUE) %>%                          # Calculate codon pair biases for wgs
#   saveRDS(paste0("S3\\data\\ext\\cds\\allflu_wgs_", set, "_compbias.rds"))

######################################################################
# Check pre-calculated and saved protein features from iFeatureOmega # Note these need to run on WSL installation of Python, and after cleaning FASTA with custom function
######################################################################

for (feat in c("2mer", "ctriad", "ctdc", "ctdt", "ctdd", "pseaac")){
  x <- list.files("S3\\data\\ext\\prot\\", pattern = set, full.names = TRUE) %>% .[grepl(feat, .)] %>%
    read.csv(.) %>%
    bind_cols(prot_df) %>% 
    filter(!(gene %in% c("M2", "NS2", "PA-X", "PB1-F2"))) %>%
    filter(!(grepl("\\|N40\\||\\|M42\\|", fastahead)))
  
  for (j in c("HA", "M1", "NA", "NP", "NS1", "PA", "PB1", "PB2")){
    x %>%
      filter(gene == j) %>%
      select(-X, -title, -UID, -subtype, -date, -protINSDC, -protaccession, -gene, -length, -label, -src, -fastahead) %>%
      relocate(gid) %>%
      saveRDS(paste0("S3\\data\\ext\\mlready\\", set, "\\allflu_prot_",feat,"_pt_",j,".rds"))
  }
  
  rm(x)
}

##############
# Stack test #
##############

# Restrict to only complete sequences
nuc_df %<>% filter(gid %in% (nuc_df %>% group_by(gid) %>% tally %>% filter(n == 8) %>% pull(gid)))
cds_df %<>% filter(gid %in% (cds_df %>% group_by(gid) %>% tally %>% filter(n == 8) %>% pull(gid)))
prot_df %<>% filter(gid %in% (prot_df %>% group_by(gid) %>% tally %>% filter(n == 8) %>% pull(gid)))

# Formalise list of complete sequences
if (length(unique(nuc_df$gid)) <= length(unique(prot_df$gid))){
  canon <- nuc_df %>% select(gid) %>% distinct
} else {
  canon <- prot_df %>% select(gid) %>% distinct
}

ext_test <- purrr::map(list.files(path = paste0("H:\\Working\\ai_for_ai\\S3\\data\\ext\\mlready\\", set), full.names = TRUE),
                       function (x) 
                         readRDS(x) %>%
                         select(-any_of(c("segment", "cds_id", "enc", "GC_content"))) %>%
                         rename_with(~paste(., gsub(".*_pt_|.rds", "", x), sep = "_"), -c(gid)) %>%
                         right_join(canon) %>%
                         arrange(gid) %>%
                         select(-gid)
) %>% 
  purrr::list_cbind(name_repair = "unique_quiet")

# Check failing n rows, and which columns with NAs
which(is.na(ext_test), arr.ind=TRUE) %>% as.data.frame() %>% pull(col) %>% unique %>% names(ext_test)[.]
ext_test %>% filter(!(complete.cases(.))) %>% nrow

# Bind in gid and remove rows with NAs
ext_test <- bind_cols(gid = canon, ext_test)
ext_test %<>% filter(complete.cases(.))

# Generate predictions over all stacks
pred <- lapply(list.files("S3\\analysis\\stacks\\", pattern = ".rds", full.names = TRUE),
               function(x)
                 
                 data.frame(set = set,
                            stack = gsub("S3\\\\analysis\\\\stacks\\\\stack_|.rds", "", x),
                            hzoon = predict(readRDS(x), newdata=ext_test, type = "prob"),
                            gid = ext_test$gid) %>%
                 left_join(v_nuc %>% select(UID, date) %>% distinct, by = c("gid" = "UID")) 
               
               
) %>%
  bind_rows

pred %>% write.csv(paste0("S3\\data\\ext\\preds_", set, "_raw.csv"), row.names=FALSE)

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

pred %>%
  filter(stack == "H5N1") %>%
  # filter(as.Date(date) > "2024-01-01") %>%
  ggplot(aes(x = as.Date(date), y = log(hzoon))) +
  geom_point(color = "#F8766d") +
  geom_hline(aes(yintercept = read.csv("S3/analysis/stack_results.csv") %>% pull(threshold) %>% log()), linetype = "dashed", color = "gray30", linewidth = 1.2, lwd = 1.2) +
  theme_bw() +
  ylab("log(p(zoonotic))") +
  xlab("Date")

# Error bars over all stacks

pred %>%
  group_by(gid, date) %>% 
  # filter(as.Date(date) > "2024-01-01") %>%
  summarise_at(vars("hzoon"), list(med = median, upper = ~quantile(., probs = 0.25), lower = ~quantile(., probs = 0.75))) %>%
  ggplot(aes(x = as.Date(date), y = log(med), ymin = log(upper), ymax = log(lower))) +
  geom_point(color = "#F8766d") +
  geom_errorbar(alpha = 0.4, linewidth = 1.2, width=0, color = "#F8766d") +
  geom_hline(aes(yintercept = read.csv("S3/analysis/stack_results.csv") %>% pull(threshold) %>% log()), linetype = "dashed", color = "gray30", linewidth = 1.2, lwd = 1.2) +
  theme_bw() +
  ylab("log(p(zoonotic))") +
  xlab("Date")

######################
# Figures and tables #
######################

to_plot <- c("khm", "dairyc", "pinnip")

allpred <- purrr::map(list.files("S3\\data\\ext\\", pattern = "preds_.*.csv", full.names = TRUE),
                      function (x) 
                        read.csv(x)) %>% 
  bind_rows

# Manually set annotation labels based on external refs
allpred %<>% 
  mutate(host = case_when(set == "pinnip" ~ "m",
                          set == "dairyc" ~ "m",
                          gid == "EPI_ISL_17394087" ~ "a",                                                                 # A/Red Tailed Hawk/ON/FAV-0473-4/2022
                          gid %in% c("EPI_ISL_19162802", "EPI_ISL_19027114") ~ "zm",                                       # dairy cattle zoonoses: A/Michigan/90/2024, A/Texas/37/2024
                          set == "khm" & gid %in% c("EPI_ISL_19033330", "EPI_ISL_327782", "EPI_ISL_321150") ~ "x",         # environmental seqs: A/Environment/Cambodia/i5T241WW2E/2024, A/environment/Cambodia/Z12CW7M3/2015, A/environment/Cambodia/Z2EP1e3W7M1/2015
                          set == "khm" & gid %in% c("EPI_ISL_18543643", "EPI_ISL_18373263", "EPI_ISL_17069010", "EPI_ISL_18366401", "EPI_ISL_18543355", "EPI_ISL_17024123") ~ "za",  # avian zoonoses: A/Cambodia/2311257/2023, A/Cambodia/NPH230776/2023, A/Cambodia/2302009/2023, A/Cambodia/2310209/2023, A/Cambodia/KSH230332/2023, A/Cambodia/NPH230032/2023
                          set == "khm" ~ "a")
  ) %>%
  mutate(set = case_when(gid %in% c("EPI_ISL_19162802", "EPI_ISL_19027114") ~ "dairyc",         # assign dairy cattle zoonoses to dairy cattle outbreak group
                         TRUE ~ set)
         
  ) %>%
  filter(set %in% to_plot)

# Combi plot

allpred %>%
  filter(stack == "H5N1") %>%   # select stack
  mutate(host = case_when(host == "za" ~ "avian, human-isolated",
                          host == "a" ~ "avian",
                          host == "zm" ~ "mammal, human-isolated",
                          host == "m" ~ "mammal",
                          host == "x" ~ "environmental")) %>%
  filter(as.Date(date) > "2023-01-01") %>%
  ggplot(aes(x = as.Date(date), y = log(hzoon), color = host, fill = host, shape = set)) +
  geom_point(size = 4, alpha = 0.33) +
  scale_color_manual(values = c("#F8766D", "#7CAE00", "black", "#00BFC4", "#C77CFF")) +
  scale_fill_manual(values = c("#F8766D", "#7CAE00", "black", "#00BFC4", "#C77CFF")) +  
  geom_hline(aes(yintercept = read.csv("S3/analysis/stack_results.csv") %>% pull(threshold) %>% log()), linetype = "dashed", color = "gray30", linewidth = 1.2) +
  theme_bw() +
  ylab("log(p(zoonotic))") +
  xlab("Date")

allpred %>%
  group_by(set, gid, date, host) %>% 
  summarise_at(vars("hzoon"), list(med = median, upper = ~quantile(., probs = 0.25), lower = ~quantile(., probs = 0.75))) %>%
  mutate(host = case_when(host == "za" ~ "avian, human-isolated",
                          host == "a" ~ "avian",
                          host == "zm" ~ "mammal, human-isolated",
                          host == "m" ~ "mammal",
                          host == "x" ~ "environmental")) %>%
  filter(as.Date(date) > "2023-01-01") %>%
  ggplot(aes(x = as.Date(date), y = log(med), ymin = log(upper), ymax = log(lower), color = host, fill = host)) +
  geom_point(size = 4, alpha = 0.3) +
  geom_errorbar(alpha = 0.2, linewidth = 1.2, width=0) +
  scale_color_manual(values = c("#F8766D", "#7CAE00", "black", "#00BFC4", "#C77CFF")) +
  scale_fill_manual(values = c("#F8766D", "#7CAE00", "black", "#00BFC4", "#C77CFF")) +  
  geom_hline(aes(yintercept = read.csv("S3/analysis/stack_results.csv") %>% pull(threshold) %>% log()), linetype = "dashed", color = "gray30", linewidth = 1.2) +
  theme_bw() +
  facet_wrap(. ~ set, scales="free_x", ncol=1) +
  ylab("log(p(zoonotic))") +
  xlab("Date")
