#######################################
# Process other sequences of interest #
#######################################

# Define seqs to be used

set <- "khm"

################################
# Read-in and process metadata #
################################

# Read-in and construct summary df, tidy up where known to be mislabelled
# Nucleotide sequences of segments
g_nuc <- process_GISAID_seq(x = readSet(file = paste0("S3\\data\\ext\\GISAID_", set, "_nuc.fasta")), label = set, type = "nuc")

# Coding sequences for proteins
g_orfs <- g_nuc %>% pull(string) %>% lapply(bind_ORF) %>% bind_rows() # Create GISAID coding sequences from transcripts
g_cds <- g_nuc %>% bind_cols(g_orfs) %>% rowwise %>% mutate(string = substr(string, start, end)) %>% as.data.frame

# Clean names of protein sequence FASTA files and resave for iFeatureOmega
# files <- list.files(path = "S3\\data\\ext\\", pattern = "prot.fasta", full.names = TRUE)
# for(i in 1:length(files)){
#  fasta_name_clean(files[i])
# }

g_prot <- process_GISAID_seq(x = readSet(file = paste0("S3\\data\\ext\\GISAID_", set, "_prot.fasta")), label = set, type = "prot")

##################################################
# Process nucleotide, coding sequences, proteins #
##################################################

nuc_df <- g_nuc

cds_df <- g_cds %>%  # Relabel GISAID genes to reflect likely ORF captured (M1 for MP, NS1 for NS), filter out NCBI genes not present in GISAID (M2, NS2, PA-X, PB1-F2)
  mutate(gene = case_when(
    gene == "MP" ~ "M1",
    gene == "NS" ~ "NS1",
    TRUE ~ gene
  )) %>%
  filter(!(gene %in% c("M2", "NS2", "PA-X", "PB1-F2"))) %>%
  select(-start, -end) %>%
  mutate(cds_id = accession) %>%
  select(-INSDC, -accession)

prot_df <- g_prot %>%
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
  names(cds) <- cds_df %>% filter(gene == unique(cds_df$gene)[i]) %>% pull(cds_id)
  
  temp_df <- cds %>%
    calc_composition_counts(codonpairs = TRUE)
  
  rm(cds)
  
  temp_df %>% 
    saveRDS(paste0("S3\\data\\ext\\cds\\allflu_cds_compcounts_pt_", set, "_", unique(cds_df$gene)[i],".rds"))
  
  temp_df %>%
    calc_composition_bias(codonpairs = FALSE) %>%                     # Do not calculate codon pair biases for individual protein cds
    select(cds_id, matches("_Bias$")) %>%
    left_join(cds_df %>% select(cds_id, gid)) %>%
    relocate(gid) %>%
    saveRDS(paste0("S3\\data\\ext\\mlready\\", set, "\\allflu_cds_compbias_pt_",unique(cds_df$gene)[i],".rds"))
  
  rm(temp_df)
  
}

gc()

# Calculate composition biases for whole genome sequences based off of all proteins - limiting to complete whole genomes only
list.files(path = "S3\\data\\ext\\cds\\", pattern = paste0("allflu_cds_compcounts_pt_", set, "_.*\\.rds"), full.names = TRUE) %>%
  map_dfr(readRDS) %>%   
  left_join(cds_df %>% select(cds_id, gid)) %>%
  group_by(gid) %>%
  summarise_if(is.numeric, sum) %>%
  ungroup %>%
  distinct %>%
  calc_composition_bias(codonpairs = TRUE) %>%                          # Calculate codon pair biases for wgs
  saveRDS(paste0("S3\\data\\ext\\cds\\allflu_wgs_", set, "_compbias.rds"))

######################################################################
# Check pre-calculated and saved protein features from iFeatureOmega #
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

stack <- readRDS("H:\\Working\\ai_for_ai\\S3\\analysis\\stacks\\stack_H7N7.rds")


ext_test <- purrr::map(list.files(path = paste0("H:\\Working\\ai_for_ai\\S3\\data\\ext\\mlready\\", set), full.names = TRUE),
                        function (x) 
                          readRDS(x) %>%
                          select(-any_of(c("segment", "cds_id", "enc", "GC_content"))) %>%
                          rename_with(~paste(., gsub(".*_pt_|.rds", "", x), sep = "_"), -c(gid)) %>%
                          arrange(gid)
) %>% 
  purrr::list_cbind(name_repair = "unique_quiet")

element <- list(main = data.frame(hzoon = predict(plr_stack, newdata=allfeats, type = "prob"), 
                                  label = allfeats$label,
                                  subtype = subtypepicked),
                irat = data.frame(hzoon = predict(plr_stack, newdata=iratfeats, type = "prob"),
                                  rawpred = predict(plr_stack, newdata=iratfeats),
                                  gid = irat_df$gid),
                coef = coef(plr_stack$ens_model$finalModel, 
                            s = plr_stack$ens_model$finalModel$lambdaOpt) %>% 
                  as.matrix %>% 
                  as.data.frame %>%
                  tibble::rownames_to_column(var = "param") %>%
                  bind_rows(data.frame(s1 = plr_stack$ens_model$finalModel$lambdaOpt, param = "lambda")) %>%
                  mutate(subtype = subtypepicked))


######################
# Figures and tables #
######################

library(magrittr)
library(dplyr)
library(tidyr)
library(ggplot2)

df <- read.csv("S3/data/irat/cdc_irat.csv") %>% filter(incomplete != "Y")
allflu_wgs_ref <- read.csv("S3\\data\\full\\allflu_wgs_ref.csv") %>%
  mutate(label = factor(case_when(label == "zoon" ~ "hzoon", label == "nz" ~ "nz"))) # Rearrange factor levels for better compatibility with model functions
holdout_zoon <- c("H7N9", "H5N1", "H9N2", "H5N6", "H10N8", "H7N3", "H3N8", "H7N7", "H7N4")
holdout_nz <- c("H4N6", "H16N3", "H4N8", "H8N4")
holdouts <- c(holdout_zoon, holdout_nz)

result_irat <- read.csv("S3/analysis/stack_irat.csv")
line <- read.csv("S3/analysis/stack_results.csv")

g1 <- result_irat %>%
  left_join(df) %>%
  mutate(host = case_when(host == "za" ~ "avian, human-isolated",
                          host == "a" ~ "avian",
                          host == "zm" ~ "mammal, human-isolated",
                          host == "m" ~ "mammal")) %>%
  ggplot(aes(x = emergence, xmin = emerge_lo, xmax = emerge_hi, y = log(med), ymin = log(upper), ymax = log(lower), color = host, fill = host, label = id)) +
  #  geom_errorbarh(alpha = 0.4, linewidth = 1.2, height = 0) +
  geom_errorbar(alpha = 0.4, linewidth = 1.2, width=0) +
  geom_point() +
  geom_text(hjust=-0.5, vjust=-0.5) +
  geom_hline(aes(yintercept = log(line$threshold)), linetype = "dashed", color = "gray30", linewidth = 1.2) +
  theme_bw() +
  ylab("log(p(zoonotic))") +
  xlab("CDC IRAT emergence score") +
  guides(label = "none") +
  theme(legend.position = c(.85,.75))

g2 <- result_irat %>%
  left_join(df) %>%
  mutate(host = case_when(host == "za" ~ "avian, human-isolated",
                          host == "a" ~ "avian",
                          host == "zm" ~ "mammal, human-isolated",
                          host == "m" ~ "mammal")) %>%
  ggplot(aes(x = impact, xmin = impact_lo, xmax = impact_hi, y = log(med), ymin = log(upper), ymax = log(lower), color = host, fill = host, label = id)) +
  #  geom_errorbarh(alpha = 0.4, linewidth = 1.2, height = 0) +
  geom_errorbar(alpha = 0.4, linewidth = 1.2, width=0) +
  geom_point() +
  geom_text(hjust=-0.5, vjust=-0.5) +
  geom_hline(aes(yintercept = log(line$threshold)), linetype = "dashed", color = "gray30", linewidth = 1.2) +
  theme_bw() +
  ylab("log(p(zoonotic))") +
  xlab("CDC IRAT impact score") +
  guides(label = "none")  +
  theme(legend.position = c(.85,.75))
