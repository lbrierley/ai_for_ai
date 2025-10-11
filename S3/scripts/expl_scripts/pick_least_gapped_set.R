library(Biostrings) ## installable from Bioconductor via the BiocManager::install function
library(magrittr)
library(stringr)
library(tidyverse)

genes <- c("PB1", "PB2", "PA", "HA", "NP", "NA", "M1", "NS1")

for (i in 1:8) {
  prot_seq <- readAAStringSet(file = paste0("S3\\data\\full\\mapping\\mafft_result\\prot\\aligned_", genes[i], ".fasta"))
  prot_seq[prot_seq %>%
             as.data.frame() %>%
             pull(x) %>%
             str_count("-") %>%
             which.min()] %>%
    writeXStringSet(file = paste0("S3\\data\\full\\mapping\\mafft_result\\prot\\leastgapped_", genes[i], ".fasta"))
  prot_seq[prot_seq %>%
             as.data.frame() %>%
             pull(x) %>%
             str_count("-") %>%
             which.min()] %>%
    DECIPHER::RemoveGaps() %>%
    writeXStringSet(file = paste0("S3\\data\\full\\mapping\\mafft_result\\prot\\leastgapped_gapless_", genes[i], ".fasta"))
  temp <- prot_seq %>%
    consensusMatrix() %>%
    msa::msaConsensusSequence(type = "upperlower", thresh = c(0, 0), ignoreGaps = TRUE) %>%
    AAStringSet() %>%
    setNames(paste0("AVIAN_CONSENSUS_ALL_INSERTS_", genes[i]))
  temp %>%
    writeXStringSet(file = paste0("S3\\data\\full\\mapping\\mafft_result\\prot\\consensus_", genes[i], ".fasta"))
}
