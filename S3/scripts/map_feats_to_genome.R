library(dplyr)
library(magrittr)
library(DECIPHER)
library(ggmsa)
library(ggplot2)
library(memes)
library(purrr)
library(stringr)
library(tidyr)

# Reference external MEME library
options(meme_bin = "C:\\Users\\lb357c\\meme-5.4.1\\src\\")

# Read in varimp

varimp <- list.files(path = "S3\\analysis\\", pattern = "varimp_perm", full.names = TRUE) %>%
  map_dfr(read.csv) %>%
  separate(var, into = c("feat", "focgene"), sep="_(?=[^_]+$)", remove=FALSE)

varimp %>% arrange(AUC_loss)
varimp %>% group_by(focgene) %>% summarise(mean = mean(AUC_loss)) %>% arrange(mean)
varimp %>% group_by(feat) %>% summarise(mean = mean(AUC_loss)) %>% arrange(mean)

# could plot these by prot..
varimp %>%
  ggplot(aes(x = focgene, y = AUC_loss)) +
  geom_boxplot(alpha = 0.9) +
  theme_bw()

# # Test with IRAT seqs
# irat <- readDNAStringSet(file = "S3\\data\\full\\GISAID_irat_nuc.fasta")

x <- readDNAStringSet(file = "S3\\data\\full\\mapping\\nuc\\zoon_clusterreps_PB2.FASTA")

## DO I NEED TO ALIGN TO A REFERENCE? JB: prob not

## JB: MAFFT is the best aligner

# CorrectFrameshifts needed, then AlignTranslation(x, readingFrame = 1)

align_x <- AlignSeqs(x)
BrowseSeqs(align_x, highlight=0)
align_prot_x <- AlignTranslation(x)
BrowseSeqs(align_prot_x, highlight=0)
ggmsa(align_x, char_width = 0.5, seq_name = T)


align_x %>% pull(seq)
vmatchPattern('TGTGG', align_x) %>% unlist() # simple k-mer matching

# more complex matching: https://www.rdocumentation.org/packages/Biostrings/versions/2.40.2/topics/lowlevel-matching
# match a bookended gap: https://www.rdocumentation.org/packages/Biostrings/versions/2.40.2/topics/matchLRPatterns
# match against a dictionary? could do em all in one go https://www.rdocumentation.org/packages/Biostrings/versions/2.40.2/topics/matchPDict


fas <- system.file("extdata", "50S_ribosomal_protein_L2.fas", package="DECIPHER")
dna <- readDNAStringSet(fas)


db <- system.file("extdata", "Influenza.sqlite", package="DECIPHER")
synteny <- FindSynteny(db, verbose=FALSE)
synteny
InfluenzaA <- AlignSynteny(synteny, db, verbose=FALSE)

# nt_sequences <- system.file("extdata", "LeaderRepeat_All.fa", package = "ggmsa")
# 
# ggmsa(nt_sequences, seq_name = T) + geom_point(aes(x = 10, y = 10), color = "red") # Can overlay stuff
