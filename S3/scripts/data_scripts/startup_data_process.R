####################################################################
#                                                                  #
#          Supporting R scripts for: Brierley et al. 2025          #
#                                                                  #
#     An AI for an AI: identifying zoonotic potential of avian     #
#          influenza viruses via genomic machine learning          #
#                                                                  #
#                    Compiled by L. Brierley                       #
#      University of Liverpool, University of Glasgow, 2025        #
#                                                                  #
####################################################################


##################
# Setup packages #
##################

rm(list=ls())

library(Biostrings) ## installable from Bioconductor via the BiocManager::install function
library(coRdon)     ## installable from Bioconductor via the BiocManager::install function
library(magrittr)
library(ORFik)      ## installable from Bioconductor via the BiocManager::install function
library(reshape2)
library(R.utils)
library(seqinr)
library(stringr)
library(tidyverse)

###############
# Set options #
###############

# Format a reference table of codons, amino acids and degeneracy values
codon_ref <- data.frame(aminoacid = Biostrings::GENETIC_CODE) %>%
  rownames_to_column("codon") %>%
  mutate(aminoacid = gsub("\\*","X",aminoacid)) %>%     # replace stop codon symbol "*" as "X"
  group_by(aminoacid) %>%
  mutate(deg = n())

# Do you want to recalculate genomic and proteomic feature sets? Note this can take several hours.
recalculate_feat_sets <- FALSE

# Do you want to recalculate sequence clustering?
recalculate_cluster <- FALSE

# Set MMseqs2 Linclust parameters to try
cluster_minseqid <- c(0.60, 0.70, 0.80)
cluster_C <- c(0.6, 0.7, 0.8)

# Define chosen MMseqs2 Linclust parameters
cluster_chosen <- c("70_7")

# Define zoonotic and non-zoonotic holdout subtypes of influenza virus
holdout_zoon <- c("H7N9", "H5N1", "H9N2", "H5N6", "H10N8", "H7N3", "H3N8", "H7N7", "H7N4")
holdout_nz <- c("H4N6", "H16N3", "H4N8", "H8N4")

###############################
# Run data processing scripts #
###############################

header(verbose, "Loading custom functions", padding=0)
source("S3\\scripts\\data_scripts\\functions.R")

header(verbose, "Extracting and processing sequence data", padding=0)
source("S3\\scripts\\data_scripts\\process_GISAID_NCBI_data.R")

if(recalculate_feat_sets == TRUE) {
  header(verbose, "Calculating genomic and proteomic feature sets", padding=0)
  source("S3\\scripts\\data_scripts\\calc_feats.R")
}

if(recalculate_cluster == TRUE) {
header(verbose, "Clustering sequences with MMseqs2", padding=0)
source("S3\\scripts\\data_scripts\\cluster_seqs.R")
}

header(verbose, "Processing chosen sequence clusters", padding=0)
source("S3\\scripts\\data_scripts\\process_clusts.R")