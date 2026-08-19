####################################################################
#                                                                  #
#          Supporting R scripts for: Brierley et al. 2026          #
#                                                                  #
#     An AI for an AI: identifying zoonotic potential of avian     #
#          influenza viruses via genomic machine learning          #
#                                                                  #
#         Compiled by L. Brierley, E. Abdallah, R. Sanders,        # 
#                      C. Belsey, L. Cattarino                     #
#   University of Liverpool, University of Glasgow, UKHSA, 2026    #
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
library(philentropy)
library(cluster)
library(readxl)
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

# Do you want to recalculate genomic and proteomic feature sets?
recalculate_feat_sets <- FALSE

# Do you want to recalculate sequence clustering?
recalculate_cluster <- FALSE

# Set MMseqs2 Linclust parameters to try
cluster_minseqid <- c(0.85, 0.9, 0.925, 0.95, 0.975, 0.99)
cluster_C <- c(0.7, 0.8, 0.9, 0.95)

###############################
# Run data processing scripts #
###############################

header(verbose, "Loading custom functions", padding=0)
source("S3\\scripts\\data_scripts\\01_functions.R")

header(verbose, "Extracting and processing sequence data", padding=0)
source("S3\\scripts\\data_scripts\\02_process_GISAID_NCBI_data.R")

if(recalculate_feat_sets == TRUE) {
  header(verbose, "Calculating genomic and proteomic feature sets", padding=0)
  source("S3\\scripts\\data_scripts\\03_calc_feats.R") # Set to run with parallelisation
}

if(recalculate_cluster == TRUE) {
header(verbose, "Clustering sequences with MMseqs2", padding=0)
source("S3\\scripts\\data_scripts\\04_cluster_seqs.R")
}

header(verbose, "Processing chosen sequence clusters", padding=0)
source("S3\\scripts\\data_scripts\\05_process_clusts.R")

header(verbose, "Assigning training and test folds from clusters", padding=0)
source("S3\\scripts\\data_scripts\\06_cluster_folds.R")