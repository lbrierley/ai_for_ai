##################
# Setup packages #
##################

rm(list=ls())

library(Biostrings)
library(coRdon)
library(janitor)
library(ggmosaic)
library(gplots)
library(kableExtra)
library(knitr)
library(lattice)
library(magrittr)
library(ORFik)
library(patchwork)
library(pbapply)
library(plotly)
library(rentrez)
library(reshape2)
library(rmarkdown)
library(R.utils)
library(seqinr)
library(stringr)
library(taxizedb)
library(tidyverse)
library(VennDiagram)
library(XML)
library(caret)

###############
# Set options #
###############

# Set global variables
#setwd("C:\\Users\\Liam\\Desktop\\CoV Genomics")

# Format a reference table of codons, amino acids and degeneracy values
codon_ref <- data.frame(aminoacid = Biostrings::GENETIC_CODE) %>%
  rownames_to_column("codon") %>%
  mutate(aminoacid = gsub("\\*","X",aminoacid)) %>%     # replace stop codon symbol "*" as "X"
  group_by(aminoacid) %>%
  mutate(deg = n())

# Do you want to load previously calculated features data instead of re-calculating?
load_prev_calcs <- TRUE

# Cluster parameters
cluster_minseqid <- c(0.60, 0.70, 0.80)
cluster_C <- c(0.6, 0.7, 0.8)

# Holdout subtypes
holdout_zoon <- c("H7N9", "H5N1", "H9N2", "H5N6", "H10N8", "H7N3", "H3N8", "H7N7", "H7N4")
holdout_nz <- c("H4N6", "H16N3", "H4N8", "H8N4")

###############
# Run scripts #
###############

# # Save data needed for ML
# save(allcov_df, cov_spikes_df, cov_wg_df, file = paste0("cov_ML_dfs_", format(Sys.time(), "%d_%m_%y"), ".RData"))
# #save(allcov_df, cov_wg_df, file = paste0("cov_ML_dfs_noframeshift_", format(Sys.time(), "%d_%m_%y"), ".RData"))
# save(cov_S1_df, file = paste0("cov_S1_dfs_", format(Sys.time(), "%d_%m_%y"), ".RData"))
# save(cov_env_df, file = paste0("cov_env_dfs_", format(Sys.time(), "%d_%m_%y"), ".RData"))
# 
# # Render lab books
# render("C:\\Users\\Liam\\Desktop\\CoV Genomics\\markdown\\data_summary.Rmd", 
#        output_file="C:\\Users\\Liam\\Desktop\\CoV Genomics\\markdown\\data_summary.html")
# 
# render("C:\\Users\\Liam\\Desktop\\CoV Genomics\\markdown\\pca_output_spikes.Rmd", 
#        output_file="C:\\Users\\Liam\\Desktop\\CoV Genomics\\markdown\\pca_output_spikes.html")
# 
# render("C:\\Users\\Liam\\Desktop\\CoV Genomics\\markdown\\pca_output_wgs.Rmd", 
#        output_file="C:\\Users\\Liam\\Desktop\\CoV Genomics\\markdown\\pca_output_wgs.html")
# 
# render("C:\\Users\\Liam\\Desktop\\CoV Genomics\\markdown\\ml_vector_output_spikes.Rmd", 
#        output_file="C:\\Users\\Liam\\Desktop\\CoV Genomics\\markdown\\ml_vector_output_spikes.html")
# 
# render("C:\\Users\\Liam\\Desktop\\CoV Genomics\\markdown\\ml_vector_output_wgs.Rmd",
#        output_file="C:\\Users\\Liam\\Desktop\\CoV Genomics\\markdown\\ml_vector_output_wgs.html")
# 
# # render("C:\\Users\\Liam\\Desktop\\CoV Genomics\\markdown\\ml_vector_output_wgs_noframeshift.Rmd",
# #        output_file="C:\\Users\\Liam\\Desktop\\CoV Genomics\\markdown\\ml_vector_output_wgs_noframeshift.html")
# 
# render("C:\\Users\\Liam\\Desktop\\CoV Genomics\\markdown\\ml_matrix_output_spikes.Rmd", 
#        output_file="C:\\Users\\Liam\\Desktop\\CoV Genomics\\markdown\\ml_matrix_output_spikes.html")
# 
# # render("C:\\Users\\Liam\\Desktop\\CoV Genomics\\markdown\\ml_matrix_output_wgs.Rmd", 
# #        output_file="C:\\Users\\Liam\\Desktop\\CoV Genomics\\markdown\\ml_matrix_output_wgs.html")