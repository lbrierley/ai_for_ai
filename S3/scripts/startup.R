##################
# Setup packages #
##################

rm(list=ls())

library(Biostrings)
library(coRdon)
#library(ggbiplot)
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


#######################################################################################################################


# Set length filters for segment nucleotide sequences and coding sequences


###############
# Run scripts #
###############

header(verbose, "Loading custom functions", padding=0)
source("S3\\scripts\\functions.R" )

header(verbose, "Extracting and processing sequence data", padding=0)


# Save data needed for ML
save(allcov_df, cov_spikes_df, cov_wg_df, file = paste0("cov_ML_dfs_", format(Sys.time(), "%d_%m_%y"), ".RData"))
#save(allcov_df, cov_wg_df, file = paste0("cov_ML_dfs_noframeshift_", format(Sys.time(), "%d_%m_%y"), ".RData"))
save(cov_S1_df, file = paste0("cov_S1_dfs_", format(Sys.time(), "%d_%m_%y"), ".RData"))
save(cov_env_df, file = paste0("cov_env_dfs_", format(Sys.time(), "%d_%m_%y"), ".RData"))

# Render lab books
render("C:\\Users\\Liam\\Desktop\\CoV Genomics\\markdown\\data_summary.Rmd", 
       output_file="C:\\Users\\Liam\\Desktop\\CoV Genomics\\markdown\\data_summary.html")

render("C:\\Users\\Liam\\Desktop\\CoV Genomics\\markdown\\pca_output_spikes.Rmd", 
       output_file="C:\\Users\\Liam\\Desktop\\CoV Genomics\\markdown\\pca_output_spikes.html")

render("C:\\Users\\Liam\\Desktop\\CoV Genomics\\markdown\\pca_output_wgs.Rmd", 
       output_file="C:\\Users\\Liam\\Desktop\\CoV Genomics\\markdown\\pca_output_wgs.html")

render("C:\\Users\\Liam\\Desktop\\CoV Genomics\\markdown\\ml_vector_output_spikes.Rmd", 
       output_file="C:\\Users\\Liam\\Desktop\\CoV Genomics\\markdown\\ml_vector_output_spikes.html")

render("C:\\Users\\Liam\\Desktop\\CoV Genomics\\markdown\\ml_vector_output_wgs.Rmd",
       output_file="C:\\Users\\Liam\\Desktop\\CoV Genomics\\markdown\\ml_vector_output_wgs.html")

# render("C:\\Users\\Liam\\Desktop\\CoV Genomics\\markdown\\ml_vector_output_wgs_noframeshift.Rmd",
#        output_file="C:\\Users\\Liam\\Desktop\\CoV Genomics\\markdown\\ml_vector_output_wgs_noframeshift.html")

render("C:\\Users\\Liam\\Desktop\\CoV Genomics\\markdown\\ml_matrix_output_spikes.Rmd", 
       output_file="C:\\Users\\Liam\\Desktop\\CoV Genomics\\markdown\\ml_matrix_output_spikes.html")

# render("C:\\Users\\Liam\\Desktop\\CoV Genomics\\markdown\\ml_matrix_output_wgs.Rmd", 
#        output_file="C:\\Users\\Liam\\Desktop\\CoV Genomics\\markdown\\ml_matrix_output_wgs.html")