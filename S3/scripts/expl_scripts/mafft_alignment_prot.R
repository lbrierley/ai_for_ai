rm(list = ls())

# Install and load required libraries
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!requireNamespace("Biostrings", quietly = TRUE))
  BiocManager::install("Biostrings", ask = FALSE)

library(Biostrings)

# Define directories
input_dir <- "S3/data/full/mapping/mafft_input/prot"
output_dir <- "S3/data/full/mapping/mafft_result/prot"

# Define protein types
protein_types <- c("HA", "M1", "NA", "NP", "NS1", "PA", "PB1", "PB2")

for (protein in protein_types) {
  # Paths
  avian_file <- file.path(input_dir, paste0("nz_clusterreps_", protein, ".FASTA"))
  human_file <- file.path(input_dir, paste0("zoon_clusterreps_", protein, ".FASTA"))
  
  # Read sequences
  avian <- readAAStringSet(avian_file)
  human <- readAAStringSet(human_file)
  
  # Rename
  names(avian) <- paste0("avian_", seq_along(avian))
  names(human) <- paste0("human_", seq_along(human))
  
  # Combine and write
  combined <- c(avian, human)
  combined_path <- file.path(output_dir, paste0("combined_", protein, ".fasta"))
  writeXStringSet(combined, filepath = combined_path)
  
  # Align with MAFFT
  aligned_path <- file.path(output_dir, paste0("aligned_", protein, ".fasta"))
  cmd <- sprintf("mafft --auto --leavegappyregion %s > %s", combined_path, aligned_path)
  system(cmd, intern = TRUE, show.output.on.console = TRUE, ignore.stdout=FALSE, ignore.stderr=FALSE, wait=FALSE)
  
  message("Completed: ", protein)
}
