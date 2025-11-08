#############################################################################################################################################
# This script generates binary sequences for PseAAC (Pseudo-amino acid composition)

# Note that only amino acid composition traits are considered (lambda correlation traits cannot be mapped to sequence as they are universal)
#############################################################################################################################################

# Load library
library(Biostrings)
library(parallel)
library(doParallel)
library(foreach)

# Define paths
project_root <- normalizePath("./S3/")
input_dir <- file.path(project_root, "data/full/mapping/mafft_result/prot")
output_base <- file.path(project_root, "data/full/mapping/binary_result/prot/pseaac")

# Set parallelisation - to disable, comment out, and change %dopar% to %do%
cl <- makePSOCKcluster(6)
registerDoParallel(cl)

# Function to generate binary vector
binary_conversion_pseaac <- function(sequences, output_folder) {
  aa_all <- unique(Biostrings::GENETIC_CODE)
  aa_all <- aa_all[aa_all != "*"] # remove wildcard
  binary_matrices <- list()
  # Loop through each amino acid
  foreach(aa = aa_all) %dopar% {
    # Only run if not already done, or if already done and failed (small file size)
    if (!file.exists(file.path(output_folder, paste0("binary_PAAC_Xc1.", aa, ".csv"))) |
        (file.exists(file.path(output_folder, paste0("binary_PAAC_Xc1.", aa, ".csv"))) &
         file.size(file.path(output_folder, paste0("binary_PAAC_Xc1.", aa, ".csv"))) < 50000)) {
      binary_matrix <- list()
      for (i in seq_along(sequences)) {
        name <- names(sequences)[i]
        seq_raw <- as.character(sequences[[i]])
        seq_raw <- gsub("\\s+", "", seq_raw)
        gapless_seq <- gsub("-", "", seq_raw)
        orig_pos <- which(strsplit(seq_raw, "")[[1]] != "-")
        binary_vector <- rep(0, nchar(seq_raw))
        aa_vector <- strsplit(gapless_seq, "")[[1]]
        aa_binary <- as.integer(aa_vector == aa)
        binary_vector[orig_pos] <- aa_binary
        binary_matrix[[name]] <- binary_vector
      }
      binary_df <- as.data.frame(do.call(rbind, binary_matrix))
      rownames(binary_df) <- names(binary_matrix)
      binary_matrices[[aa]] <- binary_df
      output_file <- file.path(output_folder, paste0("binary_PAAC_Xc1.", aa, ".csv"))
      dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
      write.csv(binary_df, output_file, row.names = TRUE)
    }
  }
  return(binary_matrices)
}

protein_types <- c("HA", "M1", "NA", "NP", "NS1", "PA", "PB1", "PB2")

for (protein in protein_types) {
  fasta_file <- file.path(input_dir, paste0("aligned_", protein, ".fasta"))
  sequences <- readAAStringSet(fasta_file)
  output_dir <- file.path(output_base, protein)
  # Make sure output directory exists
  if (!dir.exists(output_dir)) {
    stop(paste("Output directory does not exist:", output_dir))
  }
  binary_conversion_pseaac(sequences, output_dir)
  message("Completed binary encoding for: ", protein)
  
}

stopCluster(cl)
