######################################################################################################
# This script generates binary sequences for DPC (dipeptide composition)
######################################################################################################

# Load library
library(Biostrings)
library(parallel)
library(doParallel)
library(foreach)

# Define paths
project_root <- normalizePath("./S3/")
input_dir <- file.path(project_root, "data/full/mapping/mafft_result/prot")
output_base <- file.path(project_root, "data/full/mapping/binary_result/prot/dpc")

# Set parallelisation - to disable, comment out, and change %dopar% to %do%
cl <- makePSOCKcluster(8)
registerDoParallel(cl)

# Function to generate all possible k-mers
generate_kmers <- function(k) {
  amino_acids <- c(
    "A", # Alanine
    "R", # Arginine
    "N", # Asparagine
    "D", # Aspartic acid
    "C", # Cysteine
    "E", # Glutamic acid
    "Q", # Glutamine
    "G", # Glycine
    "H", # Histidine
    "I", # Isoleucine
    "L", # Leucine
    "K", # Lysine
    "M", # Methionine
    "F", # Phenylalanine
    "P", # Proline
    "S", # Serine
    "T", # Threonine
    "W", # Tryptophan
    "Y", # Tyrosine
    "V" # Valine
  )
  aa_list <- rep(list(amino_acids), k)
  # Generate all possible combinations of the bases of length k
  kmers_grid <- expand.grid(aa_list)
  kmers <- apply(kmers_grid, 1, paste0, collapse = "")
  return(kmers)
}

k <- 2 # only use k = 2 for dipeptides

# Function to generate binary vector
binary_conversion <- function(sequences, k, output_folder) {
  # Generate all possible k-mers
  all_kmers <- generate_kmers(k)
  binary_matrices <- list()
  # Loop through each kmer
  foreach(kmer = all_kmers) %dopar% {
    binary_matrix <- list() # stores the binary vectors of all sequences
    # Loop through each sequence
    for (i in seq_along(sequences)) {
      name <- names(sequences)[i] # this gets sequence name
      seq_raw <- as.character(sequences[[i]])
      seq_raw <- gsub("\\s+", "", seq_raw) # Remove whitespace
      # Remove gaps to get gapless sequence
      gapless_seq <- gsub("-", "", seq_raw)
      gapless_length <- nchar(gapless_seq)
      # Map gapless positions to original positions
      orig_pos <- which(strsplit(seq_raw, "")[[1]] != "-")
      # Create a binary vector filled with 0, it should be same length as seq_orig
      binary_vector <- rep(0, nchar(seq_raw))
      # Scan gapless sequence for k-mer matches
      for (j in 1:(gapless_length - k + 1)) {
        window <- substr(gapless_seq, j, j + k - 1)
        if (window == kmer) {
          # Map k-mer back to the gapped sequence
          orig_kmer_pos <- orig_pos[j:(j + k - 1)]
          binary_vector[orig_kmer_pos] <- 1
        }
      }
      binary_matrix[[name]] <- binary_vector
    }
    # Convert binary vector into dataframe
    binary_df <- as.data.frame(do.call(rbind, binary_matrix))
    rownames(binary_df) <- names(binary_matrix)
    binary_matrices[[kmer]] <- binary_df
    # Save dataframe as a CSV file
    output_file <- file.path(output_folder, paste0("binary_DPC_", kmer, ".csv"))
    write.csv(binary_df, output_file, row.names = TRUE)
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
  binary_conversion(sequences, k, output_dir)
  message("Completed binary encoding for: ", protein)
}

stopCluster(cl)
