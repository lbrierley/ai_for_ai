######################################################################################################
# This script generates binary sequences for nucleotide k-mers
######################################################################################################

# Load library
library(Biostrings)
library(parallel)
library(doParallel)
library(foreach)

# Set k

k <- 3

# Define paths
project_root <- normalizePath("./S3/")
input_dir <- file.path(project_root, "data/full/mapping/mafft_result/nuc")
output_base <- file.path(project_root, paste0("data/full/mapping/binary_result/nuc/", k, "mers"))

# Set parallelisation - to disable, comment out, and change %dopar% to %do%
cl <- makePSOCKcluster(6)
registerDoParallel(cl)

# Function to generate all possible k-mers
generate_kmers <- function(k) {
  bases <- c("A", "C", "G", "T")
  base_list <- rep(list(bases), k)
  # Generate all possible combinations of the bases of length k
  kmers_grid <- expand.grid(base_list)
  kmers <- apply(kmers_grid, 1, paste0, collapse = "")
  return(kmers)
}

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
    output_file <- file.path(output_folder, paste0("binary_", kmer, ".csv"))
    write.csv(binary_df, output_file, row.names = TRUE)
  }
  return(binary_matrices)
}

# Define list of protein segments
protein_types <- c("HA", "M1", "NA", "NP", "NS1", "PA", "PB1", "PB2")

# Loop through each protein segment
for (protein in protein_types) {
  # Input and output paths
  fasta_file <- file.path(input_dir, paste0("aligned_", protein, ".fasta"))
  output_dir <- file.path(output_base, protein)
  # Make sure output directory exists
  if (!dir.exists(output_dir)) {
    stop(paste("Output directory does not exist:", output_dir))
  }
  # Read aligned sequences
  sequences <- readDNAStringSet(fasta_file)
  # Run binary conversion
  binary_conversion(sequences, k, output_dir)
  message("Completed: ", protein)
}

stopCluster(cl)
