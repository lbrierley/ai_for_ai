######################################################################################################
# This script generates binary sequences for CTriad (Conjoint Triad method)
######################################################################################################


# Load library
library(Biostrings)
library(parallel)
library(doParallel)
library(foreach)

# Define paths
project_root <- normalizePath("./S3/")
input_dir <- file.path(project_root, "data/full/mapping/mafft_result/prot")
output_base <- file.path(project_root, "data/full/mapping/binary_result/prot/ctriad")

# Set parallelisation - to disable, comment out, and change %dopar% to %do%
cl <- makePSOCKcluster(8)
registerDoParallel(cl)

# Function to generate binary vector
binary_conversion_ctriad <- function(sequences, output_folder) {
  ctriad_groups <- list(
    g1 = c("A", "G", "V"),
    g2 = c("I", "L", "F", "P"),
    g3 = c("Y", "M", "T", "S"),
    g4 = c("H", "N", "Q", "W"),
    g5 = c("R", "K"),
    g6 = c("D", "E"),
    g7 = c("C")
  )
  aa_to_group <- unlist(lapply(names(ctriad_groups), function(g) {
    setNames(rep(g, length(ctriad_groups[[g]])), ctriad_groups[[g]])
  }))
  # Generate all possible combinations
  all_patterns <- generate_ctriad_patterns()
  binary_matrices <- list()
  # Loop through each combination
  foreach(pattern = all_patterns) %dopar% {
    # Skip if already done
    if (!file.exists(file.path(output_folder, paste0("binary_CTriad_", pattern, ".csv")))) {
      binary_matrix <- list()
      # Loop through each sequence
      for (i in seq_along(sequences)) {
        name <- names(sequences)[i]
        seq_raw <- as.character(sequences[[i]])
        seq_raw <- gsub("\\s+", "", seq_raw)
        gapless_seq <- gsub("-", "", seq_raw)
        gapless_length <- nchar(gapless_seq)
        orig_pos <- which(strsplit(seq_raw, "")[[1]] != "-")
        binary_vector <- rep(0, nchar(seq_raw))
        # Scan gapless sequence for triads
        for (j in 1:(gapless_length - 2)) {
          triad <- substr(gapless_seq, j, j + 2)
          triad_chars <- strsplit(triad, "")[[1]]
          if (all(triad_chars %in% names(aa_to_group))) {
            triad_group <- paste(aa_to_group[triad_chars], collapse = ".")
            if (triad_group == pattern) {
              orig_kmer_pos <- orig_pos[j:(j + 2)]
              binary_vector[orig_kmer_pos] <- 1
            }
          }
        }
        binary_matrix[[name]] <- binary_vector
      }
      # Convert binary vector into dataframe
      binary_df <- as.data.frame(do.call(rbind, binary_matrix))
      rownames(binary_df) <- names(binary_matrix)
      binary_matrices[[pattern]] <- binary_df
      # Save  dataframe as a CSV file
      output_file <- file.path(output_folder, paste0("binary_CTriad_", pattern, ".csv"))
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
  binary_conversion_ctriad(sequences, output_dir)
  message("Completed binary encoding for: ", protein)
}


stopCluster(cl)