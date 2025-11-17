######################################################################################################
# This script generates binary sequences for CTD (Composition-Transition-Distribution)
######################################################################################################

# Load library
library(Biostrings)
library(parallel)
library(doParallel)
library(foreach)

# Define paths
project_root <- normalizePath("./S3/")
input_dir <- file.path(project_root, "data/full/mapping/mafft_result/prot")
output_base <- file.path(project_root, "data/full/mapping/binary_result/prot/ctdt")

# Set parallelisation - to disable, comment out, and change %dopar% to %do%
cl <- makePSOCKcluster(6)
registerDoParallel(cl)

# define decriptor groups
descriptor_groups <- list(
  hydrophobicity_PRAM900101 = list(
    G1 = c("R", "K", "E", "D", "Q", "N"),
    G2 = c("G", "A", "S", "T", "P", "H", "Y"),
    G3 = c("C", "L", "V", "I", "M", "F", "W")
  ),
  hydrophobicity_ARGP820101 = list(
    G1 = c("Q", "S", "T", "N", "G", "D", "E"),
    G2 = c("R", "A", "H", "C", "K", "M", "V"),
    G3 = c("L", "Y", "P", "F", "I", "W")
  ),
  hydrophobicity_ZIMJ680101 = list(
    G1 = c("Q", "N", "G", "S", "W", "T", "D", "E", "R", "A"),
    G2 = c("H", "M", "C", "K", "V"),
    G3 = c("L", "P", "F", "Y", "I")
  ),
  hydrophobicity_PONP930101 = list(
    G1 = c("K", "P", "D", "E", "S", "N", "Q", "T"),
    G2 = c("G", "R", "H", "A"),
    G3 = c("Y", "M", "F", "W", "L", "C", "V", "I")
  ),
  hydrophobicity_CASG920101 = list(
    G1 = c("K", "D", "E", "Q", "P", "S", "R", "N", "T", "G"),
    G2 = c("A", "H", "Y", "M", "L", "V"),
    G3 = c("F", "I", "W", "C")
  ),
  hydrophobicity_ENGD860101 = list(
    G1 = c("R", "D", "K", "E", "N", "Q", "H", "Y", "P"),
    G2 = c("S", "G", "T", "A", "W"),
    G3 = c("C", "V", "L", "I", "M", "F")
  ),
  hydrophobicity_FASG890101 = list(
    G1 = c("K", "E", "R", "S", "Q", "D"),
    G2 = c("N", "T", "P", "G"),
    G3 = c("A", "Y", "H", "W", "V", "M", "F", "L", "I", "C")
  ),
  normwaalsvolume = list(
    G1 = c("G", "A", "S", "T", "P", "D"),
    G2 = c("N", "V", "E", "Q", "I", "L"),
    G3 = c("M", "H", "K", "F", "R", "Y", "W")
  ),
  polarity = list(
    G1 = c("L", "I", "F", "W", "C", "M", "V", "Y"),
    G2 = c("P", "A", "T", "G", "S"),
    G3 = c("H", "Q", "R", "K", "N", "E", "D")
  ),
  polarizability = list(
    G1 = c("G", "A", "S", "D", "T"),
    G2 = c("G", "P", "N", "V", "E", "Q", "I", "L"),
    G3 = c("K", "M", "H", "F", "R", "Y", "W")
  ),
  charge = list(
    G1 = c("K", "R"),
    G2 = c("A", "N", "C", "Q", "G", "H", "I", "L", "M", "F", "P", "S", "T", "W", "Y", "V"),
    G3 = c("D", "E")
  ),
  secondarystruct = list(
    G1 = c("E", "A", "L", "M", "Q", "K", "R", "H"),
    G2 = c("V", "I", "Y", "C", "W", "F", "T"),
    G3 = c("G", "N", "P", "S", "D")
  ),
  solventaccess = list(
    G1 = c("A", "L", "F", "C", "G", "I", "V", "W"),
    G2 = c("R", "K", "Q", "E", "N", "D"),
    G3 = c("M", "P", "S", "T", "H", "Y")
  )
)

# function that encodes transitions per descriptor
binary_conversion_ctd_t <- function(sequences, descriptor_name, output_folder) {
  groups <- descriptor_groups[[descriptor_name]]
  transitions <- c("1221", "1331", "2332")
  binary_matrices <- list()
  foreach(transition_name = transitions) %dopar% {
    # Only run if not already done, or if already done and failed (small file size)
    if (!file.exists(file.path(output_folder, paste0("binary_CTDT_", descriptor_name, ".Tr", transition_name, ".csv"))) |
        (file.exists(file.path(output_folder, paste0("binary_CTDT_", descriptor_name, ".Tr", transition_name, ".csv"))) &
         file.size(file.path(output_folder, paste0("binary_CTDT_", descriptor_name, ".Tr", transition_name, ".csv"))) < 50000)) {
      # identify groups from transition code string
      transitions_valid <- as.numeric(unique(unlist(strsplit(transitions[1], "")))) 
      g_a <- groups[[transitions_valid[1]]]
      g_b <- groups[[transitions_valid[2]]]
      # identify valid dipeptides for given transition code string
      transition_kmers <- c(as.vector(outer(g_a, g_b, paste0)), as.vector(outer(g_b, g_a, paste0))) 
      binary_matrix <- list()
      for (i in seq_along(sequences)) {
        name <- names(sequences)[i]
        seq_raw <- as.character(sequences[[i]])
        seq_raw <- gsub("\\s+", "", seq_raw)
        gapless_seq <- gsub("-", "", seq_raw)
        gapless_length <- nchar(gapless_seq)
        orig_pos <- which(strsplit(seq_raw, "")[[1]] != "-")
        binary_vector <- rep(0, nchar(seq_raw))
        for (j in 1:(gapless_length - 1)) {
          window <- substr(gapless_seq, j, j + 1)
          if (window %in% transition_kmers) {
            # Map k-mer back to the gapped sequence
            orig_kmer_pos <- orig_pos[j:(j + 1)]
            binary_vector[orig_kmer_pos] <- 1
          }
        }
        binary_matrix[[name]] <- binary_vector
      }
      binary_df <- as.data.frame(do.call(rbind, binary_matrix))
      rownames(binary_df) <- names(binary_matrix)
      binary_matrices[[transition_name]] <- binary_df
      output_file <- file.path(output_folder, paste0("binary_CTDT_", descriptor_name, ".Tr", transition_name, ".csv"))
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
  for (descriptor in names(descriptor_groups)) {
    output_dir <- file.path(output_base, protein)
    # Make sure output directory exists
    if (!dir.exists(output_dir)) {
      stop(paste("Output directory does not exist:", output_dir))
    }
    binary_conversion_ctd_t(sequences, descriptor, output_dir)
    message("Completed binary encoding for: ", protein, " using ", descriptor)
  }
}

stopCluster(cl)
