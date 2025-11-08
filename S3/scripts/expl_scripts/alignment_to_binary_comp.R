######################################################################################################
# This script generates binary sequences for coding sequence composition
######################################################################################################

# Load library
library(Biostrings)
library(parallel)
library(doParallel)
library(foreach)
library(stringr)

# Define paths
project_root <- normalizePath("./S3/")
input_dir <- file.path(project_root, "data/full/mapping/mafft_result/cds")
output_base <- file.path(project_root, "data/full/mapping/binary_result/cds/comp")

# Format a reference table of codons, amino acids and degeneracy values
codon_ref <- tibble::rownames_to_column(data.frame(aminoacid = Biostrings::GENETIC_CODE), "codon")
codon_ref <- dplyr::mutate(codon_ref, aminoacid = gsub("\\*","X",aminoacid))    # replace stop codon symbol "*" as "X"

# Set parallelisation - to disable, comment out, and change %dopar% to %do%
cl <- makePSOCKcluster(6)
registerDoParallel(cl)

## Generate binary vectors for different aspects of genome composition

# Nucleotide composition
binary_conversion_cds_nuc <- function(sequences, output_folder) {
  binary_matrices <- list()
  # Loop through each nucleotide
  foreach(nuc = c("A", "C", "G", "T")) %dopar% {
    # Only run if not already done, or if already done and failed (small file size)
    if (!file.exists(file.path(output_folder, paste0("binary_", nuc, "_Bias.csv"))) |
        (file.exists(file.path(output_folder, paste0("binary_", nuc, "_Bias.csv"))) &
         file.size(file.path(output_folder, paste0("binary_", nuc, "_Bias.csv"))) < 50000)) {
      binary_matrix <- list()
      for (i in seq_along(sequences)) {
        name <- names(sequences)[i]
        seq_raw <- as.character(sequences[[i]])
        seq_raw <- gsub("\\s+", "", seq_raw)
        gapless_seq <- gsub("-", "", seq_raw)
        orig_pos <- which(strsplit(seq_raw, "")[[1]] != "-")
        binary_vector <- rep(0, nchar(seq_raw))
        cds_vector <- strsplit(gapless_seq, "")[[1]]
        cds_binary <- as.integer(cds_vector == nuc)
        binary_vector[orig_pos] <- cds_binary
        binary_matrix[[name]] <- binary_vector
      }
      binary_df <- as.data.frame(do.call(rbind, binary_matrix))
      rownames(binary_df) <- names(binary_matrix)
      binary_matrices[[nuc]] <- binary_df
      output_file <- file.path(output_folder, paste0("binary_", nuc, "_Bias.csv"))
      dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
      write.csv(binary_df, output_file, row.names = TRUE)
    }
  }
  return(binary_matrices)
}

# Dinucleotide composition at position 1 (i.e. 1-2 within codons)
binary_conversion_cds_dinuc_p1 <- function(sequences, output_folder) {
  binary_matrices <- list()
  # Loop through each codon
  foreach(dinuc = c("AA","AC","AG","AT","CA","CC","CG","CT","GA","GC","GG","GT","TA","TC","TG","TT")) %dopar% {
    # Only run if not already done, or if already done and failed (small file size)
    if (!file.exists(file.path(output_folder, paste0("binary_", dinuc, "_p1_Bias.csv"))) |
        (file.exists(file.path(output_folder, paste0("binary_", dinuc, "_p1_Bias.csv"))) &
         file.size(file.path(output_folder, paste0("binary_", dinuc, "_p1_Bias.csv"))) < 50000)) {
      binary_matrix <- list()
      for (i in seq_along(sequences)) {
        name <- names(sequences)[i]
        seq_raw <- as.character(sequences[[i]])
        seq_raw <- gsub("\\s+", "", seq_raw)
        gapless_seq <- gsub("-", "", seq_raw)
        orig_pos <- which(strsplit(seq_raw, "")[[1]] != "-")
        binary_vector <- rep(0, nchar(seq_raw))
        starts <- seq(1, nchar(gapless_seq) - 2, by = 3) # define window to look at given dinucs only
        cds_match <- as.integer(stringr::str_sub(gapless_seq, starts, starts + 1) == dinuc) # compare each on a dinuc-by-dinuc basis
        cds_binary <- as.vector(rbind(cds_match, cds_match, 0)) # apply labels to dinucs at positions 1-2 and fill position 3 with zeroes
        length(cds_binary) <- length(orig_pos)
        cds_binary[is.na(cds_binary)] <- 0      # avoid warnings by padding back up to length of orig pos, essentially ignores incomplete codons at end of string
        binary_vector[orig_pos] <- cds_binary
        binary_matrix[[name]] <- binary_vector
      }
      binary_df <- as.data.frame(do.call(rbind, binary_matrix))
      rownames(binary_df) <- names(binary_matrix)
      binary_matrices[[dinuc]] <- binary_df
      output_file <- file.path(output_folder, paste0("binary_", dinuc, "_p1_Bias.csv"))
      dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
      write.csv(binary_df, output_file, row.names = TRUE)
    }
  }
  return(binary_matrices)
}

# Dinucleotide composition at position 2 (i.e. 2-3 within codons)
binary_conversion_cds_dinuc_p2 <- function(sequences, output_folder) {
  binary_matrices <- list()
  # Loop through each codon
  foreach(dinuc = c("AA","AC","AG","AT","CA","CC","CG","CT","GA","GC","GG","GT","TA","TC","TG","TT")) %dopar% {
    # Only run if not already done, or if already done and failed (small file size)
    if (!file.exists(file.path(output_folder, paste0("binary_", dinuc, "_p2_Bias.csv"))) |
        (file.exists(file.path(output_folder, paste0("binary_", dinuc, "_p2_Bias.csv"))) &
         file.size(file.path(output_folder, paste0("binary_", dinuc, "_p2_Bias.csv"))) < 50000)) {
      binary_matrix <- list()
      for (i in seq_along(sequences)) {
        name <- names(sequences)[i]
        seq_raw <- as.character(sequences[[i]])
        seq_raw <- gsub("\\s+", "", seq_raw)
        gapless_seq <- gsub("-", "", seq_raw)
        orig_pos <- which(strsplit(seq_raw, "")[[1]] != "-")
        binary_vector <- rep(0, nchar(seq_raw))
        starts <- seq(2, nchar(gapless_seq) - 1, by = 3) # define window to look at given dinucs only
        cds_match <- as.integer(stringr::str_sub(gapless_seq, starts, starts + 1) == dinuc) # compare each on a dinuc-by-dinuc basis
        cds_binary <- c(0,as.vector(rbind(cds_match, cds_match, 0))) # apply labels to dinucs at positions 2-3 and fill position 1 with zeroes
        length(cds_binary) <- length(orig_pos)
        cds_binary[is.na(cds_binary)] <- 0      # avoid warnings by padding back up to length of orig pos, essentially ignores incomplete codons at end of string
        binary_vector[orig_pos] <- cds_binary
        binary_matrix[[name]] <- binary_vector
      }
      binary_df <- as.data.frame(do.call(rbind, binary_matrix))
      rownames(binary_df) <- names(binary_matrix)
      binary_matrices[[dinuc]] <- binary_df
      output_file <- file.path(output_folder, paste0("binary_", dinuc, "_p2_Bias.csv"))
      dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
      write.csv(binary_df, output_file, row.names = TRUE)
    }
  }
  return(binary_matrices)
}

# Dinucleotide composition at position 3 (i.e. 3-1 across codons)
binary_conversion_cds_dinuc_p3 <- function(sequences, output_folder) {
  binary_matrices <- list()
  # Loop through each codon
  foreach(dinuc = c("AA","AC","AG","AT","CA","CC","CG","CT","GA","GC","GG","GT","TA","TC","TG","TT")) %dopar% {
    # Only run if not already done, or if already done and failed (small file size)
    if (!file.exists(file.path(output_folder, paste0("binary_", dinuc, "_p3_Bias.csv"))) |
        (file.exists(file.path(output_folder, paste0("binary_", dinuc, "_p3_Bias.csv"))) &
         file.size(file.path(output_folder, paste0("binary_", dinuc, "_p3_Bias.csv"))) < 50000)) {
      binary_matrix <- list()
      for (i in seq_along(sequences)) {
        name <- names(sequences)[i]
        seq_raw <- as.character(sequences[[i]])
        seq_raw <- gsub("\\s+", "", seq_raw)
        gapless_seq <- gsub("-", "", seq_raw)
        orig_pos <- which(strsplit(seq_raw, "")[[1]] != "-")
        binary_vector <- rep(0, nchar(seq_raw))
        starts <- seq(3, nchar(gapless_seq) - 2, by = 3) # define window to look at given dinucs only
        cds_match <- as.integer(stringr::str_sub(gapless_seq, starts, starts + 1) == dinuc) # compare each on a dinuc-by-dinuc basis
        cds_binary <- c(0,0,as.vector(rbind(cds_match, cds_match, 0))) # apply labels to dinucs at positions 3-1 and fill position 2 with zeroes
        length(cds_binary) <- length(orig_pos)
        cds_binary[is.na(cds_binary)] <- 0      # avoid warnings by padding back up to length of orig pos, essentially ignores incomplete codons at end of string
        binary_vector[orig_pos] <- cds_binary
        binary_matrix[[name]] <- binary_vector
      }
      binary_df <- as.data.frame(do.call(rbind, binary_matrix))
      rownames(binary_df) <- names(binary_matrix)
      binary_matrices[[dinuc]] <- binary_df
      output_file <- file.path(output_folder, paste0("binary_", dinuc, "_p3_Bias.csv"))
      dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
      write.csv(binary_df, output_file, row.names = TRUE)
    }
  }
  return(binary_matrices)
}

# Codon composition
binary_conversion_cds_codons <- function(sequences, output_folder) {
  binary_matrices <- list()
  # Loop through each codon
  foreach(codon = codon_ref$codon) %dopar% {
    # Only run if not already done, or if already done and failed (small file size)
    if (!file.exists(file.path(output_folder, paste0("binary_", codon, "_Bias.csv"))) |
        (file.exists(file.path(output_folder, paste0("binary_", codon, "_Bias.csv"))) &
         file.size(file.path(output_folder, paste0("binary_", codon, "_Bias.csv"))) < 50000)) {
      binary_matrix <- list()
      for (i in seq_along(sequences)) {
        name <- names(sequences)[i]
        seq_raw <- as.character(sequences[[i]])
        seq_raw <- gsub("\\s+", "", seq_raw)
        gapless_seq <- gsub("-", "", seq_raw)
        orig_pos <- which(strsplit(seq_raw, "")[[1]] != "-")
        binary_vector <- rep(0, nchar(seq_raw))
        starts <- seq(1, nchar(gapless_seq) - 2, by = 3) # define window to look at codons only
        cds_binary <- rep(as.integer(stringr::str_sub(gapless_seq, starts, starts + 2) == codon), each = 3) # compare each on a codon-by-codon basis
        length(cds_binary) <- length(orig_pos)
        cds_binary[is.na(cds_binary)] <- 0      # avoid warnings by padding back up to length of orig pos, essentially ignores incomplete codons at end of string
        binary_vector[orig_pos] <- cds_binary
        binary_matrix[[name]] <- binary_vector
      }
      binary_df <- as.data.frame(do.call(rbind, binary_matrix))
      rownames(binary_df) <- names(binary_matrix)
      binary_matrices[[codon]] <- binary_df
      output_file <- file.path(output_folder, paste0("binary_", codon, "_Bias.csv"))
      dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
      write.csv(binary_df, output_file, row.names = TRUE)
    }
  }
  return(binary_matrices)
}

# Amino acid composition
binary_conversion_cds_aacids <- function(sequences, output_folder) {
  binary_matrices <- list()
  # Loop through each amino acid
  foreach(aacid = codon_ref$aacid) %dopar% {
    # Only run if not already done, or if already done and failed (small file size)
    if (!file.exists(file.path(output_folder, paste0("binary_", aacid, "_aa_Bias.csv"))) |
        (file.exists(file.path(output_folder, paste0("binary_", aacid, "_aa_Bias.csv"))) &
         file.size(file.path(output_folder, paste0("binary_", aacid, "_aa_Bias.csv"))) < 50000)) {
      binary_matrix <- list()
      for (i in seq_along(sequences)) {
        name <- names(sequences)[i]
        seq_raw <- as.character(sequences[[i]])
        seq_raw <- gsub("\\s+", "", seq_raw)
        gapless_seq <- gsub("-", "", seq_raw)
        orig_pos <- which(strsplit(seq_raw, "")[[1]] != "-")
        binary_vector <- rep(0, nchar(seq_raw))
        starts <- seq(1, nchar(gapless_seq) - 2, by = 3) # define window to look at codons only
        cds_binary <- rep(as.integer(stringr::str_sub(gapless_seq, starts, starts + 2) %in% codon_ref[codon_ref$aminoacid == aacid,]$codon), each = 3) # compare each on a codon-by-codon basis, whilst cross-referencing codon table
        length(cds_binary) <- length(orig_pos)
        cds_binary[is.na(cds_binary)] <- 0      # avoid warnings by padding back up to length of orig pos, essentially ignores incomplete codons at end of string
        binary_vector[orig_pos] <- cds_binary
        binary_matrix[[name]] <- binary_vector
      }
      binary_df <- as.data.frame(do.call(rbind, binary_matrix))
      rownames(binary_df) <- names(binary_matrix)
      binary_matrices[[aacid]] <- binary_df
      output_file <- file.path(output_folder, paste0("binary_", aacid, "_aa_Bias.csv"))
      dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
      write.csv(binary_df, output_file, row.names = TRUE)
    }
  }
  return(binary_matrices)
}

protein_types <- c("HA", "M1", "NA", "NP", "NS1", "PA", "PB1", "PB2")

for (protein in protein_types) {
  fasta_file <- file.path(input_dir, paste0("aligned_", protein, ".fasta"))
  sequences <- readDNAStringSet(fasta_file)
  output_dir <- file.path(output_base, protein)
  # Make sure output directory exists
  if (!dir.exists(output_dir)) {
    stop(paste("Output directory does not exist:", output_dir))
  }
  binary_conversion_cds_nuc(sequences, output_dir)
  binary_conversion_cds_dinuc_p1(sequences, output_dir)
  binary_conversion_cds_dinuc_p2(sequences, output_dir)
  binary_conversion_cds_dinuc_p3(sequences, output_dir)
  binary_conversion_cds_codons(sequences, output_dir)
  binary_conversion_cds_aacids(sequences, output_dir)
  message("Completed binary encoding for: ", protein)
  
}

stopCluster(cl)
