#######################################################################
# This script takes binary sequences and calculates raw odds ratios
#######################################################################

rm(list = ls())

project_root <- normalizePath("./S3/")

# Input and output paths for chosen nucleotide/protein feature set
input_base <- file.path(project_root, "data/full/mapping/binary_result/nuc/3mers")
output_base <- file.path(project_root, "data/full/mapping/kmer_or_results/nuc/3mers")

# Function to compute odds ratio per-position
or_function <- function(
    binary_df,
    human_pat = "HUMAN",
    avian_pat = "AVIAN",
    haldane = TRUE, #  can be turned to False
    drop_all_zero_cols = TRUE,
    conf.level = 0.95) {
  rn <- rownames(binary_df)
  is_human <- grepl(human_pat, rn, ignore.case = TRUE)
  is_avian <- grepl(avian_pat, rn, ignore.case = TRUE)
  if (!any(is_human) || !any(is_avian)) stop("Adjust human_pat / avian_pat.") # debugging
  mat <- as.matrix(binary_df)
  storage.mode(mat) <- "numeric"
  # drop all zero columns
  if (drop_all_zero_cols) {
    keep <- colSums(mat[is_human, , drop = FALSE]) + colSums(mat[is_avian, , drop = FALSE]) > 0
    mat <- mat[, keep, drop = FALSE]
  } else {
    keep <- rep(TRUE, ncol(mat))
  }
  if (ncol(mat) == 0) {
    return(list(
      by_position = data.frame()
    ))
  }
  # per position count
  a <- colSums(mat[is_human, , drop = FALSE] == 1) # human 1s
  c0 <- colSums(mat[is_human, , drop = FALSE] == 0) # human 0s
  b <- colSums(mat[is_avian, , drop = FALSE] == 1) # avian 1s
  d <- colSums(mat[is_avian, , drop = FALSE] == 0) # avian 0s
  # using haldane methods to avoid infinity
  # However, +0.5 has large influences given the difference in avian and human sequence size, so add a scaled factor
  if (haldane) {
    human_avian_ratio <- (a + c0) / (b + d)
    human_odds <- (a + human_avian_ratio) / (c0 + human_avian_ratio)
    avian_odds <- (b + human_avian_ratio) / (d + human_avian_ratio)
  } else {
    human_odds <- ifelse(c0 == 0, Inf, a / c0)
    avian_odds <- ifelse(d == 0, Inf, b / d)
  }
  or <- human_odds / avian_odds
  # Per-position Fisher's p-values and 95% CI for OR
  pval <- numeric(length(a))
  ci_low <- numeric(length(a))
  ci_high <- numeric(length(a))
  for (j in seq_along(a)) {
    # rows: human, avian ; cols: present(1), absent(0)
    tbl <- matrix(c(a[j], c0[j], b[j], d[j]), nrow = 2, byrow = TRUE)
    ft <- fisher.test(tbl, alternative = "two.sided", conf.level = conf.level)
    pval[j] <- ft$p.value
    ci_low[j] <- ft$conf.int[1]
    ci_high[j] <- ft$conf.int[2]
  }
  # construting the dataframe
  by_position <- data.frame(
    position = which(keep),
    human_ones = a,
    human_zeros = c0,
    avian_ones = b,
    avian_zeros = d,
    human_odds = human_odds,
    avian_odds = avian_odds,
    odds_ratio = or,
    or_ci_low = ci_low,
    or_ci_high = ci_high,
    p_value = pval,
    row.names = NULL, check.names = FALSE
  )
  list(by_position = by_position)
}

# Generalized CSV reader
read_feature_csvs <- function(in_dir, pattern = "^binary_.*\\.csv$") {
  csvs <- list.files(in_dir, pattern = pattern, full.names = TRUE)
  if (!length(csvs)) {
    message("Skip (no matching CSVs): ", in_dir)
    return(NULL)
  }
  lapply(csvs, function(f) {
    feature <- sub("^binary_(.*)\\.csv$", "\\1", basename(f))
    df <- read.csv(f, row.names = 1, check.names = FALSE)
    list(feature = feature, df = df)
  })
}

# process one protein function
process_protein <- function(
    protein,
    input_base,
    output_base,
    human_pat = "HUMAN",
    avian_pat = "AVIAN",
    haldane = TRUE,
    drop_all_zero_cols = TRUE,
    conf.level = 0.95) {
  in_dir <- file.path(input_base, protein)
  out_dir <- file.path(output_base, protein)
  if (!dir.exists(in_dir)) {
    message("Skip (no input dir): ", in_dir)
    return(invisible(NULL))
  }
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  feature_list <- read_feature_csvs(in_dir, pattern = "^binary_.*\\.csv$")
  if (is.null(feature_list)) {
    return(invisible(NULL))
  }
  res_list <- lapply(feature_list, function(item) {
    r <- or_function(
      binary_df = item$df,
      human_pat = human_pat,
      avian_pat = avian_pat,
      haldane = haldane,
      drop_all_zero_cols = drop_all_zero_cols,
      conf.level = conf.level
    )
    if (!is.null(r$by_position) && nrow(r$by_position) > 0) {
      r$by_position$kmer <- item$feature
    }
    r
  })
  bypos_tbl <- do.call(rbind, lapply(res_list, `[[`, "by_position"))
  # FDR across per-position tests for this protein
  bypos_tbl$p_adj_bh <- p.adjust(bypos_tbl$p_value, method = "BH")
  # Save table
  write.csv(
    bypos_tbl[, c(
      "kmer", "position", "human_ones", "human_zeros", "avian_ones", "avian_zeros",
      "human_odds", "avian_odds", "odds_ratio", "or_ci_low", "or_ci_high",
      "p_value", "p_adj_bh"
    )],
    file.path(out_dir, "odds_ratio_by_position_with_p.csv"),
    row.names = FALSE
  )
  message("Done: ", protein)
}

# Loop all proteins
proteins <- c("HA", "M1", "NA", "NP", "NS1", "PA", "PB1", "PB2")
invisible(lapply(proteins, process_protein,
                 input_base = input_base,
                 output_base = output_base
))

################################################################################
# Uncomment this block and replace the fixed protein loop above,
# if you want to run the pipeline across all k-mer/protein composition folders
# directory should contain folders like "2mers", "3mers", "4mers", etc.
###############################################################################
# feature_base <- file.path()  # add correct path
# feature_folders <- list.dirs(feature_base, full.names = TRUE, recursive = FALSE)
#
# for (feature_path in feature_folders) {
#  feature_name <- basename(feature_path)
#  input_base  <- file.path(feature_path, "binary_result", "nuc")
#  output_base <- file.path(feature_path, "kmer_or_results", "nuc")
#
#  protein_folders <- list.dirs(input_base, full.names = FALSE, recursive = FALSE)
#
#  for (protein in protein_folders) {
#    process_protein(protein,
#                    input_base  = input_base,
#                    output_base = output_base)
#  }
# }
