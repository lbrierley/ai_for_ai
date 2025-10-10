################################################################################
# this script calculate weighted mean (log OR), and plots heatmap for proteins
################################################################################

# Load libraries
library(readr)
library(dplyr)
library(stringr)
library(tibble)
library(ggplot2)
library(ggrepel)
library(scales)
library(tidyverse)


# Define feature (all_prot for all protein features)
feat <- "allprot"

proteins <- c("HA", "M1", "NA", "NP", "NS1", "PA", "PB1", "PB2")
# Define paths
project_root <- normalizePath("./S3/")
varimp_dir <- file.path(project_root, "analysis")
if (feat == "all_prot") {
  output_path <- file.path(project_root, "data/full/mapping/kmer_or_results/prot/")
} else {
  output_path <- file.path(project_root, paste0("data/full/mapping/kmer_or_results/prot/", feat))
}

################################################################################
# Read in variable importance outside loop so do this once
# Function to read and clean varimp files
read_varimp_file <- function(path) {
  read_csv(path, show_col_types = FALSE) %>%
    select(feature = last_col(1), importance = last_col()) %>%
    mutate(
      feature = trimws(gsub('^["\\\']+|["\\\']+$', "", feature)),
      importance = as.numeric(importance)
    ) %>%
    filter(!is.na(importance) & nzchar(feature))
}

varimp_paths <- list.files(varimp_dir, pattern = "^varimp_perm_weight_pt.*\\.csv$", full.names = TRUE)
if (!length(varimp_paths)) {
  message("No varimp files found in: ", varimp_dir)
  return(invisible(NULL))
}

# combining all varimp csv and cleaning
varimp_weights <- varimp_paths %>%
  lapply(read_varimp_file) %>%
  bind_rows() %>%
  filter(!is.na(importance) & !is.na(feature)) %>%
  mutate(w = LDATS::softmax(importance)) # Softmax function for weights (new weights sum to 1)
################################################################################
# Function to process weighted log(OR) for one protein
process_weighted_logOR <- function(protein, feature_type) {
  if (feature_type == "all_prot") {
    feat_ORs <- list.files("H:\\Working\\ai_for_ai\\S3/data/full/mapping/kmer_or_results/prot",
                           pattern = "odds_ratio_by_position_with_p.csv",
                           recursive = TRUE,
                           full.names = TRUE
    ) %>%
      .[grepl(protein, .)]
    # Clean and transform OR value
    bypos <- feat_ORs %>% purrr::map_dfr(read.csv)
  } else {
    bypos_csv <- file.path(output_path, protein, "odds_ratio_by_position_with_p.csv")
    if (!file.exists(bypos_csv)) {
      message("Missing OR file for: ", protein)
      return(invisible(NULL))
    }
    bypos <- read_csv(bypos_csv, show_col_types = FALSE)
  }
  bypos <- bypos %>%
    select(kmer, position, odds_ratio) %>%
    mutate(
      position = as.integer(position),
      kmer = paste0(trimws(kmer), "_", protein),
      odds_ratio = suppressWarnings(as.numeric(odds_ratio)),
      odds_ratio = if_else(is.na(odds_ratio) | !is.finite(odds_ratio) | odds_ratio <= 0,
                           NA_real_, odds_ratio
      )
    ) %>%
    complete(kmer, position, fill = list(odds_ratio = 1)) %>% # ensure uninformative traits contribute to weighted average
    mutate(abs_logOR = abs(log(odds_ratio)))
  pos_unweighted <- bypos %>%
    group_by(position) %>%
    summarise(
      mean_abs_logOR = mean(abs_logOR, na.rm = TRUE),
      .groups = "drop"
    )
  # Join weight to OR table and compute weighted mean log OR
  bypos_w <- bypos %>%
    left_join(varimp_weights, by = c("kmer" = "feature")) %>%
    filter(is.finite(abs_logOR))
  pos_weighted <- bypos_w %>%
    group_by(position) %>%
    summarise(
      w_sum = sum(w),
      w_mean_abs_logOR = if_else(w_sum > 0, sum(w * abs_logOR) / w_sum, NA_real_),
      .groups = "drop"
    ) %>%
    filter(!is.na(w_mean_abs_logOR)) %>%
    select(-w_sum)
  if (feature_type == "all_prot") {
    # save both unweighted and weigthed results
    write_csv(pos_unweighted, file.path(paste0(output_path, "/", protein, "_position_mean_abs_logOR.csv")))
    write_csv(pos_weighted, file.path(paste0(output_path, "/", protein, "_position_weighted_mean_abs_logOR.csv")))
  } else {
    # save both unweighted and weigthed results
    write_csv(pos_unweighted, file.path(paste0(output_path, "/", protein, "/", protein, "_position_mean_abs_logOR.csv")))
    write_csv(pos_weighted, file.path(paste0(output_path, "/", protein, "/", protein, "_position_weighted_mean_abs_logOR.csv")))
  }
  message("Completed: ", protein, " [Feature: ", feature_type, "]")
}

# function for plot
plot_weighted_logOR <- function(
    protein,
    output_path,
    feature_type,
    vmin = 0.0, vmax = 3.0, center = 1.5,
    fig_width = 12, fig_height = 2, legend_height_cm = 4) {
  # Construct input CSV path
  if (feature_type == "all_prot") {
    csv_path <- file.path(paste0(output_path, "/", protein, "_position_weighted_mean_abs_logOR.csv"))
  } else {
    csv_path <- file.path(paste0(output_path, "/", protein, "/", protein, "_position_weighted_mean_abs_logOR.csv"))
  }
  if (!file.exists(csv_path)) {
    message("Missing file for: ", protein)
    return(invisible(NULL))
  }
  # Read and clean data
  df <- readr::read_csv(csv_path, show_col_types = FALSE) %>%
    dplyr::mutate(
      position = as.integer(position),
      w_mean_abs_logOR = pmin(pmax(w_mean_abs_logOR, vmin), vmax)
    ) %>%
    dplyr::arrange(position)
  # Fill missing positions with vmin
  all_pos <- seq(min(df$position, na.rm = TRUE), max(df$position, na.rm = TRUE))
  df <- full_join(tibble(position = all_pos), df, by = "position") %>%
    mutate(
      w_mean_abs_logOR = pmin(pmax(w_mean_abs_logOR, vmin), vmax),
      w_mean_abs_logOR = ifelse(is.na(w_mean_abs_logOR), vmin, w_mean_abs_logOR)
    )
  # Generate x-axis ticks
  xticks <- scales::pretty_breaks(n = 10)(range(df$position, na.rm = TRUE))
  # Build plot
  p <- ggplot(df, aes(x = position, y = 1, fill = w_mean_abs_logOR)) +
    geom_tile(height = 1) +
    scale_fill_viridis_c(
      begin = 0,
      end = 1,
      # limits = c(0, 1),
      # oob = scales::squish,
      option = "inferno",
      name = "importance index (weighted mean abs log-OR)"
    ) +
    scale_x_continuous(breaks = xticks, expand = c(0, 0)) +
    scale_y_continuous(breaks = NULL, labels = NULL, expand = c(0, 0)) +
    labs(
      x = "Residue", y = NULL,
      title = paste0(protein, " — importance [", feature_type, "]")
    ) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background = element_rect(fill = "white", colour = NA),
      legend.background = element_rect(fill = "white", colour = NA),
      legend.title = element_text(angle = 90, size = 5),
      axis.title.x = element_text(face = "bold", size = 13),
      axis.text.x = element_text(margin = margin(t = 4)),
      plot.margin = margin(6, 10, 6, 6),
      legend.position = "left"
    ) +
    guides(
      fill = guide_colorbar(
        title.position = "left",
        barheight = unit(legend_height_cm, "cm"),
        barwidth = unit(6, "pt"),
        ticks = TRUE
      )
    )
  # Save plot
  if (feature_type == "all_prot") {
    out_path <- file.path(paste0(output_path, "/weighted_mean_heatmap_", feature_type, "_", protein, ".png"))
  } else {
    out_path <- file.path(paste0(output_path, "/", protein, "/weighted_mean_heatmap_", feature_type, "_", protein, ".png"))
  }
  ggsave(filename = out_path, plot = p, width = fig_width, height = fig_height, dpi = 300, bg = "white")
  invisible(p)
}

# Loop through all proteins folders
for (protein in proteins) {
  process_weighted_logOR(protein, feat)
  plot_weighted_logOR(protein, output_path, feat)
}