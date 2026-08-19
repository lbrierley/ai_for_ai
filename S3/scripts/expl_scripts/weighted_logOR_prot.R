######################################################################################
# this script calculate weighted mean (log OR), and plots heatmap for protein seqs
######################################################################################

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
feat <- "all_prot"

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
  mutate(importance = ifelse(importance < 0, 0, importance)) %>% # Set importance to zero if negative (implies that feature for that protein was not at all informative in final model)
  mutate(w_all = importance/sum(importance)) # Rescale so global weights sum to 1

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
    # Create shell for every feature-position combination
    shell <- expand.grid(kmer = read.csv("H:\\Working\\ai_for_ai\\S3/data/full/mapping/prot_feats_ref.csv") %>% pull(x), # features
                         position = 1:max(bypos$position)) # positions
    bypos <- bypos %>% right_join(shell, by = c("kmer", "position"))
  } else {
    
    ### NOTE THIS DOES NOT INCLUDE ZERO-FEATURES AND ZERO-POSITIONS YET - COULD DO BY LABELLING THEM BY FEATURE_TYPE AND SELECTING THEM BASED ON FEATURE_TYPE DYNAMICALLY
    
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
    complete(kmer, position, fill = list(odds_ratio = 1)) %>% # ensure positions where no odds ratio difference still contribute to weighted average
    mutate(abs_logOR = abs(log(odds_ratio))) %>%
    filter(!is.na(position))
  pos_unweighted <- bypos %>%
    group_by(position) %>%
    summarise(
      mean_abs_logOR = mean(abs_logOR, na.rm = TRUE),
      .groups = "drop"
    )
  # Join weight to OR table and compute weighted mean log OR
  bypos_w <- bypos %>%
    left_join(varimp_weights, by = c("kmer" = "feature")) %>%
    mutate(importance = replace_na(importance, 0),
           w_all = replace_na(w_all, 0)) %>%  # Set importance/global weight to zero if NA (implies that feature set for that protein was not retained in final model)
    # WEIGHTING FUNCTION GOES HERE
    mutate(w_loc =  importance/sum(importance)) # Rescale so global weights sum to 1
  pos_weighted <- bypos_w %>%
    group_by(position) %>%
    summarise(
      w_loc_mean_abs_logOR = sum(w_loc * abs_logOR), 
      w_all_mean_abs_logOR = sum(w_all * abs_logOR),
      .groups = "drop"
    ) %>%
    filter(!is.na(w_loc_mean_abs_logOR))
  if (feature_type == "all_prot") {
    # save both unweighted and weigthed results
    write_csv(pos_unweighted, file.path(paste0(output_path, "/", protein, "_position_mean_abs_logOR.csv")))
    write_csv(pos_weighted, file.path(paste0(output_path, "/", protein, "_position_weighted_mean_abs_logOR.csv")))
  } else {
    # save both unweighted and weigthed results
    write_csv(pos_unweighted, file.path(paste0(output_path, "/", protein, "/", protein, "_position_mean_abs_logOR.csv")))
    write_csv(pos_weighted, file.path(paste0(output_path, "/", protein, "/", protein, "_position_weighted_mean_abs_logOR.csv")))
  }
  message("Completed weight calc: ", protein, " [Feature: ", feature_type, "]")
}

# function for plot
plot_weighted_logOR <- function(
    protein,
    output_path,
    feature_type,
    vmin = 0.0,
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
      plot_value = w_all_mean_abs_logOR,
    ) %>%
    dplyr::arrange(position)
  if (nrow(df) == 0) {
    message(paste("Feature", feature_type, "for protein", protein, "not in final stacks, passing plot"))
    return(invisible(NULL))
  }
  # Fill missing positions with vmin
  all_pos <- seq(min(df$position, na.rm = TRUE), max(df$position, na.rm = TRUE))
  df <- full_join(tibble(position = all_pos), df, by = "position") %>%
    mutate(
      plot_value = ifelse(is.na(plot_value), vmin, plot_value)
    )
  # Generate x-axis ticks
  xticks <- scales::pretty_breaks(n = 10)(range(df$position, na.rm = TRUE))
  # Build plot
  p <- ggplot(df, aes(x = position, y = 1, fill = plot_value)) +
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

# function for plot all proteins
plot_weighted_logOR_allprot <- function(
    output_path,
    feature_type,
    vmin = 0.0) {
  prot_dfs <- list()
  # Construct input list
  for (i in 1:length(proteins)){
    # Construct input CSV path
    if (feature_type == "all_prot") {
      prot_dfs[[i]] <- read.csv(file.path(paste0(output_path, "/", proteins[i], "_position_weighted_mean_abs_logOR.csv"))) %>%
        mutate(prot = proteins[i])
    } else {
      csv_path <- file.path(paste0(output_path, "/", proteins[i], "/", proteins[i], "_position_weighted_mean_abs_logOR.csv")) %>%
        mutate(prot = proteins[i])
    }
    
  }
  df <- bind_rows(prot_dfs) %>%
    dplyr::mutate(
      position = as.integer(position),
      plot_value = w_all_mean_abs_logOR,
    ) %>%
    dplyr::arrange(position)
  # Fill missing positions with vmin
  all_pos <- seq(min(df$position, na.rm = TRUE), max(df$position, na.rm = TRUE))
  df <- full_join(tibble(position = all_pos), df, by = "position") %>%
    mutate(
      plot_value = ifelse(is.na(plot_value), vmin, plot_value)
    )
  # Build plot
  p <- ggplot(df, aes(x = position, y = 1, fill = plot_value)) +
    geom_tile(height = 1) +
    scale_fill_viridis_c(
      begin = 0,
      end = 1,
      # limits = c(0, 1),
      # oob = scales::squish,
      option = "inferno",
      name = "importance index (weighted mean abs log-OR)"
    ) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(breaks = NULL, labels = NULL, expand = c(0, 0)) +
    labs(
      x = "Residue", y = NULL,
    ) +
    facet_wrap(vars(prot), scales = "free", nrow = 8) +
    theme_minimal(base_size = 8) +
    theme(
      panel.grid = element_blank(),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background = element_rect(fill = "white", colour = NA),
      legend.background = element_rect(fill = "white", colour = NA),
      legend.title = element_text(angle = 90, size = 6),
      axis.title.x = element_text(face = "bold"),
      plot.margin = margin(6, 10, 6, 6),
      legend.position = "left"
    ) +
    guides(
      fill = guide_colorbar(
        title.position = "left",
        barheight = unit(10, "cm"),
        barwidth = unit(6, "pt"),
        ticks = TRUE
      )
    )
  # Save plot
  if (feature_type == "all_prot") {
    out_path <- file.path(paste0(output_path, "/weighted_mean_heatmap_", feature_type, "_all.png"))
  } else {
    out_path <- file.path(paste0(output_path, "/", protein, "/weighted_mean_heatmap_", feature_type, "_all.png"))
  }
  ggsave(filename = out_path, plot = p, width = 5, height = 5.5, dpi = 300, bg = "white")
  invisible(p)
}

# function for plot all proteins (unweighted, for reference only)
plot_unweighted_logOR_allprot <- function(
    output_path,
    feature_type,
    vmin = 0.0) {
  prot_dfs <- list()
  # Construct input list
  for (i in 1:length(proteins)){
    # Construct input CSV path
    if (feature_type == "all_prot") {
      prot_dfs[[i]] <- read.csv(file.path(paste0(output_path, "/", proteins[i], "_position_mean_abs_logOR.csv"))) %>%
        mutate(prot = proteins[i])
    } else {
      csv_path <- file.path(paste0(output_path, "/", proteins[i], "/", proteins[i], "_position_mean_abs_logOR.csv")) %>%
        mutate(prot = proteins[i])
    }
    
  }
  df <- bind_rows(prot_dfs) %>%
    dplyr::mutate(
      position = as.integer(position),
      plot_value = mean_abs_logOR,
    ) %>%
    dplyr::arrange(position)
  # Fill missing positions with vmin
  all_pos <- seq(min(df$position, na.rm = TRUE), max(df$position, na.rm = TRUE))
  df <- full_join(tibble(position = all_pos), df, by = "position") %>%
    mutate(
      plot_value = ifelse(is.na(plot_value), vmin, plot_value)
    )
  # Build plot
  p <- ggplot(df, aes(x = position, y = 1, fill = plot_value)) +
    geom_tile(height = 1) +
    scale_fill_viridis_c(
      begin = 0,
      end = 1,
      # limits = c(0, 1),
      # oob = scales::squish,
      option = "inferno",
      name = "importance index (mean abs log-OR)"
    ) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(breaks = NULL, labels = NULL, expand = c(0, 0)) +
    labs(
      x = "Residue", y = NULL,
    ) +
    facet_wrap(vars(prot), scales = "free", nrow = 8) +
    theme_minimal(base_size = 8) +
    theme(
      panel.grid = element_blank(),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background = element_rect(fill = "white", colour = NA),
      legend.background = element_rect(fill = "white", colour = NA),
      legend.title = element_text(angle = 90, size = 6),
      axis.title.x = element_text(face = "bold"),
      plot.margin = margin(6, 10, 6, 6),
      legend.position = "left"
    ) +
    guides(
      fill = guide_colorbar(
        title.position = "left",
        barheight = unit(10, "cm"),
        barwidth = unit(6, "pt"),
        ticks = TRUE
      )
    )
  # Save plot
  if (feature_type == "all_prot") {
    out_path <- file.path(paste0(output_path, "/unweighted_mean_heatmap_", feature_type, "_all.png"))
  } else {
    out_path <- file.path(paste0(output_path, "/", protein, "/unweighted_mean_heatmap_", feature_type, "_all.png"))
  }
  ggsave(filename = out_path, plot = p, width = 5, height = 5.5, dpi = 300, bg = "white")
  invisible(p)
}

# Loop through all proteins folders
for (protein in proteins) {
  process_weighted_logOR(protein, feat)
  plot_weighted_logOR(protein, output_path, feat)
}

# Plot all in one
plot_weighted_logOR_allprot(output_path, feat)
plot_unweighted_logOR_allprot(output_path, feat)
