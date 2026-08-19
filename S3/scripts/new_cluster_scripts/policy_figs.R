###############################################
# Overall figures across all model algorithms #
###############################################

#################
# Load packages #
#################

library(dplyr)
library(tidyr)
library(forcats)
library(magrittr)
library(purrr)
library(ggplot2)
library(ggdist)
library(tidytext)

############################
# Load IDs and set options #
############################

results_date <- "17_07_26"
results_path <- "//filepor10/DOP$/X037_AVI_GeneticMarkers/S3/analysis/"

allflu_wgs_ref <- read.csv("//filepor10/DOP$/X037_AVI_GeneticMarkers/S3/data/Segmented/Raw_Sequences/allflu_wgs_df.csv", na.strings = "")
allflu_nuc_ref <- read.csv("//filepor10/DOP$/X037_AVI_GeneticMarkers/S3/data/Segmented/Raw_Sequences/allflu_nuc_df.csv", na.strings = "")
allflu_cds_ref <- read.csv("//filepor10/DOP$/X037_AVI_GeneticMarkers/S3/data/Segmented/Raw_Sequences/allflu_cds_df.csv", na.strings = "")

cbbPalette_ordered <- c("#D55E00", "#E69F00", "#F0E442", "#009E73", "#56B4E9", "#0072B2", "#CC79A7", "#999999")

cbbPalette_12 <- c("#D55E00", "#E69F00", "#F0E442", "#009E73", "#56B4E9", "white", "#0072B2","#DDDDDD", "#999999", "black", "#CC79A7", "#5D3A9B")

##########################################
# Performance figures: individual models #
##########################################

# Performance metrics on holdout sets

results_list <- list.files(path = results_path, pattern = paste0("test_results_", results_date), full.names = TRUE)
all_res <- map_df(results_list, read.csv, na.strings = "")

# Fig 2
# AUC heatmap of best model methods

fig_results_heat_AUC <- all_res %>%
  group_by(method, featset, focgene) %>%
  summarise(AUC = mean(AUC)) %>% 
  ungroup %>%
  mutate(featset = case_when(
    featset == "cds_compbias" ~ "nuc: composition",
    featset == "nuc_2mer" ~ "nuc: 2-mers",
    featset == "nuc_3mer" ~ "nuc: 3-mers",
    featset == "nuc_4mer" ~ "nuc: 4-mers",
    featset == "nuc_5mer" ~ "nuc: 5-mers",
    featset == "nuc_6mer" ~ "nuc: 6-mers",
    featset == "prot_2mer" ~ "prot: 2-mers",
    featset == "prot_ctdc" ~ "prot: CTD-C",
    featset == "prot_ctdt" ~ "prot: CTD-T",
    featset == "prot_ctdd" ~ "prot: CTD-D",
    featset == "prot_ctriad" ~ "prot: CTriad",
    featset == "prot_pseaac" ~ "prot: PseAAC"),
    featset = as.factor(featset),
    method = case_when(
      method == "glmnet" ~ "PLR",
      method == "rf" ~ "RF",
      method == "svm" ~ "RSVM",
      method == "svmlin" ~ "SVM",
      method == "xgb" ~ "XGB",
    ),
    focgene = factor(focgene, levels = c("PB2", "PB1", "PA", "HA", "NP", "NA", "M1", "NS1")),
  ) %>%
  gather(metric, value, -method, -featset, -focgene) %>%
  filter(metric %in% c("AUC")) %>%
  mutate(value = as.numeric(value)) %>%
  ggplot(aes(x = focgene, y = factor(featset, levels = rev(levels(featset))), fill = value)) + 
  geom_tile(color="white") +
  scale_fill_distiller("AUROC", palette = "RdBu", limits = c(0,1)) +
  facet_wrap(~ method, nrow=1) +
  theme_bw() +
  xlab("Influenza virus gene/protein") +
  ylab("Feature set")

ggsave("//filepor10/DOP$/X037_AVI_GeneticMarkers/S3/figures_tables/all_results_heat_AUC_segmentwise.png", plot = fig_results_heat_AUC, width = 15, height = 3)

# Supplementary Fig 2
# AUC heatmap of all model methods

fig_results_heat_AUC_one <- all_res %>%
  group_by(method, featset, focgene) %>%
  summarise(AUC = mean(AUC)) %>% 
  ungroup %>%
  group_by(featset, focgene) %>%
  slice_max(AUC, with_ties = FALSE) %>%
  ungroup %>%
  mutate(featset = case_when(
    featset == "cds_compbias" ~ "nuc: composition",
    featset == "nuc_2mer" ~ "nuc: 2-mers",
    featset == "nuc_3mer" ~ "nuc: 3-mers",
    featset == "nuc_4mer" ~ "nuc: 4-mers",
    featset == "nuc_5mer" ~ "nuc: 5-mers",
    featset == "nuc_6mer" ~ "nuc: 6-mers",
    featset == "prot_2mer" ~ "prot: 2-mers",
    featset == "prot_ctdc" ~ "prot: CTD-C",
    featset == "prot_ctdt" ~ "prot: CTD-T",
    featset == "prot_ctdd" ~ "prot: CTD-D",
    featset == "prot_ctriad" ~ "prot: CTriad",
    featset == "prot_pseaac" ~ "prot: PseAAC"),
    featset = as.factor(featset),
    method = case_when(
      method == "glmnet" ~ "PLR",
      method == "rf" ~ "RF",
      method == "svm" ~ "RSVM",
      method == "svmlin" ~ "SVM",
      method == "xgb" ~ "XGB",
    ),
    focgene = factor(focgene, levels = c("PB2", "PB1", "PA", "HA", "NP", "NA", "M1", "NS1"))
  ) %>%
  gather(metric, value, -method, -featset, -focgene) %>%
  filter(metric %in% c("AUC")) %>%
  mutate(value = as.numeric(value)) %>%
  ggplot(aes(x = focgene, y = factor(featset, levels = rev(levels(featset))), fill = value)) + 
  geom_tile(color="white") +
  geom_point(aes(pch = method), color = "black", size = 2) +
  scale_fill_distiller("AUROC", palette = "RdBu", limits = c(0,1)) +
  scale_shape_manual("algorithm", values = c(19, 15, 17, 18, 8)) + 
  theme_bw() +
  xlab("Influenza virus gene/protein") +
  ylab("Feature set")

ggsave("//filepor10/DOP$/X037_AVI_GeneticMarkers/S3/figures_tables/all_results_heat_AUC_segmentwise_one.png", plot = fig_results_heat_AUC_one, width = 6, height = 3.5)

# Test plot Policy Piece

fig_results_bar_AUC <- all_res %>%
  group_by(method, featset, focgene) %>%
  summarise(AUC = mean(AUC)) %>% 
  ungroup %>%
  mutate(featset = case_when(
    featset == "cds_compbias" ~ "nuc: composition",
    featset == "nuc_2mer" ~ "nuc: 2-mers",
    featset == "nuc_3mer" ~ "nuc: 3-mers",
    featset == "nuc_4mer" ~ "nuc: 4-mers",
    featset == "nuc_5mer" ~ "nuc: 5-mers",
    featset == "nuc_6mer" ~ "nuc: 6-mers",
    featset == "prot_2mer" ~ "prot: 2-mers",
    featset == "prot_ctdc" ~ "prot: CTD-C",
    featset == "prot_ctdt" ~ "prot: CTD-T",
    featset == "prot_ctdd" ~ "prot: CTD-D",
    featset == "prot_ctriad" ~ "prot: CTriad",
    featset == "prot_pseaac" ~ "prot: PseAAC"),
    featset = as.factor(featset),
    method = case_when(
      method == "glmnet" ~ "PLR",
      method == "rf" ~ "RF",
      method == "svm" ~ "RSVM",
      method == "svmlin" ~ "SVM",
      method == "xgb" ~ "XGB",
    ),
    focgene = factor(focgene, levels = c("PB2", "PB1", "PA", "HA", "NP", "NA", "M1", "NS1")),
    label = paste0(method, ", ", featset),
  ) %>%
  gather(metric, value, -method, -featset, -focgene, -label) %>%
  filter(metric %in% c("AUC")) %>%
  mutate(value = as.numeric(value)) %>%
  mutate(focgene = fct_reorder(focgene, value, .fun = max, .desc = TRUE)) %>%
  group_by(focgene) %>%                            
  slice_max(order_by = value, n = 5) %>%  # Only plot top 5 models for each segment
  ggplot(aes(y = reorder_within(label, value, focgene), x = value, fill = featset)) + 
  geom_bar(alpha = 0.9, color = "black", stat = "identity", show.legend = TRUE) +
  stat_identity(aes(pch = method, x = value + 0.05), color = "black", size = 4) +
  scale_shape_manual(values = c(19, 15, 17, 18, 8), guide = "none") + 
  scale_fill_manual(values = cbbPalette_12) +
  scale_y_reordered() +
  facet_wrap(~ focgene, ncol=2, scales = "free_y") +
  theme_bw() +
  guides(fill = guide_legend(override.aes = list(shape = NA))) +
  ylab("Feature set x Model") +
  xlab("AUC")

ggsave("//filepor10/DOP$/X037_AVI_GeneticMarkers/S3/figures_tables/all_results_bar_AUC_segmentwise.png", plot = fig_results_bar_AUC, width = 15, height = 7)

# Identify best models per segment

best_mods <- all_res %>%
  group_by(method, featset, focgene) %>%
  summarise(AUC = mean(AUC), threshold = mean(threshold)) %>% 
  ungroup %>%
  group_by(focgene) %>%
  slice_max(AUC, with_ties = FALSE) %>%
  ungroup

all_pred_raw <- pmap_dfr(
  best_mods[, c("method", "featset", "focgene", "AUC", "threshold")],
  function(method, featset, focgene, AUC, threshold) {
    read.csv(paste0(results_path, "pred_raw/test_pred_raw_", method, "_", featset, "_", focgene, "_", results_date, ".csv"), na.strings = "") %>% 
      mutate(focgene = focgene, AUC = AUC, threshold)
  }
) %>%
  select(focgene, AUC, threshold, gid, hzoon, nz, label, pred) %>%
  left_join(allflu_cds_ref %>% select(gid, subtype) %>% distinct) %>%
  mutate(focgene = fct_reorder(focgene, AUC, .fun = max, .desc = TRUE)) %>%
  mutate(label = as.character(label)) %>%
  mutate(label = case_when(label == "hzoon" ~ "Zoonotic", 
                           label == "nz" ~ "Avian") %>% as.factor)

fig_results_const_raw <- all_pred_raw %>% 
  ggplot(aes(x = label, y = hzoon, fill = label)) +
  geom_dots(data = . %>% filter(label == "Avian"), binwidth = unit(1.5, "mm"),color = "gray30", dotsize = 2, alpha = 0.6, side = "both", overflow = "compress") +
  geom_dots(data = . %>% filter(label == "Zoonotic"), binwidth = unit(1.5, "mm"), color = "gray30", dotsize = 2, alpha = 0.6, side = "both", overflow = "compress") +
  geom_hline(aes(yintercept = threshold), linetype = "longdash", color = "gray10",linewidth = 1.2, lwd = 1.2) +
  theme_bw() +
  scale_fill_discrete(direction = -1) +
  scale_color_discrete(direction = -1) +
  guides(fill = "none") +
  facet_wrap(~ focgene, ncol=2, scales = "free_y") +
  ylab("p(zoonotic)") +
  xlab(NULL) +
  guides(color = "none")

ggsave("//filepor10/DOP$/X037_AVI_GeneticMarkers/S3/figures_tables/all_results_const_raw.png", plot = fig_results_const_raw, width = 10, height = 6)

all_pred_raw %>% select(focgene, gid, subtype, label, hzoon) %>% arrange(focgene, desc(hzoon)) %>% write.csv("//filepor10/DOP$/X037_AVI_GeneticMarkers/S3/figures_tables/table_all_pred_raw_best_const_models.csv")