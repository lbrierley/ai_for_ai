###############################################
# Overall figures across all model algorithms #
###############################################

#################
# Load packages #
#################

library(dplyr)
library(tidyr)
library(magrittr)
library(ggplot2)
library(ggrepel)
library(rentrez)
library(patchwork)

############################
# Load IDs and set options #
############################

allflu_wgs_ref <- read.csv("S3\\data\\full\\allflu_wgs_ref.csv")

cbbPalette_ordered <- c("#D55E00", "#E69F00", "#F0E442", "#009E73", "#56B4E9", "#0072B2", "#CC79A7", "#999999")

cluster_chosen <- "70_7"

################
# Data figures #
################

# Supplementary Fig S1
# Sequence representation over time

S1 <- allflu_wgs_ref %>%
  mutate(label = case_when(
    label == "nz" ~ "avian",
    label == "zoon" ~ "zoonotic",
  )) %>%  
  filter(!is.na(date) & date > as.Date("1990-01-01")) %>%
  add_count(subtype, name = "sub_n") %>%
  filter(sub_n > 500) %>%
  ggplot(aes(x = as.Date(date), fill = subtype)) +
  geom_histogram(position = "stack", binwidth=365) +
  scale_fill_manual(values = rev(c(RColorBrewer::brewer.pal(10, "Paired")))) +
  scale_x_date(limits = c(as.Date("1990-01-01"), as.Date("2022-12-31")), date_labels =  "%Y") +
  facet_grid(rows = vars(label), cols = vars(src), scales = "free_y") +
  theme_bw() +
  xlab("Date") +
  ylab("Frequency")

ggsave("S3\\figures_tables\\time_dist_wgs.png", plot = S1, width = 14, height = 6)

##########################################
# Performance figures: individual models #
##########################################

# Performance metrics on holdout sets

results_rf <- read.csv(paste0("S3\\analysis\\results_", "14_02_24", ".csv"), na.strings = "NaN") %>% filter(cluster_set == cluster_chosen) %>% mutate(method = "rf") 
results_plr <- read.csv(paste0("S3\\analysis\\results_", "15_02_24", ".csv"), na.strings = "NaN") %>% filter(cluster_set == cluster_chosen) %>% mutate(method = "glmnet")
results_xgb <- read.csv(paste0("S3\\analysis\\results_", "16_02_24", ".csv"), na.strings = "NaN") %>% filter(cluster_set == cluster_chosen) %>% mutate(method = "xgb")
results_svmlin <- read.csv(paste0("S3\\analysis\\results_", "17_02_24", ".csv"), na.strings = "NaN") %>% filter(cluster_set == cluster_chosen) %>% mutate(method = "svmlin")
results_svmrad <- read.csv(paste0("S3\\analysis\\results_", "18_02_24", ".csv"), na.strings = "NaN") %>% filter(cluster_set == cluster_chosen) %>% mutate(method = "svm")

all_res <- bind_rows(results_rf,
                     results_svmlin,
                     results_svmrad,
                     results_xgb,
                     results_plr) 

# Save for reference
all_res %>% write.csv("S3\\analysis\\results_all_methods.csv")

# Fig 2
# AUC heatmap of best model methods

fig_results_heat_AUC_70_7 <- all_res %>%
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

ggsave("S3\\figures_tables\\all_results_heat_AUC_70_7.png", plot = fig_results_heat_AUC_70_7, width = 15, height = 3)

# Supplementary Fig 2
# AUC heatmap of all model methods

fig_results_heat_AUC_70_7_one <- all_res %>%
  group_by(method, featset, focgene) %>%
  summarise(AUC = mean(AUC)) %>% 
  ungroup %>%
  group_by(featset, focgene) %>%
  slice_max(AUC) %>%
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

ggsave("S3\\figures_tables\\all_results_heat_AUC_70_7_one.png", plot = fig_results_heat_AUC_70_7_one, width = 6, height = 3.5)

# Supplementary Fig 3
# F1 heatmap of all model methods

fig_results_heat_F1_70_7 <- all_res %>%
  group_by(method, featset, focgene) %>%
  summarise(F1 = mean(F1)) %>%
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
  filter(metric %in% c("F1")) %>%
  mutate(value = as.numeric(value)) %>%
  ggplot(aes(x = focgene, y = factor(featset, levels = rev(levels(featset))), fill = value)) +
  geom_tile(color="white") +
  scale_fill_viridis_c("F1") +
  facet_wrap(~ method, nrow=1) +
  theme_bw() +
  # theme(plot.background = element_rect(fill = "#F2F6F9", color = "#F2F6F9"),
  #       legend.key = element_rect(fill = "#F2F6F9", color = "#F2F6F9")) +
  xlab("Influenza virus gene/protein") +
  ylab("Feature set")

ggsave("S3\\figures_tables\\all_results_heat_F1_70_7.png", plot = fig_results_heat_F1_70_7, width = 15, height = 3)

#######################################
# Performance figures: stacked models #
#######################################

# Read in stacked model coefficients (on individual models)

stacked_coef <- read.csv("S3/analysis/stack_weight_coef.csv") %>%
  select(-X) %>%
  filter(param != "(Intercept)" & param != "lambda")

# Arrange individual models by count of retention

mod_count <- stacked_coef %>%
  group_by(param) %>%
  tally() %>%
  arrange(-n) 

# Arrange individual models by coefficient magnitude

mod_mag <- stacked_coef %>% 
  tidyr::expand(param, subtype) %>% 
  left_join(stacked_coef) %>%
  mutate(s1 = replace_na(s1, 0)) %>%
  group_by(param) %>% 
  summarise_at(vars("s1"), list(mean = mean)) %>%
  arrange(-abs(mean))

# Save as a table

mod_mag %>% left_join(mod_count) %>%
  separate_wider_delim(param, delim = "_", names = c("method", "feattype", "feat", "focgene")) %>%
  write.csv("S3/analysis/stack_weight_coef_table.csv")


# Load predictions for all test set sequences

stacked_raw <- read.csv("S3/analysis/stack_weight_subtypeacc_raw.csv") %>%
  bind_cols(allflu_wgs_ref %>% filter((subtype %in% holdout_zoon & label == "zoon")|(subtype %in% holdout_nz)) %>% arrange(subtype) %>% select(-X, -label, -subtype)) %>%
  left_join(allflu_wgs_ref %>%
              filter(subtype %in% holdout_zoon & label == "zoon"|subtype %in% holdout_nz) %>%
              dplyr::count(subtype)) %>%
  mutate(orig_subtype = subtype) %>%
  mutate(subtype = case_when(n < 10 ~ "rare subtypes (< 10)",
                             TRUE ~ subtype),
         subtype = factor(subtype, levels = c("H4N6", "H4N8", "H8N4", "H16N3", "H5N1", "H5N6", "H7N9", "H9N2", "rare subtypes (< 10)")))

# List all avian test set sequences that were above threshold for zoonotic prediction

zoon_risk <- stacked_raw %>% 
  filter(label == "nz" & pred == "hzoon") %>% 
  arrange(-hzoon) %>%
  as.data.frame


# Fig 3
# Plot all test set sequence predictions

fig_results_stack_raw_p01 <- stacked_raw %>% 
  ggplot(aes(x = subtype, y = hzoon, colour = label, size = orig_subtype, pch = orig_subtype)) +
  geom_jitter(alpha = 0.4, position = position_jitter(width = 0.2, seed = 1649)) +
  geom_hline(aes(yintercept = read.csv("S3/analysis/stack_weight_results.csv") %>% pull(threshold)), linetype = "dashed", color = "gray30", linewidth = 1.2, lwd = 1.2) +
  scale_shape_manual("subtype", 
                     values = c("H10N8" = 12, 
                                "H16N3" = 19, 
                                "H3N8" = 15, 
                                "H4N6" = 19,
                                "H4N8" = 19,
                                "H5N1" = 19, 
                                "H5N6" = 19,
                                "H7N3" = 17, 
                                "H7N4" = 8, 
                                "H7N7" = 4, 
                                "H7N9" = 19, 
                                "H8N4" = 19, 
                                "H9N2" = 19),
                     breaks = c("H10N8", "H3N8", "H7N3", "H7N4", "H7N7")) + 
  scale_size_manual("subtype", 
                    values = c("H10N8" = 2.5, 
                               "H16N3" = 1.5, 
                               "H3N8" = 2.5, 
                               "H4N6" = 1.5,
                               "H4N8" = 1.5,
                               "H5N1" = 1.5, 
                               "H5N6" = 1.5,
                               "H7N3" = 2.5, 
                               "H7N4" = 2.5, 
                               "H7N7" = 2.5, 
                               "H7N9" = 1.5, 
                               "H8N4" = 1.5, 
                               "H9N2" = 1.5),
                    breaks = c("H10N8", "H3N8", "H7N3", "H7N4", "H7N7")) + 
  theme_bw() +
  guides(colour = "none") +
  theme(legend.title=element_blank(),
        legend.margin=margin(t = 0, unit='cm'),
        legend.position = "inside",
        legend.position.inside = c(.93,.87),
        legend.key.size = unit(0.4, 'cm')) +
  ylab("p(zoonotic)") +
  xlab("Subtype") +
  guides(color = "none")

ggsave(paste0("S3\\figures_tables\\fig_results_stack_weight_raw_p01.png"), plot = fig_results_stack_raw_p01, width = 10, height = 5.5)


# Log-scaled version of Fig 3

fig_results_stack_raw <- stacked_raw %>% 
  ggplot(aes(x = subtype, y = log(hzoon), colour = label, size = orig_subtype, pch = orig_subtype)) +
  geom_jitter(alpha = 0.4, position = position_jitter(width = 0.2, seed = 1649)) +
  geom_hline(aes(yintercept = read.csv("S3/analysis/stack_weight_results.csv") %>% pull(threshold) %>% log), linetype = "dashed", color = "gray30", linewidth = 1.2, lwd = 1.2) +
  scale_shape_manual("subtype", 
                     values = c("H10N8" = 12, 
                                "H16N3" = 19, 
                                "H3N8" = 15, 
                                "H4N6" = 19,
                                "H4N8" = 19,
                                "H5N1" = 19, 
                                "H5N6" = 19,
                                "H7N3" = 17, 
                                "H7N4" = 8, 
                                "H7N7" = 4, 
                                "H7N9" = 19, 
                                "H8N4" = 19, 
                                "H9N2" = 19),
                     breaks = c("H10N8", "H3N8", "H7N3", "H7N4", "H7N7")) + 
  scale_size_manual("subtype", 
                    values = c("H10N8" = 2.5, 
                               "H16N3" = 1.5, 
                               "H3N8" = 2.5, 
                               "H4N6" = 1.5,
                               "H4N8" = 1.5,
                               "H5N1" = 1.5, 
                               "H5N6" = 1.5,
                               "H7N3" = 2.5, 
                               "H7N4" = 2.5, 
                               "H7N7" = 2.5, 
                               "H7N9" = 1.5, 
                               "H8N4" = 1.5, 
                               "H9N2" = 1.5),
                    breaks = c("H10N8", "H3N8", "H7N3", "H7N4", "H7N7")) + 
  theme_bw() +
  guides(colour = "none") +
  theme(legend.title=element_blank(),
        legend.margin=margin(t = 0, unit='cm'),
        legend.position = "inside",
        legend.position.inside = c(.93,.27),
        legend.key.size = unit(0.4, 'cm')) +
  ylab("log(p(zoonotic))") +
  xlab("Subtype") +
  guides(color = "none")

ggsave(paste0("S3\\figures_tables\\fig_results_stack_weight_raw.png"), plot = fig_results_stack_raw, width = 10, height = 5.5)


###########################################
# Permutation variable importance figures #
###########################################

# Read through all permutation variable importance csv files in case they were batch-processed

varimp <- list.files(path = "S3\\analysis\\", pattern = "varimp_perm", full.names = TRUE) %>%
  purrr::map_dfr(read.csv) %>%
  select(-X) %>%
  separate(var, into = c("feat", "focgene"), sep="_(?=[^_]+$)", remove=FALSE) %>%
  mutate(featset = case_when(
    grepl("_Bias$", feat) ~ "nuc: composition",
    grepl(paste0("^", paste0(rep("[A|C|G|T]", 2), collapse = ""), "$"), feat) ~ "nuc: 2-mers",
    grepl(paste0("^", paste0(rep("[A|C|G|T]", 3), collapse = ""), "$"), feat) ~ "nuc: 3-mers",
    grepl(paste0("^", paste0(rep("[A|C|G|T]", 4), collapse = ""), "$"), feat) ~ "nuc: 4-mers",
    grepl(paste0("^", paste0(rep("[A|C|G|T]", 5), collapse = ""), "$"), feat) ~ "nuc: 5-mers",
    grepl(paste0("^", paste0(rep("[A|C|G|T]", 6), collapse = ""), "$"), feat) ~ "nuc: 6-mers",
    grepl("^DPC", feat) ~ "prot: 2-mers",
    grepl("^CTDC", feat) ~ "prot: CTD-C",
    grepl("^CTDT", feat) ~ "prot: CTD-T",
    grepl("^CTDD", feat) ~ "prot: CTD-D",
    grepl("^PAAC", feat) ~ "prot: PseAAC",
    grepl("^CTriad", feat) ~ "prot: CTriad")
  )

# Fig 4
# Permutation variable importance boxplots

fig_varimp_ordered <- varimp %>% 
  filter(AUC_loss > 0) %>%
  mutate(label = paste0(focgene, ", ", featset),
         label = forcats::fct_reorder(label, AUC_loss, .fun = mean, .desc = TRUE)) %>%
  ggplot(aes(x = label, y = AUC_loss, fill = focgene)) +
  geom_boxplot(alpha = 0.9) +
  theme_bw() +
  scale_fill_manual(values = cbbPalette_ordered) +
  scale_y_log10() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title=element_blank()) +
  xlab("Influenza virus gene/protein") +
  ylab("AUROC loss")

ggsave("S3\\figures_tables\\permut_varimp_ordered.png", plot = fig_varimp_ordered, width = 14, height = 5)


# Identify individual models that contributed to at least five holdout subtype stack models
pmc_df <- mod_count %>%
  separate_wider_delim(param, delim = "_", names = c("method", "feattype", "feat", "focgene"), cols_remove = FALSE) %>%
  mutate(focgene = as.factor(focgene)) %>%
  filter(n >= 5) %>%
  mutate(method = case_when(
    method == "glmnet" ~ "PLR",
    method == "rf" ~ "RF",
    method == "svm" ~ "RSVM",
    method == "svmlin" ~ "SVM",
    method == "xgb" ~ "XGB",
  )) %>%
  mutate(label = paste0(method, ", ", feattype, "_", feat),
         label = gsub("cds_compbias", "nuc: composition", label),
         label = gsub("nuc_2mer", "nuc: 2-mers", label),
         label = gsub("nuc_3mer", "nuc: 3-mers", label),
         label = gsub("nuc_4mer", "nuc: 4-mers", label),
         label = gsub("nuc_5mer", "nuc: 5-mers", label),
         label = gsub("nuc_6mer", "nuc: 6-mers", label),
         label = gsub("prot_2mer", "prot: 2-mers", label),
         label = gsub("prot_ctdc", "prot: CTD-C", label),
         label = gsub("prot_ctdt", "prot: CTD-T", label),
         label = gsub("prot_ctdd", "prot: CTD-D", label),
         label = gsub("prot_ctriad", "prot: CTriad", label),
         label = gsub("prot_pseaac", "prot: PseAAC", label)) %>%
  mutate(param = forcats::fct_reorder(param, n, .desc = FALSE))


# Supp Fig 4
# Bar counts of holdout subtype stack models each individual model contributed to
fig_modcount <- pmc_df %>%
  ggplot(aes(x = param, 
             y = n, 
             fill = focgene)) +
  geom_bar(alpha = 0.9, stat = "identity", show.legend = TRUE) +
  theme_bw() +
  scale_fill_manual(values = cbbPalette_ordered, drop = FALSE) +
  scale_x_discrete(labels=pmc_df$label) +
  scale_y_continuous(expand = c(0, 0.1), breaks = c(0:11)) +
  ylab("N stacks") +
  xlab("Machine learning model") +
  coord_flip() +
  theme(legend.title=element_blank()) +
  guides(fill = guide_legend(nrow=1))

# Permutation variable importance boxplots grouped by feature set
fig_varimp <- varimp %>% 
  filter(AUC_loss > 0) %>%
  ggplot(aes(x = focgene, y = AUC_loss, fill = focgene)) +
  geom_boxplot(alpha = 0.9) +
  theme_bw() +
  scale_fill_manual(values = cbbPalette_ordered) +
  scale_y_log10() +
  ylab("AUROC loss") +
  xlab("Influenza virus gene/protein") +
  facet_wrap(~ featset, nrow=4) +
  guides(fill = "none")

# Combine into multipanel

fig_varimp_combi <- fig_modcount + fig_varimp +
  plot_annotation(tag_levels = 'A') +
  plot_layout(widths = c(1,3), guides = 'collect') &
  theme(legend.position = "bottom")

ggsave("S3\\figures_tables\\model_count_permut_varimp.png", plot = fig_varimp_combi, width = 12, height = 6.5)