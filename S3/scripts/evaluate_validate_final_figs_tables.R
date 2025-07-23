###############################################
# Overall figures across all model algorithms #
###############################################

################################
# Load packages, reference IDs #
################################

library(dplyr)
library(tidyr)
library(magrittr)
library(ggplot2)
library(ggrepel)
library(rentrez)

allflu_wgs_ref <- read.csv("S3\\data\\full\\allflu_wgs_ref.csv")

#############
# Data figs #
#############

# Supplementary Fig S1

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

#############################################
# Performance figs on individual algorithms #
#############################################

# Performance metrics on holdout sets

results_rf <- read.csv(paste0("S3\\analysis\\results_", "14_02_24", ".csv"), na.strings = "NaN") %>% filter(cluster_set == "70_7") %>% mutate(method = "rf") 
results_plr <- read.csv(paste0("S3\\analysis\\results_", "15_02_24", ".csv"), na.strings = "NaN") %>% filter(cluster_set == "70_7") %>% mutate(method = "glmnet")
results_xgb <- read.csv(paste0("S3\\analysis\\results_", "16_02_24", ".csv"), na.strings = "NaN") %>% filter(cluster_set == "70_7") %>% mutate(method = "xgb")
results_svmlin <- read.csv(paste0("S3\\analysis\\results_", "17_02_24", ".csv"), na.strings = "NaN") %>% filter(cluster_set == "70_7") %>% mutate(method = "svmlin")
results_svmrad <- read.csv(paste0("S3\\analysis\\results_", "18_02_24", ".csv"), na.strings = "NaN") %>% filter(cluster_set == "70_7") %>% mutate(method = "svm")

all_res <- bind_rows(results_rf,
                     results_svmlin,
                     results_svmrad,
                     results_xgb,
                     results_plr) 

#all_res %>% write.csv("S3\\analysis\\results_all_methods.csv")

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
    featset == "prot_ctriad" ~ "prot: c.triad",
    featset == "prot_pseaac" ~ "prot: pseudo-aac"),
    focgene = factor(focgene, levels = c("PB2", "PB1", "PA", "HA", "NP", "NA", "M1", "NS1"))
  ) %>%
  gather(metric, value, -method, -featset, -focgene) %>%
  filter(metric %in% c("AUC")) %>%
  mutate(value = as.numeric(value)) %>%
  ggplot(aes(x = focgene, y = featset, fill = value)) + 
  geom_tile(color="white") +
  scale_fill_distiller("AUC", palette = "RdBu", limits = c(0,1)) +
  #  scale_fill_viridis_c("AUC") +
  facet_wrap(~ method, nrow=1) +
  theme_bw() +
  # theme(plot.background = element_rect(fill = "#F2F6F9", color = "#F2F6F9"),      # poster colours
  #       legend.key = element_rect(fill = "#F2F6F9", color = "#F2F6F9")) +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank())

ggsave("S3\\figures_tables\\all_results_heat_AUC_70_7.png", plot = fig_results_heat_AUC_70_7, width = 15, height = 3)

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
    featset == "prot_ctriad" ~ "prot: c.triad",
    featset == "prot_pseaac" ~ "prot: pseudo-aac"),
    method = case_when(
      method == "glmnet" ~ "LASSO",
      method == "rf" ~ "RF",
      method == "svm" ~ "SVM (lin)",
      method == "svmlin" ~ "SVM (rad)",
      method == "xgb" ~ "XGBoost",
    ),
    focgene = factor(focgene, levels = c("PB2", "PB1", "PA", "HA", "NP", "NA", "M1", "NS1"))
  ) %>%
  gather(metric, value, -method, -featset, -focgene) %>%
  filter(metric %in% c("AUC")) %>%
  mutate(value = as.numeric(value)) %>%
  ggplot(aes(x = focgene, y = featset, fill = value)) + 
  geom_tile(color="white") +
  geom_point(aes(pch = method), color = "black", size = 2) +
  scale_fill_distiller("AUROC", palette = "RdBu", limits = c(0,1)) +
  scale_shape_manual("algorithm", values = c(19, 15, 17, 18, 8)) + 
  #  scale_fill_viridis_c("AUC") +
  theme_bw() +
  # theme(plot.background = element_rect(fill = "#F2F6F9", color = "#F2F6F9"),      # poster colours
  #       legend.background = element_rect(fill = "#F2F6F9", color = "#F2F6F9"),
  #       axis.title.x = element_text(size=10),
  #       axis.title.y = element_text(size=10)) +
  xlab("Protein model") +
  ylab("Feature set")

ggsave("S3\\figures_tables\\all_results_heat_AUC_70_7_one.png", plot = fig_results_heat_AUC_70_7_one, width = 6, height = 3.5)



# fig_results_heat_F1_70_7 <- all_res %>%
#   group_by(method, featset, focgene) %>%
#   summarise(F1 = mean(F1)) %>% 
#   ungroup %>%
#   mutate(featset = case_when(
#     featset == "cds_compbias" ~ "nuc: composition",
#     featset == "nuc_2mer" ~ "nuc: 2-mers",
#     featset == "nuc_3mer" ~ "nuc: 3-mers",
#     featset == "nuc_4mer" ~ "nuc: 4-mers",
#     featset == "nuc_5mer" ~ "nuc: 5-mers",
#     featset == "nuc_6mer" ~ "nuc: 6-mers",
#     featset == "prot_2mer" ~ "prot: 2-mers",
#     featset == "prot_ctdc" ~ "prot: CTD-C",
#     featset == "prot_ctdt" ~ "prot: CTD-T",
#     featset == "prot_ctdd" ~ "prot: CTD-D",
#     featset == "prot_ctriad" ~ "prot: c.triad",
#     featset == "prot_pseaac" ~ "prot: pseudo-aac"),
#     focgene = factor(focgene, levels = c("PB2", "PB1", "PA", "HA", "NP", "NA", "M1", "NS1"))
#   ) %>%
#   gather(metric, value, -method, -featset, -focgene) %>%
#   filter(metric %in% c("F1")) %>%
#   mutate(value = as.numeric(value)) %>%
#   ggplot(aes(x = focgene, y = featset, fill = value)) + 
#   geom_tile(color="white") +
#   scale_fill_viridis_c("F1") +
#   facet_wrap(~ method, nrow=1) +
#   theme_bw() +
#   theme(plot.background = element_rect(fill = "#F2F6F9", color = "#F2F6F9"),
#         legend.key = element_rect(fill = "#F2F6F9", color = "#F2F6F9")) +
#   theme(axis.title.x=element_blank(),
#         axis.title.y=element_blank())
# 
# ggsave("S3\\figures_tables\\all_results_heat_F1_70_7.png", plot = fig_results_heat_F1_70_7, width = 15, height = 3)

# fig_results_heat_Sensitivity_70_7 <- all_res %>%
#   group_by(method, featset, focgene) %>%
#   summarise(Sensitivity = mean(Sensitivity)) %>% 
#   ungroup %>%
#   mutate(featset = case_when(
#     featset == "cds_compbias" ~ "nuc: composition",
#     featset == "nuc_2mer" ~ "nuc: 2-mers",
#     featset == "nuc_3mer" ~ "nuc: 3-mers",
#     featset == "nuc_4mer" ~ "nuc: 4-mers",
#     featset == "nuc_5mer" ~ "nuc: 5-mers",
#     featset == "nuc_6mer" ~ "nuc: 6-mers",
#     featset == "prot_2mer" ~ "prot: 2-mers",
#     featset == "prot_ctdc" ~ "prot: CTD-C",
#     featset == "prot_ctdt" ~ "prot: CTD-T",
#     featset == "prot_ctdd" ~ "prot: CTD-D",
#     featset == "prot_ctriad" ~ "prot: c.triad",
#     featset == "prot_pseaac" ~ "prot: pseudo-aac"),
#     focgene = factor(focgene, levels = c("PB2", "PB1", "PA", "HA", "NP", "NA", "M1", "NS1"))
#   ) %>%
#   gather(metric, value, -method, -featset, -focgene) %>%
#   filter(metric %in% c("Sensitivity")) %>%
#   mutate(value = as.numeric(value)) %>%
#   ggplot(aes(x = focgene, y = featset, fill = value)) + 
#   geom_tile(color="white") +
#   scale_fill_viridis_c("Sensitivity") +
#   facet_wrap(~ method, nrow=1) +
#   theme_bw() +
#   theme(plot.background = element_rect(fill = "#F2F6F9", color = "#F2F6F9"),
#         legend.key = element_rect(fill = "#F2F6F9", color = "#F2F6F9")) +
#   theme(axis.title.x=element_blank(),
#         axis.title.y=element_blank())
# 
# ggsave("S3\\figures_tables\\all_results_heat_Sensitivity_70_7.png", plot = fig_results_heat_Sensitivity_70_7, width = 15, height = 3)


######################################
# Performance figs on stacked models #
######################################

# What models are being selected?

stacked_coef <- read.csv("S3/analysis/stack_weight_coef.csv") %>%
  select(-X) %>%
  filter(param != "(Intercept)" & param != "lambda")

# By count of retention

mod_count <- stacked_coef %>%
  group_by(param) %>%
  tally() %>%
  arrange(-n) 

mod_count %>% 
  dplyr::slice(1:20)

# By magnitude (≈ variable importance)

mod_mag <- stacked_coef %>% 
  tidyr::expand(param, subtype) %>% 
  left_join(stacked_coef) %>%
  mutate(s1 = replace_na(s1, 0)) %>%
  group_by(param) %>% 
  summarise_at(vars("s1"), list(mean = mean)) %>%
  arrange(-abs(mean))

mod_mag %>% 
  dplyr::slice(1:20)

mod_mag %>% left_join(mod_count) %>% dplyr::slice(1:20)

mod_mag %>% left_join(mod_count) %>%
  separate_wider_delim(param, delim = "_", names = c("method", "feattype", "feat", "focgene")) %>%
  write.csv("S3/analysis/stack_weight_coef_table.csv")


# Predictions for all test set sequences

stacked_raw <- read.csv("S3/analysis/stack_weight_subtypeacc_raw.csv") %>%
  bind_cols(allflu_wgs_ref %>% filter((subtype %in% holdout_zoon & label == "zoon")|(subtype %in% holdout_nz)) %>% arrange(subtype) %>% select(-X, -label, -subtype)) %>%
  left_join(allflu_wgs_ref %>%
              filter(subtype %in% holdout_zoon & label == "zoon"|subtype %in% holdout_nz) %>%
              dplyr::count(subtype)) %>%
  mutate(orig_subtype = subtype) %>%
  mutate(subtype = case_when(n < 10 ~ "rare subtypes (< 10)",
                             TRUE ~ subtype),
         subtype = factor(subtype, levels = c("H4N6", "H4N8", "H8N4", "H16N3", "H5N1", "H5N6", "H7N9", "H9N2", "rare subtypes (< 10)")))

# All avian sequences above threshold for zoonotic prediction

zoon_risk <- stacked_raw %>% 
  filter(label == "nz" & pred == "hzoon") %>% 
  arrange(-hzoon) %>%
  as.data.frame

# Search GenBank sequences for location

# metadata_extraction <- lapply(zoon_risk$title, function(x){
#   g <- tryCatch(entrez_search(db = "nuccore", term = x, retmax = 10000)$ids[1] %>%
#                   entrez_summary(db = "nuccore") %>% 
#                   purrr::flatten() %>%
#                   .$subname,
#                 error=function(e) NULL)
#   Sys.sleep(5)
#   return(g)
# }
# )

#### MANUALLY RESOLVE PIECES FOR NOW BUT IN FUTURE ALSO PULL $subtype AND USE TO DYNAMICALLY ASSIGN COLUMNS TO METADATA
metadata_data <- data.frame(metadata_data = unlist(metadata_extraction)) %>%
  mutate(metadata_data = gsub("CEIRS#.*#\\|","",metadata_data)) %>%
  separate_wider_delim(metadata_data, delim = "|", 
                       too_few = "debug",
                       too_many = "debug",
                       names = stringr::str_split_1("strain|serotype|host|lab_host|country|segment|collection_date", pattern = "\\|"))

# TIDY THIS UP LATER - ASSIGNS ABOUT 80% RIGHT BUT STILL NEEDS MANUAL FIXING FOR REMAINING 20%
# USING COPYPASTE AND ALT + : TO ONLY SELECT THOSE IN FILTER ETC
extract_last_capital_segment <- function(input_string) {
  segments <- strsplit(input_string, "\\|")[[1]]
  capital_segments <- segments[grepl("^[A-Z]", segments)]
  if (length(capital_segments) > 0) {
    return(tail(capital_segments, 1))
  } else {
    return(NA)  # or "" if you prefer
  }
}

metadata_countries <- data.frame(metadata_data = unlist(metadata_extraction)) %>%
  rowwise %>%
  mutate(metadata_data = extract_last_capital_segment(metadata_data)) %>%
  rename(clean_country = metadata_data)

# bind_cols(metadata_data, metadata_countries) %>% write.csv("S3\\data\\full\\metadata_risk_preds_nz.csv")

metadata_clean <- read.csv("S3\\data\\full\\metadata_risk_preds_nz.csv")

zoon_risk %>% 
  left_join(metadata_clean %>% mutate(strain = gsub(" ", "_", strain)),
            by = c("title" = "strain")) %>%
  rename(p_zoon = hzoon, prediction = pred) %>%
  select(p_zoon, label, subtype, prediction,	title, src,	date,	clean_country) %>%
  rename(country = clean_country) %>%
  distinct() %>%
  write.csv("S3\\data\\full\\zoon_risk_preds_nz.csv") 

# Plot

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
  # theme(plot.background = element_rect(fill = "#F2F6F9", color = "#F2F6F9"),      # poster colours
  #       legend.background = element_rect(fill = "#F2F6F9", color = "#F2F6F9"),
  #       axis.title.x = element_text(size=10),
  #       axis.title.y = element_text(size=10)) +
  ylab("log(p(zoonotic))") +
  xlab("Subtype") +
  guides(color = "none")

ggsave(paste0("S3\\figures_tables\\fig_results_stack_weight_raw.png"), plot = fig_results_stack_raw, width = 10, height = 5.5)

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
  # theme(plot.background = element_rect(fill = "#F2F6F9", color = "#F2F6F9"),      # poster colours
  #       legend.background = element_rect(fill = "#F2F6F9", color = "#F2F6F9"),
  #       axis.title.x = element_text(size=10),
  #       axis.title.y = element_text(size=10)) +
  ylab("p(zoonotic)") +
  xlab("Subtype") +
  guides(color = "none")

ggsave(paste0("S3\\figures_tables\\fig_results_stack_weight_raw_p01.png"), plot = fig_results_stack_raw_p01, width = 10, height = 5.5)
# 
# ggsave(paste0("S3\\figures_tables\\fig_results_stack_weight_raw_p01_poster.png"), plot = fig_results_stack_raw_p01, width = 10, height = 2)

#######################################################
# Permutation variable importance from stacked models #
#######################################################

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
    grepl("^PAAC", feat) ~ "prot: pseudo-aac",
    grepl("^CTriad", feat) ~ "prot: c.triad")
  )

varimp %>% arrange(AUC_loss)
varimp %>% group_by(focgene) %>% summarise(mean = mean(AUC_loss)) %>% arrange(mean)
varimp %>% group_by(feat) %>% summarise(mean = mean(AUC_loss)) %>% arrange(mean)

varimp %>%
  ggplot(aes(x = focgene, y = AUC_loss)) +
  geom_boxplot(alpha = 0.9) +
  theme_bw()

### THIS ONLY SHOWS THE VARIMP PERMUTATIONS TESTED, I.E. IN THE ACTUAL STACKS.
### WE COULD ADD IN ALL THE OTHER VARIMP VARIABLES NOT INCLUDED IN THE STACKS AS ZERO?
### SO FILL OUT EVERY COMBO OF FEAT - PROT AND LEFT JOIN THE VARIMP VALUES

# Box per featset*protein combo?

varimp %>% 
  ggplot(aes(x = focgene, y = AUC_loss, fill = featset)) +
  geom_jitter(alpha = 0.7, shape = 21, size = 3.5, position = position_jitter(width = 0.2, seed = 1615)) +
  theme_bw() +
  scale_fill_manual(values =   c("#E69F00", "#F0E442", "#B2DF8A","#56B4E9","#CC79A7", "#33A02C", "#0072B2", "#D55E00","#6A3D9A", "#EEEEEE", "#666666", "#000000")) +
  # theme(plot.background = element_rect(fill = "#F2F6F9", color = "#F2F6F9"),      # poster colours
  #       legend.background = element_rect(fill = "#F2F6F9", color = "#F2F6F9"),
  #       axis.title.x = element_text(size=10),
  #       axis.title.y = element_text(size=10)) +
  ylab("AUROC loss") +
  xlab("Protein")

# Plot avg per feat per protein?


#####################################
# Performance figs on external data #
#####################################

# Predictions for external sequence sets - H5N1 2.3.4.4b cattle outbreak

pred_dairy <- read.csv("S3\\data\\ext\\preds_weight_dairyc_raw.csv")
pred_misc <- read.csv("S3\\data\\ext\\preds_weight_misc_raw.csv") %>% filter(gid %in% c("EPI_ISL_19162802", "EPI_ISL_19027114"))

fig_results_dairyc <- bind_rows(pred_dairy, pred_misc) %>%
  filter(stack == "H5N1") %>%                                                     # ONLY USE THE H5N1 STACK
  filter(as.Date(date) > "2024-01-01") %>%
  ggplot(aes(x = as.Date(date), y = log(hzoon), colour = set)) +
  geom_point(alpha = 0.8) +
  geom_hline(aes(yintercept = read.csv("S3/analysis/stack_weight_results.csv") %>% pull(threshold) %>% log()), linetype = "dashed", color = "gray30", linewidth = 1.2, lwd = 1.2) +
  scale_color_manual(values = c("#C77CFF", "#F8766D")) +
  #  scale_y_continuous(limits = c(-5, -3)) +
  theme_bw() +
  theme(plot.background = element_rect(fill = "#F2F6F9", color = "#F2F6F9"),      # poster colours
        legend.background = element_rect(fill = "#F2F6F9", color = "#F2F6F9"),
        axis.title.x = element_text(size=10),
        axis.title.y = element_text(size=10)) +
  ylab("log(p(zoonotic))") +
  xlab("Month") +
  guides(color = "none")

ggsave(paste0("S3\\figures_tables\\fig_results_dairyc_poster_weight.png"), plot = fig_results_dairyc, width = 4.5, height = 2.25)


# Predictions for external sequence sets - CDC IRAT

pred_weight <- read.csv(paste0("S3\\data\\irat\\preds_weight_", set, "_raw.csv"))
irat_df <- read.csv("S3/data/irat/cdc_irat.csv", fileEncoding="UTF-8-BOM") %>% filter(incomplete != "Y")

fig_results_cdc_error_weight <- pred_weight %>%
  left_join(irat_df) %>%
  group_by(gid, id, host, emergence) %>% 
  summarise_at(vars("hzoon"), list(med = median, upper = ~quantile(., probs = 0.25), lower = ~quantile(., probs = 0.75))) %>%
  mutate(host = case_when(host == "za" ~ "avian, human-isolated",
                          host == "a" ~ "avian",
                          host == "zm" ~ "mammal, human-isolated",
                          host == "m" ~ "mammal")) %>%
  ggplot(aes(x = emergence, y = med, ymin = upper, ymax = lower, color = host, fill = host, label = id)) +
  geom_errorbar(alpha = 0.4, lwd = 1.2, linewidth = 1.2, width=0) +
  geom_point() +
  geom_hline(aes(yintercept = read.csv("S3/analysis/stack_weight_results.csv") %>% pull(threshold)), linetype = "dashed", color = "gray30", linewidth = 1.2, lwd = 1.2) +
  geom_text(hjust=-0.5, vjust=-0.2, show.legend  = FALSE) +
  scale_y_continuous(limits = c(0, 1), expand = c(0,0.1)) +
  scale_x_continuous(limits = c(2.7, 7.7), expand = c(0,0)) +
  ylab("p(zoonotic))") +
  xlab("CDC IRAT emergence score") +
  guides(label = "none") +
  theme_bw() +
  theme(plot.background = element_rect(fill = "#F2F6F9", color = "#F2F6F9"),
        axis.title.x = element_text(size=10),
        axis.title.y = element_text(size=10),
        legend.title=element_blank(),
        legend.margin=margin(t = 0, unit='cm'),
        legend.position = c(.82,.82),
        legend.key.size = unit(0.4, 'cm'))

ggsave(paste0("S3\\figures_tables\\fig_results_cdc_poster_weight.png"), plot = fig_results_cdc_error_weight, width = 6.5*2.5/3, height = 2.5)
