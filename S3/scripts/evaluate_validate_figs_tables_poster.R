#################
# Load packages #
#################

library(dplyr)
library(tidyr)
library(magrittr)
library(ggplot2)

results_date <- "16_02_24"

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
  theme(plot.background = element_rect(fill = "#F2F6F9", color = "#F2F6F9"),
        legend.key = element_rect(fill = "#F2F6F9", color = "#F2F6F9")) +
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
  theme(plot.background = element_rect(fill = "#F2F6F9", color = "#F2F6F9"),
        legend.background = element_rect(fill = "#F2F6F9", color = "#F2F6F9"),
        axis.title.x = element_text(size=10),
        axis.title.y = element_text(size=10)) +
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



# Predictions for individual sequences

stacked_raw <- read.csv("S3/analysis/stack_subtypeacc_raw.csv") %>%
  bind_cols(allflu_wgs_ref %>% filter((subtype %in% holdout_zoon & label == "zoon")|(subtype %in% holdout_nz)) %>% arrange(subtype) %>% select(-X, -label, -subtype)) %>%
  left_join(allflu_wgs_ref %>%
              filter(subtype %in% holdout_zoon & label == "zoon"|subtype %in% holdout_nz) %>%
              count(subtype)) %>%
  mutate(subtype = case_when(n < 10 ~ "rare subtypes (< 10)",
                             TRUE ~ subtype),
         subtype = factor(subtype, levels = c("H5N1", "H5N6", "H7N9", "H9N2", "rare subtypes (< 10)", "H4N6", "H4N8", "H8N4", "H16N3")))

stacked_raw %>% filter(label == "nz") %>% arrange(-hzoon) %>% head

fig_results_stack_raw <- stacked_raw %>%
  ggplot(aes(x = subtype, y = log(hzoon), colour = label)) +
  geom_jitter(alpha = 0.4, width = 0.2) +
  geom_hline(aes(yintercept = read.csv("S3/analysis/stack_results.csv") %>% pull(threshold) %>% log()), linetype = "dashed", color = "gray30", linewidth = 1.2) +
  theme_bw() +
  theme(plot.background = element_rect(fill = "#F2F6F9", color = "#F2F6F9"),
        legend.background = element_rect(fill = "#F2F6F9", color = "#F2F6F9"),
        axis.title.x = element_text(size=10),
        axis.title.y = element_text(size=10)) +
  ylab("log(p(zoonotic))") +
  xlab("Subtype") +
  guides(color = "none")

ggsave(paste0("S3\\figures_tables\\fig_results_stack_rawposter_",results_date,".png"), plot = fig_results_stack_raw, width = 10, height = 2)


# Predictions for external set - cattle

pred_dairy <- read.csv("S3\\data\\ext\\preds_dairyc_raw.csv")
pred_misc <- read.csv("S3\\data\\ext\\preds_misc_raw.csv") %>% filter(gid %in% c("EPI_ISL_19162802", "EPI_ISL_19027114"))

fig_results_dairyc <- bind_rows(pred_dairy, pred_misc) %>%
  filter(stack == "H5N1") %>%                 # ONLY USE THE H5N1 STACK
  filter(as.Date(date) > "2024-01-01") %>%
  ggplot(aes(x = as.Date(date), y = log(hzoon), colour = set)) +
  geom_point(alpha = 0.8) +
  geom_hline(aes(yintercept = read.csv("S3/analysis/stack_results.csv") %>% pull(threshold) %>% log()), linetype = "dashed", color = "gray30", linewidth = 1.2, lwd = 1.2) +
  scale_color_manual(values = c("#C77CFF", "#F8766D")) +
#  scale_y_continuous(limits = c(-5, -3)) +
  theme_bw() +
  theme(plot.background = element_rect(fill = "#F2F6F9", color = "#F2F6F9"),
        legend.background = element_rect(fill = "#F2F6F9", color = "#F2F6F9"),
        axis.title.x = element_text(size=10),
        axis.title.y = element_text(size=10)) +
  ylab("log(p(zoonotic))") +
  xlab("Month") +
  guides(color = "none")

ggsave(paste0("S3\\figures_tables\\fig_results_dairyc_poster_",results_date,".png"), plot = fig_results_dairyc, width = 4.5, height = 2.25)


pred_dairy <- read.csv("S3\\data\\ext\\preds_weight_dairyc_raw.csv")
pred_misc <- read.csv("S3\\data\\ext\\preds_weight_misc_raw.csv") %>% filter(gid %in% c("EPI_ISL_19162802", "EPI_ISL_19027114"))

fig_results_dairyc <- bind_rows(pred_dairy, pred_misc) %>%
  filter(stack == "H5N1") %>%                 # ONLY USE THE H5N1 STACK
  filter(as.Date(date) > "2024-01-01") %>%
  ggplot(aes(x = as.Date(date), y = hzoon, colour = set)) +
  geom_point(alpha = 0.8) +
  geom_hline(aes(yintercept = read.csv("S3/analysis/stack_results.csv") %>% pull(threshold)), linetype = "dashed", color = "gray30", linewidth = 1.2, lwd = 1.2) +
  scale_color_manual(values = c("#C77CFF", "#F8766D")) +
  #  scale_y_continuous(limits = c(-5, -3)) +
  theme_bw() +
  theme(plot.background = element_rect(fill = "#F2F6F9", color = "#F2F6F9"),
        legend.background = element_rect(fill = "#F2F6F9", color = "#F2F6F9"),
        axis.title.x = element_text(size=10),
        axis.title.y = element_text(size=10)) +
  ylab("p(zoonotic)") +
  xlab("Month") +
  guides(color = "none")

ggsave(paste0("S3\\figures_tables\\fig_results_dairyc_poster_weight_",results_date,".png"), plot = fig_results_dairyc, width = 4.5, height = 2.25)


# Predictions for external set - CDC IRAT

pred <- read.csv(paste0("S3\\data\\irat\\preds_", set, "_raw.csv"))
pred_weight <- read.csv(paste0("S3\\data\\irat\\preds_weight_", set, "_raw.csv"))
irat_df <- read.csv("S3/data/irat/cdc_irat.csv", fileEncoding="UTF-8-BOM") %>% filter(incomplete != "Y")

fig_results_cdc_error <- pred %>%
  left_join(irat_df) %>%
  group_by(gid, id, host, emergence) %>% 
  summarise_at(vars("hzoon"), list(med = median, upper = ~quantile(., probs = 0.25), lower = ~quantile(., probs = 0.75))) %>%
  mutate(host = case_when(host == "za" ~ "avian, human-isolated",
                          host == "a" ~ "avian",
                          host == "zm" ~ "mammal, human-isolated",
                          host == "m" ~ "mammal")) %>%
  ggplot(aes(x = emergence, y = log(med), ymin = log(upper), ymax = log(lower), color = host, fill = host, label = id)) +
  geom_errorbar(alpha = 0.4, lwd = 1.2, linewidth = 1.2, width=0) +
  geom_point() +
  geom_hline(aes(yintercept = read.csv("S3/analysis/stack_results.csv") %>% pull(threshold) %>% log()), linetype = "dashed", color = "gray30", linewidth = 1.2, lwd = 1.2) +
  geom_text(hjust=-0.5, vjust=-0.2, show.legend  = FALSE) +
  scale_x_continuous(limits = c(2.7, 7.7), expand = c(0,0)) +
  ylab("log(p(zoonotic))") +
  xlab("CDC IRAT emergence score") +
  guides(label = "none") +
  theme_bw() +
  theme(plot.background = element_rect(fill = "#F2F6F9", color = "#F2F6F9"),
        axis.title.x = element_text(size=10),
        axis.title.y = element_text(size=10),
        legend.title=element_blank(),
        legend.margin=margin(t = 0, unit='cm'),
        legend.position = c(.14,.87),
        legend.key.size = unit(0.4, 'cm'))

ggsave(paste0("S3\\figures_tables\\fig_results_cdc_poster.png"), plot = fig_results_cdc_error, width = 6.5, height = 3)

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
  geom_hline(aes(yintercept = read.csv("S3/analysis/stack_results.csv") %>% pull(threshold)), linetype = "dashed", color = "gray30", linewidth = 1.2, lwd = 1.2) +
  geom_text(hjust=-0.5, vjust=-0.2, show.legend  = FALSE) +
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
        legend.position = c(.14,.87),
        legend.key.size = unit(0.4, 'cm'))

ggsave(paste0("S3\\figures_tables\\fig_results_cdc_poster_weight.png"), plot = fig_results_cdc_error_weight, width = 6.5, height = 3)





bestmodels_rf <- results_rf %>% group_by(focgene) %>% slice_max(AUC) %>% slice_max(F1)
bestmodels_svm <- results_svm %>% group_by(focgene) %>% slice_max(AUC) %>% slice_max(F1)
bestmodels_plr <- results_plr %>% group_by(focgene) %>% slice_max(AUC) %>% slice_max(F1)

allflu_wgs_ref <- read.csv("S3\\data\\full\\allflu_wgs_ref.csv")
holdout_zoon <- c("H7N9", "H5N1", "H9N2", "H5N6", "H10N8", "H7N3", "H3N8", "H7N7", "H7N4")
holdout_nz <- c("H4N6", "H16N3", "H4N8", "H8N4")
holdouts <- c(holdout_zoon, holdout_nz)

results_subtype_rf <- read.csv(paste0("S3\\analysis\\resultsbysubtype_", "03_11_23", ".csv"), na.strings = "NaN") %>% 
  left_join(allflu_wgs_ref %>%
              filter(subtype %in% holdout_zoon & label == "zoon"|subtype %in% holdout_nz) %>%
              count(subtype)) %>% 
  mutate(method = "RF") 

results_subtype_svm <- read.csv(paste0("S3\\analysis\\resultsbysubtype_", "15_11_23", ".csv"), na.strings = "NaN") %>% 
  left_join(allflu_wgs_ref %>%
              filter(subtype %in% holdout_zoon & label == "zoon"|subtype %in% holdout_nz) %>%
              count(subtype)) %>% 
  mutate(method = "SVM") 

results_subtype_plr <- read.csv(paste0("S3\\analysis\\resultsbysubtype_", "13_11_23", ".csv"), na.strings = "NaN") %>% 
  left_join(allflu_wgs_ref %>%
              filter(subtype %in% holdout_zoon & label == "zoon"|subtype %in% holdout_nz) %>%
              count(subtype)) %>% 
  mutate(method = "PLR") 

order <- bind_rows(results_subtype_rf,
                   results_subtype_svm,
                   results_subtype_plr) %>%
  mutate(subtype = case_when(n < 10 ~ "other",
                             TRUE ~ subtype),
         subtype = factor(subtype, levels = c("H5N1", "H5N6", "H7N9", "H9N2", "other", "H4N6", "H4N8", "H8N4", "H16N3")),
         method = factor(method, levels = c("RF", "PLR", "SVM"))) %>%
  group_by(subtype) %>%
  summarise(accuracy = mean(accuracy)) %>%
  arrange(-accuracy) %>%
  pull(subtype) %>%
  as.character

# Boxplots over clusters: best performing feature set only
plot_dat  <- bind_rows(results_subtype_rf,
                       results_subtype_svm,
                       results_subtype_plr) %>%
  mutate(subtype = case_when(n < 10 ~ "other",
                             TRUE ~ subtype),
         subtype = factor(subtype, levels = order),
         method = factor(method, levels = c("RF", "PLR", "SVM"))) %>%
  group_by(cluster_set, featset, focgene, subtype, method) %>%
  summarise(accuracy = mean(accuracy)) %>%
  mutate(label = case_when(subtype %in% holdout_nz ~ "avian",
                           TRUE ~ "zoon"))

# plot_dat$subtype_cols <- ifelse(plot_dat$subtype %in% holdout_zoon, "bold", "plain")

fig_results_subtype_boxes_all <- plot_dat %>%
  ggplot(aes(x = subtype, y = accuracy, fill = method)) +
  geom_boxplot(alpha = 0.9, outlier.shape=NA) +
  facet_wrap(vars(focgene), nrow = 4) +
  scale_fill_manual(values = c("#E69F00", "#009E73", "#CC79A7")) +
  theme_bw() +
  #  theme(axis.text.x = element_text(face = subtype_cols)) +
  theme(plot.background = element_rect(fill = "#F2F6F9", color = "#F2F6F9")) +
  theme(legend.key = element_rect(fill = "transparent", colour = "transparent")) +
  ylab("Accuracy") +
  xlab("Subtype")

ggsave("E:\\Working\\Talks\\2023\\11 Epidemics\\boxes_allmeths.png", plot = fig_results_subtype_boxes_all, width = 14, height = 6)



# Boxplots over clusters: best performing feature set only
# order_rf <- results_subtype_rf %>%
#   mutate(subtype = case_when(n < 10 ~ "other",
#                              TRUE ~ subtype),
#          subtype = factor(subtype, levels = c("H5N1", "H5N6", "H7N9", "H9N2", "other", "H4N6", "H4N8", "H8N4", "H16N3")),
#          method = factor(method, levels = c("RF", "PLR"))) %>%
#   group_by(subtype, label) %>%
#   summarise(accuracy = mean(accuracy)) %>%
#   arrange(-accuracy) %>%
#   pull(subtype) %>%
#   as.character

plot_dat_rf  <- results_subtype_rf %>%
  mutate(subtype = case_when(n < 10 ~ "other",
                             TRUE ~ subtype),
         subtype = factor(subtype, levels = c("H4N8", "H4N6", "H8N4", "H16N3", "H5N6","H9N2", "H7N9", "H5N1",  "other")),
         method = factor(method, levels = c("RF", "PLR")),
         focgene = factor(focgene, levels = c("PB2", "PB1", "PA", "HA", "NP", "NA", "M1", "NS1"))) %>%
  group_by(cluster_set, featset, focgene, subtype, method) %>%
  summarise(accuracy = mean(accuracy)) %>%
  mutate(label = case_when(subtype %in% holdout_nz ~ "avian",
                           TRUE ~ "zoon"))

fig_results_subtype_boxes_rf <- plot_dat_rf %>%
  ggplot(aes(x = subtype, y = accuracy, fill = label)) +
  geom_boxplot(alpha = 0.9, outlier.shape=NA) +
  facet_wrap(vars(focgene), nrow = 2) +
  theme_bw() +
  theme(plot.background = element_rect(fill = "#F2F6F9", color = "#F2F6F9"),
        legend.position = "none",
        axis.title.x=element_blank()) +
  ylab("Accuracy (proportion correctly predicted)")

ggsave("E:\\Working\\Talks\\2023\\11 Epidemics\\boxes.png", plot = fig_results_subtype_boxes_rf, width = 15, height = 4)




# Points: best performing model (feature set and cluster)
# IF CHANGING TO SVM, CHANGE THE DATE TOO
best_files <- bestmodels_rf %>% 
  filter(focgene == "PB2") %>%
  mutate(file = paste0("S3\\analysis\\subtyperaw\\subtypeacc_raw_", cluster_set, "_", featset, "_", focgene, "_", "03_11_23", ".csv")) %>%
  pull(file)

threshold <- bestmodels_rf %>% 
  filter(focgene == "PB2") %>%
  pull(threshold)

bestmodels_raw <- best_files %>% 
  purrr::map_dfr(function(x) read.csv(x, na.strings = "NaN")) %>% 
  bind_rows %>%
  mutate(subtype = factor(subtype, levels = c("H4N8", "H4N6", "H8N4", "H16N3", "H5N6","H9N2", "H7N9", "H5N1",  "other")),
         label = case_when(label == "nz" ~ "avian",
                           TRUE ~ "zoon")) %>%
  filter(focgene == "PB2" & subtype %in% c("H4N8", "H4N6", "H8N4", "H16N3"))

fig_results_subtype_PB2 <- bestmodels_raw %>%
  filter(zoon < 0.4) %>%
  ggplot(aes(x = subtype, y = zoon)) +
  geom_jitter(alpha = 0.4, width = 0.2, fill = "#F8766D", color = "#F8766D", shape=21) +
  geom_jitter(aes(x = subtype, y = zoon, fill = label), 
              data = bestmodels_raw %>% filter(zoon >= 0.4), alpha = 0.5, width = 0.2, colour = "black", stroke = 1.5,  shape=21) +
  geom_hline(aes(yintercept = threshold), linetype = "dashed", colour = "gray30", linewidth = 1.2) +
  facet_wrap(vars(focgene), nrow = 4) +
  theme_bw() +
  theme(plot.background = element_rect(fill = "#F2F6F9", color = "#F2F6F9"),
        legend.position = "none",
        axis.title.x=element_blank()) +
  ylab("p(zoonotic)")

ggsave("E:\\Working\\Talks\\2023\\11 Epidemics\\raw.png", plot = fig_results_subtype_PB2, width = 3.2, height = 2)
