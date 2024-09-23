#################
# Load packages #
#################

library(dplyr)
library(tidyr)
library(magrittr)
library(ggplot2)


# Set options
results_date <- "16_02_24"
method = "xgb"

# Grid search validation results

gridsearch <- read.csv(paste0("S3\\analysis\\gridsearch_", results_date, ".csv"))

if (method == "rf"){
  
  fig_grid_acc <- gridsearch %>%
    group_by(cluster_set, featset, focgene, min.node.size) %>%
    summarise_at(vars("Accuracy"), list(med = median, min = min, max = max)) %>%
    as.data.frame %>%
    mutate(min.node.size = factor(min.node.size)) %>%
    ggplot(aes(x = min.node.size, y = med, ymin = min, ymax = max, color = featset)) +
    geom_errorbar(width=0, position=position_dodge(0.6)) +
    geom_point(position=position_dodge(0.6)) +
    facet_grid(focgene ~ cluster_set) +
    theme_bw() +
    ylab("Accuracy") +
    xlab("Min node size")
  
  fig_grid_kappa <- gridsearch %>%
    group_by(cluster_set, featset, focgene, min.node.size) %>%
    summarise_at(vars("Kappa"), list(med = median, min = min, max = max)) %>%
    as.data.frame %>%
    mutate(min.node.size = factor(min.node.size)) %>%
    ggplot(aes(x = min.node.size, y = med, ymin = min, ymax = max, color = featset)) +
    geom_errorbar(width=0, position=position_dodge(0.6)) +
    geom_point(position=position_dodge(0.6)) +
    facet_grid(focgene ~ cluster_set) +
    theme_bw() +
    ylab("Kappa") +
    xlab("Min node size")
  
} else if (method == "glmnet"){
  
  fig_grid_acc <- gridsearch %>%
    group_by(alpha, featset, focgene, lambda) %>%
    summarise_at(vars("Accuracy"), list(med = median, min = min, max = max)) %>%
    as.data.frame %>%
    mutate(lambda = factor(lambda)) %>%
    ggplot(aes(x = lambda, y = med, ymin = min, ymax = max, color = featset)) +
    geom_errorbar(width=0, position=position_dodge(0.6)) +
    geom_point(position=position_dodge(0.6)) +
    facet_grid(focgene ~ alpha) +
    theme_bw() +
    ylab("Accuracy") +
    xlab("Params")
  
  fig_grid_kappa <- gridsearch %>%
    group_by(alpha, featset, focgene, lambda) %>%
    summarise_at(vars("Kappa"), list(med = median, min = min, max = max)) %>%
    as.data.frame %>%
    mutate(lambda = factor(lambda)) %>%
    ggplot(aes(x = lambda, y = med, ymin = min, ymax = max, color = featset)) +
    geom_errorbar(width=0, position=position_dodge(0.6)) +
    geom_point(position=position_dodge(0.6)) +
    facet_grid(focgene ~ alpha) +
    theme_bw() +
    ylab("Kappa") +
    xlab("Params")
  
} else if (method == "svm"){
  
  fig_grid_acc <- gridsearch %>%
    mutate(sigma_C = paste(sigma, C, sep = "_")) %>%
    group_by(Weight, featset, focgene, sigma_C) %>%
    summarise_at(vars("Accuracy"), list(med = median, min = min, max = max)) %>%
    as.data.frame %>%
    mutate(sigma_C = factor(sigma_C)) %>%
    ggplot(aes(x = sigma_C, y = med, ymin = min, ymax = max, color = featset)) +
    geom_errorbar(width=0, position=position_dodge(0.6)) +
    geom_point(position=position_dodge(0.6)) +
    facet_grid(focgene ~ Weight) +
    theme_bw() +
    ylab("Accuracy") +
    xlab("Params")
  
  fig_grid_kappa <- gridsearch %>%
    mutate(sigma_C = paste(sigma, C, sep = "_")) %>%
    group_by(Weight, featset, focgene, sigma_C) %>%
    summarise_at(vars("Kappa"), list(med = median, min = min, max = max)) %>%
    as.data.frame %>%
    mutate(sigma_C = factor(sigma_C)) %>%
    ggplot(aes(x = sigma_C, y = med, ymin = min, ymax = max, color = featset)) +
    geom_errorbar(width=0, position=position_dodge(0.6)) +
    geom_point(position=position_dodge(0.6)) +
    facet_grid(focgene ~ Weight) +
    theme_bw() +
    ylab("Kappa") +
    xlab("Params")
  
} else if (method == "svmlin"){
  
  fig_grid_acc <- gridsearch %>%
    group_by(weight, featset, focgene, cost) %>%
    summarise_at(vars("Accuracy"), list(med = median, min = min, max = max)) %>%
    as.data.frame %>%
    mutate(cost = factor(cost)) %>%
    ggplot(aes(x = cost, y = med, ymin = min, ymax = max, color = featset)) +
    geom_errorbar(width=0, position=position_dodge(0.6)) +
    geom_point(position=position_dodge(0.6)) +
    facet_grid(focgene ~ weight) +
    theme_bw() +
    ylab("Accuracy") +
    xlab("Cost")
  
  fig_grid_kappa <- gridsearch %>%
    group_by(weight, featset, focgene, cost) %>%
    summarise_at(vars("Kappa"), list(med = median, min = min, max = max)) %>%
    as.data.frame %>%
    mutate(cost = factor(cost)) %>%
    ggplot(aes(x = cost, y = med, ymin = min, ymax = max, color = featset)) +
    geom_errorbar(width=0, position=position_dodge(0.6)) +
    geom_point(position=position_dodge(0.6)) +
    facet_grid(focgene ~ weight) +
    theme_bw() +
    ylab("Kappa") +
    xlab("Cost")
  
} else if (method == "xgb"){
  
  fig_grid_acc <- gridsearch %>%
    mutate(depth_weight = paste(max_depth, min_child_weight, sep = "_")) %>%
    group_by(eta, featset, focgene, depth_weight) %>%
    summarise_at(vars("Accuracy"), list(med = median, min = min, max = max)) %>%
    as.data.frame %>%
    mutate(depth_weight = factor(depth_weight)) %>%
    ggplot(aes(x = depth_weight, y = med, ymin = min, ymax = max, color = featset)) +
    geom_errorbar(width=0, position=position_dodge(0.6)) +
    geom_point(position=position_dodge(0.6)) +
    facet_grid(focgene ~ eta) +
    theme_bw() +
    ylab("Accuracy") +
    xlab("Params")
  
  fig_grid_kappa <- gridsearch %>%
    mutate(depth_weight = paste(max_depth, min_child_weight, sep = "_")) %>%
    group_by(eta, featset, focgene, depth_weight) %>%
    summarise_at(vars("Kappa"), list(med = median, min = min, max = max)) %>%
    as.data.frame %>%
    mutate(depth_weight = factor(depth_weight)) %>%
    ggplot(aes(x = depth_weight, y = med, ymin = min, ymax = max, color = featset)) +
    geom_errorbar(width=0, position=position_dodge(0.6)) +
    geom_point(position=position_dodge(0.6)) +
    facet_grid(focgene ~ eta) +
    theme_bw() +
    ylab("Kappa") +
    xlab("Params")
  
}

ggsave(paste0("S3\\figures_tables\\grid_acc_",results_date,".png"), plot = fig_grid_acc, width = 20, height = 10)
ggsave(paste0("S3\\figures_tables\\grid_kappa_",results_date,".png"), plot = fig_grid_kappa, width = 20, height = 10)

# Performance metrics on holdout sets

results <- read.csv(paste0("S3\\analysis\\results_", results_date, ".csv"), na.strings = "NaN")

# fig_results_line <- results %>%
#   gather(metric, value, -cluster_set, -featset, -focgene) %>%
#   filter(metric %in% c("AUC", "F1", "Sensitivity")) %>%
#   mutate(value = as.numeric(value)) %>%
#   ggplot(aes(x = cluster_set, y = value, color = featset, group = featset)) + 
#   geom_point(size = 0.5) +
#   geom_line() +
#   facet_grid(metric ~ focgene, scales = "free") +
#   theme_bw() +
#   theme(axis.title.y=element_blank())
# 
# ggsave(paste0("S3\\figures_tables\\fig_results_line_",results_date,".png"), plot = fig_results_line, width = 22, height = 6)

# fig_results_heat_AUC <- results %>%
#   gather(metric, value, -cluster_set, -featset, -focgene) %>%
#   filter(metric %in% c("AUC")) %>%
#   mutate(value = as.numeric(value)) %>%
#   ggplot(aes(x = featset, y = cluster_set, fill = value)) + 
#   geom_tile(color="white") +
#   scale_fill_distiller("AUC", palette = "RdBu", limits = c(0,1)) +
#   #  scale_fill_viridis_c("AUC") +
#   facet_wrap(~ focgene, nrow=2) +
#   theme_bw() +
#   theme(axis.text.x=element_text(angle = 90, hjust = 0))

# fig_results_heat_F1 <- results %>%
#   gather(metric, value, -cluster_set, -featset, -focgene) %>%
#   filter(metric %in% c("F1")) %>%
#   mutate(value = as.numeric(value)) %>%
#   ggplot(aes(x = featset, y = cluster_set, fill = value)) + 
#   geom_tile(color="white") +
#   scale_fill_viridis_c("F1", limits = c(0,1)) +
#   facet_wrap(~ focgene, nrow=2) +
#   theme_bw() +
#   theme(axis.text.x=element_text(angle = 90, hjust = 0))

# fig_results_heat_sens <- results %>%
#   gather(metric, value, -cluster_set, -featset, -focgene) %>%
#   filter(metric %in% c("Sensitivity")) %>%
#   mutate(value = as.numeric(value)) %>%
#   ggplot(aes(x = featset, y = cluster_set, fill = value)) + 
#   geom_tile(color="white") +
#   scale_fill_viridis_c("Sensitivity", limits = c(0,1)) +
#   facet_wrap(~ focgene, nrow=2) +
#   theme_bw() +
#   theme(axis.text.x=element_text(angle = 90, hjust = 0))

# ggsave(paste0("S3\\figures_tables\\fig_results_heat_AUC_",results_date,".png"), plot = fig_results_heat_AUC, width = 12, height = 4)
# ggsave(paste0("S3\\figures_tables\\fig_results_heat_F1_",results_date,".png"), plot = fig_results_heat_F1, width = 12, height = 4)
# ggsave(paste0("S3\\figures_tables\\fig_results_heat_sens_",results_date,".png"), plot = fig_results_heat_sens, width = 12, height = 4)

bestmodels <- read.csv(paste0("S3\\analysis\\results_", results_date, ".csv"), na.strings = "NaN") %>% 
  filter(cluster_set == "70_7") %>%
  group_by(focgene) %>% 
  slice_max(AUC) %>% 
  slice_max(F1)   # If multiple models tie in AUC, pick the one with highest F1 score

allflu_wgs_ref <- read.csv("S3\\data\\full\\allflu_wgs_ref.csv")
holdout_zoon <- c("H7N9", "H5N1", "H9N2", "H5N6", "H10N8", "H7N3", "H3N8", "H7N7", "H7N4")
holdout_nz <- c("H4N6", "H16N3", "H4N8", "H8N4")
holdouts <- c(holdout_zoon, holdout_nz)

results_subtype <- read.csv(paste0("S3\\analysis\\resultsbysubtype_", results_date, ".csv"), na.strings = "NaN") %>% 
  filter(cluster_set == "70_7") %>%
  left_join(allflu_wgs_ref %>%
              filter(subtype %in% holdout_zoon & label == "zoon"|subtype %in% holdout_nz) %>%
              count(subtype))

# # Boxplots over clusters: best performing feature set only
# fig_results_subtype_boxes <- bestmodels %>% 
#   select(featset, focgene) %>%
#   left_join(results_subtype) %>%
#   mutate(subtype = case_when(n < 10 ~ "rare subtypes (< 10)",
#                              TRUE ~ subtype),
#          subtype = factor(subtype, levels = c("H5N1", "H5N6", "H7N9", "H9N2", "rare subtypes (< 10)", "H4N6", "H4N8", "H8N4", "H16N3"))) %>%
#   group_by(cluster_set, featset, focgene, subtype) %>%
#   summarise(accuracy = mean(accuracy)) %>%
#   mutate(label = case_when(subtype %in% holdout_nz ~ "avian",
#                            TRUE ~ "zoon")) %>%
#   ggplot(aes(x = focgene, y = accuracy, fill = label)) +
#   geom_boxplot(alpha = 0.9) +
#   facet_wrap(vars(subtype), nrow = 2) +
#   theme_bw() +
#   ylab("Accuracy") +
#   xlab("Protein model")
# 
# ggsave(paste0("S3\\figures_tables\\fig_results_subtype_boxes_",results_date,".png"), plot = fig_results_subtype_boxes, width = 12, height = 4)
# 
# fig_results_subtype_lines <- bestmodels %>% 
#   select(featset, focgene) %>%
#   left_join(results_subtype) %>%
#   mutate(subtype = case_when(n < 10 ~ "rare subtypes (< 10)",
#                              TRUE ~ subtype),
#          subtype = factor(subtype, levels = c("H5N1", "H5N6", "H7N9", "H9N2", "rare subtypes (< 10)", "H4N6", "H4N8", "H8N4", "H16N3"))) %>%
#   group_by(cluster_set, featset, focgene, subtype) %>%
#   summarise(accuracy = mean(accuracy)) %>%
#   ggplot(aes(x = cluster_set, y = accuracy, color = featset, group = featset)) + 
#   geom_point(size = 2, alpha = 0.8) +
#   facet_grid(focgene ~ subtype, scales = "free") +
#   geom_line() +
#   theme_bw() +
#   theme(axis.title.y=element_blank())
# 
# ggsave(paste0("S3\\figures_tables\\fig_results_subtype_lines_",results_date,".png"), plot = fig_results_subtype_lines, width = 24, height = 11)

# bestmodels %>% 
#   select(featset, focgene) %>%
#   left_join(results_subtype) %>%
#   mutate(subtype = case_when(n < 10 ~ "rare subtypes (< 10)",
#                              TRUE ~ subtype),
#          subtype = factor(subtype, levels = c("H5N1", "H5N6", "H7N9", "H9N2", "rare subtypes (< 10)", "H4N6", "H4N8", "H8N4", "H16N3"))) %>%
#   group_by(cluster_set, featset, focgene, subtype) %>%
#   summarise(accuracy = mean(accuracy)) %>%
#   ungroup %>% 
#   mutate(set = paste0(subtype, featset, focgene)) %>%
#   ggplot(aes(x = cluster_set, y = accuracy)) + 
#   geom_boxplot(alpha = 0.9) +
#   theme_bw() +
#   guides(color = FALSE) +
#   theme(axis.title.y=element_blank())

# # Bars: best performing model (feature set and cluster)
# fig_results_subtype_bars <- bestmodels %>% 
#   select(cluster_set, featset, focgene) %>%
#   left_join(results_subtype) %>%
#   mutate(subtype = case_when(n < 10 ~ "rare subtypes (< 10)",
#                              TRUE ~ subtype),
#          subtype = factor(subtype, levels = c("H5N1", "H5N6", "H7N9", "H9N2", "rare subtypes (< 10)", "H4N6", "H4N8", "H8N4", "H16N3"))) %>%
#   group_by(cluster_set, featset, focgene, subtype) %>%
#   summarise(accuracy = mean(accuracy)) %>%
#   mutate(label = case_when(subtype %in% holdout_nz ~ "avian",
#                            TRUE ~ "zoon")) %>%
#   ggplot(aes(x = focgene, y = accuracy, fill = label)) +
#   geom_bar(stat = "identity") +
#   facet_wrap(vars(subtype), nrow = 2) +
#   theme_bw() +
#   ylab("Accuracy") +
#   xlab("Protein model")

# Points: best performing model (feature set and cluster)
best_files <- bestmodels %>% 
  mutate(file = paste0("S3\\analysis\\subtyperaw\\subtypeacc_raw_", cluster_set, "_", featset, "_", focgene, "_", results_date, ".csv")) %>%
  pull(file)

bestmodels_raw <- best_files %>% purrr::map_dfr(function(x) read.csv(x, na.strings = "NaN")) %>% bind_rows

# # Find top predicted zoonotic capability among avian viruses
# top_zoon <- bestmodels_raw %>% filter(focgene %in% c("PB2", "NP") & label == "nz" & zoon > 0.4) %>% arrange(-zoon)
# 
# test_set_df <- lapply(sort(holdouts), function(x)   
#   allflu_wgs_ref %>% 
#     filter(subtype == x) %>%
#     filter(subtype %in% holdout_zoon & label == "zoon"|subtype %in% holdout_nz)
# ) %>% bind_rows()
# 
# top_zoon %<>% bind_cols(test_set_df %>% slice(top_zoon$X) %>% select(gid, title, src, date)) %>% select(-nz, -pred, -X)
# top_zoon
# 
# for (x in unique(top_zoon$subtype)){
#   allflu_wgs_ref %>% filter(subtype == x) %>% .[(top_zoon %>% filter(subtype == x) %>% pull(X)),] %>% print()
# }

fig_results_subtype_raw <- bestmodels_raw %>%
  mutate(focgene = factor(focgene),
         subtype = factor(subtype, levels = c("H3N8", "H5N1", "H5N6", "H7N3", "H7N4", "H7N7", "H7N9", "H9N2", "H10N8",
                                              "H4N6", "H4N8", "H8N4", "H16N3")),
         label = case_when(label == "nz" ~ "avian",
                           TRUE ~ "zoon")) %>%
  ggplot(aes(x = subtype, y = zoon, colour = label)) +
  geom_jitter(alpha = 0.4, width = 0.2) +
  geom_text(data = bestmodels %>% mutate(focgene = factor(focgene)), aes(label = featset, y = 0.95, x = 12.5, colour = "black"), size = 2.5, colour = "black") +
  geom_hline(data = bestmodels %>% mutate(focgene = factor(focgene)), aes(yintercept = read.csv("S3/analysis/stack_results.csv") %>% pull(threshold) %>% log()), linetype = "dashed", colour = "gray30", linewidth = 1.2) +
  facet_wrap(vars(focgene), nrow = 4) +
  theme_bw() +
  ylab("p(zoonotic)") +
  xlab("Subtype") +
  guides(color = "none")

ggsave(paste0("S3\\figures_tables\\fig_results_subtype_raw_",results_date,".png"), plot = fig_results_subtype_raw, width = 12, height = 8)


# results %>% group_by(focgene) %>% slice_max(Accuracy) %>% as.data.frame
# results %>% group_by(focgene) %>% slice_max(AUC) %>% as.data.frame
# # Performance density plots
# bind_cols(label = test_used %>% bind_rows %>% pull(label) %>% droplevels,
# zoon = predict_prob_test %>% pull(zoon)) %>%
# ggplot(aes(x = zoon, fill = label)) +
# geom_density(alpha = 0.4) +
# facet_wrap(~ label, scales = "free") +
# theme_bw()
# # Performance by date
# bind_cols(cluster_rep = test_used %>% bind_rows %>% pull(cluster_rep),
# zoon = predict_prob_test %>% pull(zoon)) %>%
# inner_join(read.csv("S3\\data\\full\\cluster_rep_labels.csv")) %>%
# ggplot(aes(x = as.Date(date), y = zoon, color = label)) +
# geom_point(alpha = 0.4) +
# theme_bw()
# # Performance by subtype
# bind_cols(cluster_rep = test_used %>% bind_rows %>% pull(cluster_rep),
# zoon = predict_prob_test %>% pull(zoon)) %>%
# inner_join(read.csv("S3\\data\\full\\cluster_rep_labels.csv")) %>%
# add_count(subtype) %>%
# filter(n > 25) %>%
# ggplot(aes(x = as.Date(date), y = zoon, color = label)) +
# geom_point(alpha = 0.4) +
# facet_wrap(~ subtype) +
# theme_bw() 



#################
# Stacked model #
#################

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
  ylab("log(p(zoonotic))") +
  xlab("Subtype") +
  guides(color = "none")

ggsave(paste0("S3\\figures_tables\\fig_results_stack_raw_",results_date,".png"), plot = fig_results_stack_raw, width = 10, height = 5.5)

# What models are being selected?

stacked_coef <- read.csv("S3/analysis/stack_coef.csv") %>%
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
  expand(param, subtype) %>% 
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
  write.csv("S3/analysis/stack_coef_table.csv")


# Variable importance (RF for now)

for (featset in list.files(path = paste0("S3\\data\\full\\mlready")) %>% gsub("allflu_|_pt.*.rds", "", .) %>% unique()){
  
  varimp <- read.csv(file=paste0("S3\\analysis\\varimp_", results_date, "_", featset, ".csv"),  na.strings = "NaN")
  varimp_order <- varimp %>%
    group_by(name) %>%
    summarise(mean = mean(relGini)) %>%
    arrange(-mean) %>%
    pull(name)
  
  # fig_varimp_errors <- varimp %>%
  #   filter(name %in% varimp_order[1:20]) %>%
  #   group_by(name, focgene) %>%
  #   summarise_at(vars("relGini"), list(med = median, lower = ~quantile(., probs = 0.25), upper = ~quantile(., probs = 0.75))) %>%
  #   as.data.frame %>%
  #   ggplot(aes(x = name, y = med, ymin = lower, ymax = upper, color = focgene)) +
  #   geom_errorbar(width=0, position=position_dodge(0.6)) +
  #   geom_point(position=position_dodge(0.6)) +
  #   scale_x_discrete(limits = varimp_order[1:20]) +
  #   theme_bw() +
  #   facet_wrap(~ focgene, nrow = 8) +
  #   ylab("Relative Gini decrease") +
  #   xlab("Feature")
  
  fig_varimp_heat <- varimp %>%
    filter(name %in% varimp_order[1:100]) %>%
    group_by(name, focgene) %>%
    summarise_at(vars("relGini"), list(mean = mean)) %>%
    as.data.frame %>%
    ggplot(aes(x = focgene, y = name, fill = mean)) +
    geom_tile(color=NA) +
    scale_fill_viridis_c("MRGD") +
    scale_x_discrete(limits = rev(c("HA", "M1", "NA", "NP", "NS1", "PA", "PB1", "PB2"))) +
    scale_y_discrete(limits = varimp_order[1:100]) +
    theme_bw() +
    coord_flip() +
    theme(panel.background = element_rect(fill="#440154"),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.text.x=element_text(angle = 90, hjust = 0)) +
    xlab("Gene") +
    ylab("Feature")
  
  ggsave(paste0("S3\\figures_tables\\fig_varimp_heat_", featset, "_", results_date, ".png"), plot = fig_varimp_heat, width = 18, height = 2.5)
  
  # varimp %>%
  #   group_by(name, focgene) %>%
  #   summarise_at(vars("relGini"), list(mean = mean)) %>%
  #   pivot_wider(names_from = focgene, values_from = mean) %>%
  #   as.data.frame %>%
  #   pivot_longer(HA:PB2, names_to="x_var", values_to="x_val") %>% 
  #   pivot_longer(HA:PB2, names_to="y_var", values_to="y_val") %>% 
  #   nest(data=c(x_val, y_val)) %>%
  #   mutate(cor = map(data, ~cor(.x$x_val, .x$y_val))) %>% 
  #   unnest(cor)
  
  cor_mat <- varimp %>%
    group_by(name, focgene) %>%
    summarise_at(vars("relGini"), list(mean = mean)) %>%
    pivot_wider(names_from = focgene, values_from = mean) %>%
    as.data.frame %>%
    select(-name) %>%
    cor(use="complete.obs")
  
  cor_mat[upper.tri(cor_mat)]<- NA
  cor_mat[diag(cor_mat)]<- NA
  
  fig_varimp_cor <- cor_mat %>%
    as.data.frame %>%
    tibble::rownames_to_column("focgene_x") %>%
    pivot_longer(cols = c("HA", "M1", "NA", "NP", "NS1", "PA", "PB1", "PB2"),
                 names_to = "focgene_y",
                 values_to = "cor") %>%
    as.data.frame %>%
    mutate(cor = na_if(cor, 1)) %>%
    ggplot(aes(x = focgene_x, y = focgene_y, fill = cor, label = round(cor,2))) +
    geom_tile(color=NA) +
    geom_text(color = "black", size = 3) +
    scale_fill_distiller("rho", palette = "RdBu", limits = c(-1,1)) +
    scale_x_discrete(limits = rev(c("HA", "M1", "NA", "NP", "NS1", "PA", "PB1", "PB2"))) +
    scale_y_discrete(limits = rev(c("HA", "M1", "NA", "NP", "NS1", "PA", "PB1", "PB2"))) +
    theme_bw() +
    theme(axis.title.x=element_blank(),
          axis.title.y=element_blank())
  
  ggsave(paste0("S3\\figures_tables\\fig_varimp_cor_", featset, "_", results_date, ".png"), plot = fig_varimp_cor, width = 4, height = 3)
  
}







# # Test plot
# varimp_order[1:50,] %>% ggplot(aes(x=name, y=mean)) +
#   geom_bar(stat="identity") +
#   scale_x_discrete(limits = (varimp_order %>% dplyr::slice(1:50) %>% arrange(mean) %>% pull(name))) +
#   coord_flip()

