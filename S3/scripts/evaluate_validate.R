#################
# Load packages #
#################

library(caret)
library(e1071)
library(matrixStats)
library(magrittr)
library(pROC)
library(randomForest)
library(janitor)
library(dplyr)
library(tidyr)
library(readr)
library(purrr)
library(forcats)
library(stringr)
library(tibble)
library(ranger)
# library(patchwork)
# library(aricode)


## Stray varimp order code

varimp_order_cds <- lapply(rf_list, function(x)
  varImp(x)$importance %>% rownames_to_column("name") %>% mutate(Overall = Overall/100) %>% rename(relGini = Overall)
) %>% 
  bind_rows() %>%
  group_by(name) %>% 
  summarise(mean = mean(relGini)) %>% 
  arrange(-mean) %>%
  pull(name)






# Read in list of most important variables for use in PD
# vars_to_pd<- readRDS("varimp_order.rds")[1:4]
# 
# # Load ML models
# load("S3\\data\\full\\listresults_ml_vector_cds_20_05_23_notune.RData")
# 
# # # Parameter optimisation
# # gridsearch <- lapply(rf_list, function(x)
# #   x$results
# # ) %>% bind_rows() %>% mutate(min.node.size = factor(min.node.size, levels = c("5", "12.5", "20")),
# #                              mtry = factor(mtry, levels = c("5", "12.5", "20")))
# 
# # Select holdout sets (each cluster in turn)
# validate_used <- lapply(unique(model_df$gid), function(x)
#   model_df %>% filter(gid == x)
# )
# 
# # Create validations
# predict_class_test <- Map(function(model, newdata)
#   
#   if (nrow(newdata)>0){
#     levels(model_df$outcome)[predict(model, newdata=newdata, type="raw")]
#   },
#   model = rf_list,
#   newdata = validate_used
#   
# ) %>% unlist %>% as.factor
# 
# predict_prob_test <- Map(function(model, newdata)
#   
#   if (nrow(newdata)>0){
#     predict(model, newdata=newdata, type="prob")
#   },
#   model = rf_list,
#   newdata = validate_used
# ) %>% bind_rows
# 
# matrix_test <- confusionMatrix(predict_class_test, validate_used %>% bind_rows %>% pull(outcome) %>% droplevels)
# 
# # Calculate variable importance
# varimp <- lapply(rf_list, function(x)
#   varImp(x)$importance %>% rownames_to_column("name") %>% mutate(Overall = Overall/100) %>% rename(relGini = Overall)
# ) %>% bind_rows() %>%
#   mutate(name =  str_replace_all(name, c("_Bias" = "", 
#                                          "T" = "U", 
#                                          "_p" = " \\(p",
#                                          "1$" = "1\\)",
#                                          "2$" = "2\\)",
#                                          "3$" = "3\\)")))
# 
# varimp_order <- varimp %>%                               
#   group_by(name) %>% 
#   summarise(mean = mean(relGini)) %>% 
#   mutate(name = fct_reorder(name, -mean)) %>% 
#   pull(name) %>% 
#   levels

# load("S3\\data\\full\\listresults_pd_spike_03_10_20.RData")

load("S3\\data\\full\\results_bundle_cds_21_05_23.RData")

cluster_cds_wgs_df <- get(load(file = "S3\\data\\full\\cluster_cds_wgs_df_16_5_23.RData"))




################################
# Save/plot tables and figures #
################################

# # Descriptive stats - mean, SD over each CORONAVIRUS (not sequence)
# enc_per_spp <- model_df_predownsample %>% 
#   filter(accessionversion %in% model_df$accessionversion) %>% 
#   group_by(genus, accessionversion) %>% 
#   summarise(spp_mean = mean(enc)) %>%
#   ungroup
# 
# bind_rows(enc_per_spp %>%
#             group_by(genus) %>%
#             summarise(mean = mean(spp_mean),
#                       sd = sd(spp_mean)), 
#           data.frame(genus = "Total", mean = mean(enc_per_spp$spp_mean), sd = sd(enc_per_spp$spp_mean))) %>%
#   write.csv("figures//Desc_ENC_spike.csv")

confusionMatrix(predict_prob_test %>% mutate(thresh = as.factor(ifelse(zoon > 0.3, "zoon", "nz"))) %>% pull(thresh), 
                validate_used %>% bind_rows %>% pull(outcome), positive="zoon")



# Table 1
matrix_test$byClass %>%
  round(., 3) %>%
  t() %>%
  cbind(., AUC = multiclass.roc(response = validate_used %>% 
                                  bind_rows() %>% 
                                  pull(outcome), predictor = predict_prob_test) %>% 
          .$auc %>% 
          as.numeric() %>% 
          round(3)) %>%
  write.csv("S3\\figures_tables\\Table_1.csv")

# Figure 2, variable importance
F2 <- ggplot(varimp %>%
                     mutate(name = factor(name, levels = varimp_order),
                            type = case_when(grepl("^[A|C|G|U]$", name) ~ "nucleotide biases",
                                             grepl("^[A|C|G|U][A|C|G|U] ", name) ~ "dinucleotide biases",
                                             grepl("^[A|C|G|U][A|C|G|U][A|C|G|U]$", name) ~ "codon biases")) %>%
                     filter(name %in% varimp_order[1:30]), 
                   aes(x=name, y=relGini, fill=type)) +
  stat_summary(geom = "bar", fun.y = "mean") +
  stat_summary(geom = "errorbar", fun.data=mean_sdl, fun.args=list(mult=1), lwd=0.25, width=0) +
  coord_flip(ylim=c(0,1.02), expand=FALSE) + 
  scale_y_continuous(breaks=seq(0,1,0.1)) +
  theme_bw(base_size = 11) + xlab("Genomic feature") + ylab("Mean rel. Gini decrease") +
  theme(legend.justification=c(1,1), legend.position=c(.98,.98), panel.spacing = unit(1.1, "lines"), legend.title = element_blank()) +
  scale_fill_manual(values=c("#F0E442","#CC79A7","#56B4E9"))

# # Figure 3, partial dependence
# for (i in 1:length(vars_to_pd_spike)){
#   
#   plot_df <- lapply(list_PD, function(x) x[[vars_to_pd_spike[i]]]) %>% bind_rows %>% 
#     rename_with(., ~ str_replace_all(., c("_Bias" = "", 
#                                           "T" = "U", 
#                                           "_p" = " \\(p",
#                                           "1$" = "1\\)",
#                                           "2$" = "2\\)",
#                                           "3$" = "3\\)",
#                                           "aves" = "bird", 
#                                           "camel" = "camelid", 
#                                           "yangbat" = "yangochiroptera", 
#                                           "yinbat" = "yinpterochiroptera")))
#   
#   plot_df_summ <- plot_df %>% melt(id.vars = names(plot_df)[1]) %>% group_by(!!sym(names(plot_df)[1]), variable) %>%
#     summarise(median = median(value), 
#               lower = quantile(value, probs = .025),
#               upper = quantile(value, probs = .975))
#   
#   #plot_df[,2:ncol(plot_df)] <- (plot_df[,2:ncol(plot_df)]/(1-plot_df[,2:ncol(plot_df)])) %>% mapply('/', ., (vector_probs/(1-vector_probs))) # calculate odds ratios
#   
#   assign(paste0("F3_spike_spike_", i),
#          ggplot(plot_df_summ, 
#                 aes(x = !!sym(names(plot_df)[1]))) +  
#            geom_ribbon(aes(fill=variable, ymin = lower, ymax = upper), alpha = 0.2) +
#            geom_line(aes(colour = variable, y = median), lwd = 1.5, alpha = 0.8) +
#            scale_x_continuous(expand=c(0,0)) +
#            scale_color_manual(values=cbbPalette_ordered_bw, name = "host") +
#            scale_fill_manual(values=cbbPalette_ordered_bw, name = "host") +
#            ylab("Probability") +
#            guides(fill = guide_legend(title = "Prediction"), colour = guide_legend(title = "Prediction")) +
#            #facet_wrap(~ variable) +
#            theme_bw())
# }
# 
# for (i in 1:length(vars_to_pd_wg)){
#   
#   plot_df <- lapply(list_PD, function(x) x[[vars_to_pd_wg[i]]]) %>% bind_rows %>% 
#     rename_with(., ~ str_replace_all(., c("_Bias" = "", 
#                                           "T" = "U", 
#                                           "_p" = " \\(p",
#                                           "1$" = "1\\)",
#                                           "2$" = "2\\)",
#                                           "3$" = "3\\)",
#                                           "aves" = "bird", 
#                                           "camel" = "camelid", 
#                                           "yangbat" = "yangochiroptera", 
#                                           "yinbat" = "yinpterochiroptera")))
#   
#   plot_df_summ <- plot_df %>% melt(id.vars = names(plot_df)[1]) %>% group_by(!!sym(names(plot_df)[1]), variable) %>%
#     summarise(median = median(value), 
#               lower = quantile(value, probs = .025),
#               upper = quantile(value, probs = .975))
#   
#   #plot_df[,2:ncol(plot_df)] <- (plot_df[,2:ncol(plot_df)]/(1-plot_df[,2:ncol(plot_df)])) %>% mapply('/', ., (vector_probs/(1-vector_probs))) # calculate odds ratios
#   
#   assign(paste0("F3_spike_wg_", i),
#          ggplot(plot_df_summ, 
#                 aes(x = !!sym(names(plot_df)[1]))) +  
#            geom_ribbon(aes(fill=variable, ymin = lower, ymax = upper), alpha = 0.2) +
#            geom_line(aes(colour = variable, y = median), lwd = 1.5, alpha = 0.8) +
#            scale_x_continuous(expand=c(0,0)) +
#            scale_color_manual(values=cbbPalette_ordered_bw, name = "host") +
#            scale_fill_manual(values=cbbPalette_ordered_bw, name = "host") +
#            ylab("Probability") +
#            guides(fill = guide_legend(title = "Prediction"), colour = guide_legend(title = "Prediction")) +
#            #facet_wrap(~ variable) +
#            theme_bw())
# }


valid_df_raw <- data.frame(outcome = validate_used %>% bind_rows() %>% pull(outcome),
                           gid = validate_used %>% bind_rows() %>% pull(gid),
                           subtype = validate_used %>% bind_rows() %>% pull(subtype),
                           predict_prob_test)

valid_df_ordered <- rbind(valid_df_raw %>% filter(outcome == "nz") %>% arrange(-nz),
                                valid_df_raw %>% filter(outcome == "zoon") %>% arrange(-nz)
) %>%
  group_by(outcome) %>% 
  mutate(index = row_number()) %>%
  ungroup %>%
  melt(id.vars = c("outcome", "gid", "subtype", "index")) %>%
  mutate_at(vars(outcome, variable), funs(str_replace_all(.,
                                                                c("nz" = "non-zoonotic", "^zoon$" = "zoonotic"))))

F4 <- valid_df_ordered %>% filter(variable == "zoonotic") %>%
  ggplot(aes(x = index, y = value, fill = variable)) + 
  geom_bar(position="stack", stat="identity", width=1.1) +
  scale_fill_manual(values = c("orange")) +
  geom_hline(yintercept=0.3, alpha = 0.4, color="grey20", lty="dashed", size=1) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  ylab("Probability zoonotic") +
  theme_bw() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank()) +
  guides(fill = FALSE) +
  facet_wrap(~ fct_rev(outcome), scales = "free_x", nrow=2,
             labeller = as_labeller(
               c(`zoonotic` = "zoonotic clusters (n = 88)",
                 `non-zoonotic` = "non-zoonotic clusters (n = 3958)")))

F4 %>% ggsave(filename = "figures//ViBiom_Figure.png", width = 5, height = 3)



valid_df_ordered_spike_genus <- rbind(valid_df_raw %>% filter(host_category == "aves") %>% arrange(genus, -aves),
                                      valid_df_raw %>% filter(host_category == "camel") %>% arrange(genus, -camel),
                                      valid_df_raw %>% filter(host_category == "carnivore") %>% arrange(genus, -carnivore),
                                      valid_df_raw %>% filter(host_category == "human") %>% arrange(genus, -human),
                                      valid_df_raw %>% filter(host_category == "rodent") %>% arrange(genus, -rodent),
                                      valid_df_raw %>% filter(host_category == "swine") %>% arrange(genus, -swine),
                                      valid_df_raw %>% filter(host_category == "yangbat") %>% arrange(genus, -yangbat),
                                      valid_df_raw %>% filter(host_category == "yinbat") %>% arrange(genus, -yinbat)) %>%
  mutate(accessionversion = factor(accessionversion, levels=accessionversion)) %>%
  group_by(host_category) %>% 
  mutate(index = row_number()) %>%
  ungroup %>%
  melt(id.vars = c("accessionversion", "host_category", "childtaxa_name", "index", "genus")) %>%
  mutate_at(vars(host_category, variable), funs(str_replace_all(.,
                                                                c("aves" = "bird", "camel" = "camelid", "yangbat" = "yangochiroptera", "yinbat" = "yinpterochiroptera"))))


known_combs <- valid_df_ordered_spike_genus %>% distinct(host_category, genus)
spike_genus_annots <- list()
for (i in 1:nrow(known_combs)){
  spike_genus_annots[[i]] <- data.frame(host_category = known_combs$host_category[i],
                                        genus = known_combs$genus[i],
                                        xmin = valid_df_ordered_spike_genus %>% 
                                          filter(host_category == known_combs$host_category[i] & genus == known_combs$genus[i]) %>% 
                                          pull(index) %>% min - 0.5,
                                        xmax = valid_df_ordered_spike_genus %>% 
                                          filter(host_category == known_combs$host_category[i] & genus == known_combs$genus[i]) %>% 
                                          pull(index) %>% max + 0.5,
                                        xmean = valid_df_ordered_spike_genus %>% 
                                          filter(host_category == known_combs$host_category[i] & genus == known_combs$genus[i]) %>% 
                                          pull(index) %>% mean
  ) %>%
    mutate(label = case_when(genus == "Alphacoronavirus" ~ "alpha",
                             genus == "Betacoronavirus" ~ "beta",
                             genus == "Gammacoronavirus" ~ "gamma",
                             genus == "Deltacoronavirus" ~ "delta",
                             genus == "unclassified" ~ "U")
    )
  
  
}

spike_genus_annots <- bind_rows(spike_genus_annots)

F4_spike_genus_d <-valid_df_ordered_spike_genus %>%
  ggplot(aes(x = accessionversion, y = value, fill = variable)) + 
  geom_bar(position="stack", stat="identity", width=1) +
  geom_rect(data=spike_genus_annots, 
            aes(x = NULL, y = NULL, fill = NULL, 
                xmin = xmin, xmax = xmax, ymin = 0, ymax = 1), colour = "black", fill = NA, lwd=0.5) +
  geom_text(data=spike_genus_annots, size = 3,
            aes(x = xmean, y = 0.05, label = label, fill = NULL)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = cbbPalette_ordered_bw) +
  ylab("Probability") +
  theme_bw() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank()) +
  guides(fill = guide_legend(title = "Prediction")) +
  facet_wrap(~ host_category, scales = "free_x", nrow=2)

valid_df_ordered_spike_spp <- rbind(rearrange_to_plot_fig_4_spp("aves"),
                                    rearrange_to_plot_fig_4_spp("camel"),
                                    rearrange_to_plot_fig_4_spp("carnivore"),
                                    rearrange_to_plot_fig_4_spp("human"),
                                    rearrange_to_plot_fig_4_spp("rodent"),
                                    rearrange_to_plot_fig_4_spp("swine"),
                                    rearrange_to_plot_fig_4_spp("yangbat"),
                                    rearrange_to_plot_fig_4_spp("yinbat")) %>%
  mutate(accessionversion = factor(accessionversion, levels=accessionversion)) %>%
  group_by(host_category) %>% 
  mutate(index = row_number()) %>%
  ungroup %>%
  select(-genus) %>%
  melt(id.vars = c("accessionversion", "host_category", "childtaxa_name", "index")) %>%
  mutate_at(vars(host_category, variable), funs(str_replace_all(.,
                                                                c("aves" = "bird", "camel" = "camelid", "yangbat" = "yangochiroptera", "yinbat" = "yinpterochiroptera"))))


known_combs <- valid_df_ordered_spike_spp %>% distinct(host_category, childtaxa_name)
spike_spp_annots <- list()
for (i in 1:nrow(known_combs)){
  spike_spp_annots[[i]] <- data.frame(host_category = known_combs$host_category[i],
                                      childtaxa_name = known_combs$childtaxa_name[i],
                                      xmin = valid_df_ordered_spike_spp %>% 
                                        filter(host_category == known_combs$host_category[i] & childtaxa_name == known_combs$childtaxa_name[i]) %>% 
                                        pull(index) %>% min - 0.5,
                                      xmax = valid_df_ordered_spike_spp %>% 
                                        filter(host_category == known_combs$host_category[i] & childtaxa_name == known_combs$childtaxa_name[i]) %>% 
                                        pull(index) %>% max + 0.5,
                                      xmean = valid_df_ordered_spike_spp %>% 
                                        filter(host_category == known_combs$host_category[i] & childtaxa_name == known_combs$childtaxa_name[i]) %>% 
                                        pull(index) %>% mean
  )
}

spike_spp_annots <- bind_rows(spike_spp_annots) %>% 
  mutate(childtaxa_name = factor(childtaxa_name, levels = unique(childtaxa_name))) %>%
  group_by(childtaxa_name) %>% 
  mutate(label = cur_group_id()) %>%
  ungroup

F4_spike_spp_d <- valid_df_ordered_spike_spp %>%
  ggplot(aes(x = accessionversion, y = value, fill = variable)) + 
  geom_bar(position="stack", stat="identity", width=1) +
  geom_rect(data=spike_spp_annots, 
            aes(x = NULL, y = NULL, fill = NULL, 
                xmin = xmin, xmax = xmax, ymin = 0, ymax = 1), colour = "black", fill = NA, lwd=0.3) +
  geom_text(data=spike_spp_annots, size = 1.4, angle=90,
            aes(x = xmean, y = -0.05, label = label, fill = NULL)) +
  scale_y_continuous(expand = c(0, 0.05)) +
  scale_fill_manual(values = cbbPalette_ordered_bw) +
  ylab("Probability") +
  theme_bw() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  guides(fill = guide_legend(title = "Prediction")) +
  facet_wrap(~ host_category, scales = "free_x", nrow=2)

# Figure 5, prediction for human epidemic viruses
F5_spike <- rbind(predict_prob_epidemic %>% filter(flag == "MERS-CoV") %>% arrange(-camel),
                  predict_prob_epidemic %>% filter(flag == "SARS-CoV") %>% arrange(-rodent),
                  predict_prob_epidemic %>% filter(flag == "SARS-CoV-2") %>% arrange(-yinbat)) %>%
  mutate(accessionversion = factor(accessionversion, levels=accessionversion)) %>%
  melt(id.vars = c("accessionversion", "taxid", "flag")) %>%
  mutate_at(vars(variable), funs(str_replace_all(.,
                                                 c("aves" = "bird", "camel" = "camelid", "yangbat" = "yangochiroptera", "yinbat" = "yinpterochiroptera")))) %>%
  
  ggplot(aes(x = accessionversion, y = value, fill = variable)) + 
  geom_bar(position="stack", stat="identity", width=1) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = cbbPalette_ordered_bw) +
  ylab("Probability") +
  theme_bw() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank()) +
  guides(fill = guide_legend(title = "Prediction")) +
  facet_grid(~ flag, scales = "free_x", space="free_x")

# predict_pang %>% arrange(-yinbat) %>%
#   mutate(accessionversion = factor(accessionversion, levels=accessionversion)) %>%
#   melt(id.vars = c("accessionversion", "taxid")) %>%
#   ggplot(aes(x = accessionversion, y = value, fill = variable)) + 
#   geom_bar(position="stack", stat="identity", width=1) +
#   scale_y_continuous(expand = c(0, 0)) +
#   scale_fill_manual(values = cbbPalette_ordered_bw) +
#   ylab("Probability") +
#   theme_bw() +
#   theme(axis.title.x=element_blank(),
#         axis.text.x=element_blank()) +
#   guides(fill = guide_legend(title = "Prediction"))


# Full data table, featuring misclassifications
validate_used %>%
  bind_rows() %>%
  add_column(prediction = predict_class_test) %>%
  cbind(round(predict_prob_test,5)) %>%
  left_join(df_tax_matcher, by = "accessionversion") %>%
  select(genus, childtaxa_name, accessionversion, name_txt, outcome, prediction, aves, camel, carnivore, human, rodent, swine, yangbat, yinbat) %>%
  rename(virus = childtaxa_name, accession = accessionversion, metadata_host = name_txt, host_category = outcome, bird = aves, camelid = camel, yangochiroptera = yangbat, yinpterochiroptera = yinbat) %>%
  mutate_at(vars(host_category, prediction), funs(str_replace_all(.,
                                                                  c("aves" = "bird", "camel" = "camelid", "yangbat" = "yangochiroptera", "yinbat" = "yinpterochiroptera")))) %>%
  write.csv("figures\\Supp_Data_spike.csv")

validate_epidemic %>%
  cbind(round(predict_prob_epidemic %>% select(aves:yinbat),5)) %>%
  select(genus, childtaxa_name, accessionversion, aves, camel, carnivore, human, rodent, swine, yangbat, yinbat) %>%
  rename(virus = childtaxa_name, accession = accessionversion, bird = aves, camelid = camel, yangochiroptera = yangbat, yinpterochiroptera = yinbat) %>%
  arrange(virus) %>%
  write.csv("figures\\Supp_Data_spike_epidemic.csv")

# Table S1
model_df %>% 
  group_by(outcome) %>% 
  summarise(n_sequences = n(), n_species = n_distinct(taxid)) %>% 
  write.csv("figures\\S1_Table_spike.csv")

# Numerical confusion matrices
matrix_test$table %>%  write.csv("figures\\conf_matrix_spike.csv")

# Figure S1 - dinuc biases versus position (boxplots or heatmap)
S1_spike <- model_df %>% 
  rename_with(., ~ str_replace_all(., c("_Bias" = "", 
                                        "T" = "U", 
                                        "_p" = " \\(p",
                                        "1$" = "1\\)",
                                        "2$" = "2\\)",
                                        "3$" = "3\\)"))) %>%
  select(matches("^[A|C|G|U][A|C|G|U] ")) %>% 
  melt %>%
  mutate(dinuc = substr(variable, 1, 2),
         position = substr(variable, 6, 6),
         position = case_when(position == "1" ~ "1-2",
                              position == "2" ~ "2-3",
                              position == "3" ~ "3-1")) %>%
  ggplot(aes(y = value, x = position, fill = position)) +
  geom_boxplot(outlier.shape=NA, width = 0.9) +
  facet_wrap(. ~ variable, scales="free_y") +
  theme_bw() +
  theme(axis.text.x=element_blank()) +
  xlab("") +
  ylab("Dinucleotide bias") +
  geom_hline(yintercept=1, alpha = 0.3, color="grey50", lty="dashed", size=1) +
  facet_wrap(~ dinuc)

# Figure S2 - DISCONTINUED

# Figure S3 - grid search results
S3_spike <- ggplot(gridsearch, aes(y = Accuracy, x = min.node.size, fill = mtry)) +
  geom_boxplot() +
  theme_bw() +
  ylab("Cross-validation (inner loop) accuracy")

####################
# Combined figures #
####################

# F2_spike + (F2_wg + guides(fill = "none")) + F2_comp +
#   plot_annotation(tag_levels = "A") +
#   plot_layout(ncol = 3) +
#   ggsave("figures//Figure_2.png", width = 14, height = 5)

F2_comp +
  ggsave("figures//Figure_2.png", width = 6, height = 6)

S_PD_a <- F3_spike_spike_1 + 
  (F3_spike_spike_2 + guides(colour = "none", fill = "none")) + 
  (F3_spike_spike_3 + guides(colour = "none", fill = "none")) + 
  (F3_spike_spike_4 + guides(colour = "none", fill = "none")) + 
  (F3_wg_spike_1 + guides(colour = "none", fill = "none")) + 
  (F3_wg_spike_2 + guides(colour = "none", fill = "none")) + 
  (F3_wg_spike_3 + guides(colour = "none", fill = "none")) + 
  (F3_wg_spike_4 + guides(colour = "none", fill = "none")) + 
  plot_annotation(tag_levels = "A") +
  plot_layout(ncol = 4, guides = "collect") &
  theme(legend.position = "bottom")

S_PD_a + ggsave("figures//Figure_S_PD_a.png", width = 16, height = 8)

S_PD_b <-   F3_wg_wg_1 +
  (F3_wg_wg_2 + guides(colour = "none", fill = "none")) +
  (F3_wg_wg_3 + guides(colour = "none", fill = "none")) +
  (F3_wg_wg_4 + guides(colour = "none", fill = "none")) +
  (F3_spike_wg_1 + guides(colour = "none", fill = "none")) + 
  (F3_spike_wg_2 + guides(colour = "none", fill = "none")) + 
  (F3_spike_wg_3 + guides(colour = "none", fill = "none")) + 
  (F3_spike_wg_4 + guides(colour = "none", fill = "none")) + 
  plot_annotation(tag_levels = "A") +
  plot_layout(ncol = 4, guides = "collect") &
  theme(legend.position = "bottom")

S_PD_b + ggsave("figures//Figure_S_PD_b.png", width = 16, height = 8)

F4 <- F4_spike_d + (F4_wg_d + guides(colour = "none", fill = "none")) + 
  plot_annotation(tag_levels = "A") + 
  plot_layout(nrow = 2) & 
  theme(legend.position = "right") 

F4 + ggsave("figures\\Figure_4_type_d.png", width = 15, height = 9)

F4_genus <- F4_spike_genus_d + (F4_wg_genus_d + guides(colour = "none", fill = "none")) + 
  plot_annotation(tag_levels = "A") + 
  plot_layout(nrow = 2) & 
  theme(legend.position = "right") 

F4_genus + ggsave("figures\\Figure_4_type_d_genus.png", width = 15, height = 9)

F4_spp <- F4_spike_spp_d + (F4_wg_spp_d + guides(colour = "none", fill = "none")) + 
  plot_annotation(tag_levels = "A") + 
  plot_layout(nrow = 2) & 
  theme(legend.position = "right") 

F4_spp + ggsave("figures\\Figure_4_type_d_spp.png", width = 17, height = 9)

spike_spp_annots %>% 
  select(label, childtaxa_name) %>% 
  distinct() %>% 
  rename(virus = childtaxa_name) %>% 
  write.table("figures//Figure_4_labels.txt", sep = "   ", row.names = FALSE)

F5 <- F5_spike + (F5_wg + guides(colour = "none", fill = "none")) + 
  plot_annotation(tag_levels = "A") + 
  plot_layout(nrow = 2) & 
  theme(legend.position = "right") 

F5 + ggsave("figures\\Figure_5.png", width = 10.5, height = 5)

S1 <- S1_spike + (S1_wg + guides(fill = "none")) +
  plot_annotation(tag_levels = "A") +
  plot_layout(ncol = 2, guides = "collect") &
  theme(legend.position = "bottom")

S1 + ggsave("figures//S2_Figure.png", width = 10, height = 6)

S3 <- S3_spike + (S3_wg + guides(fill = "none")) +
  plot_annotation(tag_levels = "A") +
  plot_layout(ncol = 2, guides = "collect") &
  theme(legend.position = "bottom")

S3 + ggsave("figures//S3_Figure.png", width = 10, height = 4)