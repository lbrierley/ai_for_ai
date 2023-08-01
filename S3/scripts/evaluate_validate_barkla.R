#################
# Load packages #
#################

library(caret)
library(e1071)
library(matrixStats)
# library(patchwork)
library(pROC)
library(randomForest)
# library(aricode)

# Read in list of most important variables for use in PD
# vars_to_pd<- readRDS("varimp_order.rds")[1:4]

# Load ML models
load("listresults_ml_vector_cds_20_05_23_notune.RData")

# # Parameter optimisation
# gridsearch <- lapply(rf_list, function(x)
#   x$results
# ) %>% bind_rows() %>% mutate(min.node.size = factor(min.node.size, levels = c("5", "12.5", "20")),
#                              mtry = factor(mtry, levels = c("5", "12.5", "20")))

# Select holdout sets (each cluster in turn)
validate_used <- lapply(unique(model_df$gid), function(x)
  model_df %>% filter(gid == x)
)

# Create validations
predict_class_test <- Map(function(model, newdata)
  
  if (nrow(newdata)>0){
    levels(model_df$outcome)[predict(model, newdata=newdata, type="raw")]
  },
  model = rf_list,
  newdata = validate_used
  
) %>% unlist %>% as.factor

predict_prob_test <- Map(function(model, newdata)
  
  if (nrow(newdata)>0){
    predict(model, newdata=newdata, type="prob")
  },
  model = rf_list,
  newdata = validate_used
) %>% bind_rows

matrix_test <- confusionMatrix(predict_class_test, validate_used %>% bind_rows %>% pull(outcome) %>% droplevels)

# Calculate variable importance
varimp <- lapply(rf_list, function(x)
  varImp(x)$importance %>% rownames_to_column("name") %>% mutate(Overall = Overall/100) %>% rename(relGini = Overall)
) %>% bind_rows() %>%
  mutate(name =  str_replace_all(name, c("_Bias" = "", 
                                         "T" = "U", 
                                         "_p" = " \\(p",
                                         "1$" = "1\\)",
                                         "2$" = "2\\)",
                                         "3$" = "3\\)")))

varimp_order <- varimp %>%                               
  group_by(name) %>% 
  summarise(mean = mean(relGini)) %>% 
  mutate(name = fct_reorder(name, -mean)) %>% 
  pull(name) %>% 
  levels

save(validate_used, matrix_test, varimp, varimp_order, predict_class_test, predict_prob_test,
     file=paste0("results_bundle_cds_", format(Sys.time(), "%d_%m_%y"), ".RData"))