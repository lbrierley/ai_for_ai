library(R.utils)

# Load packages and set overall parameters
header(verbose, "Loading packages and defining variables", padding = 1)
source("S3/scripts/startup.R", echo = TRUE)

# Defined functions for processing sequence data
header(verbose, "Loading custom functions", padding = 1)
source("S3/scripts/functions.R", echo = TRUE)

# # Process all training sequence data
# header(verbose, "Extracting and processing sequence data", padding = 1)
# source("S3/scripts/process_GISAID_NCBI_data.R", echo = TRUE)

# protein_feat_extract.py
# protein_params.json
# Generate protein features

# Define same folds in 5-fold cross validation for use across all ML algorithms
header(verbose, "Creating 5-fold cross validation indices", padding = 1)
source("S3/scripts/create_training_fold_indices.R", echo = TRUE)

# Construct ML models for each feature set-gene combination. XGboost split into 5 scripts as longer run time.
header(verbose, "Training a Lasso Regularised Linear Model", padding = 1)
source("S3/scripts/build_glmnet_vectorised_barkla.R", echo = TRUE)
header(verbose, "Training a Random Forest Model", padding = 1)
source("S3/scripts/build_rf_vectorised_barkla.R", echo = TRUE)
header(verbose, "Training a Support Vector Machine Model with Non-Linear Kernel", padding = 1)
source("S3/scripts/build_svm_vectorised_barkla.R", echo = TRUE)
header(verbose, "Training a Support Vector Machine Model with Linear Kernel", padding = 1)
source("S3/scripts/build_svmlin_vectorised_barkla.R", echo = TRUE)
header(verbose, "Training an eXtreme Gradient Boosting Classification Tree Model - Part 1", padding = 1)
source("S3/scripts/build_xgb_vectorised_barkla.R", echo = TRUE)
header(verbose, "Training an eXtreme Gradient Boosting Classification Tree Model - Part 2", padding = 1)
source("S3/scripts/build_xgb_vectorised_barkla2.R", echo = TRUE)
header(verbose, "Training an eXtreme Gradient Boosting Classification Tree Model - Part 3", padding = 1)
source("S3/scripts/build_xgb_vectorised_barkla3.R", echo = TRUE)
header(verbose, "Training an eXtreme Gradient Boosting Classification Tree Model - Part 4", padding = 1)
source("S3/scripts/build_xgb_vectorised_barkla4.R", echo = TRUE)

# # Stack 96 best models (best algorithm for each feature set-gene combination) into a meta-learner or 'stack' models, trained on the outputted probabilities from those 96 models
# # Include class weighting in the stack model construction (i.e. upweight zoonotic as the rarer label in training data)
# source("S3/scripts/stack_weight_models_barkla.R", echo = TRUE)
# source("S3/scripts/stack_weight_varimp_barkla.R", echo = TRUE)
# 
# # Calculate model predictions on holdout sets and calculate metrics of performance
# source("S3/scripts/evaluate_validate_barkla.R", echo = TRUE)
# source("S3/scripts/evaluate_validate_local.R", echo = TRUE)
# 
# # Generate figures and tables
# source("S3/scripts/evaluate_validate_figs_tables.R", echo = TRUE)
# source("S3/scripts/evaluate_validate_figs_tables_poster.R", echo = TRUE)
