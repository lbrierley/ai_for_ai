# Load packages and set overall parameters
header(verbose, "Loading packages and define variables", padding=0)
source("S3/scripts/startup.R", echo = TRUE)

# Defined functions for processing sequence data
header(verbose, "Loading custom functions", padding=0)
source("S3/scripts/functions.R", echo = TRUE)

# # Process all training sequence data
# header(verbose, "Extracting and processing sequence data", padding=0)
# source("S3/scripts/process_GISAID_NCBI_data.R", echo = TRUE)

# protein_feat_extract.py
# protein_params.json
# Generate protein features

# Construct ML models for each feature set-gene combination. XGboost split into 5 scripts as longer run time.
header(verbose, "Train ML models", padding=0)
source("S3/scripts/build_glmnet_vectorised_barkla.R", echo = TRUE)
source("S3/scripts/build_rf_vectorised_barkla.R", echo = TRUE)
source("S3/scripts/build_svm_vectorised_barkla.R", echo = TRUE)
# source("S3/scripts/build_svmlin_vectorised_barkla.R", echo = TRUE)
# source("S3/scripts/build_xgb_vectorised_barkla.R", echo = TRUE)
# source("S3/scripts/build_xgb_vectorised_barkla2.R", echo = TRUE)
# source("S3/scripts/build_xgb_vectorised_barkla3.R", echo = TRUE)
# source("S3/scripts/build_xgb_vectorised_barkla4.R", echo = TRUE)
# 
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
