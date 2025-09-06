# Code repository associated with Brierley et al., 'An AI for AI..'

This repository contains supporting data, code, and result files associated with <b>Brierley et al. (2025)</b>, <i>&quot;An AI for an AI: identifying zoonotic potential of avian influenza viruses via genomic machine learning.&quot;</i>

This work takes a large set of avian influenza virus genome sequences from both avian and human hosts (indicative of zoonotic spillover), clusters them, calculates various genomic and proteomic feature representations, and trains supervised machine learning models to identify zoonotic sequences.

Additional large data files containing the genome sequence feature sets and final trained stack ensemble models are available at [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17068424.svg)](https://doi.org/10.5281/zenodo.17068424). These files should nest within the folders herein when downloaded.


## scripts
### data_scripts/

- `00_startup_data_process.R` sets options for sequence data processing, and calls the remaining R scripts in turn.
- `01_functions.R` defines custom functions for data processing and cleaning
- `02_process_GISAID_NCBI_data.R` filters and cleans sequence data and corresponding metadata.
- `03_calc_feats.R` calculates features from sequences:
  - overlapping k-mers of segment nucleotide sequence (2 ≤ k ≤ 6);
  - genome composition of coding sequences (nucleotide bias, dinucleotide bias, Relative Synonymous Codon Usage, amino acid bias)
  - processes pre-calculated protein feature sequences (see below)
- `04_cluster_seqs.R` calls a separate installation of [MMseqs2](https://github.com/soedinglab/MMseqs2) to cluster whole genome sequences by shared identity
- `05_process_clusts.R` formats clustering outputs from MMseqs2 and reselects cluster representatives for mixed-label clusters.
- `protein_feat_extract.py` calls [iFeatureOmega-CLI](https://github.com/Superzchen/iFeatureOmega-CLI) to calculate functional physicochemistry measures of protein sequences (Conjoint Triad, Composition-Transition-Distribution, Pseudo-Amino Acid Composition) according to parameters of `protein_params.json`.

### ml_scripts/

- `01a-01f` R scripts construct individual machine learning models (12 feature sets * 8 gene/proteins * 13 holdout sets) to predict zoonotic status using five different binary classification algorithms, which are parallelised by default.
- The exception is `01e_create_training_fold_indices.R` which extracts the previously defined folds for cross-validation during parameter optimisation, to supply to XGBoost separately in `01f_build_xgb.R`
- `02_evaluate_validate.R` loads all individual models and evaluates performance on holdout subtypes.
- `03_stack_weighted_model.R` constructs a LASSO logistic regression `stack`, or meta-learner model, using inputs from the best individual models as new features, before evaluating performance on holdout subtypes
- `04_stack_weighted_varimp.R` calculates variable importance by permuting each raw feature used within models within the stack one-by-one; <b>note this can require lengthy computation and is best split into a batch process along the vector `varnames`</b>
- `05_figs_tables.R` generates the figures in folder "figures_tables"

## data

`fold_indices_list.rds` defines 5-fold cross-validation data folds for consistency between model runs
`allflu_wgs_ref.csv` defines ID, source, host label, and date for all sequences considered for analysis.
`holout_clusters/` contains cluster members and cluster representative IDs selected when sequences were clustered with different parameter sets and excluding different holdout subtypes.

## analysis/

Analytical result outputs describing:
- grid search optimisation of individual model hyperparameters
- performance of individual models on held out test subtypes
- performance of stack ensemble models on held out test subtypes
- permutation variable importance of genomic and proteomic features

## figures_tables

Contains figures as used in manuscript.
