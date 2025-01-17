# Modelling Avian Influenza Zoonotic Risk

> A machine learning model to identify risk of zoonosis from Avian Influenza Viruses using genomic data

## Aim

The aim of the project is to identify Avian Influenza strains with zoonotic potential - which are able to infect humans - using viral genomic sequence characteristics and machine learning. The project allows to monitor the risk of spillover of Avian Influenza from avian and mammal species to humans.

The model used in this project has been developed by [Dr Liam Brierley](https://www.gla.ac.uk/schools/infectionimmunity/staff/liambrierley/) at the University of Glasgow. Dr Brierley has an Academic Honorary Contract with UKHSA. This is a [link](https://vibelab.co.uk/) to Dr Brierley lab's page.

This project contributes to 

1. Strengthening pandemic preparedness and H5 monitoring thus ensuring the health and well being of the population.
2. Building key capabilities in AI and machine learning across the agency.

Project contact: [lorenzo.cattarino@ukhsa.gov.uk](mailto:lorenzo.cattarino@ukhsa.gov.uk)

## Data 

The project uses publicly available viral DNA sequences from NCBI GenBank and GISAID. These viral sequences are derived from biological samples taken from human or animal hosts. The sequences are uploaded by individuals (usually researchers from academic institutions or government and non-government organizations) to platforms such as  NCBI GenBank and GISAID where they are available to anyone.

Different sequence types are used: 

* entire nucleotide of segment (coding and non-coding parts)
* coding sequence
* protein sequence

Raw sequences are available as FASTA files.

The genomic data (raw sequences and processed outputs) are saved in the project's shared drive (`\\filepor10\DOP$\X037_AVI_GeneticMarkers`) and also in the project directory on the HPC (`/data/projects/zoonosis-risk-ai/zoonosis-risk-ai-modelling`).

## Analysis pipeline 

The code pipeline in this repository consists of three main steps:

1. Pre-process the raw sequences (`process_GISAID_NCBI_data.R`). This step includes: 
    - Quality check 
    - Calculate the nucleotide/coding sequence features for the ML models (e.g., k-mers)
    - Cluster raw sequences (to reduce overfitting as some lineages are oversampled and to reduce phylogenetic bias as some sequences are really homogenous within subtype or lineage). To make sure models generalise well to unseen data, clusters are created by holding out (from the full dataset) an entire virus subtype each time. The sequences from the held-out subtype are used for testing the ML models.
    - Save processed data as labels and features:
        - Labels: contain all the sequences representative of each cluster for each hold-out subtype 
        - Features: 12 sets of features for each of the 8 genes considered (total of 96 tables) used for training the ML models. 
    Features are stored separately from holdout cluster labels to avoid saving copies of features columns for the same sequences across holdout clusters (improve efficiency of data storge). **PLEASE NOTE: Labels and features tables were shared by University of Glasgow and the `process_GISAID_NCBI_data.R` script was never run.** 

2. Calculate the protein features for the ML models (`protein_feat_extract.py `) using a well-documented existing package (iFeatureOmega) that's only available for Python 
3. Train Machine Learning models using 5 different ML algorithms (e.g., `build_glmnet_vectorised_barkla.R`) :

    - Penalised LASSO Regression
    - Random Forests
    - Linear Support Vector Machines
    - Radial Support Vector Machines
    - Boosted Classification Trees
    
For each ML algorithm, one separate model is trained for each of 13 hold out subtypes, 12 feature sets and 8 genes.

## Computing requirements and system access

To be able to run the code in this repository, you are strongly advised to use a High Performance Computer (HPC). The requirements for HPC and data access are:  

1. Follow the instructions in this [page](https://confluence.collab.test-and-trace.nhs.uk/display/DEDT/HPC) to request access to the HPC.
2. Once you are logged into the HPC (Zeus), navigate to the project directory `/data/projects/zoonosis-risk-ai` and GitHub repository `zoonosis-risk-ai-modelling` within.
3. Make sure you are using on the HPC the latest version of the code by pulling from the remote GitHub repo. To be able to work with git from the HPC (i.e., pull code developed locally) you might need to configure a Secure Shell Protocol (SSH) connection. Follow the instructions [here](https://confluence.collab.test-and-trace.nhs.uk/display/DEDT/HPC#HPC-Cloningtherepo) to do that.
4. The script for running the whole code routine is located in `S3/batch_scripts/run_full_analysis.sbatch`. This script contains also the specific resources required on the HPC for running the routine. Follow the instruction in this [page](https://confluence.collab.test-and-trace.nhs.uk/display/DEDT/HPC) to familiarize yourself with HPC resources and to submit the script to the HPC. The `S3/batch_scripts/run_full_analysis.sbatch` script calls the conda environment (*zoonosis-risk-ml*) used for this project. The environment has already been set up using the `S3/environment.yml` file which contains the necessary code dependencies (packages).
5. Request access to the project shared drive through this [LAMa agreement application](http://datascience.phe.gov.uk/Lama/SpecialProjects/SpecialProjectDetail?pro=X037&ver=1.0). The data used in this project have already been copied from this location to the HPC (so you do not need to do it) but it is advised that you request access to the space for future work and development. 

### Packages
The main R packages used and their purpose are as follows (must be installed through [Bioconductor](https://cran.r-project.org/web/packages/BiocManager/vignettes/BiocManager.html) or the Bioconda channel when running on the HPC):

- [coRdon](https://www.bioconductor.org/packages/release/bioc/vignettes/coRdon/inst/doc/coRdon.html) - for some simple shortcuts to calculating codon frequency, and also a simple function to read in nucleotide sequences (followed by process_NCBI_seq and process_GISAID_seq, which are my own functions to use regexp and text tools in tidyr to reconstruct the metadata from the FASTA header)
- [Biostrings](https://bioconductor.org/packages/release/bioc/html/Biostrings.html) - functions to calculate k-mer frequency, and additional functions for reading/converting nucleotide and protein sequences
- [ORFik](https://github.com/Roleren/ORFik) - finds open reading frames in a nucleotide sequence, used to define the coding sequences of GISAID sequences which are not retained when you download as FASTA (just the whole segment sequence)

## Repository ground rules

1.  No data, output or blob files (.CSV, .DOC, .XLS, .PDF, .HTML, .PNG, .JPEG or .SVG) should ever enter the repository
2.  Ensure the pre-commit hooks are working before trying to push any commits to the repository (Note: you may need a local installation of python, with R, Python and Git added to the local user environment variables or PATH)
3.  If a new folder is added where data could potentially be stored, ensure the `.gitignore` is updated to reflect this
4.  Any raw or intermediate data files created should be stored in the `\data` folder. No CSV/RDS/JSON/favourite_data_format in any folder aside `\data`, unless it is small and **open source**.
5.  All outputs (e.g images containing graphs etc) should be stored into the `\outputs` folder.
6.  All new `.R` files should be stored in `\R`
7. Pushing to *main* requires the code to be reviewed by a peer before entering *main*

## Links
[Project Confluence page](https://confluence.collab.test-and-trace.nhs.uk/display/DEDT/Avian+Influenza+Zoonotic+Risk+Modelling+ML)
