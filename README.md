# Modelling Avian Influenza Zoonotic Risk

> A machine learning model to identify risk of zoonosis from Avian Influenza Viruses using genomic data

## Aim

The aim of the project is to identify Avian Influenza strains with zoonotic potential - which are able to infect humans - using viral genomic sequence characteristics and machine learning. The project allows to monitor the risk of spillover of Avian Influenza from avian and mammal species to humans.

The model used in this project has been developed by [Dr Liam Brierley](https://www.gla.ac.uk/schools/infectionimmunity/staff/liambrierley/) at the University of Glasgow. Dr Brierley has an Academic Honorary Contract with UKHSA. This is a [link](https://vibelab.co.uk/) to Dr Brierley lab's page.

This project contributes to 
1. Strengthening pandemic preparedness and H5 monitoring thus ensuring the health and well being of the population.
2. Building key capabilities in AI and machine learning across the agency.

Project contact: [lorenzo.cattarino@ukhsa.gov.uk](mailto:lorenzo.cattarino@ukhsa.gov.uk)

## Computing requirements and system access

To be able to run the code in this repository, you are strongly advised to use a High Performance Computer (HPC). The requirements for HPC and data access are:  

1. Follow the instructions in this [page](https://confluence.collab.test-and-trace.nhs.uk/display/DEDT/HPC) to request access to the HPC.
2. Once you are logged into the HPC (Zeus), navigate to the project space in `/data/projects/zoonosis-risk-ai`
3. The script for running the whole code routine is located in `S3/batch_scripts/run_full_analysis.sbatch`. This script contains also the specific resources required on the HPC for running the job. Follow the instruction in this [page](https://confluence.collab.test-and-trace.nhs.uk/display/DEDT/HPC) to familiarize yourself with HPC resources and to submit the script to the HPC.
4. Request access to the project shared drive through this [LAMa agreement application](http://datascience.phe.gov.uk/Lama/SpecialProjects/SpecialProjectDetail?pro=X037&ver=1.0). The data have already been copied from this location to the HPC (so you do not need to do it) but it is advised that you request access to the space for future work and development. 

## Data 

The project uses publicly available viral DNA sequences from NCBI GenBank and GISAID. These viral sequences are derived from biological samples taken from human or animal hosts. The sequences are uploaded by individuals (usually researchers from academic institutions or government and non-government organizations) to platforms such as  NCBI GenBank and GISAID where they are available to anyone.

The viral genome sequence data (raw sequences and processed outputs) are saved in the project's shared drive (*\\filepor10\DOP$\X037_AVI_GeneticMarkers*)
 
## Guidance/Ground rules

1.  No data, output or blob files (.CSV, .DOC, .XLS, .PDF, .HTML, .PNG, .JPEG or .SVG) should ever enter the repository
2.  Ensure the pre-commit hooks are working before trying to push any commits to the repository (Note: you may need a local installation of python, with R, Python and Git added to the local user environment variables or PATH)
3.  If a new folder is added where data could potentially be stored, ensure the `.gitignore` is updated to reflect this
4.  Any raw or intermediate data files created should be stored in the `\data` folder. No CSV/RDS/JSON/favourite_data_format in any folder aside `\data`, unless it is small and **open source** (such as the ONS population data), in which case it can be stored in `\data-raw`.
5.  All outputs (e.g images containing graphs etc) should be stored into the `\outputs` folder.
6.  All new `.R` files should be stored in `\R`
7. Pushing to *main* requires the code to be reviewed by a peer before entering *main*

## Links
[Project Confluence page](https://confluence.collab.test-and-trace.nhs.uk/display/DEDT/Avian+Influenza+Zoonotic+Risk+Modelling+ML)
