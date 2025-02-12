# Modelling Avian Influenza Zoonotic Risk

> A machine learning model to identify risk of zoonosis from Avian Influenza Viruses using genomic data

## Aim

The aim of the project is to identify Avian Influenza strains with zoonotic potential - which are able to infect humans - using viral genomic sequence characteristics and machine learning. The project allows to monitor the risk of spillover of Avian Influenza from avian and mammal species to humans.

The model used in this project has been developed by [Dr Liam Brierley](https://www.gla.ac.uk/schools/infectionimmunity/staff/liambrierley/) at the University of Glasgow. Dr Brierley has an Academic Honorary Contract with UKHSA. This is a [link](https://vibelab.co.uk/) to Dr Brierley lab's page.

This project contributes to 

1. Strengthening pandemic preparedness and H5 monitoring thus ensuring the health and well being of the population.
2. Building key capabilities in AI and machine learning across the agency.

Project contact: [lorenzo.cattarino@ukhsa.gov.uk](mailto:lorenzo.cattarino@ukhsa.gov.uk)

## Access requirements 

1. We use the UKHSA High Performance Computer (HPC) to run the code in this repository. Please email the HPC team (HPC@PHEcloud.onmicrosoft.com) to request an HPC account.
2. Install Putty and WinSCP from the company portal. 
3. As part of this project, the shared drive `\\filepor10\DOP$\X037_AVI_GeneticMarkers` is used. Please request access to the project shared drive through this [LAMa agreement application](http://datascience.phe.gov.uk/Lama/SpecialProjects/SpecialProjectDetail?pro=X037&ver=1.0). The data used in this project have already been copied from this location to the High Performance Computer (so you do not need to do it) but it is advised that you request access to the space for future work and development.

## Data 

The project uses viral genome sequences from NCBI GenBank and GISAID. Data (both raw sequences and processed outputs) are saved in the project shared drive and also in the project directory on the Porton Down HPC (`/data/projects/zoonosis-risk-ai/zoonosis-risk-ai-modelling/data`).

**03/02/2025 - NOTE: Due a water leak issue in Porton Down, the project directory has been migrated to the HPC cluster in Colindale (path: `/hpscol02/tenant1/zoonosis-risk-ai`)**    

## How to run the code
Log into the HPC using Putty. Please see this [HPC guidance](https://confluence.collab.test-and-trace.nhs.uk/display/DEDT/High+Performance+Computer) for more details. 

Navigate to the project GitHub repository in the project directory and update the main branch using Putty, with 

```
cd /data/projects/zoonosis-risk-ai/zoonosis-risk-ai-modelling
git pull
```

To be able to work with git from the HPC (i.e., pull code developed locally) you might need to configure a Secure Shell Protocol (SSH) connection. Follow the instructions [here](https://confluence.collab.test-and-trace.nhs.uk/display/DEDT/HPC#HPC-Cloningtherepo) to do that.

The full routine can be run using the `run_full_analysis.sbatch` batch script. The script contains the specific resources required on the HPC for running the full routine. Please see this [HPC guidance](https://confluence.collab.test-and-trace.nhs.uk/display/DEDT/HPC) to familiarize yourself with HPC resources and to submit the script to the HPC. It also calls the conda environment (*zoonosis-risk-ml*) used for this project. The environment has already been set up using the `S3/environment.yml` file which contains the necessary code dependencies (packages). **Please note this script currently can take up to 2 days to run.** You can run the whole routine from Putty with: 

```
sbatch S3/batch_scripts/run_full_analysis.sbatch
```

Sometimes you might want to run only part of the routine. To do that, you can edit locally the `run_full_analysis.sbatch` file by commenting off any lines corresponding to the scripts you do not want to re-run. Then you commit the change. After synchronizing the Git tree on the HPC with the remote tree (`git pull`), you can resubmit the batch job.

## Links
[Project Confluence page](https://confluence.collab.test-and-trace.nhs.uk/display/DEDT/Avian+Influenza+Zoonotic+Risk+Modelling+ML)
