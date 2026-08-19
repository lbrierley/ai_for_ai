#packages 
library(foreach)
library(parallel)
library(tidyr)
library(doParallel)

#Setting variables 
cluster_minseqid <- c(0.85, 0.90, 0.925, 0.95, 0.975, 0.99)
cluster_C <- c(0.7, 0.8, 0.90, 0.95)
combos <- crossing(seg = 1:8, minseqid = cluster_minseqid, C = cluster_C)

#Set parallelisation
# cl <- makePSOCKcluster(10)
# registerDoParallel(cl)
# clusterSetRNGStream(cl, 1429)

#Running code through alll 8 segments 
foreach(seg = 1:8) %:%
  foreach(minseqid_param = cluster_minseqid) %:%
  foreach(C_param = cluster_C) %do% {   #change to %dopar% for parrallel runs
    
    #Setting input and output folders 
    out_tsv <- paste0(
      "/hpscol02/tenant1/zoonosis-risk-ai/zoonosis-risk-ai-modelling/Clustering_RS_2026_03_25/Clustering_Results/clust_",
      seg, "_", minseqid_param * 100, "_", C_param * 10, "_cluster.tsv"
    )
    
    input_fasta <- paste0(
      "/hpscol02/tenant1/zoonosis-risk-ai/zoonosis-risk-ai-modelling/Clustering_RS_2026_03_25/Clustering_FASTAs/allflu_cds_seg_",
      seg, "_full.FASTA"
    )
    
    out_prefix <- paste0(
      "/hpscol02/tenant1/zoonosis-risk-ai/zoonosis-risk-ai-modelling/Clustering_RS_2026_03_25/Clustering_Results/clust_",
      seg, "_", minseqid_param * 100, "_", C_param * 10
    )
    
    #Running MMSeq2
    if (!file.exists(out_tsv)) {
      system2(
        "mmseqs",
        args = c(
          "easy-linclust",
          input_fasta,
          out_prefix,
          "/hpscol02/tenant1/zoonosis-risk-ai/zoonosis-risk-ai-modelling/Clustering_RS_2026_03_25",
          "--min-seq-id", as.character(minseqid_param),
          "-c", as.character(C_param),
          "--cov-mode", "0"
        ),
        stdout = TRUE,
        stderr = TRUE
      )
      
      #Deletes FASTAs created during clustering 
      file.remove(list.files(
        path = "/hpscol02/tenant1/zoonosis-risk-ai/zoonosis-risk-ai-modelling/Clustering_RS_2026_03_25/Clustering_Results",
        pattern = "\\.fasta$",
        full.names = TRUE
      ))
    }
  }

# stopCluster(cl)