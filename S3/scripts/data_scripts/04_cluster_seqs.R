###########################################
# Cluster to collapse redundant sequences #
###########################################
##########################################################################################
# Calculate clustering of each segment's set of coding sequences based on k-mer overlaps #
##########################################################################################

# Read in reference of all sequences used
allflu_cds_ref <- read.csv("S3\\data\\segmentwise\\Raw_Sequences\\allflu_cds_df.csv", na.strings = "")

# Parameter sets to test
combos <- crossing(seg = 1:8, minseqid = cluster_minseqid, C = cluster_C)

# Call MMseqs2 14-7e284 to cluster on similarity score - replace with your own MMseqs2 path if needed
foreach(seg = 1:8) %:%
  foreach(minseqid_param = cluster_minseqid) %:%
  foreach(C_param = cluster_C) %do% {
    
    if(!file.exists(paste0("S3\\data\\segmentwise\\Clustering\\clust_", seg, "_", minseqid_param*100, "_", C_param*10, "_cluster.tsv"))){
      system(paste0("S3\\scripts\\mmseqs-win64\\mmseqs.bat easy-linclust S3\\data\\segmentwise\\allflu_cds_seg_", seg, "_full.FASTA S3\\data\\segmentwise\\Clustering\\clust_", seg, "_", minseqid_param*100, "_", C_param*10, " tmp --min-seq-id ", minseqid_param ," -c ", C_param, " --cov-mode 0"),
             intern = TRUE, show.output.on.console = TRUE, ignore.stdout=FALSE, ignore.stderr=FALSE, wait=FALSE)
      
      file.remove(list.files(path = "S3\\data\\segmentwise\\Clustering\\", pattern = "*.fasta", full.names = TRUE)) # Keep only the .tsv output from MMseqs2
    }
  }