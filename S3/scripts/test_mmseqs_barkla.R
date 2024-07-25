# Cluster parameters
minseqid_param <- 0.70
C_param <- 0.6

system(paste0("/users/lbrier/mmseqs-win64/mmseqs.bat easy-linclust /users/lbrier/allflu_nuc_wgs.FASTA /users/lbrier/wgs_", minseqid_param*100, "_", C_param*10, " tmp --min-seq-id ", minseqid_param ," -c ", C_param, " --cov-mode 0"),
       intern = TRUE, show.output.on.console = TRUE, ignore.stdout=FALSE, ignore.stderr=FALSE, wait=FALSE)