read_RDS <- function(x, focgene, n_col) {
  
  df <- readRDS(x) %>% 
    select(-any_of(c("segment", "cds_id", "enc", "GC_content"))) %>%
    rename_with(~paste(., focgene, sep = "_"), -c(gid))
  
  if (!is.null(n_col)) {
    
    # select n_col features at random
    # first column is the sequence ID
    n_col <- min(n_col, ncol(df)-1)
    idx <- sample(2:ncol(df), n_col)
    
    df <- df[, c(1, idx)]
    
  }
  
  df
  
}

read_featset <- function(focgene, n_col=NULL) {
  
  feature_sets_files <- list.files(
    path = "S3/data/full/mlready", 
    pattern = paste0("(nuc_2mer|prot_2mer|ctdc|ctdt|ctdd|pseaac|ctriad).+", focgene), 
    full.names = TRUE)
  
  feature_sets <- map(feature_sets_files, read_RDS, focgene = focgene, n_col = n_col) %>%
    reduce(left_join, by = "gid")
  
}
