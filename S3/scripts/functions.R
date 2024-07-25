####################
# Define functions #
####################

process_NCBI_seq <- function(x, label, type){
  df <- data.frame(title = x %>% names(), 
                   length = x %>% width()) 
  
  if (type == "nuc"){
    df %<>% tidyr::separate(title, sep = "\\|", into = c("title", "subtype", "date", "accession", "gene", "segment"), extra = "drop")
  } else if (type == "cds"){
    df %<>% tidyr::separate(title, sep = "\\|", into = c("title", "subtype", "date", "cds_id", "gene", "segment", "gb", "accession"), extra = "drop")
  } else if (type == "prot"){
    df %<>% tidyr::separate(title, sep = "\\|", into = c("title", "subtype", "date", "prot_id", "gene", "segment",  "gb", "accession"), extra = "drop")  
  } else {
    stop("invalid type (must be 'nuc', 'cds', or 'prot')")
  }
  
  df %<>%
    mutate(
      accession = gsub("\\:.*","",accession),
      gene = str_sub(str_sub(gene, 4, -1), 1, -2),
      label = label,
      src = "NCBI",
      subtype = case_when(                    # fix mislabelled records where known
        subtype == "H3N6,H3" ~ "H3N6",
        subtype == "H6N1,H6" ~ "H6N1",
        subtype == "H01N2" ~ "H1N2",
        title == "A/mallard/Korea/M198/2013" ~ "H4N2",
        title == "A/Duck/China/FJ2253/2014" ~ "H6N6",
        TRUE ~ toupper(subtype)
      ),
      title = case_when(
        title == "Influenza A virus" ~ paste0(title, subtype),
        title %in% c(
          "A/mallard/Italy/3401/2005",
          "A/Avian/Viet Nam/Egg/2014",
          "A/Avian/Viet Nam/Egg/2017",
          "A/chicken/China/embryonated chicken eggs/2014",
          "A/chicken/MUWRP-Uganda/853/2018",
          "A/wild duck/South Korea/KNU2018-26/2020"
        ) ~ paste0(title, date),
        TRUE ~ title),
      title = gsub(" ", "_", gsub("\\.[0-9]\\.A", "", title)),
      gid = paste0(title, str_sub(accession, 0, -8)), # Create genome id that can be used for both NCBI and GISAID data
      fastahead = x %>% names() %>% gsub(" ", "_", .),
      string = x %>% as.character(use.names=FALSE),
    ) %>% 
    select(-any_of("gb")) %>%
    relocate(fastahead, string, .after = last_col())
  
  if (type == "cds"|type == "prot"){
    df %<>%
      mutate(gene = case_when(
        segment == 1 & grepl("\\|PB2\\|", fastahead) ~ "PB2",                  # Attempt to assign based on fasta header first
        segment == 2 & grepl("\\|PB1-F2\\|", fastahead)  ~ "PB1-F2",                              
        segment == 2 & grepl("\\|PB1\\|", fastahead) ~ "PB1",  
        segment == 3 & grepl("\\|PA-X\\|", fastahead) ~ "PA-X",         
        segment == 3 & grepl("\\|PA\\|", fastahead) ~ "PA",                                                  
        segment == 4 & grepl("\\|HA\\|", fastahead) ~ "HA",
        segment == 5 & grepl("\\|NP\\|", fastahead) ~ "NP",
        segment == 6 & grepl("\\|NA\\|", fastahead) ~ "NA",
        segment == 7 & grepl("\\|M2\\|", fastahead) ~ "M2",           
        segment == 7 & grepl("\\|M1\\|", fastahead) ~ "M1",                                                   
        segment == 8 & grepl("\\|NS2\\|", fastahead) ~ "NS2",          
        segment == 8 & grepl("\\|NS1\\|", fastahead) ~ "NS1" ,      
        segment == 1 ~ "PB2",
        segment == 2 & !!type == "cds" & length < 500  ~ "PB1-F2",               # Else distinguish PB1-F2 based on size
        segment == 2 & !!type == "prot" & length < 166  ~ "PB1-F2",              # Else distinguish PB1-F2 based on size        
        segment == 2 ~ "PB1",                                                  # Else set to PB1
        segment == 3 & grepl("(gb|.*\\:.*\\,).*)",fastahead) ~ "PA-X",         # Distinguish PA-X based on joining ORFs
        segment == 3 ~ "PA",                                                   # Else set to PA
        segment == 4 ~ "HA",
        segment == 5 ~ "NP",
        segment == 6 ~ "NA",
        segment == 7 & grepl("(gb|.*\\:.*\\,).*)",fastahead) ~ "M2",           # Distinguish M2 based on joining ORFs
        segment == 7 ~ "M1",                                                   # Else set to M1
        segment == 8 & grepl("(gb|.*\\:.*\\,).*)",fastahead) ~ "NS2",          # Distinguish NS2 based on joining ORFs
        segment == 8 ~ "NS1"                                                   # Else set to NS1
      )) 
  }
  
  if (type == "nuc"|type == "cds"){
    df %<>% filter(!(grepl("\\|N40\\||\\|M42\\|", fastahead)))
    
  }
  
  return(df)
}

process_GISAID_seq <- function(x, label, type){
  df <- data.frame(title = x %>% names(), 
                   length = x %>% width()) 
  
  if (type == "nuc"){
    df %<>% tidyr::separate(title, sep = "\\|", into = c("title", "UID", "subtype", "null", "date", "INSDC", "accession", "title2", "gene", "segment"), extra = "drop")
  } else if (type == "prot"){
    df %<>% tidyr::separate(title, sep = "\\|", into = c("title", "UID", "subtype", "null", "date", "protINSDC", "protaccession", "gene"), extra = "drop")  
  } else {
    stop("invalid type (must be 'nuc' or 'prot')")
  }
  
  if (type == "nuc"){
    df %<>%
      mutate(
        accession = ifelse(accession %in% "", NA, accession),
        segment = case_when(
          gene == "PB2" ~ "1",
          gene == "PB1" ~ "2",
          gene == "PA" ~ "3",
          gene == "HA" ~ "4",
          gene == "NP" ~ "5",
          gene == "NA" ~ "6",
          gene == "MP" ~ "7",
          gene == "NS" ~ "8"
        ),
        subtype = str_sub(subtype, 5, -1),
        label = label,
        src = "GISAID",
        title = gsub(" ", "_", title),
        gid = UID,
        fastahead = x %>% names() %>% gsub(" ", "_", .),
        string = x %>% as.character(use.names=FALSE)
      ) %>%
      select(-any_of(c("null", "title2")))
  }
  
  if (type == "prot"){
    df %<>%
      mutate(
        subtype = str_sub(subtype, 5, -1),
        label = label,
        src = "GISAID",
        title = gsub(" ", "_", title),
        gid = UID,
        fastahead = x %>% names() %>% gsub(" ", "_", .),
      ) %>%
      select(-any_of(c("null")))
  }
  
  return(df)
}

bind_ORF <- function(x) {
  
  df <- findORFs(x, startCodon="ATG") %>%
    as.data.frame() 
  
  if(nrow(df) == 0){
    
    df <- data.frame(start = 0, end = 0) # Set start/end to zero if no valid open reading frames found
    
  } else {
    
    df <- df %>% 
      arrange(-width) %>%
      slice(1) %>%
      select(start, end)
  }   
  
  return(df)
}

calc_composition_counts <- function(x, codonpairs = FALSE){
  
  df <- cbind(data.frame(
    cds_id = x %>% names(),
    enc = x %>% codonTable() %>% ENC(stop.rm = FALSE), # Calculate Effective Number of Codons (including STOP codons)
    x %>% letterFrequency("GC", as.prob = TRUE) * 100 %>% as.vector(), # Calculate % GC content
    x %>% letterFrequency(c("A", "C", "G", "T")), # Nucleotide counts
    x %>% dinucleotideFrequency(), # Dinucleotide counts
    x %>%
      DNAStringSet(start = 1) %>%
      dinucleotideFrequency(step = 3) %>%
      as.data.frame() %>%
      rename_all(., ~ paste0(., "_p1")), # Dinucleotide counts between positions 1-2 only
    x %>%
      DNAStringSet(start = 2) %>%
      dinucleotideFrequency(step = 3) %>%
      as.data.frame() %>%
      rename_all(., ~ paste0(., "_p2")), # Dinucleotide counts between positions 2-3 only
    x %>%
      DNAStringSet(start = 3) %>%
      dinucleotideFrequency(step = 3) %>%
      as.data.frame() %>%
      rename_all(., ~ paste0(., "_p3")), # Dinucleotide counts between positions 3-1 only
    x %>%
      codonTable() %>%
      codonCounts() # Codon counts
  )) %>% rename_at(vars(G.C), ~"GC_content")
  
  if(codonpairs == TRUE){
    df <- cbind(df, x %>% oligonucleotideFrequency(6, step=3))
  }
  
  return(df)
  
}

calc_kmer_counts <- function(x, k, overlap = TRUE, rescale = TRUE){
  
  df <- bind_cols(data.frame(
    gid = x %>% pull(gid),
    segment = x %>% pull(segment),
    x %>% pull(string) %>% DNAStringSet() %>% oligonucleotideFrequency(k, step=ifelse(overlap == TRUE, 1, k))
  )) 
  
  if(rescale == TRUE){
    df %<>%
      mutate(across(-c(1,2))/rowSums(across(-c(1,2))))
  }
  
  return(df)
  
}

calc_composition_bias <- function(df, codonpairs = FALSE){
  
  # Calculate total nucleotides
  df %<>% mutate(length = A+C+G+T)
  
  # Calculate amino acid frequencies
  for (i in 1:length(unique(codon_ref$aminoacid))) {
    df[, paste0(unique(codon_ref$aminoacid)[i], "_aa")] <- 
      df %>% select(subset(codon_ref, aminoacid == unique(codon_ref$aminoacid)[i])$codon) %>% rowSums()
  }
  
  # Calculate total dinucleotides (pos 1-2, pos 2-3, pos 3-1), total codons
  df %<>% mutate(n_dinucs = (select(., matches("^[A|C|G|T][A|C|G|T]$")) %>% rowSums),
                 n_dinucs_p1 = (select(., matches("^[A|C|G|T][A|C|G|T]_p1$")) %>% rowSums),
                 n_dinucs_p2 = (select(., matches("^[A|C|G|T][A|C|G|T]_p2$")) %>% rowSums),
                 n_dinucs_p3 = (select(., matches("^[A|C|G|T][A|C|G|T]_p3$")) %>% rowSums),
                 n_codons = (select(., matches("^[A|C|G|T][A|C|G|T][A|C|G|T]$")) %>% rowSums))
  
  # Calculate nucleotide biases
  df %<>% mutate_at(vars(matches("^[A|C|G|T]$")), .funs = list(Bias = ~./length))
  
  # Calculate dinucleotide biases
  for (i in 1:ncol(df)) {
    
    focalcol <- colnames(df)[i]
    
    if (grepl("^[A|C|G|T][A|C|G|T]$",focalcol)){
      df[, paste0(focalcol, "_Bias")] <- 
        (df[,i]/df$n_dinucs)/
        (df[,substr(focalcol,1,1)]/df$length * df[,substr(focalcol,2,2)] / df$length)
    }
  }
  
  # Calculate dinucleotide biases separately for positions 1-2, 2-3, 3-1
  # But first, need to calculate nucleotide frequencies for these positions - slow but easiest way of doing the calculation
  
  for (nuc in c("A","C","G","T")){
    df[paste0(nuc,"_p1")] <- apply(df %>% select(matches("^[A|C|G|T][A|C|G|T]_p1$")),     # Select only columns describing dinucleotides at position 1-2
                                   1, function(x) x %*% (df %>% select(matches("^[A|C|G|T][A|C|G|T]_p1$")) %>% names %>% str_count(nuc)))   # Use dot product to calculate number of respective nucleotides at position 1-2
    
    df[paste0(nuc,"_p2")] <- apply(df %>% select(matches("^[A|C|G|T][A|C|G|T]_p2$")),     # Select only columns describing dinucleotides at position 2-3
                                   1, function(x) x %*% (df %>% select(matches("^[A|C|G|T][A|C|G|T]_p2$")) %>% names %>% str_count(nuc)))   # Use dot product to calculate number of respective nucleotides at position 2-3
    
    df[paste0(nuc,"_p3")] <- apply(df %>% select(matches("^[A|C|G|T][A|C|G|T]_p3$")),     # Select only columns describing dinucleotides at position 3-1
                                   1, function(x) x %*% (df %>% select(matches("^[A|C|G|T][A|C|G|T]_p3$")) %>% names %>% str_count(nuc)))   # Use dot product to calculate number of respective nucleotides at position 3-1
  }
  
  df$n_nucs_p1 <- df %>% select(matches("^[A|C|G|T]_p1$")) %>% rowSums
  df$n_nucs_p2 <- df %>% select(matches("^[A|C|G|T]_p2$")) %>% rowSums
  df$n_nucs_p3 <- df %>% select(matches("^[A|C|G|T]_p3$")) %>% rowSums
  
  for (i in 1:ncol(df)) {
    
    focalcol <- colnames(df)[i]
    
    if (grepl("^[A|C|G|T][A|C|G|T]_p1$",focalcol)){
      df[, paste0(focalcol, "_Bias")] <- 
        (df[,i]/df$n_dinucs_p1)/
        (df[,paste0(substr(focalcol,1,1),"_p1")]/df$n_nucs_p1 * df[,paste0(substr(focalcol,2,2),"_p1")] / df$n_nucs_p1)
    } else if (grepl("^[A|C|G|T][A|C|G|T]_p2$",focalcol)){
      df[, paste0(focalcol, "_Bias")] <- 
        (df[,i]/df$n_dinucs_p2)/
        (df[,paste0(substr(focalcol,1,1),"_p2")]/df$n_nucs_p2 * df[,paste0(substr(focalcol,2,2),"_p2")] / df$n_nucs_p2)
    } else   if (grepl("^[A|C|G|T][A|C|G|T]_p3$",focalcol)){
      df[, paste0(focalcol, "_Bias")] <- 
        (df[,i]/df$n_dinucs_p3)/
        (df[,paste0(substr(focalcol,1,1),"_p3")]/df$n_nucs_p3 * df[,paste0(substr(focalcol,2,2),"_p3")] / df$n_nucs_p3)
    }
  }
  
  # Calculate Relative Synonymous Codon Usage
  for (i in 1:ncol(df)) {
    
    focalcol <- colnames(df)[i]
    
    if (grepl("^[A|C|G|T][A|C|G|T][A|C|G|T]$",focalcol)){
      df[, paste0(focalcol, "_Bias")] <- 
        df[,i]*subset(codon_ref, codon == focalcol)$deg/
        (df %>%
           select(subset(codon_ref, aminoacid == subset(codon_ref, codon == focalcol)$aminoacid)$codon) %>%
           rowSums())
    }
  }
  
  
  # Calculate amino acid biases - denominator uses total amino acids, including stop codons, so can just use total codons
  for (i in 1:length(unique(codon_ref$aminoacid))) {
    df[, paste0(unique(codon_ref$aminoacid)[i], "_aa_Bias")] <- 
      df[, paste0(unique(codon_ref$aminoacid)[i], "_aa")]/df$n_codons
  }
  
  if(codonpairs == TRUE){
    
    # Calc codon pair bias as log (freq codon pair)/((freq codon A*freq codon B/freq aacid A*freq aacid B) * freq aacid pair) following Coleman et al. 2008
    
    # Calculate codon pair biases
    for (i in 1:nrow(codon_ref)) {
      for (j in 1:nrow(codon_ref)) {
        
        # Calculate frequency of pairs of corresponding amino acids for codons i and j
        aminoacidpaircounts <- df %>%
          select(do.call(paste0,
                         expand.grid(codon_ref %>% subset(aminoacid == codon_ref$aminoacid[i]) %>% .$codon,
                                     codon_ref %>% subset(aminoacid == codon_ref$aminoacid[j]) %>% .$codon))) %>% rowSums()
        
        df[, paste0(codon_ref$codon[i], "_", codon_ref$codon[j], "_Bias")] <-
          log(
            df[, paste0(codon_ref$codon[i], codon_ref$codon[j])]/
              (((df[, codon_ref$codon[i]]*df[, codon_ref$codon[j]])/
                  (df[, paste0(codon_ref$aminoacid[i], "_aa")]*df[, paste0(codon_ref$aminoacid[j], "_aa")]))*
                 aminoacidpaircounts))
        
        # If frequency of codon pair = 0 but frequency of amino acid pair != 0 specifiy codon pair bias as NA as a marker to replace later
        df[which(aminoacidpaircounts != 0 & df[, paste0(codon_ref$codon[i], codon_ref$codon[j])] == 0), paste0(codon_ref$codon[i], "_", codon_ref$codon[j], "_Bias")] <- NA
        
      }
    }
    
    # Below calculations are easily changeable!
    # If frequency of amino acid = 0 replace NaN with 1 (always equivalent to mean RSCU)
    df %<>% mutate_at(vars(matches("^[A|C|G|T][A|C|G|T][A|C|G|T]_Bias$")), ~ifelse(is.nan(.), 1, .))
    
    # Work out column of mean codon pair bias per virus across all non-NaN and non-NA values
    df %<>% mutate(mean_CPB = rowMeans(select(., (matches("^[A|C|G|T][A|C|G|T][A|C|G|T]_[A|C|G|T][A|C|G|T][A|C|G|T]_Bias$"))), na.rm=TRUE))
    
    # Do it excluding pairs involving stop codons
    df %<>% mutate(mean_CPB_nostop = rowMeans(select(., matches("^[A|C|G|T][A|C|G|T][A|C|G|T]_[A|C|G|T][A|C|G|T][A|C|G|T]_Bias$"), -matches("TGA|TAG|TAA")), na.rm=TRUE))
    df %<>% mutate(mean_CPB_nostop1 = rowMeans(select(., matches("^[A|C|G|T][A|C|G|T][A|C|G|T]_[A|C|G|T][A|C|G|T][A|C|G|T]_Bias$"), -matches("^TGA|^TAG|^TAA")), na.rm=TRUE))
    df %<>% mutate(mean_CPB_nostop2 = rowMeans(select(., matches("^[A|C|G|T][A|C|G|T][A|C|G|T]_[A|C|G|T][A|C|G|T][A|C|G|T]_Bias$"), -matches("TGA.Bias|TAG.Bias|TAA.Bias")), na.rm=TRUE))
    
    # Following Babayan et al. rules: if frequency of amino acid pair = 0 replace NaN with mean codon pair bias
    df %<>% mutate_at(vars(matches("^[A|C|G|T][A|C|G|T][A|C|G|T]_[A|C|G|T][A|C|G|T][A|C|G|T]_Bias$")), ~ifelse(is.nan(.), mean_CPB, .))
    
    # Following Babayan et al. rules: if frequency of codon pair = 0 but frequency of amino acid pair != 0 replace NA with -9999 to indicate extreme underrepresentation
    df %<>% mutate_at(vars(matches("^[A|C|G|T][A|C|G|T][A|C|G|T]_[A|C|G|T][A|C|G|T][A|C|G|T]_Bias$")),~ifelse(is.na(.), -9999, .))
  }
  
  return(df)
}

fasta_name_clean <- function(file){
  fasta <- readSet(file = file)
  names(fasta) <- names(fasta) %>% 
    gsub(" ", "_", .) %>% # Replace spaces with underscores, else iFeaturesOmega will not parse FASTA properly
    gsub(">", "", .)      # Remove special characters (that should not be present anyway), else iFeaturesOmega will not parse FASTA properly
  writeXStringSet(fasta, file = file)
  rm(fasta)
}

####################
# Define functions - old and unused, if you use them put them in the above
####################

# Adapted from https://stackoverflow.com/a/27626007
batch <- function(x, n) 
{mapply(function(a, b) (x[a:b]), seq.int(from=1, to=length(x), by=n), pmin(seq.int(from=1, to=length(x), by=n)+(n-1), length(x)), SIMPLIFY=FALSE)}


# Set up accessible colour blindness palette
cbbPalette <- c("#E69F00", "#F0E442", "#56B4E9", "#009E73", "#D55E00")
cbbPalette_full <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cbbPalette_bw <- c("#FFFFFF", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cbbPalette_ordered <- c("#D55E00", "#E69F00", "#F0E442", "#009E73", "#56B4E9", "#0072B2", "#CC79A7", "#999999")
cbbPalette_ordered_bw <- c("#D55E00", "#E69F00", "#F0E442", "#009E73", "#56B4E9", "#0072B2", "#CC79A7", "#b9b9b9")

genomic_pca <- function(df, vars, outcome, choices = 1:2){
  
  df_name <- deparse(substitute(df))
  
  df %<>% as.data.frame
  
  # Merge in relevant outcome column if it doesn't already exist
  if (!(outcome %in% names(df))){
    df %<>% left_join(allcov_df %>% select(childtaxa_id, !! sym(outcome)),
                      by = c("taxid" = "childtaxa_id"))
  }
  
  # Create relevant principal components analysis
  
  if (vars == "dinucs"){
    pca <- df %>% select(matches("^[A|C|G|T][A|C|G|T]_p[1|2|3]_Bias$")) %>% prcomp
  } else if (vars == "codons"){
    pca <- df %>% select(matches("^[A|C|G|T][A|C|G|T][A|C|G|T]_Bias$")) %>% prcomp
  } else if (vars == "codons_nostop"){
    pca <- df %>% select(matches("^[A|C|G|T][A|C|G|T][A|C|G|T]_Bias$")) %>% select(-c(TAG_Bias, TAA_Bias, TGA_Bias)) %>% prcomp
  } else if (vars == "aa"){
    pca <- df %>% select(matches("^.*_aa_Bias$")) %>% prcomp
  } else {
    stop("no valid variables selected")
  }
  
  # Write summary
  sink(paste0("figs\\pcasumm_",df_name,"_",vars,".txt"))
  summary(pca)
  sink()
  
  # Screeplot
  ggscreeplot(pca) +
    geom_hline(yintercept=1/length(pca$sdev), alpha = 0.5, color="dodgerblue", lty="dashed", size=1.5) +
    theme_bw() +
    ggsave(paste0("figs\\scree_",df_name,"_",vars,".png"), width = 8, height = 5)
  
  # Biplot
  g <- ggbiplot(pca,
                choices = choices,
                groups = df[, outcome],
                ellipse = TRUE,
                alpha = 0.4,
                varname.abbrev=TRUE) +
    geom_point(alpha=0, aes(fill= df[, outcome], label=df$childtaxa_name)) +
    theme(legend.position='none') +
    theme_bw()
  ggplotly(g) %>% hide_legend()
  
}

ml_extractor <- function(f, type="full"){
  
  if (f == "AUC"){
    y <- function(x) {x$AUC}
  }
  
  if (f == "acc"){
    y <- function(x) {x$matrix_test$overall["Accuracy"]}
  }
  
  if (f == "TSS"){
    y <- function(x) {as.numeric(x$matrix_test$byClass["Sensitivity"] + x$matrix_test$byClass["Specificity"] - 1)}
  }
  
  if (type == "full"){
    return(data.frame(method = c(rep("RF",nloops),rep("LASSO LR",nloops),rep("GBM ADA",nloops)),
                      metric = c(lapply(rf_list, function(x) x %>% y) %>% unlist,
                                 lapply(lr_list, function(x) x %>% y) %>% unlist,
                                 lapply(gbm_list, function(x) x %>% y) %>% unlist)))
  }
  
  if (type == "mean"){
    return(data.frame(method = c("RF","LASSO LR","GBM ADA"),
                      metric = c(lapply(rf_list, function(x) x %>% y) %>% unlist %>% mean(),
                                 lapply(lr_list, function(x) x %>% y) %>% unlist %>% mean(),
                                 lapply(gbm_list, function(x) x %>% y) %>% unlist %>% mean())))
  }
} 

# Micro and macro F1 scores from https://www.datascienceblog.net/post/machine-learning/performance-measures-multi-class-problems/
get.conf.stats <- function(cm) {
  out <- vector("list", length(cm))
  for (i in seq_along(cm)) {
    x <- cm[[i]]
    tp <- x$table[x$positive, x$positive] 
    fp <- sum(x$table[x$positive, colnames(x$table) != x$positive])
    fn <- sum(x$table[colnames(x$table) != x$positive, x$positive])
    # TNs are not well-defined for one-vs-all approach
    elem <- c(tp = tp, fp = fp, fn = fn)
    out[[i]] <- elem
  }
  df <- do.call(rbind, out)
  rownames(df) <- unlist(lapply(cm, function(x) x$positive))
  return(as.data.frame(df))
}

get.micro.f1 <- function(cm) {
  cm.summary <- get.conf.stats(cm)
  tp <- sum(cm.summary$tp)
  fn <- sum(cm.summary$fn)
  fp <- sum(cm.summary$fp)
  pr <- tp / (tp + fp)
  re <- tp / (tp + fn)
  f1 <- 2 * ((pr * re) / (pr + re))
  return(f1)
}

get.macro.f1 <- function(cm) {
  c <- cm[[1]]$byClass # a single matrix is sufficient
  re <- sum(c[, "Recall"]) / nrow(c)
  pr <- sum(c[, "Precision"]) / nrow(c)
  f1 <- 2 * ((re * pr) / (re + pr))
  return(f1)
}

# Convenience function for setting up data to plot species-labelled prediction figure
rearrange_to_plot_fig_4_spp <- function(host) {
  order <- valid_df_raw %>% 
    filter(host_category == host) %>% 
    group_by(childtaxa_name) %>% 
    summarise(correct = mean(!!sym(host))) %>% 
    arrange(-correct) %>% 
    pull(childtaxa_name)
  valid_df_raw %>% 
    filter(host_category == host) %>% 
    arrange(factor(childtaxa_name, levels = order), -!!sym(host)) %>%
    return
}

# Apply models to predict host in complete start to finish approach from accession number(s)
predict_host_from_accession <- function(accession, rf_list, use_spike) {
  
  Seq_FASTA(accession) %>% write(file = "data\\temp.fasta")
  
  seq <- readDNAStringSet("data\\temp.fasta")
  
  # Construct summary df for individual CDS
  cds_df <- data.frame(
    title = seq %>% names(),
    accessionversion = seq %>% names() %>% str_match("lcl\\|(.*?)_cds_") %>% .[, 2] %>% as.character(),
    gene = seq %>% names() %>% str_match("\\[gene=(.*?)\\]") %>% .[, 2],
    protein = seq %>% names() %>% str_match("\\[protein=(.*?)\\]") %>% .[, 2],
    length = seq %>% width()
  )
  
  if(use_spike == FALSE & frameshift_correct == TRUE) {
    
    # Extract and append CDS location information by text extraction from title
    cds_df %<>% cbind(cds_df %>%
                        pull(title) %>%
                        str_match("\\[location=(.*?)\\]") %>%
                        .[, 2] %>%
                        gsub("join\\(|\\)|>|<", "", .) %>%
                        gsub(",", "..", .) %>%
                        str_split(., "\\..") %>%
                        lapply(., function(x) {
                          x %>%
                            t() %>%
                            as.data.frame()
                        }) %>%
                        bind_rows() %>%
                        mutate_all(., as.numeric) %>%
                        rename_all(., ~ gsub("V", "loc_", .)))
    
    # Correct frameshift/slippage causing double accounting of elements of some ORF1ab CDS by specifying new location to start calculating metrics at
    cds_df %<>% left_join(
      cds_df %>% filter( # Only apply correction where needed
        accessionversion %in% (cds_df %>% filter(grepl("join", title)) %>% pull(accessionversion))) %>% # Filter to sequences containing a CDS that join overlapping components
        group_by(accessionversion, loc_1) %>% arrange(accessionversion, loc_1) %>% filter(n() > 1) %>% # Filter to sequences with two different CDS that both start at the same location
        mutate(start_new = ifelse(is.na(loc_3), lag(loc_2) - (loc_1 - 1) + 1, 1)) %>% # For the shorter sequence of ORF1a, set new start location based on slippage location of ORF1b
        ungroup() %>%
        select(title, start_new),
      by = "title"
    ) %>%
      replace_na(list(start_new = 1))
    
    cds_df %<>% mutate(nested_flag = ifelse(start_new >= length, "nested", NA))
    
    if ("nested" %in% cds_df$nested){
      
      # If ORF1a nested entirely within ORF1b, then just drop ORF1a
      seq %<>% .[-which(cds_df$nested_flag == "nested")]
      cds_df %<>% filter(is.na(nested_flag))
      
    }
    
    # Apply new CDS start locations
    seq %<>% subseq(start = cds_df$start_new)
    
  }
  
  # Assign each cds as whole S, S1, S2 or other
  cds_df %<>% mutate(seqtype = case_when(
    grepl("^s1| s1|^s-1| s-1|subunit 1$|spike 1|spike protein 1|spike glycoprotein 1", protein, ignore.case = TRUE) ~ "S1",
    grepl("^s2| s2|^s-2| s-2|subunit 2$|spike 2|spike protein 2|spike glycoprotein 2", protein, ignore.case = TRUE) ~ "S2",
    grepl("spike|^surface|surface gly|s gly|s prot|peplom|protein S$|^S$", protein, ignore.case = TRUE) ~ "S"
  )) %>% replace_na(list(seqtype = "other"))
  
  # Override with gene information where available
  cds_df %<>% mutate(seqtype = case_when(
    grepl("^s1| s1|^s-1| s-1|subunit 1$|spike 1|spike protein 1|spike glycoprotein 1", gene, ignore.case = TRUE) ~ "S1",
    grepl("^s2| s2|^s-2| s-2|subunit 2$|spike 2|spike protein 2|spike glycoprotein 2", gene, ignore.case = TRUE) ~ "S2",
    TRUE ~ seqtype
  ))
  
  # Count nucs, dinucs, codons
  df <- data.frame(accessionversion = cds_df$accessionversion,
                   seqtype = cds_df$seqtype,
                   enc = seq %>% codonTable() %>% ENC(stop.rm = FALSE), # Calculate Effective Number of Codons (including STOP codons)
                   seq %>% letterFrequency("GC", as.prob = TRUE) * 100 %>% as.vector(), # Calculate % GC content
                   seq %>% letterFrequency(c("A", "C", "G", "T")), # Nucleotide counts
                   seq %>% dinucleotideFrequency(), # Dinucleotide counts
                   seq %>%
                     DNAStringSet(start = 1) %>%
                     dinucleotideFrequency(step = 3) %>%
                     as.data.frame() %>%
                     rename_all(., ~ paste0(., "_p1")), # Dinucleotide counts between positions 1-2 only
                   seq %>%
                     DNAStringSet(start = 2) %>%
                     dinucleotideFrequency(step = 3) %>%
                     as.data.frame() %>%
                     rename_all(., ~ paste0(., "_p2")), # Dinucleotide counts between positions 2-3 only
                   seq %>%
                     DNAStringSet(start = 3) %>%
                     dinucleotideFrequency(step = 3) %>%
                     as.data.frame() %>%
                     rename_all(., ~ paste0(., "_p3")), # Dinucleotide counts between positions 3-1 only
                   seq %>%
                     codonTable() %>%
                     codonCounts()) # Codon counts
  
  # If doing spike based predictions, select spike sequence, otherwise sum counts over whole genome
  if (use_spike == TRUE){
    df %<>% filter(seqtype == "S")
  } else{
    df %<>%
      group_by(accessionversion) %>%
      summarise_if(is.numeric, sum) %>%
      mutate(GC_content = 100 * (G + C) / (A + C + G + T), CDS_length = A + G + C + T)
  }
  
  # Calculate composition biases
  df %<>% calc_composition_bias
  
  # Strip out any remaining sequences that have missing values that have slipped through (usually truncated sequences with no stop codon)
  if (df %>% filter(!complete.cases(.)) %>% nrow > 0){
    df %>% filter(!complete.cases(.)) %>% pull(accessionversion) %>% paste(collapse = ", ") %>% paste0("Removed accessions due to missing data: ", .) %>% print
    df %<>% filter(complete.cases(.))
  }
  
  # Produce rf predictions
  output <- lapply(rf_list, function(model)
    
    data.frame(accessionversion = df$accessionversion, predict(model, newdata = df, type="prob"))
    
  ) %>% 
    bind_rows %>%
    group_by(accessionversion) %>%
    summarise_if(is.numeric, mean)
  
  return(as.data.frame(output))
}
