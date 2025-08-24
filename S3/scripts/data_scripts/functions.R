####################
# Define functions #
####################

# Create data from specific NCBI Flu fasta header and clean metadata
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

# Create data from specific GISAID fasta header and clean metadata
process_GISAID_seq <- function(x, label, type){
  df <- data.frame(title = x %>% names(), 
                   length = x %>% width()) 
  
  if (type == "nuc"){
    df %<>% separate_wider_regex(title, patterns = c(title = ".*?", "\\|",
                                                     UID = ".*?", "\\|",
                                                     subtype = ".*?", "\\|",
                                                     null = ".*?", "\\|",
                                                     date = ".*?", "\\|",
                                                     INSDC = ".*?", "\\|",
                                                     accession = ".*?", "\\|",
                                                     title2 = ".*", "\\|",
                                                     gene = ".*?", "\\|",
                                                     segment = ".*?"))
  } else if (type == "prot"){
    df %<>% tidyr::separate_wider_delim(title, delim = "|", names = c("title", "UID", "subtype", "null", "date", "protINSDC", "protaccession", "gene"), too_many = "drop")  
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

# Create data from specific NCBI Virus fasta header and clean metadata
process_NV_seq <- function(x, label, type){
  df <- data.frame(title = x %>% names(), 
                   length = x %>% width()) 
  
  if (type == "nuc"){
    df %<>% tidyr::separate(title, sep = "\\|", into = c("title", "subtype", "date", "accession", "segment"), extra = "drop")
  } else if (type == "cds"){
    df %<>% tidyr::separate(title, sep = "\\|", into = c("title", "subtype", "date", "accession", "segment"), extra = "drop")
  } else if (type == "prot"){
    df %<>% tidyr::separate(title, sep = "\\|", into = c("title", "subtype", "date", "accession", "segment"), extra = "drop")  
  } else {
    stop("invalid type (must be 'nuc', 'cds', or 'prot')")
  }
  
  df %<>%
    mutate(
      accession = gsub("\\:.*","",accession),
      gene = case_when(
        segment == 1 ~ "PB2",                  # Assign based on segment field as should be pre-filtered on NCBI Virus to genes of interest
        segment == 2 ~ "PB1",
        segment == 3 ~ "PA",
        segment == 4 ~ "HA",
        segment == 5 ~ "NP",
        segment == 6 ~ "NA",
        segment == 7 ~ "M1",
        segment == 8 ~ "NS1"),
      label = label,
      src = "NCBI",
      title = gsub("Influenza A virus \\(", "", title),
      title = gsub(") .*", "", title),
      title = gsub("_.*|\\/[1-2][0-9][0-9][0-9].*", "", title),
      gid = title,
      fastahead = x %>% names() %>% gsub(" ", "_", .),
      string = x %>% as.character(use.names=FALSE),
    ) %>%
    relocate(fastahead, string, .after = last_col())
  
  return(df)
  
}

# Bind ORFs into a data frame, and return the longest ORF
bind_ORF <- function(x) {
  
  df <- findORFs(x, startCodon="ATG") %>%
    as.data.frame() 
  
  if(nrow(df) == 0){
    
    df <- data.frame(start = 0, end = 0) # Set start/end to zero if no valid open reading frames found
    
  } else {
    
    df <- df %>% 
      arrange(-width) %>%
      slice(1) %>%  # Select and keep only the longest ORF
      select(start, end)
  }   
  
  return(df)
}

# Produce data frame of nucleotide, dinucleotide, codon, and optionally, codon pair frequency
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

# Produce data frame of k-mer frequency
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

# Produce nucleotide, dinucleotide, codon, and optionally, codon pair frequency into composition bias features (e.g., adjusting for expectation based on underlying nucleotide composition)
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

# Replace incompatible characters in nucleotide fasta headers for iFeaturesOmega
nuc_fasta_name_clean <- function(filename){
  fasta <- readSet(file = filename)
  names(fasta) <- names(fasta) %>% 
    gsub(" ", "_", .) %>% # Replace spaces with underscores, else iFeaturesOmega will not parse FASTA properly
    gsub(">", "", .)      # Remove special characters (that should not be present anyway), else iFeaturesOmega will not parse FASTA properly
  writeXStringSet(fasta, file = filename)
  rm(fasta)
}

# Replace incompatible characters in protein fasta headers for iFeaturesOmega
prot_fasta_name_clean <- function(filename){
  fasta <- readAAStringSet(file = filename)
  names(fasta) <- names(fasta) %>% 
    gsub(" ", "_", .) %>% # Replace spaces with underscores, else iFeaturesOmega will not parse FASTA properly
    gsub(">", "", .)      # Remove special characters (that should not be present anyway), else iFeaturesOmega will not parse FASTA properly
  writeXStringSet(fasta, file = filename)
  rm(fasta)
}