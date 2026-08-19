# Load libraries
library(tidyverse)


list.files("H:\\Working\\ai_for_ai\\S3/data/full/mapping/binary_result/prot",
           recursive = TRUE,
           full.names = TRUE
) %>% 
  gsub(".*/binary_", "", .) %>%
  unique() %>%
  gsub(".csv", "", .) %>%
  write.csv("H:\\Working\\ai_for_ai\\S3/data/full/mapping/prot_feats_ref.csv")


list.files("H:\\Working\\ai_for_ai\\S3/data/full/mapping/binary_result/cds",
           recursive = TRUE,
           full.names = TRUE
) %>% 
  gsub(".*/binary_", "", .) %>%
  unique() %>%
  gsub(".csv", "", .) %>%
  write.csv("H:\\Working\\ai_for_ai\\S3/data/full/mapping/cds_feats_ref.csv")

list.files("H:\\Working\\ai_for_ai\\S3/data/full/mapping/binary_result/nuc",
           recursive = TRUE,
           full.names = TRUE
) %>% 
  gsub(".*/binary_", "", .) %>%
  unique() %>%
  gsub(".csv", "", .) %>%
  write.csv("H:\\Working\\ai_for_ai\\S3/data/full/mapping/nuc_feats_ref.csv")

# OR plots

list <- list()
proteins <- c("HA", "M1", "NA", "NP", "NS1", "PA", "PB1", "PB2")
for(i in 1:8){
  feat_ORs <- list.files("H:\\Working\\ai_for_ai\\S3/data/full/mapping/kmer_or_results/prot",
                         pattern = "odds_ratio_by_position_with_p.csv",
                         recursive = TRUE,
                         full.names = TRUE
  ) %>%
    .[grepl(proteins[i], .)]
  # Clean and transform OR value
  bypos <- feat_ORs %>% purrr::map_dfr(read.csv)
  # Create shell for every feature-position combination
  shell <- expand.grid(kmer = read.csv("H:\\Working\\ai_for_ai\\S3/data/full/mapping/prot_feats_ref.csv") %>% pull(x), # features
                       position = 1:max(bypos$position)) # positions
  list[[i]] <- bypos %>% right_join(shell, by = c("kmer", "position")) %>%
    select(kmer, position, odds_ratio) %>%
    mutate(
      position = as.integer(position),
      odds_ratio = suppressWarnings(as.numeric(odds_ratio)),
      odds_ratio = if_else(is.na(odds_ratio) | !is.finite(odds_ratio) | odds_ratio <= 0,
                           NA_real_, odds_ratio
      )
    ) %>%
    complete(kmer, position, fill = list(odds_ratio = 1)) %>% # ensure positions where no odds ratio difference still contribute to weighted average
    mutate(abs_logOR = abs(log(odds_ratio))) %>%
    filter(!is.na(position)) %>%
    mutate(focgene = proteins[i])
}

all_ORs <- bind_rows(list)

all_ORs %>%
  filter(abs_logOR != 0) %>%
  ggplot(aes(x = abs_logOR , color = focgene)) +
  geom_density(lwd=2)

all_ORs %>%
  group_by(position, focgene) %>%
  summarise(
    mean_abs_logOR = mean(abs_logOR, na.rm = TRUE),
  ) %>%
  ggplot(aes(x = mean_abs_logOR , color = focgene)) +
  geom_density(lwd=2)



list <- list()
proteins <- c("HA", "M1", "NA", "NP", "NS1", "PA", "PB1", "PB2")
for(i in 1:8){
  feat_ORs <- list.files("H:\\Working\\ai_for_ai\\S3/data/full/mapping/kmer_or_results/nuc",
                         pattern = "odds_ratio_by_position_with_p.csv",
                         recursive = TRUE,
                         full.names = TRUE
  ) %>%
    .[grepl(proteins[i], .)]
  # Clean and transform OR value
  bypos <- feat_ORs %>% purrr::map_dfr(read.csv)
  # Create shell for every feature-position combination
  shell <- expand.grid(kmer = read.csv("H:\\Working\\ai_for_ai\\S3/data/full/mapping/nuc_feats_ref.csv") %>% pull(x), # features
                       position = 1:max(bypos$position)) # positions
  list[[i]] <- bypos %>% right_join(shell, by = c("kmer", "position")) %>%
    select(kmer, position, odds_ratio) %>%
    mutate(
      position = as.integer(position),
      odds_ratio = suppressWarnings(as.numeric(odds_ratio)),
      odds_ratio = if_else(is.na(odds_ratio) | !is.finite(odds_ratio) | odds_ratio <= 0,
                           NA_real_, odds_ratio
      )
    ) %>%
    complete(kmer, position, fill = list(odds_ratio = 1)) %>% # ensure positions where no odds ratio difference still contribute to weighted average
    mutate(abs_logOR = abs(log(odds_ratio))) %>%
    filter(!is.na(position)) %>%
    mutate(focgene = proteins[i])
}

all_ORs <- bind_rows(list)

all_ORs %>%
  filter(abs_logOR != 0) %>%
  ggplot(aes(x = abs_logOR , color = focgene)) +
  geom_density(lwd=2)

all_ORs %>%
  group_by(position, focgene) %>%
  summarise(
    mean_abs_logOR = mean(abs_logOR, na.rm = TRUE),
  ) %>%
  ggplot(aes(x = mean_abs_logOR , color = focgene)) +
  geom_density(lwd=2)

