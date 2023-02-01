library(tidyverse)
library(vroom)
library(magrittr)
virion <- vroom("C:\\Users\\Liam\\Documents\\GitHub\\virion\\Virion\\Virion.csv.gz")

virion %>% filter(Virus == "influenza a virus" & HostClass == "aves") %>% filter(!is.na(NCBIAccession)) %>% pull(NCBIAccession) %>% unique %>% length

 # 22,433 sequences

virion %>% filter(Database != "GenBank" & Virus == "influenza a virus" & HostClass == "aves") %>% filter(!is.na(NCBIAccession)) %>% nrow

# 3679 not in GenBank as Database ID, all from EID2, which may be listing proteins separately compared to main GenBank entries above