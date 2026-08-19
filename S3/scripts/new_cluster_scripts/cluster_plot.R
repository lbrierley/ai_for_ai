library(tidyverse)
library(ggplot2)
library(ggsankey)
library(ggalluvial)
library(infotheo)
library(aricode)

# Read in reference of all sequences used
allflu_wgs_ref <- read.csv("S3\\data\\segmentwise\\Raw_Sequences\\allflu_wgs_df.csv", na.strings = "")

# Read in cluster codes
clust_codes <- read.csv("S3\\data\\segmentwise\\cds_wgs_clust_codes.csv") %>%
  left_join(allflu_wgs_ref %>% select(gid, subtype, label))

clust_codes %>% group_by(code) %>% tally() %>% with(., table(n))
clust_codes %>% group_by(code) %>% tally() %>% with(., table(n)) %>% barplot
clust_codes %>% group_by(code) %>% tally() %>% arrange(-n)

# Calculate mutual information between each pair of segments' cluster codes

mutinf <- clust_codes %>% select(X1:X8) %>% NMI()
mutinf[upper.tri(mutinf, diag = FALSE)] <- NA

mutinf %>% write.csv("S3\\data\\segmentwise\\cds_wgs_clust_code_mi.csv")

nmi <- expand.grid(C1 = names(select(clust_codes, X1:X8)),
                   C2 = names(select(clust_codes, X1:X8)),
                   stringsAsFactors = FALSE) %>%
  rowwise() %>%
  mutate(nmi = NMI(clust_codes[[C1]], clust_codes[[C2]])) %>%
  ungroup() %>%
  pivot_wider(names_from = C2, values_from = nmi) %>%
  column_to_rownames("C1")

nmi[upper.tri(nmi, diag = FALSE)] <- NA
nmi %>% write.csv("S3\\data\\segmentwise\\cds_wgs_clust_code_nmi.csv")


# Sankey version - displays in order consistently between all segments but can't colour by subtype

S_all <- clust_codes %>% 
  filter(label == "zoon") %>%
  make_long(X4, X6, X1, X2, X3, X5, X7, X8) %>%
  ggplot(aes(x = x, 
             next_x = next_x, 
             node = node, 
             next_node = next_node,
             fill = factor(node),
             label = node),) +
  geom_sankey(na.rm = TRUE, flow.alpha = 0.4, node.color = 1, color = "black") +
  geom_sankey_text(aes(label = after_stat(paste0(node))), size = 3.5, color = "white") +
  scale_x_discrete(limits = c("X4", "X6", "X1", "X2", "X3", "X5", "X7", "X8"), 
                   labels = c("HA", "NA", "PB2", "PB1", "PA",  "NP", "M1", "NS1"),
                   expand = c(.2, .05)) +
  theme_sankey() +
  theme(legend.position = "none",
        axis.title.x = element_blank())

ggsave("S3\\figures_tables\\sankey_segmentwise.png", plot = S_all, width = 30, height = 15)
  
# Alluvial version version - colour by subtype but can't display in order consistently between all segments
# Not sure this is correct, some are moving from 1 to 1 but not sticking together

A_all <- clust_codes %>% 
  filter(label == "zoon") %>%
  group_by(X1, X2, X3, X4, X5, X6, X7, X8, subtype) %>% 
  tally %>%
  ungroup %>%
  mutate(across(X1:X8, as.factor)) %>%
  ggplot(aes(axis1 = X4, axis2 = X6, axis3 = X1, axis4 = X2, axis5 = X3, axis6 = X5, axis7 = X7, axis8 = X8,
             y = n)) +
  scale_x_discrete(limits = c("X4", "X6", "X1", "X2", "X3", "X5", "X7", "X8"), 
                   labels = c("HA", "NA", "PB2", "PB1", "PA",  "NP", "M1", "NS1"),
                   expand = c(.2, .05)) +
  geom_alluvium(aes(fill = subtype)) +
  geom_stratum() +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  theme_bw()


ggsave("S3\\figures_tables\\alluvial_all_segmentwise.png", plot = A_all, width = 20, height = 10)


A_HA_NS1 <- clust_codes %>% 
  filter(label == "zoon") %>%
  group_by(X4, X8, subtype) %>% 
  tally %>%
  ungroup %>%
  mutate(across(X4:X8, as.factor)) %>%
  ggplot(aes(axis1 = X4, axis2 = X8,
             y = n)) +
  scale_x_discrete(limits = c("X4", "X8"), 
                   labels = c("HA", "NS1"),
                   expand = c(.2, .05)) +
  geom_alluvium(aes(fill = subtype)) +
  geom_stratum() +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  theme_bw()


ggsave("S3\\figures_tables\\alluvial_HA_NS1_segmentwise.png", plot = A_HA_NS1, width = 4, height = 10)

A_HA_M1 <- clust_codes %>% 
  filter(label == "zoon") %>%
  group_by(X4, X7, subtype) %>% 
  tally %>%
  ungroup %>%
  mutate(across(X4:X7, as.factor)) %>%
  ggplot(aes(axis1 = X4, axis2 = X7,
             y = n)) +
  scale_x_discrete(limits = c("X4", "X7"), 
                   labels = c("HA", "M1"),
                   expand = c(.2, .05)) +
  geom_alluvium(aes(fill = subtype)) +
  geom_stratum() +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  theme_bw()


ggsave("S3\\figures_tables\\alluvial_HA_M1_segmentwise.png", plot = A_HA_M1, width = 4, height = 10)

A_HA_PB2 <- clust_codes %>% 
  filter(label == "zoon") %>%
  group_by(X4, X1, subtype) %>% 
  tally %>%
  ungroup %>%
  mutate(across(X4:X1, as.factor)) %>%
  ggplot(aes(axis1 = X4, axis2 = X1,
             y = n)) +
  scale_x_discrete(limits = c("X4", "X1"), 
                   labels = c("HA", "PB2"),
                   expand = c(.2, .05)) +
  geom_alluvium(aes(fill = subtype)) +
  geom_stratum() +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  theme_bw()


ggsave("S3\\figures_tables\\alluvial_HA_PB2_segmentwise.png", plot = A_HA_PB2, width = 4, height = 10)

A_HA_PB1 <- clust_codes %>% 
  filter(label == "zoon") %>%
  group_by(X4, X2, subtype) %>% 
  tally %>%
  ungroup %>%
  mutate(across(X4:X2, as.factor)) %>%
  ggplot(aes(axis1 = X4, axis2 = X2,
             y = n)) +
  scale_x_discrete(limits = c("X4", "X2"), 
                   labels = c("HA", "PB1"),
                   expand = c(.2, .05)) +
  geom_alluvium(aes(fill = subtype)) +
  geom_stratum() +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  theme_bw()


ggsave("S3\\figures_tables\\alluvial_HA_PB1_segmentwise.png", plot = A_HA_PB1, width = 4, height = 10)

A_HA_PA <- clust_codes %>% 
  filter(label == "zoon") %>%
  group_by(X4, X3, subtype) %>% 
  tally %>%
  ungroup %>%
  mutate(across(X4:X3, as.factor)) %>%
  ggplot(aes(axis1 = X4, axis2 = X3,
             y = n)) +
  scale_x_discrete(limits = c("X4", "X3"), 
                   labels = c("HA", "PA"),
                   expand = c(.2, .05)) +
  geom_alluvium(aes(fill = subtype)) +
  geom_stratum() +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  theme_bw()


ggsave("S3\\figures_tables\\alluvial_HA_PA_segmentwise.png", plot = A_HA_PA, width = 4, height = 10)

A_HA_NP <- clust_codes %>% 
  filter(label == "zoon") %>%
  group_by(X4, X5, subtype) %>% 
  tally %>%
  ungroup %>%
  mutate(across(X4:X5, as.factor)) %>%
  ggplot(aes(axis1 = X4, axis2 = X5,
             y = n)) +
  scale_x_discrete(limits = c("X4", "X5"), 
                   labels = c("HA", "NP"),
                   expand = c(.2, .05)) +
  geom_alluvium(aes(fill = subtype)) +
  geom_stratum() +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  theme_bw()


ggsave("S3\\figures_tables\\alluvial_HA_NP_segmentwise.png", plot = A_HA_NP, width = 4, height = 10)

A_HA_NA <- clust_codes %>% 
  filter(label == "zoon") %>%
  group_by(X4, X6, subtype) %>% 
  tally %>%
  ungroup %>%
  mutate(across(X4:X6, as.factor)) %>%
  ggplot(aes(axis1 = X4, axis2 = X6,
             y = n)) +
  scale_x_discrete(limits = c("X4", "X6"), 
                   labels = c("HA", "NA"),
                   expand = c(.2, .05)) +
  geom_alluvium(aes(fill = subtype)) +
  geom_stratum() +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  theme_bw()


ggsave("S3\\figures_tables\\alluvial_HA_NA_segmentwise.png", plot = A_HA_NA, width = 4, height = 10)

