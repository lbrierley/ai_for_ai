
install.packages("BiocManager", repos = "https://cloud.r-project.org")

bioconductor_pkgs <- c("Biostrings", "coRdon", "ORFik")

BiocManager::install(bioconductor_pkgs)

cran_pkgs <- c(
  "ggmosaic",
  "gplots",
  "kableExtra",
  "knitr",
  "lattice",
  "magrittr",
  "patchwork",
  "pbapply",
  "plotly",
  "rentrez",
  "reshape2",
  "rmarkdown",
  "R.utils",
  "seqinr",
  "stringr",
  "tidyverse",
  "VennDiagram",
  "XML"
)

install.packages(cran_pkgs)

remotes::install_github("ropensci/taxizedb")
