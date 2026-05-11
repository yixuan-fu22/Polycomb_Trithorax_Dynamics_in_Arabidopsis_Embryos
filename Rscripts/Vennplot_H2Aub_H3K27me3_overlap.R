# Overlap of H2Aub marked genes in embryos and these with a inscreased H3K27me3 in seedlings

#renv
library(renv)
#renv::init()
renv::activate()
renv::snapshot()

suppressMessages({
  library(stringr)
  library(dplyr)
  library(tidyverse)
  library(tidyr)
  library(ggplot2)
  #
  library(GenomicRanges)
  library(ChIPseeker)
  #library(ggupset)
  # Venn Plot
  library(VennDiagram)
})















