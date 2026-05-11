#renv
library(renv)
# renv::init(bioconductor = "3.22")
renv::activate()
# renv::restore()

suppressMessages({
  # utilities
  library(stringr)
  library(dplyr)
  library(tidyr)
  library(tidyverse)
  
  # differential analysis, annnotation, and enrichment
  library(ChIPseeker)
  library(DiffBind)
  library(AnnotationDbi)
  library(TxDb.Athaliana.BioMart.plantsmart25)
  library(org.At.tair.db)
  library(clusterProfiler)
  
  # plots
  library(ggplot2)
  library(ggpubr)
  # library(ggupset)
})

install.packages("stringr")
install.packages("dplyr")
install.packages("tidyr")
install.packages("tidyverse")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ChIPseeker")
BiocManager::install("DiffBind")
BiocManager::install("AnnotationDbi")
BiocManager::install("TxDb.Athaliana.BioMart.plantsmart25")
BiocManager::install("org.At.tair.db")
BiocManager::install("clusterProfiler")

install.packages("ggplot2")
install.packages("ggpubr")

renv::snapshot()



library(ChIPpeakAnno)
#
# motif analysis
library(BSgenome.Athaliana.TAIR.TAIR9)
# library(TFBSTools)

# library(motifStack)

BiocManager::install("ChIPpeakAnno")
# Error: package 'TFMPvalue' is not available
BiocManager::install("BSgenome.Athaliana.TAIR.TAIR9")
# BiocManager::install("motifStack")
# BiocManager::install("TFBSTools")

BiocManager::install("TxDb.Athaliana.BioMart.plantsmart25")

install.packages("ggbeeswarm")
install.packages("ggpubr")


install.packages("corrplot")








