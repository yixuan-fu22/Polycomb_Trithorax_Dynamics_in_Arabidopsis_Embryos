#
# This script is an example of diffBind differential analysis for H3K4me3 between two stages
# and getting count matrix from it

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

# specify the histone modification to analyze
histoneModification <- "H3K4ME3" # use H3K4me3 first, and process the another two histone modifications with goto functions

# specific comparison and make folder
stage1 <- "WS"
stage2 <- "MG" # choose from three combinations: WS vs MG, WS vs SDL, and MG vs SDL 
              # Avoid inputing all stages at once, as seedling samples may screw normalization
tag <- paste0(stage1, "_vs_", stage2)
resDir <- paste0("../results/diffBind_H3K4me3_", tag)
if(!file.exists(resDir)){dir.create(resDir)}
qcDir <- file.path(resDir, "qcPlots")
if(!file.exists(qcDir)){dir.create(qcDir)} # make a directory to store plots for quality control that are not intended to be used for publication, incl. the defaudt ones in the package

# read samples and select tissues
samples <- read.csv('./sampleSheet_diffBind.csv', sep = ",")

samples_filtered <- samples %>% filter(Factor == histoneModification) %>% 
  filter(Tissue %in% c(stage1, stage2)) %>%
  filter(Condition %in% c("wt"))

#### Differential Analysis by DiffBind (DEseq2) ####

# make objects
dbObj <- dba(sampleSheet=samples_filtered, scoreCol=5)
dbObj$config$cores <- 8

# count reads
dbObj <- dba.count(dbObj, bUseSummarizeOverlaps=TRUE, 
                   bParallel = TRUE,
                   bRemoveDuplicates = TRUE)
dbObj

# Plot PCA and Coorelation for QC
png(filename = file.path(qcDir, paste0(tag, "_PCA_general.png")), height = 1000, width = 1000)
dba.plotPCA(dbObj)
dev.off()
png(filename = file.path(qcDir, paste0(tag, "_COORELATION.png")), height = 1500, width = 1500)
plot(dbObj)
dev.off()

# set contrast
dbObj <- dba.contrast(dbObj, categories=DBA_TISSUE, minMembers = 2)

# analyze
dbObj <- dba.analyze(dbObj, method=DBA_DESEQ2 , bBlacklist = FALSE, bGreylist = FALSE, bParallel=TRUE)
dba.show(dbObj, bContrasts=T)

# after-comparison plots
# PCA plots
png(filename = file.path(qcDir, paste0(tag, "_PCA_DESEQ2.png")), height = 1000, width = 1000)
dba.plotPCA(dbObj, contrast=1, method=DBA_DESEQ2)
dev.off()
# MA plots
png(filename = file.path(qcDir, paste0(tag, "_MA_DESEQ2.png")), height = 1000, width = 1000)
dba.plotMA(dbObj, method=DBA_DESEQ2)
dev.off()
# XY scatter plot 
png(filename = file.path(qcDir, paste0(tag, "_XYScatter_ALLMETHOD.png")), height = 1000, width = 1000)
dba.plotMA(dbObj, bXY=TRUE)
dev.off()
# Volcano plots
png(filename = file.path(qcDir, paste0(tag, "VOLCANO_DESEQ2.png")), height = 1000, width = 1000)
dba.plotVolcano(dbObj, method = DBA_DESEQ2)
dev.off()
# Box plots
png(filename = file.path(qcDir, paste0(tag, "BOX_DESEQ2.png")), height = 1000, width = 1000)
dba.plotBox(dbObj, method = DBA_DESEQ2)
dev.off()

# set a threshold for log2 FC
LFCthres <- 0.5

# save results
res_deseq <- dba.report(dbObj, contrast = 1, th = 1)
out <- as.data.frame(res_deseq)
summary_table <- dba.show(dbObj, bContrasts=T)
write.table(x = summary_table, 
            file = file.path(resDir, "diffBind_summary.tsv"), sep = "\t", quote = FALSE)
res_deseq <- dba.report(dbObj, contrast = 1, th = 1)
out <- as.data.frame(res_deseq)
write.table(x = out,
              file = file.path(resDir, paste0("res_deseq_", tag, ".tsv")),
              sep = "\t",
              quote = FALSE)

## Extract count matrix
resDir_countMatrix <- c(file.path(resDir, "countMatrix"))
if(!file.exists(resDir_countMatrix)){dir.create(resDir_countMatrix)}
dbObj <- dba(sampleSheet=samples_filtered, scoreCol = 5)
dbObj$config$cores <- 8 
# DESeq2 normalized
dbObj_DESEQ2normalized <- dba.count(dbObj, bUseSummarizeOverlaps=TRUE, bParallel = TRUE,
                                    bRemoveDuplicates = TRUE, score = DBA_SCORE_NORMALIZED)
countmatrix_DESEQ2normalized <- dba.peakset(dbObj_DESEQ2normalized, bRetrieve = TRUE, DataType = DBA_DATA_FRAME)
write.table(x = countmatrix_DESEQ2normalized,
            file = file.path(resDir_countMatrix, paste0("countmatrix_DESEQ2normalized.tsv")),
            sep = "\t",
            quote = FALSE)
# RPKM
dbObj_RPKM <- dba.count(dbObj, bUseSummarizeOverlaps=TRUE, bParallel = TRUE,
                        bRemoveDuplicates = TRUE, score = DBA_SCORE_RPKM)
countmatrix_RPKM <- dba.peakset(dbObj_RPKM, bRetrieve = TRUE, DataType = DBA_DATA_FRAME)
write.table(x = countmatrix_RPKM,
            file = file.path(resDir_countMatrix, paste0("countmatrix_RPKM.tsv")),
            sep = "\t",
            quote = FALSE)



































