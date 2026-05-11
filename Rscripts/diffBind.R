# This script does a differential analysis on H3K4me3 or H2AUb in MG and WS embryos
# The results can be piped into ChIPseeker for annotation and then enrichment analysis

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
  library(DiffBind) # 3.12.0
  library(GenomicRanges)
  #library(edgeR)
})

histoneModification <- "H2AUB"
stage1 <- "WS"
stage2 <- "SDL"
resDir <- paste0("../results/diffBind_ChIPseeker_", stage1, "vs", stage2, "_" ,
                 histoneModification, "_newDataset")
if(!file.exists(resDir)){dir.create(resDir)}

# read samples and select tissues
samples <- read.csv('./sampleSheet_diffBind.csv', sep = ",")

samples_filtered <- samples %>% filter(Factor == histoneModification) %>% 
  filter(Tissue %in% c(stage1, stage2)) %>%
  filter(Condition %in% c("wt"))


####

dbObj <- dba(sampleSheet=samples_filtered, scoreCol=5)

dbObj$config$cores <- 8

# count reads
dbObj <- dba.count(dbObj, bUseSummarizeOverlaps=TRUE, 
                   bParallel = TRUE,
                   bRemoveDuplicates = TRUE)
dbObj

# plot PCA
png(filename = file.path(resDir,"PCA_general.png"), height = 1000, width = 1000)
dba.plotPCA(dbObj)
dev.off()

# Plot correlation
png(filename = file.path(resDir,"COORELATION.png"), height = 1500, width = 1500)
plot(dbObj)
dev.off()

# set contrast
dbObj <- dba.contrast(dbObj, categories=DBA_TISSUE, minMembers = 2)

# analyze
dbObj <- dba.analyze(dbObj, method=DBA_DESEQ2 , bBlacklist = FALSE, bGreylist = FALSE, bParallel=TRUE)
dba.show(dbObj, bContrasts=T) # Unable to normalize datset with edge

# PCA plots
png(filename = file.path(resDir,"PCA_DESEQ2.png"), height = 1000, width = 1000)
dba.plotPCA(dbObj, contrast=1, method=DBA_DESEQ2)
dev.off()

png(filename = file.path(resDir,"PCA_EDGER.png"), height = 1000, width = 1000)
dba.plotPCA(dbObj, contrast=1, method=DBA_EDGER)
dev.off()

# MA plots
png(filename = file.path(resDir,"MA_DESEQ2.png"), height = 1000, width = 1000)
dba.plotMA(dbObj, method=DBA_DESEQ2)
dev.off()

png(filename = file.path(resDir,"MA_EDGER.png"), height = 1000, width = 1000)
dba.plotMA(dbObj, method=DBA_EDGER)
dev.off()

# XY scatter plot 
png(filename = file.path(resDir,"XYScatter_ALLMETHOD.png"), height = 1000, width = 1000)
dba.plotMA(dbObj, bXY=TRUE)
dev.off()

# Volcano plots
png(filename = file.path(resDir,"VOLCANO_DESEQ2.png"), height = 1000, width = 1000)
dba.plotVolcano(dbObj, method = DBA_DESEQ2)
dev.off()

png(filename = file.path(resDir,"VOLCANO_EDGER.png"), height = 1000, width = 1000)
dba.plotVolcano(dbObj, method = DBA_EDGER)
dev.off()

# Box plots
png(filename = file.path(resDir,"BOX_DESEQ2.png"), height = 1000, width = 1000)
dba.plotBox(dbObj, method = DBA_DESEQ2)
dev.off()

png(filename = file.path(resDir,"BOX_EDGER.png"), height = 1000, width = 1000)
dba.plotBox(dbObj, method = DBA_EDGER)
dev.off()

# save results
res_deseq <- dba.report(dbObj, contrast = 1, th = 1)
out <- as.data.frame(res_deseq)

# Write down results
summary_table <- dba.show(dbObj, bContrasts=T)
write.table(x = summary_table, 
            file = file.path(resDir, "diffBind_summary.tsv"), sep = "\t", quote = FALSE)
for (i in 1:(length(summary_table)-1)) {
  res_deseq <- dba.report(dbObj, contrast = i, th = 1)
  out <- as.data.frame(res_deseq)
  write.table(x = out,
              file = file.path(resDir, paste0("res_deseq_", i, ".tsv")),
              sep = "\t",
              quote = FALSE)
}


# Extract count matrix
# For count matrix I will have it for H2AUB in embryos, not including seelings
# because the seedling datasets seem too far off from the embryo ones
dbObj <- dba(sampleSheet=samples_filtered, scoreCol = 5)
dbObj$config$cores <- 8 

# raw
dbObj_raw <- dba.count(dbObj, bUseSummarizeOverlaps=TRUE, bParallel = TRUE,
                       bRemoveDuplicates = TRUE, score = DBA_SCORE_READS)
countmatrix_raw <- dba.peakset(dbObj_raw, bRetrieve = TRUE, DataType = DBA_DATA_FRAME)
write.table(x = countmatrix_raw,
            file = file.path(resDir, paste0("countmatrix_raw.tsv")),
            sep = "\t",
            quote = FALSE)

# DESeq2 normalized
dbObj_DESEQ2normalized <- dba.count(dbObj, bUseSummarizeOverlaps=TRUE, bParallel = TRUE,
                                    bRemoveDuplicates = TRUE, score = DBA_SCORE_NORMALIZED) # Here I actually don't know how to set the normalization method, but the default seems to be DESeq2 (RLE)
countmatrix_DESEQ2normalized <- dba.peakset(dbObj_DESEQ2normalized, bRetrieve = TRUE, DataType = DBA_DATA_FRAME)
write.table(x = countmatrix_DESEQ2normalized,
            file = file.path(resDir, paste0("countmatrix_DESEQ2normalized.tsv")),
            sep = "\t",
            quote = FALSE)

# RPKM
dbObj_RPKM <- dba.count(dbObj, bUseSummarizeOverlaps=TRUE, bParallel = TRUE,
                        bRemoveDuplicates = TRUE, score = DBA_SCORE_RPKM)
countmatrix_RPKM <- dba.peakset(dbObj_RPKM, bRetrieve = TRUE, DataType = DBA_DATA_FRAME)
write.table(x = countmatrix_RPKM,
            file = file.path(resDir, paste0("countmatrix_RPKM.tsv")),
            sep = "\t",
            quote = FALSE)



















