# plot GO terms

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
  # plots
  library(ggplot2)
  # library(ggpubr)
  # library(ggupset)
})

suppressMessages({
  library(clusterProfiler)
  library(org.At.tair.db)
})

resDir <- "../results/goTerms_substractedPeaks"
if (!file.exists(resDir)){dir.create(resDir)}

geneList_all <- keys(org.At.tair.db, keytype = "TAIR") # all genes in database


samplefiles <- list.files("../results/ChIPseeker_substractedPeaks", pattern= ".tsv", full.names=T)
samplenames <- str_remove(basename(samplefiles),".consensusPeaks.withAnnotation.tsv")
names(samplefiles) <- samplenames

peak_annotation_df <- vector(mode = "list", length = length(samplefiles))
for(i in 1:length(samplefiles)){
  peak_annotation_df[[i]] <- read.table(samplefiles[[i]], header = TRUE, sep = "\t", row.names = 1) %>% as.data.frame()
}
names(peak_annotation_df) <- samplenames

# read in gene list
geneList <- vector(mode = "list", length = length(samplenames))
for(i in 1:length(samplenames)){
  geneList_temp <- peak_annotation_df[[i]]$ASSOCIATED_GENES %>% str_split(", ") %>% unlist %>% unique
  geneList[[i]] <- geneList_temp[!is.na(geneList_temp)]
}
names(geneList) <- samplenames

egos <- enrichGO(gene = geneList$H2AUB_MG_substract_WS, 
                 keyType = "TAIR",
                 OrgDb = org.At.tair.db, 
                 ont = "BP", 
                 pAdjustMethod = "BH",
                 qvalueCutoff = 0.05, 
                 readable = TRUE,
                 pool = TRUE)

godotplots <- dotplot(egos)

cluster_summary <- data.frame(egos)

ggsave(file.path(resDir, paste0("H2AUB_MG_substract_WS", ".png")), godotplots, width = 6, height = 2+nrow(godotplots$data)/5)
write.csv(cluster_summary, file.path(resDir, paste0("H2AUB_MG_substract_WS", ".go.tsv")), quote = F)

egos <- enrichGO(gene = geneList$H2AUB_WS_substract_MG, 
                 keyType = "TAIR",
                 OrgDb = org.At.tair.db, 
                 ont = "BP", 
                 pAdjustMethod = "BH",
                 qvalueCutoff = 0.05, 
                 readable = TRUE,
                 pool = TRUE)

godotplots <- dotplot(egos)

cluster_summary <- data.frame(egos)

ggsave(file.path(resDir, paste0("H2AUB_WS_substract_MG", ".png")), godotplots, width = 6, height = 2+nrow(godotplots$data)/5)
write.csv(cluster_summary, file.path(resDir, paste0("H2AUB_WS_substract_MG", ".go.tsv")), quote = F)




















