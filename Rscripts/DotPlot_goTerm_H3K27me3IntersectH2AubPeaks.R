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

resDir <- "../results/goTerms_H3K27me3IntersectH3K4me3"
if (!file.exists(resDir)){dir.create(resDir)}

geneList_all <- keys(org.At.tair.db, keytype = "TAIR") # all genes in database

# read in the annotation from ChIPseeker

samplefiles <- list.files("../results/ChIPseeker_H3K27me3_intersect_H3K4me3_peaks", pattern= ".tsv", full.names=T)
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


# 
id="WS_H3K27ME3_OVERLAPWITH_H3K4ME3"

egos <- enrichGO(gene = geneList$WS_H3K27ME3_OVERLAPWITH_H3K4ME3, 
                      keyType = "TAIR",
                      OrgDb = org.At.tair.db, 
                      ont = "BP", 
                      pAdjustMethod = "BH",
                      qvalueCutoff = 0.05, 
                      readable = TRUE,
                      pool = TRUE)

dotplot(egos)

if(length(egos$ID) > 0){
  godotplots <- dotplot(egos, font.size = 10, showCategory=10, title = id) +
    scale_y_discrete(labels=function(x) str_wrap(x, width=25))
  cluster_summary <- data.frame(egos)
  write.csv(cluster_summary, file.path(resDir, paste0("goTable_","clusterprofiler", id, ".go.tsv")), quote = F)
  #
  ggsave(file.path(resDir, paste0("clusterProfilerDotPlots", id,".png")), godotplots, width = 6, height = 2+nrow(godotplots$data)/5)
  #
}else{
  file.create(file.path(resDir_goTable, paste0(id,".noEnrichmentFound.txt")))
}

#

#










# overlap mg H2Aub genes with sdl H3K27me3 up genes
# peaks with increased H3K27me3 from em to sdl
deseq_res_2 <- file.path(paste0("../results/diffBind_ChIPseeker_", "MGvsSDL_H3K27ME3_newDataset"), 
                         "DESeq2Res_withAnnotations.tsv") %>%
  read.table(header = TRUE, fill = TRUE, sep = "\t", quote = "")
stage <- c()
stage[1] <- names(deseq_res_2)[7] %>% str_remove("Conc_")
stage[2] <- names(deseq_res_2)[8] %>% str_remove("Conc_")
FDR_th = 0.05

peakList <- deseq_res_2 %>% filter(FDR < FDR_th) %>%
  filter(.data[[paste0("Conc_", "SDL")]] > .data[[paste0("Conc_", "MG")]])

geneList_sdlK27me3up <- peakList$ASSOCIATED_GENES %>% str_split(", ") %>% unlist %>% unique
geneList_sdlK27me3up <- geneList_sdlK27me3up[!is.na(geneList_sdlK27me3up)]

all_H3K27me3_genes <- deseq_res_2$ASSOCIATED_GENES %>% str_split(", ") %>% unlist %>% unique

# read in gene list
samplefiles <- list.files("../results/ChIPseeker_H2AubIntersectH3K27me3Peaks", pattern= ".tsv", full.names=T)
samplenames <- str_remove(basename(samplefiles),".consensusPeaks.withAnnotation.tsv")
names(samplefiles) <- samplenames

peak_annotation_df <- vector(mode = "list", length = length(samplefiles))
for(i in 1:length(samplefiles)){
  peak_annotation_df[[i]] <- read.table(samplefiles[[i]], header = TRUE, sep = "\t", row.names = 1) %>% as.data.frame()
}
names(peak_annotation_df) <- samplenames

geneList <- vector(mode = "list", length = length(samplenames))
for(i in 1:length(samplenames)){
  geneList_temp <- peak_annotation_df[[i]]$ASSOCIATED_GENES %>% str_split(", ") %>% unlist %>% unique
  geneList[[i]] <- geneList_temp[!is.na(geneList_temp)]
}
names(geneList) <- samplenames
















