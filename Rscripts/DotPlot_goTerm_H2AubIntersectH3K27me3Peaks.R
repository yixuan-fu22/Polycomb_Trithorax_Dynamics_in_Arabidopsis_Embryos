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

resDir <- "../results/goTerms_H2AubIntersectH3K27me3Peaks"
if (!file.exists(resDir)){dir.create(resDir)}

geneList_all <- keys(org.At.tair.db, keytype = "TAIR") # all genes in database

# read in the annotation from ChIPseeker

samplefiles <- list.files("../results/ChIPseeker_H2AubIntersectH3K27me3Peaks", pattern= ".tsv", full.names=T)
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


# mg
mg_all_H2Aub_genes <- c(geneList$MG_H2AUB_OVERLAPWITH_H3K27ME3, 
                         geneList$MG_H2AUB_WITHOUT_H3K27ME3) %>% unique()

egos_with <- enrichGO(gene = geneList$MG_H2AUB_OVERLAPWITH_H3K27ME3, 
                      keyType = "TAIR",
                      OrgDb = org.At.tair.db, 
                      ont = "BP", 
                      pAdjustMethod = "BH",
                      qvalueCutoff = 0.05, 
                      readable = TRUE,
                      pool = TRUE)

dotplot(egos_with)

egos_without <- enrichGO(gene = geneList$MG_H2AUB_WITHOUT_H3K27ME3, 
                         keyType = "TAIR",
                         OrgDb = org.At.tair.db, 
                         ont = "BP", 
                         pAdjustMethod = "BH",
                         qvalueCutoff = 0.05, 
                         readable = TRUE,
                         pool = TRUE,
                         universe = MG_all_H2Aub_genes)

dotplot(egos_without)

egos_all <- enrichGO(gene = sdl_all_H2Aub_genes, 
                     keyType = "TAIR",
                     OrgDb = org.At.tair.db, 
                     ont = "BP", 
                     pAdjustMethod = "BH",
                     qvalueCutoff = 0.05, 
                     readable = TRUE,
                     pool = TRUE)

dotplot(egos_all)

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

mg_all_H2Aub_genes <- c(geneList$MG_H2AUB_OVERLAPWITH_H3K27ME3, geneList$MG_H2AUB_WITHOUT_H3K27ME3) %>% unique
ws_all_H2Aub_genes <- c(geneList$WS_H2AUB_OVERLAPWITH_H3K27ME3, geneList$WS_H2AUB_WITHOUT_H3K27ME3) %>% unique



# against all H3K27me3 genes 
id = "mgH2AubGenes_overlap_sdlH3K27me3UpGenes"
title = "MG H2AUb Genes Overlap H3K27me3 Up Genes in SDL Compared To MG"
all_H3K27me3_genes <- deseq_res_2$ASSOCIATED_GENES %>% str_split(", ") %>% unlist %>% unique
egos <- enrichGO(gene = geneList_sdlK27me3up[(geneList_sdlK27me3up %in% mg_all_H2Aub_genes)], 
                 keyType = "TAIR",
                 OrgDb = org.At.tair.db, 
                 ont = "BP", 
                 pAdjustMethod = "BH",
                 qvalueCutoff = 0.05, 
                 readable = TRUE,
                 pool = TRUE,
                 universe = all_H3K27me3_genes)
godotplots <- dotplot(egos, font.size = 10, showCategory=10, 
                      title = title) +
  scale_y_discrete(labels=function(x) str_wrap(x, width=25))
cluster_summary <- data.frame(egos)
ggsave(file.path(resDir, paste0(id, ".png")), godotplots, width = 6, height = 2+nrow(godotplots$data)/5)
write.csv(cluster_summary, file.path(resDir, paste0(id, ".go.tsv")), quote = F)



# against all H2Aub genes
egos <- enrichGO(gene = mg_all_H2Aub_genes[(mg_all_H2Aub_genes %in% geneList_sdlK27me3up)], 
                 keyType = "TAIR",
                 OrgDb = org.At.tair.db, 
                 ont = "BP", 
                 pAdjustMethod = "BH",
                 qvalueCutoff = 0.05, 
                 readable = TRUE,
                 pool = TRUE,
                 universe = mg_all_H2Aub_genes)
dotplot(egos)
data.frame(egos) %>% view
















