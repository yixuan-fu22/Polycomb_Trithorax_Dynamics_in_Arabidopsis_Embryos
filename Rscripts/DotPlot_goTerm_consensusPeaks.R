# plot GO terms for all genes marked by a histone mark at one stage

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

peakFolder <- paste0("../results/ChIPseeker_consensusPeaks")

# make result dirs
GOdir <- paste0("../results/goTerms_consensusPeaks")
#GOdir = ChIPseekerFolder
if(!file.exists(GOdir)){dir.create(GOdir)}

resDir_dot_clusterProfiler <- GOdir
if(!file.exists(resDir_dot_clusterProfiler)){dir.create(resDir_dot_clusterProfiler)}

resDir_goTable <- GOdir

#
ids = c("WS_K4ME3", "WS_H2AUB", "WS_K27ME3",
  "MG_K4ME3", "MG_H2AUB", "MG_K27ME3",
  "SDL_K4ME3", "SDL_H2AUB_REPMERGED", "SDL14D_K27ME3")

tissue <- unlist(str_split(id, "_"))[1]
hist <- unlist(str_split(id, "_"))[2]

for(id in ids){

peakfileName <- paste0(id, ".consensusPeaks.withAnnotation.tsv")
peaks_df <- read.table(file = file.path(peakFolder, 
                                                peakfileName), sep = "\t", quote = "")



allGeneList <- peaks_df$ASSOCIATED_GENES %>% strsplit(", ") %>% unlist
allGeneList <- allGeneList[!is.na(allGeneList)] %>% unique()

egos <- enrichGO(gene = allGeneList, 
                 keyType = "TAIR",
                 OrgDb = org.At.tair.db, 
                 ont = "all", 
                 pAdjustMethod = "BH",
                 qvalueCutoff = 0.05, 
                 readable = TRUE,
                 pool = TRUE)

if(length(egos$ID) > 0){
  # unsimpliefied
  godotplots <- dotplot(egos, font.size = 10, showCategory=10, title = paste0(id, "consensusPeaks.go.png")) +
    scale_y_discrete(labels=function(x) str_wrap(x, width=25))
  cluster_summary <- data.frame(egos)
  write.csv(cluster_summary, file.path(resDir_goTable, paste0("goTable_","clusterprofiler", id, ".go.tsv")), quote = F)
  #
  ggsave(file.path(resDir_dot_clusterProfiler, paste0("clusterProfilerDotPlots", id,".png")), godotplots, width = 6, height = 2+nrow(godotplots$data)/5)
  #
  # simplify
  egos_simplified <- simplify(egos,
                     cutoff = 0.7,
                     by = "p.adjust",
                     select_fun = min,
                     measure = "Wang",
                     semData = NULL
    )
  godotplots_simplified <- dotplot(egos_simplified, font.size = 10, showCategory=10, title = paste0(id, ".consensusPeaks.go.simplified.png")) +
    scale_y_discrete(labels=function(x) str_wrap(x, width=25))
  cluster_summary <- data.frame(egos_simplified)
  write.csv(cluster_summary, file.path(resDir_goTable, paste0("goTable_","clusterprofiler.", id, ".go.simplified.tsv")), quote = F)
  #
  ggsave(file.path(resDir_dot_clusterProfiler, paste0("clusterProfilerDotPlots", id,".simplified.png")), godotplots_simplified, width = 6, height = 2+nrow(godotplots$data)/5)
  #
  
}else{
  file.create(file.path(resDir_goTable, paste0(id,".noEnrichmentFound.txt")))
}


dotplot(egos)
dotplot(egos_simplified)

}


