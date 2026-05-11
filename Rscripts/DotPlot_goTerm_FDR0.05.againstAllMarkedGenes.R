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

# To have cleaner results, start from annotated DEseq results
suppressMessages({
  library(clusterProfiler)
  library(org.At.tair.db)
})

# start

comparison <- "WSvsSDL_K4ME3"

{
  unlist(str_split(comparison, "_"))[2]
  hist = unlist(str_split(comparison, "_"))[2]
  ChIPseekerFolder <- paste0("../results/diffBind_ChIPseeker_", comparison)
  
  # make result dirs
  GOdir <- paste0("../results/goTerms_FDR0.01")
  #GOdir = ChIPseekerFolder
  if(!file.exists(GOdir)){dir.create(GOdir)}
  
  resDir_dot_clusterProfiler <- GOdir
  if(!file.exists(resDir_dot_clusterProfiler)){dir.create(resDir_dot_clusterProfiler)}
  
  resDir_goTable <- GOdir
  if(!file.exists(resDir_goTable)){dir.create(resDir_goTable)}
  
  DEseq2_res_annotated <- file.path(ChIPseekerFolder,
                                    "DESeq2Res_withAnnotations.tsv") %>%
    read.table(header = TRUE, fill = TRUE, sep = "\t", quote = "")
  
  stage <- c()
  stage[1] <- names(DEseq2_res_annotated)[7] %>% str_remove("Conc_")
  stage[2] <- names(DEseq2_res_annotated)[8] %>% str_remove("Conc_")
  
  # Extract all genes with the histone modification as background for GO testing
  allGeneList <- DEseq2_res_annotated$ASSOCIATED_GENES %>% strsplit(", ") %>% unlist
  allGeneList <- allGeneList[!is.na(allGeneList)] %>% unique()
  
  # Extract DE peaks
  FDR_th=0.01
  Fold_th=0
  
  geneList_upInStage <- vector("list", length = 2)
  peaknumbers <- c()
  genenumbers <- c()
  for(i in 1:2){
    j = 3 - i # opposite of i
    peakList <- DEseq2_res_annotated %>% filter(FDR < FDR_th) %>% filter(abs(Fold) > Fold_th) %>%
      filter(.data[[paste0("Conc_", stage[i])]] > .data[[paste0("Conc_", stage[j])]]) # Fold > 0 means MG > WS in the context of WSvsMG
    geneList_temp <- peakList$ASSOCIATED_GENES %>% strsplit(", ") %>% unlist
    geneList_upInStage[[i]] <- geneList_temp[!is.na(geneList_temp)]
    peaknumbers[i] <- dim(peakList)[1]
    genenumbers[i] <- length(geneList_upInStage[[i]])
  }
  deSummary <- data.frame(V1=c("comparison", "FDR", paste0("peaks_Up_in_", stage[1]), paste0("peaks_Up_in_", stage[2]),
                               paste0("genes_Up_in_", stage[1]), paste0("genes_Up_in_", stage[2])),
                          V2=c(comparison, FDR_th, peaknumbers[1], peaknumbers[2], genenumbers[1], genenumbers[2]))
  
  write.table(deSummary, file = file.path(GOdir, paste0("deSummary_", comparison, "_FDR_", FDR_th, ".tsv")),
              sep = "\t", quote = F, col.names = F, row.names = F) 
  
  # enrichment
  for(i in 1:2){
    j = 3 - i
    id <- paste0("Up_in_", stage[i], ";comparison", stage[i], "vs", stage[j], ";hist", hist, ".againstAllGenesMarked", ".FDR", FDR_th, ".BP") # font size 10, wrap 25; font size 8, wrap 40
    egos <- enrichGO(gene = geneList_upInStage[[i]], 
                     keyType = "TAIR",
                     OrgDb = org.At.tair.db, 
                     ont = "BP", 
                     pAdjustMethod = "BH",
                     qvalueCutoff = 0.05, 
                     readable = TRUE,
                     pool = TRUE,
                     universe = allGeneList)
    if(length(egos$ID) > 0){
      # simplify
     # egos <- simplify(egos,
    #                   cutoff = 0.7,
    #                   by = "p.adjust",
    #                   select_fun = min,
    #                   measure = "Wang",
    #                   semData = NULL
    #  )
      #
      godotplots <- dotplot(egos, font.size = 10, showCategory=10, title = paste0("Up in ", stage[i], "; comparison", stage[i], "vs", stage[j], "_", hist)) +
        scale_y_discrete(labels=function(x) str_wrap(x, width=25))
      cluster_summary <- data.frame(egos)
      write.csv(cluster_summary, file.path(resDir_goTable, paste0("goTable_","clusterprofiler", id, ".go.tsv")), quote = F)
      #
      ggsave(file.path(resDir_dot_clusterProfiler, paste0("clusterProfilerDotPlots", id,".png")), godotplots, width = 6, height = 2+nrow(godotplots$data)/5)
      #
    }else{
      file.create(file.path(resDir_goTable, paste0(id,".noEnrichmentFound.txt")))
    }
  }
  
  dotplot(egos)
  
}



'
# code for self ploting
# read GO term list
goTxt_path <- file.path("..", "results", "diffBind_ChIPseeker_WSvsMG_K4ME3",
                        paste0("geneList_DEs_L2FC_minus", ".go.clusterprofiler_summary.txt"))

goTable <- read.csv(goTxt_path) # these can be read well with read.csv default settings

# make dotplots
# count = size of dots, color = p.value, x axis = ratio (enriched/total)
goTable <- goTable %>% mutate(GeneRatio_numeric = (goTable$GeneRatio %>% str_split("/")) %>% 
                                sapply(function(x) as.numeric(x[1]) / as.numeric(x[2])))

d <- ggplot() +
  theme_bw() +
  geom_point(data = goTable, aes(x = GeneRatio_numeric, y = reorder(reorder(Description, -p.adjust), GeneRatio_numeric), 
                                 size = Count, color = p.adjust)) + # sort y axis first according to gene ratio, then p value 
  scale_color_gradient(high = "red", low = "blue", limits = range(exp(-16), 0.05), trans = "reverse") # high = smaller p value # exp(-16): to avoid having 0 on the color scale; if a p value small than e-16, it would appear grey

d
'
