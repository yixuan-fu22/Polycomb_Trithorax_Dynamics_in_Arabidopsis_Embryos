#
# To Plot GO terms made by clusterProfiler
# While the clusterProfiler programs make a dot plot, I will need to plot myself to make it adjustable for panels
# And, for some omparisons with a lot of GO term coming out, I need to reduce redandency
# For better visualization, I could use bar plot (like fig2 from "H2A monoubiquitination in Arabidopsisthaliana is generally independent of LHP1and PRC2 activity")
# or simply tables (like figures from "Canonical and Noncanonical Actions of Arabidopsis HistoneDeacetylases in Ribosomal RNA Processing")

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
  # plots
  library(ggplot2)
  # library(ggpubr)
  # library(ggupset)
})



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


# To have cleaner results, start from annotated DEseq results
suppressMessages({
  library(clusterProfiler)
  library(org.At.tair.db)
})

# start

comparison <- "WSvsSDL_H2AUB"
unlist(str_split(comparison, "_"))[2]
hist = unlist(str_split(comparison, "_"))[2]
ChIPseekerFolder <- paste0("../results/diffBind_ChIPseeker_", comparison)

# make result dirs
resDir_dot_clusterProfiler <- ChIPseekerFolder
if(!file.exists(resDir_dot_clusterProfiler)){dir.create(resDir_dot_clusterProfiler)}

resDir_goTable <- ChIPseekerFolder
if(!file.exists(resDir_goTable)){dir.create(resDir_goTable)}

DEseq2_res_annotated <- file.path(ChIPseekerFolder,
                                  "DESeq2Res_withAnnotations.tsv") %>%
  read.table(header = TRUE, fill = TRUE, sep = "\t", quote = "")

stage <- c()
stage[1] <- names(DEseq2_res_annotated)[7] %>% str_remove("Conc_")
stage[2] <- names(DEseq2_res_annotated)[8] %>% str_remove("Conc_")

# Extract DE peaks
geneList_upInStage <- vector("list", length = 2)
for(i in 1:2){
  j = 3 - i # opposite of i
  peakList <- DEseq2_res_annotated %>% filter(FDR < 0.05) %>% 
    filter(.data[[paste0("Conc_", stage[i])]] > .data[[paste0("Conc_", stage[j])]]) # Fold > 0 means MG > WS in the context of WSvsMG
  geneList_temp <- peakList$ASSOCIATED_GENES %>% strsplit(", ") %>% unlist
  geneList_upInStage[[i]] <- geneList_temp[!is.na(geneList_temp)]
}
  # enrichment
for(i in 1:2){
  j = 3 - i
  id <- paste0("Up_in_", stage[i], ";comparison", stage[i], "vs", stage[j], ";hist", hist, "fontSize10")
  egos <- enrichGO(gene = geneList_upInStage[[i]], 
                   keyType = "TAIR",
                   OrgDb = org.At.tair.db, 
                   ont = "all", 
                   pAdjustMethod = "BH",
                   qvalueCutoff = 0.05, 
                   readable = TRUE,
                   pool = TRUE)
  if(length(egos$ID) > 0){
    # simplify
    egos <- simplify(egos,
             cutoff = 0.7,
             by = "p.adjust",
             select_fun = min,
             measure = "Wang",
             semData = NULL
    )
    #
    godotplots <- dotplot(egos, font.size = 10, showCategory=10, title = paste0("Up in ", stage[i], "; comparison", stage[i], "vs", stage[j], "_", hist)) +
      scale_y_discrete(labels=function(x) str_wrap(x, width=25))
    cluster_summary <- data.frame(egos)
    write.csv(cluster_summary, file.path(resDir_goTable, paste0("goTable_","clusterprofiler", id, ".go.tsv")), quote = F)
    #
    ggsave(file.path(resDir_dot_clusterProfiler, paste0("clusterProfilerDotPlots", id,".simplfied.png")), godotplots, width = 6, height = 2+nrow(godotplots$data)/5)
    #
  }else{
    file.create(file.path(resDir_goTable, paste0(id,".noEnrichmentFound.txt")))
  }
}

dotplot(egos)

