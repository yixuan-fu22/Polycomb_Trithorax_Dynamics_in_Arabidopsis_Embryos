# To compare the distribution of peaks of a same histone modification between different stages

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

peakFolder <- paste0("../results/ChIPseeker_consensusPeaks")

resDir <- paste0("../results/Vennplots_consensusPeaks")

resDir_pairwise <- paste0("../results/Vennplots_consensusPeaks_pairwise")

ids <- c("WS_H2AUB", "MG_H2AUB", "SDL_H2AUB_REPMERGED")

for(ids in list(
  c("WS_H2AUB", "MG_H2AUB", "SDL_H2AUB_REPMERGED"),
  c("WS_K4ME3", "MG_K4ME3", "SDL_K4ME3"), 
  c("WS_K27ME3", "MG_K27ME3", "SDL_H3K27ME3_newDataset")
)){
  print(ids)
}

for(ids in list(
  c("WS_H2AUB", "MG_H2AUB", "SDL_H2AUB_REPMERGED"),
  c("WS_K4ME3", "MG_K4ME3", "SDL_K4ME3"), 
  c("WS_K27ME3", "MG_K27ME3", "SDL_H3K27ME3_newDataset")
)){

hist <- unlist(str_split(ids[[1]], "_"))[2]

peakfileNames <- paste0(ids, ".consensusPeaks.withAnnotation.tsv")

# read files and convert as gr ranges
peaks_dfs <-  vector(mode = "list", length = length(ids))

for(i in 1:length(ids)){
  peaks_dfs[[i]] <- read.table(file = file.path(peakFolder, 
                                        peakfileNames[i]), sep = "\t", quote = "")
}

# make gr objects
consensus_peaks <- vector(mode = "list", length = length(ids))
for(i in 1:length(ids)){
  consensus_peaks[[i]] <- GRanges(seqnames = peaks_dfs[[i]]$CHR, 
                      IRanges(start = peaks_dfs[[i]]$START, end = peaks_dfs[[i]]$END),
                      strand = "*") 
}
names(consensus_peaks) <- ids


# GRange counts
A <- consensus_peaks[[1]]
B <- consensus_peaks[[2]]
C <- consensus_peaks[[3]]

A_total = length(A)
B_total = length(B)
C_total = length(C)

A_B_C    = length(intersect(intersect(A, B), C))
A_B   = sum(countOverlaps(intersect(A, B), C) == 0)
B_C = sum(countOverlaps(intersect(B, C), A) == 0)
A_C   = sum(countOverlaps(intersect(A, C), B) == 0)
A_only = sum(as.integer(countOverlaps(A, B) > 0) + as.integer(countOverlaps(A, C) > 0) == 0)
B_only = sum(as.integer(countOverlaps(B, C) > 0) + as.integer(countOverlaps(B, A) > 0) == 0)
C_only = sum(as.integer(countOverlaps(C, A) > 0) + as.integer(countOverlaps(C, B) > 0) == 0)
counts <- c(A_B_C = A_B_C, A_B = A_B, A_C = A_C, B_C = B_C, 
            A_only = A_only, B_only = B_only, C_only = C_only)

venn.plot <- draw.triple.venn(
  area1 = A_only + A_B + A_C + A_B_C,
  area2 = B_only + A_B + B_C + A_B_C,
  area3 = C_only + A_C + B_C + A_B_C,
  n12   = A_B + A_B_C,
  n13   = A_C + A_B_C,
  n23   = B_C + A_B_C,
  n123  = A_B_C,
  category = paste0(ids, "\nTotal:", c(A_total, B_total, C_total)),
  cat.dist = c(0.06, 0.06, 0.06),
  fill = c("#CC79A7", "#56B4E9", "#F0E442"),
  alpha = 0.3,
  scaled = TRUE,
  fontfamily = "Arial",
  cat.fontfamily = "Arial"
)

png(file.path(resDir, paste0("vennplot_GRangesCounts_", hist, ".png" )))
grid.draw(venn.plot)
dev.off()

# plot pairwise venn
color_platte <- c("#CC79A7", "#56B4E9", "#F0E442")
for(i in 1:length(ids)){
  for(j in 1:i){
    if(i!=j){
      I = consensus_peaks[[i]]
      J = consensus_peaks[[j]]
      I_total = length(I)
      J_total = length(J)
      I_only <-  sum(countOverlaps(I, J) == 0)
      J_only <- sum(countOverlaps(J,I) == 0)
      I_J <- length(intersect(I,J))
      I_list <- c(paste0("I_only_", 1:I_only), paste0("I_J", 1:I_J))
      J_list <- c(paste0("J_only_", 1:J_only), paste0("I_J", 1:I_J))
      venn.diagram(list(I_list, J_list), 
                   category.names = paste0(c(ids[i], ids[j]), "\n", c(I_total, J_total)),
                   fill = c(color_platte[i], color_platte[j]),
                   filename = file.path(resDir_pairwise, paste0("vennPairwise_GRangesCounts_",
                                                       hist, "_", 
                                                       i+j-2, ".png" )),
                   fontfamily = "Arial",
                   cat.fontfamily = "Arial",
                   disable.logging = TRUE,
                   hyper.test = F, cat.dist = c(0.06, 0.06), margin = 0.1 )
    }
    
    
  }
}

}














