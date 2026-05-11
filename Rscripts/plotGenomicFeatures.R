# To plot the genomic feature distribution with ChIPseeker of peaks at different stages

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
  # plot genomic feature distribution
  library(ChIPpeakAnno)
  library(TxDb.Athaliana.BioMart.plantsmart25)
  #library(ggupset)
})


for(ids in list(
  c("WS_H2AUB", "MG_H2AUB", "SDL_H2AUB_REPMERGED"),
  c("WS_K4ME3", "MG_K4ME3", "SDL_K4ME3"), 
  c("WS_K27ME3", "MG_K27ME3", "SDL14D_K27ME3")
)){


peakFolder <- paste0("../results/ChIPseeker_consensusPeaks")

peakfileNames <- paste0(ids, ".consensusPeaks.withAnnotation.tsv")
#

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

hist <- unlist(str_split(ids[[1]], "_"))[2]

# Exon / Intron and other features # requires ChIPpeakAnno
# gene feature pie
geneFeatureDir <- file.path("../results", "peaksOnGeneFeatures")
if (!file.exists(geneFeatureDir)){dir.create(geneFeatureDir)}

for(i in 1:length(ids)){
  png(filename = file.path(geneFeatureDir,paste0(names(consensus_peaks)[i],".pie.png")), width = 500, height = 500)
  genomicElementDistribution(consensus_peaks[[i]],
                             promoterRegion = c(upstream = 500, downstream = 50),
                             TxDb = TxDb.Athaliana.BioMart.plantsmart25)
  dev.off()
}

# ChIPpeakanno bar plot
consensus_peaks_grlist <- GRangesList(consensus_peaks[[1]],
                                      consensus_peaks[[2]],
                                      consensus_peaks[[3]])

names(consensus_peaks_grlist) <- names(consensus_peaks)
names(consensus_peaks_grlist) <- paste0(length(ids):1, "_", names(consensus_peaks))

png(filename = file.path(geneFeatureDir, paste0("peaksOnGeneFeatures.ChIPpeakanno", hist , ".bar.png")), width = 750, height = 750)
genomicElementDistribution(consensus_peaks_grlist,
                           promoterRegion = c(upstream = 500, downstream = 50),
                           TxDb = TxDb.Athaliana.BioMart.plantsmart25)
dev.off()

}



# Peak Distribution around genomic features
# distance to TSS
peakDistributionDir <- file.path("../results", "peakDisAroundTSS")
if (!file.exists(peakDistributionDir)){dir.create(peakDistributionDir)}

annotation_data <- transcripts(TxDb.Athaliana.BioMart.plantsmart25)
for(i in 1:length(ids)){
  png(filename = file.path(peakDistributionDir,paste0(names(ids)[i],".png")), width = 1000, height = 1000)
  binOverFeature(consensus_peaks[[i]], 
                 featureSite = "FeatureStart",
                 nbins = 20,
                 annotationData = annotation_data,
                 xlab = "peak distance from TSS (bp)", 
                 ylab = "peak count", 
                 main = "Distribution of consensus peak numbers around TSS")
  dev.off()
}

