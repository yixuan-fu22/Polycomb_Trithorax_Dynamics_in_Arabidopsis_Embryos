# To annotate consensus peaks marked by a certain histone modification

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
  library(ChIPseeker)
  library(AnnotationDbi)
  library(TxDb.Athaliana.BioMart.plantsmart25)
  library(org.At.tair.db)
  #library(ggupset)
})

peakFolder <- paste0("../data_fromOtherProjects/bedtools_consensusPeak")

# make a result dir
resDir <- paste0("../results/ChIPseeker_consensusPeaks")
if(!file.exists(resDir)){dir.create(resDir)}

##############  Annotate Peaks #################

id = "SDL_H3K27ME3_newDataset"

{

# read peaks files
peakfileName <- paste0(id, ".consensusPeaks.bed")

peakdf <- read.table(file = file.path(peakFolder, 
                                      peakfileName))

names(peakdf) <- c("CHR", "START", "END", "SCORE")
peakdf <- peakdf %>% filter(!(CHR %in% c("ChrC", "ChrM")))


# make a gr object for peaks
peaks_gr <- GRanges(seqnames = peakdf$CHR, 
                    IRanges(start = peakdf$START, end = peakdf$END),
                    strand = "*") 

# import txdb
txdb <- TxDb.Athaliana.BioMart.plantsmart25
seqlevels(txdb) # check chromosome names
seqlevels(txdb) <- c("Chr1","Chr2","Chr3","Chr4","Chr5","ChrM","ChrC")

# annotate peaks
peakAnno <- annotatePeak(peaks_gr, TxDb = txdb, tssRegion=c(-500, 0),
                         genomicAnnotationPriority = c("Exon", "Promoter", "5UTR", "3UTR", "Intron",
                                                       "Downstream", "Intergenic"),
                         addFlankGeneInfo = TRUE, flankDistance = 500, verbose=TRUE)
annotations <- peakAnno@anno %>% as.data.frame()

# 
# import the table for mapping TAIR IDs to gene names
araport_split_tairID_and_symbol <- read.csv(file.path("..", "reference_data", "GFF_araport", "araport_split_tairID_and_symbol.tsv"),
                                            sep = "\t")

# use flank_geneIds as the annotation
expansion <- annotations %>% separate_longer_delim(cols = flank_geneIds, delim = ";")
expension_withnames <- left_join(expansion, araport_split_tairID_and_symbol, join_by(flank_geneIds == geneID), relationship = "many-to-many")
expension_withnames_collapsed <- expension_withnames %>% unique %>% group_by(across(!c(flank_geneIds, symbol))) %>% summarize_all(paste, collapse=", ")
annotations_with_names <- data.frame(CHR=expension_withnames_collapsed$seqnames, 
                                     START=expension_withnames_collapsed$start, 
                                     END=expension_withnames_collapsed$end,
                                     ASSOCIATED_GENES=expension_withnames_collapsed$flank_geneIds,
                                     GENE_NAMES=expension_withnames_collapsed$symbol)

# join the annotated peaks back to the count matrix
peaks_annotated <- annotations_with_names

# write down the count matrix with annotations
write.table(x = peaks_annotated,
            file = file.path(resDir, paste0(id, ".consensusPeaks", ".withAnnotation", 
                                            ".tsv")),
            sep = "\t",
            quote = FALSE)

}






