# The following codes annotate the results of diffBind
# it first annotates the results of the counting matrix which will be used for plotting
# then it takes the significantly differential peak and test for GO term enrichment

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
  library(clusterProfiler)
  #
  #volcano plot
  #library(gridExtra)
  #library(ggimage)
  #library(ggrepel)
})

comparison <- "WSvsSDL_H2AUB_newDataset"
diffBindFolder <- paste0("../results/diffBind_ChIPseeker_", comparison)

# make a result dir
resDir <- paste0("../results/diffBind_ChIPseeker_", comparison)
if(!file.exists(resDir)){dir.create(resDir)}

##############  Annotate Counting Matrix #################
  
# specify the count matrix to use
normalization_method <- "DESEQ2normalized" # select from "raw" "DESEQ2normalized" "RPKM"

# read count matrix
countmatrix <- read.table(file = file.path(diffBindFolder, 
                                           paste0("countmatrix_", normalization_method, ".tsv")))
countmatrix <- countmatrix %>% filter(!(CHR %in% c("ChrC", "ChrM")))


# make a gr object for peaks
peaks_gr <- GRanges(seqnames = countmatrix$CHR, 
                    IRanges(start = countmatrix$START, end = countmatrix$END),
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
countmatrix_with_annotations <- left_join(countmatrix, annotations_with_names)

# write down the count matrix with annotations
write.table(x = countmatrix_with_annotations,
            file = file.path(resDir, paste0("countmatrix_with_annotations_", 
                                            normalization_method, 
                                            ".tsv")),
            sep = "\t",
            quote = FALSE)

############## Export DE Peaks and Test for GO Enrichment #################
# read DESeq2 results
DEseq2File <- file.path(diffBindFolder, "res_deseq_1.tsv")
DESeq2Res <- read.table(DEseq2File, header = TRUE, fill = TRUE, sep = "\t") 

annotation_only <- countmatrix_with_annotations[c(1,2,3, (ncol(countmatrix_with_annotations)-1):ncol(countmatrix_with_annotations))]

DESeq2Res_withAnnotation <- left_join(DESeq2Res, annotation_only, by = join_by("seqnames" == "CHR",
                                                   "start" == "START",
                                                   "end" == "END"))

# writing down annotated DESeq2 results
write.table(DESeq2Res_withAnnotation, file = file.path(resDir, "DESeq2Res_withAnnotations.tsv"),
            sep = "\t", quote = F)

# Extract DE peaks
peakList_L2FCplus <- DESeq2Res_withAnnotation %>% filter(FDR < 0.05) %>% filter(Fold > 0) # Fold > 0 means MG > WS in the context of WSvsMG
geneList_DEpeaksL2FCplus <- peakList_L2FCplus$ASSOCIATED_GENES %>% strsplit(", ") %>% unlist
geneList_DEpeaksL2FCplus <- geneList_DEpeaksL2FCplus[!(geneList_DEpeaksL2FCplus == "NA")]

peakList_L2FCminus <- DESeq2Res_withAnnotation %>% filter(FDR < 0.05) %>% filter(Fold < 0)
geneList_DEpeaksL2FCminus <- peakList_L2FCminus$ASSOCIATED_GENES %>% strsplit(", ") %>% unlist
geneList_DEpeaksL2FCminus <- geneList_DEpeaksL2FCminus[!(geneList_DEpeaksL2FCminus == "NA")]

write.table(geneList_DEpeaksL2FCplus, file = file.path(resDir, "geneList.DEpeaks.L2FCplus.txt"),
            quote = F, col.names = F, row.names = F)
write.table(geneList_DEpeaksL2FCminus, file = file.path(resDir, "geneList.DEpeaks.L2FCminus.txt"),
            quote = F, col.names = F, row.names = F) # write down the genes

# Test for enrichment for positive FC genes
egos <- enrichGO(gene = geneList_DEpeaksL2FCplus, 
                 keyType = "TAIR",
                 OrgDb = org.At.tair.db, 
                 ont = "all", 
                 pAdjustMethod = "BH",
                 qvalueCutoff = 0.05, 
                 readable = TRUE,
                 pool = TRUE)
length(egos$ID) > 0
if(length(egos$ID) > 0){
  godotplots <- dotplot(egos, showCategory=50, title = "geneList_DEs L2FC plus")
  cluster_summary <- data.frame(egos)
  write.csv(cluster_summary, file.path(resDir, paste0("geneList_DEs_L2FC_plus",".go.clusterprofiler_summary.txt")))
  #
  ggsave(file.path(resDir, paste0("geneList_DEs_L2FC_plus",".go.clusterprofiler_summary.png")), godotplots, width = 7, height = 4+nrow(godotplots$data)/1.5)
}else{
  file.create(file.path(resDir, paste0("geneList_DEs_L2FC_plus",".noEnrichmentFound.txt")))
}

# Testing enrichment for genes with negative FC
egos <- enrichGO(gene = geneList_DEpeaksL2FCminus, 
                 keyType = "TAIR",
                 OrgDb = org.At.tair.db, 
                 ont = "all", 
                 pAdjustMethod = "BH",
                 qvalueCutoff = 0.05, 
                 readable = TRUE,
                 pool = TRUE)
length(egos$ID) > 0
if(length(egos$ID) > 0){
  godotplots <- dotplot(egos, showCategory=50, title = "geneList_DEs L2FC minus")
  cluster_summary <- data.frame(egos)
  write.csv(cluster_summary, file.path(resDir, paste0("geneList_DEs_L2FC_minus",".go.clusterprofiler_summary.txt")))
  #
  ggsave(file.path(resDir, paste0("geneList_DEs_LOG2FC_minus",".go.clusterprofiler_summary.png")), godotplots, width = 7, height = 4+nrow(godotplots$data)/1.5)
}else{
  file.create(file.path(resDir, paste0("geneList_DEs_LOG2FC_minus",".noEnrichmentFound.txt")))
}











































