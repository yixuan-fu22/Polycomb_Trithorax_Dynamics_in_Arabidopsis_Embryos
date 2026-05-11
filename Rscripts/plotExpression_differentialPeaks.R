# This script is to plot the expression of genes associated with differential peaks between two stages

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
  library(ggbeeswarm)
  library(ggpubr)
})

############## Preprocessing TPM matrix ##############
# Get Expression Matrix
# read in the transcriptome data
{
tpmfile <- "../data_fromOtherProjects/tpm_fromRNAseqOfEmbryosSeedlingsAndFlowers/embryos_sdl_flw.deseq2normalized.tpm.csv"

readtpmfile <- read.table(tpmfile, sep = ",") %>% t()

tpm_matrix <- readtpmfile[-1,]
colnames(tpm_matrix) <- readtpmfile[1,]
colnames(tpm_matrix)[1] <- "gene"
tpm_matrix <- tpm_matrix %>% as.data.frame
#view(tpm_matrix)

# calculate log mean expression at each stage

# Late torpedo stage from Hoffmann et al 2017: SRR8054371, SRR8054372, SRR8054373
# MG stage: SRR8054377, SRR8054378, SRR8054379
# SDL1: ERR10163199, ERR10163200
# SDL2: SRR2500947, SRR1931614

#
tpm_matrix <- tpm_matrix %>%
  filter(str_starts(gene, "AT"))

# calculate log mean expression
# sdl
tpm_matrix$ERR10163199 <- as.numeric(tpm_matrix$ERR10163199)
tpm_matrix$ERR10163200 <- as.numeric(tpm_matrix$ERR10163200)
tpm_matrix$SRR2500947 <- as.numeric(tpm_matrix$SRR2500947)
tpm_matrix$SRR1931614 <- as.numeric(tpm_matrix$SRR1931614)
tpm_matrix$mean_expression_sdl <- rowMeans(tpm_matrix[, c("ERR10163199", "ERR10163200",
                                                          "SRR2500947", "SRR1931614")])

tpm_matrix$logmean_expression_sdl <- log(tpm_matrix$mean_expression_sdl+1)

# ws
tpm_matrix$SRR8054371 <- as.numeric(tpm_matrix$SRR8054371)
tpm_matrix$SRR8054372 <- as.numeric(tpm_matrix$SRR8054372)
tpm_matrix$SRR8054373 <- as.numeric(tpm_matrix$SRR8054373)

tpm_matrix$mean_expression_ws <- rowMeans(tpm_matrix[, c("SRR8054371", "SRR8054372", "SRR8054373")])

tpm_matrix$logmean_expression_ws <- log(tpm_matrix$mean_expression_ws+1)

# mg
tpm_matrix$SRR8054377 <- as.numeric(tpm_matrix$SRR8054377)
tpm_matrix$SRR8054378 <- as.numeric(tpm_matrix$SRR8054378)
tpm_matrix$SRR8054379 <- as.numeric(tpm_matrix$SRR8054379)

tpm_matrix$mean_expression_mg <- rowMeans(tpm_matrix[, c("SRR8054377", "SRR8054378", "SRR8054379")])

tpm_matrix$logmean_expression_mg <- log(tpm_matrix$mean_expression_mg+1)
}
############## ##############


FDR_th = 0.05

comparison <- "WSvsMG_H2AUB" # WSvsMG WSvsMG, H2AUB, H3K27ME3, K4ME3

{
# make result dirs
resDir <- paste0("../results/plotExpression_diffPeaks_FDR", FDR_th)
if(!file.exists(resDir)){dir.create(resDir)}

# start from DEseq results
unlist(str_split(comparison, "_"))[2]
hist = unlist(str_split(comparison, "_"))[2]
ChIPseekerFolder <- paste0("../results/diffBind_ChIPseeker_", comparison)


DEseq2_res_annotated <- file.path(ChIPseekerFolder,
                                  "DESeq2Res_withAnnotations.tsv") %>%
  read.table(header = TRUE, fill = TRUE, sep = "\t", quote = "")

stage <- c()
stage[1] <- names(DEseq2_res_annotated)[7] %>% str_remove("Conc_")
stage[2] <- names(DEseq2_res_annotated)[8] %>% str_remove("Conc_")

# Extract all genes with the histone modification
allGeneList <- DEseq2_res_annotated$ASSOCIATED_GENES %>% strsplit(", ") %>% unlist
allGeneList <- allGeneList[!is.na(allGeneList)] %>% unique()


geneList_upInStage <- vector("list", length = 2)
peaknumbers <- c()
genenumbers <- c()
for(i in 1:2){
  j = 3 - i # opposite of i
  peakList <- DEseq2_res_annotated %>% filter(FDR < FDR_th) %>%
    filter(.data[[paste0("Conc_", stage[i])]] > .data[[paste0("Conc_", stage[j])]]) # Fold > 0 means MG > WS in the context of WSvsMG
  geneList_temp <- peakList$ASSOCIATED_GENES %>% strsplit(", ") %>% unlist
  geneList_upInStage[[i]] <- geneList_temp[!is.na(geneList_temp)]
  peaknumbers[i] <- dim(peakList)[1]
  genenumbers[i] <- length(geneList_upInStage[[i]])
}

deSummary <- data.frame(V1=c("comparison", "FDR", paste0("peaks_Up_in_", stage[1]), paste0("peaks_Up_in_", stage[2]),
                             paste0("genes_Up_in_", stage[1]), paste0("genes_Up_in_", stage[2])),
                        V2=c(comparison, FDR_th, peaknumbers[1], peaknumbers[2], genenumbers[1], genenumbers[2]))

write.table(deSummary, file = file.path(resDir, paste0("deSummary_", comparison, "_FDR_", FDR_th, ".tsv")),
            sep = "\t", quote = F, col.names = F, row.names = F) 

# Non-significant genes
genes_NS <- allGeneList[!(allGeneList %in% c(geneList_upInStage[[1]], geneList_upInStage[[2]]))]

# Plot Expression
tpm_matrix <- mutate(tpm_matrix, 
                     UP_IN_STAGE_1=ifelse(gene %in% geneList_upInStage[[1]], 1, 0),
                     UP_IN_STAGE_2=ifelse(gene %in% geneList_upInStage[[2]], 1, 0),
                     NS_GENES=ifelse(gene %in% genes_NS, 1, 0))

column_name_stage1 <- paste0("logmean_expression_", str_to_lower(stage[[1]]))
column_name_stage2 <- paste0("logmean_expression_", str_to_lower(stage[[2]]))

# plot genes up in stage 1
df_long_stage1 <- tpm_matrix %>% filter(UP_IN_STAGE_1 == 1) %>%
  dplyr::select(.data[[column_name_stage1]], .data[[column_name_stage2]]) %>%
  pivot_longer(
    cols = everything(),
    names_to = "condition",
    values_to = "expression"
  )

p_val_1 <- wilcox.test(
  x = df_long_stage1$expression[df_long_stage2$condition == column_name_stage1],
  y = df_long_stage1$expression[df_long_stage2$condition == column_name_stage2],
  paired = TRUE,
  exact = FALSE
)$p.value
p_txt <- paste0("p (wilcox.test, paired = T) = ", signif(p_val_1, 3))

g1 <- ggplot(df_long_stage1, aes(x = condition, y = expression, fill = condition)) +
 # geom_quasirandom(color="grey", alpha=0.5)+
  geom_boxplot(outlier.shape = NA, fill = "white") +
      geom_jitter(color="#363636", size=0.5, alpha=0.5) +
  #  scale_fill_manual(values = c("#9c71a6", "#f5ec73"))+
  theme_classic() +
  theme(
    legend.position="none",
    plot.title = element_text(size=10),
    panel.background=element_rect(fill = "white"),
    plot.background=element_rect(fill = "white"),
    text = element_text(family = "Arial", size = 8)
  ) +
  stat_compare_means(
    paired = TRUE,
    method = "wilcox.test",
    label = "p.format"
  )+
  scale_x_discrete(limits = c(column_name_stage1, column_name_stage2))+
  xlab("hist")+
  ggtitle(paste(hist, " Up In ", stage[[1]])) +
  #  ylim(0,2)+
  ylab("log(TPM+1)") +
  annotate(
    "text",
    x = -Inf, y = Inf,
    label = p_txt,
    hjust = -0.05, vjust = 1.05,
    size = 4
  )
  

ggsave(filename = file.path(resDir,paste0(hist, "_", stage[[1]], "vs", stage[[2]], 
                                          "_up_in_", stage[[1]], ".FDR", FDR_th,
                                          ".expression", ".png")),
       plot = g1, width = 12, height = 12, units = "cm")



# plot genes up in stage 2
df_long_stage2 <- tpm_matrix %>% filter(UP_IN_STAGE_2 == 1) %>%
  dplyr::select(.data[[column_name_stage1]], .data[[column_name_stage2]]) %>%
  pivot_longer(
    cols = everything(),
    names_to = "condition",
    values_to = "expression"
  )

p_val_2 <- wilcox.test(
  x = df_long_stage2$expression[df_long_stage2$condition == column_name_stage1],
  y = df_long_stage2$expression[df_long_stage2$condition == column_name_stage2],
  paired = TRUE,
  exact = FALSE
)$p.value
p_txt <- paste0("p (wilcox.test, paired = T) = ", signif(p_val_2, 3))

g2 <- ggplot(df_long_stage2, aes(x = condition, y = expression, fill = condition)) +
  # geom_quasirandom(color="grey", alpha=0.5)+
  geom_boxplot(outlier.shape = NA, fill = "white") +
  geom_jitter(color="#363636", size=0.5, alpha=0.5) +
  #  scale_fill_manual(values = c("#9c71a6", "#f5ec73"))+
  theme_classic() +
  theme(
    legend.position="none",
    plot.title = element_text(size=10),
    panel.background=element_rect(fill = "white"),
    plot.background=element_rect(fill = "white"),
    text = element_text(family = "Arial", size = 8)
  ) +
  xlab("hist")+
  ggtitle(paste(hist, " Up In ", stage[[2]])) +
  scale_x_discrete(limits = c(column_name_stage2, column_name_stage1))+
  #  ylim(0,2)+
  ylab("log(TPM+1)") +
  annotate(
    "text",
    x = -Inf, y = Inf,
    label = p_txt,
    hjust = -0.05, vjust = 1.05,
    size = 4
  )


ggsave(filename = file.path(resDir,paste0(hist, "_", stage[[1]], "vs", stage[[2]], 
                                          "_up_in_", stage[[2]], ".FDR", FDR_th, 
                                          ".expression", ".png")),
       plot = g2, width = 12, height = 12, units = "cm")
}















