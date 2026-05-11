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

comparison <- "WSvsSDL_H3K27ME3_newDataset" # WSvsMG WSvsMG, H2AUB, H3K27ME3, K4ME3

# processing
{
  # make result dirs
  resDir <- paste0("../results/plotExpressionHeatmap_diffPeaks_FDR", FDR_th)
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
}

# plot heatmap
{
heatmap_df <- tpm_matrix %>%
  filter(UP_IN_STAGE_1 == 1 | UP_IN_STAGE_2 == 1) %>%
  dplyr::select(
    gene,
    all_of(c(column_name_stage1, column_name_stage2)),
    UP_IN_STAGE_1,
    UP_IN_STAGE_2
  ) %>%
  pivot_longer(
    cols = all_of(c(column_name_stage1, column_name_stage2)),
    names_to = "condition",
    values_to = "expression"
  ) %>%
  mutate(
    category =  paste0(hist, " up in ", case_when(
      UP_IN_STAGE_1 == 1 ~ stage[[1]],
      UP_IN_STAGE_2 == 1 ~ stage[[2]]
    ) )
  )

gene_order <- heatmap_df %>%
  group_by(gene) %>%
  summarize(mean_expr = mean(expression, na.rm = TRUE)) %>%
  arrange(mean_expr) %>%
  pull(gene)

heatmap_df <- heatmap_df %>%
  mutate(gene = factor(gene, levels = gene_order))

g_heat <- ggplot(
  heatmap_df,
  aes(x = condition, y = gene, fill = expression)
) +
  geom_tile() +
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = median(heatmap_df$expression, na.rm = TRUE),
    name = "log(TPM+1)"
  ) +
  facet_grid(
    rows = vars(category),
    scales = "free_y",
    space = "free_y") +
  theme_void() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(),
    axis.title = element_blank(),
    strip.background = element_rect(fill = "white", color = NA),
    strip.text = element_text(size = 9),
    text = element_text(family = "Arial", size = 8)
  ) +
  theme(
    legend.position = "bottom"
  ) +
  ggtitle(paste(hist, " up genes up in", stage[[1]], "or", stage[[2]], ", Expression, FDR = ", FDR_th))

g_heat



ggsave(filename = file.path(resDir,paste0(hist, "_", stage[[1]], "vs", stage[[2]], 
                                          ".FDR", FDR_th, 
                                          ".expression", ".png")),
       plot = g_heat, width = 12, height = 12, units = "cm")

}


# Heatmap for all three stages

# plot heatmap
{
  heatmap_df <- tpm_matrix %>%
    filter(UP_IN_STAGE_1 == 1 | UP_IN_STAGE_2 == 1) %>%
    dplyr::select(
      gene,
      c(logmean_expression_ws, logmean_expression_mg, logmean_expression_sdl),
      UP_IN_STAGE_1,
      UP_IN_STAGE_2
    ) %>%
    pivot_longer(
      cols = c(logmean_expression_ws, logmean_expression_mg, logmean_expression_sdl),
      names_to = "condition",
      values_to = "expression"
    ) %>%
    mutate(
      category =  paste0(hist, " up in ", case_when(
        UP_IN_STAGE_1 == 1 ~ stage[[1]],
        UP_IN_STAGE_2 == 1 ~ stage[[2]]
      ) )
    )
  
  gene_order <- heatmap_df %>%
    group_by(gene) %>%
    summarize(mean_expr = mean(expression, na.rm = TRUE)) %>%
    arrange(mean_expr) %>%
    pull(gene)
  
  heatmap_df <- heatmap_df %>%
    mutate(
      condition = factor(
        condition,
        levels = c("logmean_expression_ws",
                   "logmean_expression_mg",
                   "logmean_expression_sdl")
      )
    )
  
  
  heatmap_df <- heatmap_df %>%
    mutate(gene = factor(gene, levels = gene_order))
  
  g_heat <- ggplot(
    heatmap_df,
    aes(x = condition, y = gene, fill = expression)
  ) +
    geom_tile() +
    scale_fill_gradient2(
      low = "blue",
      mid = "white",
      high = "red",
      midpoint = median(heatmap_df$expression, na.rm = TRUE),
      name = "log(TPM+1)"
    ) +
    facet_grid(
      rows = vars(category),
      scales = "free_y",
      space = "free_y") +
    theme_void() +
    theme(
      axis.text.y = element_text(),
      axis.ticks.y = element_blank(),
      axis.text.x = element_text(),
      axis.title = element_blank(),
      strip.background = element_rect(fill = "white", color = NA),
      strip.text = element_text(size = 9),
      text = element_text(family = "Arial", size = 8)
    ) +
    theme(
      legend.position = "bottom"
    ) +
    ggtitle(paste(hist, " up genes up in", stage[[1]], "or", stage[[2]], ", Expression, FDR = ", FDR_th))
  
  g_heat
  
  
  
  ggsave(filename = file.path(resDir,paste0(hist, "_", stage[[1]], "vs", stage[[2]], 
                                            ".FDR", FDR_th, 
                                            ".expression.allStages", ".png")),
         plot = g_heat, width = 12, height = 12, units = "cm")
  
}


{
# plot a box plot instead
  g <- ggplot(
    heatmap_df %>% filter(UP_IN_STAGE_1 == 1),
    aes(x = condition, y = expression)
  ) +
    geom_quasirandom(color="grey", alpha=0.5)+
    geom_boxplot(outlier.shape = NA, alpha = 0.3) +
    #    geom_jitter(color="#363636", size=0.05, alpha=0.3) +
    #  scale_fill_manual(values = c("#9c71a6", "#f5ec73"))+
    theme_minimal() +
    theme(
      legend.position="none",
      plot.title = element_text(size=10),
      panel.background=element_rect(fill = "white"),
      plot.background=element_rect(fill = "white"),
      text = element_text(family = "Arial", size = 8)
    ) 
  
  my_comparisons <- list(
    c("logmean_expression_ws", "logmean_expression_mg"),
    c("logmean_expression_mg", "logmean_expression_sdl")
  )
  
  g <- g + stat_compare_means(
    comparisons = my_comparisons,
    method = "wilcox.test",
    label = "p.format"
  )
  
  g
}

{
  peakfile <- file.path("../results/ChIPseeker_consensusPeaks", 
                        "SDL_H3K27ME3_newDataset.consensusPeaks.withAnnotation.tsv")
  peaks_df <- read.table(file = peakfile, sep = "\t", quote = "")
  allGeneList <- peaks_df$ASSOCIATED_GENES %>% strsplit(", ") %>% unlist
  allGeneList <- allGeneList[!is.na(allGeneList)] %>% unique()
  
  geneList_inquiry <- (heatmap_df %>% filter(category == "H2AUB up in MG"))$gene %>% unique
  (geneList_inquiry %in% allGeneList) %>% sum # 27
  
  peakfile <- file.path("../results/ChIPseeker_consensusPeaks", 
                        "MG_K27ME3.consensusPeaks.withAnnotation.tsv")
  peaks_df <- read.table(file = peakfile, sep = "\t", quote = "")
  allGeneList <- peaks_df$ASSOCIATED_GENES %>% strsplit(", ") %>% unlist
  allGeneList <- allGeneList[!is.na(allGeneList)] %>% unique()
  (geneList_inquiry %in% allGeneList) %>% sum # 24
  
  deseq_res_2 <- file.path(paste0("../results/diffBind_ChIPseeker_", "MGvsSDL_H3K27ME3"), 
                "DESeq2Res_withAnnotations.tsv") %>%
  read.table(header = TRUE, fill = TRUE, sep = "\t", quote = "")
  
  deseq_res_2_filterred <- separate_longer_delim(deseq_res_2, ASSOCIATED_GENES, ", ") %>% 
    filter (ASSOCIATED_GENES %in% geneList_inquiry) %>% group_by(across(1:3))
  
  box_df <- deseq_res_2_filterred %>%
    pivot_longer(
      cols = c(Conc_MG, Conc_SDL),
      names_to = "condition",
      values_to = "value"
    )
  
  b <- ggplot(box_df, aes(x = condition, y = value)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.4) +
    geom_jitter(width = 0.15, alpha = 0.4, size = 0.6) +
    theme_classic() +
    labs(
      x = NULL,
      y = "Conc H3K27ME3"
    )
  
  b <- b + stat_compare_means(
    method = "wilcox.test",
    label = "p.value"
  )
  
  b
}
















