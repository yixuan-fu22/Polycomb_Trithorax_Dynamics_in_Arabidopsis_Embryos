# Plot the expression of H2Aub genes overlapped and not-overlapped between WS and MG embryos

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
  library(org.At.tair.db)
})

resDir <- "../results/ChIPseeker_substractedPeaks"
if (!file.exists(resDir)){dir.create(resDir)}

geneList_all <- keys(org.At.tair.db, keytype = "TAIR") # all genes in database

# read in the annotation from ChIPseeker

samplefiles <- list.files("../results/ChIPseeker_substractedPeaks", pattern= ".tsv", full.names=T)
samplenames <- str_remove(basename(samplefiles),".consensusPeaks.withAnnotation.tsv")
names(samplefiles) <- samplenames

peak_annotation_df <- vector(mode = "list", length = length(samplefiles))
for(i in 1:length(samplefiles)){
  peak_annotation_df[[i]] <- read.table(samplefiles[[i]], header = TRUE, sep = "\t", row.names = 1) %>% as.data.frame()
}
names(peak_annotation_df) <- samplenames

# read in the transcriptome data
{
  tpmfile <- "../data_fromOtherProjects/tpm_fromRNAseqOfEmbryosSeedlingsAndFlowers/embryos_sdl_flw.deseq2normalized.tpm.csv"
  
  readtpmfile <- read.table(tpmfile, sep = ",") %>% t()
  
  tpm_matrix <- readtpmfile[-1,]
  colnames(tpm_matrix) <- readtpmfile[1,]
  colnames(tpm_matrix)[1] <- "gene"
  tpm_matrix <- tpm_matrix %>% as.data.frame
  
  # read in gene list
  geneList <- vector(mode = "list", length = length(samplenames))
  for(i in 1:length(samplenames)){
    geneList_temp <- peak_annotation_df[[i]]$ASSOCIATED_GENES %>% str_split(", ") %>% unlist %>% unique
    geneList[[i]] <- geneList_temp[!is.na(geneList_temp)]
  }
  names(geneList) <- samplenames
  
  # calculate log mean expression at each stage
  
  # Late torpedo stage from Hoffmann et al 2017: SRR8054371, SRR8054372, SRR8054373
  # MG stage: SRR8054377, SRR8054378, SRR8054379
  # SDL1: ERR10163199, ERR10163200
  # SDL2: SRR2500947, SRR1931614
  
  # 
  tpm_matrix <- mutate(tpm_matrix, 
                       H2AUB_MG_SUBSTRACT_WS=ifelse(gene %in% geneList$H2AUB_MG_substract_WS, 1, 0),
                       H2AUB_WS_SUBSTRACT_MG=ifelse(gene %in% geneList$H2AUB_WS_substract_MG, 1, 0))
  # 1=marked, 0=not marked
  
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

# heatplot
{
  heatmap_df <- tpm_matrix %>%
    filter(H2AUB_MG_SUBSTRACT_WS == 1 | H2AUB_WS_SUBSTRACT_MG == 1) %>%
    dplyr::select(
      gene,
      c(logmean_expression_ws, logmean_expression_mg, logmean_expression_sdl),
      H2AUB_MG_SUBSTRACT_WS,
      H2AUB_WS_SUBSTRACT_MG
    ) %>%
    pivot_longer(
      cols = c(logmean_expression_ws, logmean_expression_mg),
      names_to = "condition",
      values_to = "expression"
    ) %>%
    mutate(
      category =  paste0(case_when(
        H2AUB_MG_SUBSTRACT_WS == 1 ~ "Peaks Only In MG",
        H2AUB_WS_SUBSTRACT_MG == 1 ~ "Peaks Only in WS"
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
    )
  
  g_heat

  ggsave(filename = file.path(resDir,paste0("heatplot", "PeaksOnlyInMGorWS",
                                            ".expression", ".png")),
         plot = g_heat, width = 12, height = 12, units = "cm")
}

# boxplot
{
  g <- ggplot(
    heatmap_df %>% filter(H2AUB_MG_SUBSTRACT_WS == 1),
    aes(x = condition, y = expression, fill = condition)
  ) +
    geom_quasirandom(color="grey", alpha=0.5)+
    geom_boxplot(outlier.shape = NA, alpha = 0.3) +
    #    geom_jitter(color="#363636", size=0.05, alpha=0.3) +
    #  scale_fill_manual(values = c("#9c71a6", "#f5ec73"))+
    theme_minimal() +
    ggtitle("H2Aub Genes In MG Only")+
    theme(
      legend.position="none",
      plot.title = element_text(size=10),
      text = element_text(family = "Arial", size = 8)
    )
  
  my_comparisons <- list(
    c("logmean_expression_ws", "logmean_expression_mg")
  )
  
  g <- g + stat_compare_means(
    comparisons = my_comparisons,
    method = "wilcox.test",
    label = "p.format"
  )
  
  g
  
  ggsave(filename = file.path(resDir,paste0("H2AubGenesInMGonly",
                                            ".expression", ".png")),
         plot = g, width = 12, height = 12, units = "cm")

  }

# boxplot another direction
{
  g <- ggplot(
    heatmap_df %>% filter(H2AUB_WS_SUBSTRACT_MG == 1),
    aes(x = condition, y = expression, fill = condition)
  ) +
    geom_quasirandom(color="grey", alpha=0.5)+
    geom_boxplot(outlier.shape = NA, alpha = 0.3) +
    #    geom_jitter(color="#363636", size=0.05, alpha=0.3) +
    #  scale_fill_manual(values = c("#9c71a6", "#f5ec73"))+
    theme_minimal() +
    theme(
      legend.position="none",
      plot.title = element_text(size=10),
      text = element_text(family = "Arial", size = 8)
    )
  
  my_comparisons <- list(
    c("logmean_expression_ws", "logmean_expression_mg")
  )
  
  g <- g + stat_compare_means(
    comparisons = my_comparisons,
    method = "wilcox.test",
    label = "p.format"
  )
  
  g
  
  ggsave(filename = file.path(resDir,paste0("H2AubGenesInWSonly",
                                            ".expression", ".png")),
         plot = g, width = 12, height = 12, units = "cm")
  
}





















































