# Plot the expression of H3K27me3 peaks overlapped with H3K4me3

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


resDir <- "../results/plotExpression_H3K27me3IntersectH3K4me3Peaks"
if (!file.exists(resDir)){dir.create(resDir)}

geneList_all <- keys(org.At.tair.db, keytype = "TAIR") # all genes in database

# read in the annotation from ChIPseeker

samplefiles <- list.files("../results/ChIPseeker_H3K27me3_intersect_H3K4me3_peaks", pattern= ".tsv", full.names=T)
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
  view(tpm_matrix)
  
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
                       WS_H3K27ME3_OVERLAPWITH_H3K4ME3=ifelse(gene %in% geneList$WS_H3K27ME3_OVERLAPWITH_H3K4ME3, 1, 0),
                       MG_H3K27ME3_OVERLAPWITH_H3K4ME3=ifelse(gene %in% geneList$MG_H3K27ME3_OVERLAPWITH_H3K4ME3, 1, 0),
                       SDL_H3K27ME3_OVERLAPWITH_H3K4ME3=ifelse(gene %in% geneList$SDL_H3K27ME3_OVERLAPWITH_H3K4ME3, 1, 0))
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

# boxplot with scatters
boxplotDir <- resDir
if (!file.exists(boxplotDir)){dir.create(boxplotDir)}

expression_to_plot = "logmean_expression_ws"
peakset_to_plot <- "WS_H3K27ME3_OVERLAPWITH_H3K4ME3"

templist <- c(rep("logmean_expression_mg", 1), rep("logmean_expression_sdl", 1), 
              rep("logmean_expression_ws", 1))

for (i in 1:length(samplenames)) {
  expression_to_plot = templist[i]
  peakset_to_plot = samplenames[i]
  
  tpm_matrix <- tpm_matrix %>% mutate(expression = .data[[expression_to_plot]])
  
  p <- tpm_matrix %>% filter(expression < quantile(tpm_matrix$expression, probs = 0.99)) %>%
    ggplot(aes(x = as.factor(.data[[peakset_to_plot]]), y = expression, fill = as.factor(.data[[peakset_to_plot]]))) +
    geom_quasirandom(color="grey", alpha=0.5)+
    geom_boxplot(outlier.shape = NA) +
    #    geom_jitter(color="#363636", size=0.05, alpha=0.3) +
    #  scale_fill_manual(values = c("#9c71a6", "#f5ec73"))+
    theme_minimal() +
    theme(
      legend.position="none",
      plot.title = element_text(size=10),
      panel.background=element_rect(fill = "white"),
      plot.background=element_rect(fill = "white"),
      text = element_text(family = "Arial", size = 8)
    ) +
    ggtitle(peakset_to_plot) +
    xlab("")+
    #  ylim(0,2)+
    ylab("log(TPM+1)")
  
  my_comparisons <- list(c("0", "1"))
  p <- p+
    stat_compare_means(comparisons = my_comparisons,label = "p.format", method = "wilcox.test")
  
  p
  
  ggsave(filename = file.path(boxplotDir,paste0(peakset_to_plot, "_vs_", expression_to_plot, ".png")),
         plot = p, width = 8, height = 8, units = "cm")
}

# Plot to compare at each stage
# ws
{
  EXP_WITH_H3K27ME3 = tpm_matrix$logmean_expression_ws[tpm_matrix$WS_H2AUB_OVERLAPWITH_H3K27ME3 == 1]
  EXP_WITHOUT_H3K27ME3 = tpm_matrix$logmean_expression_ws[tpm_matrix$WS_H2AUB_WITHOUT_H3K27ME3 == 1]
  EXP_ALL_GENES = tpm_matrix$logmean_expression_ws
  df_exp <- bind_rows(
    data.frame(group = "WITH_H3K27ME3",  expression = EXP_WITH_H3K27ME3),
    data.frame(group = "WITHOUT_H3K27ME3",  expression = EXP_WITHOUT_H3K27ME3),
    data.frame(group = "ALL_GENES", expression = EXP_ALL_GENES)
  )
  df_exp$group <- factor(
    df_exp$group,
    levels = c("WITH_H3K27ME3", "WITHOUT_H3K27ME3", "ALL_GENES")
  )
  g <- ggplot(df_exp, aes(x = group, y = expression, fill = group)) +
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
    ) +
    xlab("hist")+
    ggtitle("H2Aub Marked Genes' Expression in WS, With or Without H3K27me3") +
    #  ylim(0,2)+
    ylab("log(TPM+1)")
  my_comparisons <- list(
    c("WITH_H3K27ME3", "ALL_GENES"),
    c("WITHOUT_H3K27ME3", "ALL_GENES"),
    c("WITH_H3K27ME3", "WITHOUT_H3K27ME3")
  )
  # to get p value manually # example
  g <- g +
    stat_compare_means(comparisons = my_comparisons,label = "p", method = "wilcox.test")
  g
  ggsave(filename = file.path(boxplotDir,paste0("WS", "_H2AubGenes_withAndWithoutH3K27me3", ".png")),
         plot = g, width = 12, height = 12, units = "cm")
}


# plot heatmap for all three stages
# ws
{
  id="WS_H3K27ME3_OVERLAPWITH_H3K4ME3"
  heatmap_df <- tpm_matrix %>%
    filter(WS_H3K27ME3_OVERLAPWITH_H3K4ME3 == 1) %>%
    dplyr::select(
      gene,
      c(logmean_expression_ws, logmean_expression_mg, logmean_expression_sdl),
      WS_H3K27ME3_OVERLAPWITH_H3K4ME3
    ) %>%
    pivot_longer(
      cols = c(logmean_expression_ws, logmean_expression_mg, logmean_expression_sdl),
      names_to = "condition",
      values_to = "expression"
    ) %>%
    mutate(
      category =  paste0(case_when(
        WS_H3K27ME3_OVERLAPWITH_H3K4ME3 == 1 ~ "WTIH H3K27ME3 and H3K4ME3 at WS",
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
  
  ggsave(filename = file.path(resDir,paste0(id, 
                                            ".expression",".allStages", ".png")),
         plot = g_heat, width = 12, height = 12, units = "cm")
}

# mg
{
  id="MG_H3K27ME3_OVERLAPWITH_H3K4ME3"
  heatmap_df <- tpm_matrix %>%
    filter(MG_H3K27ME3_OVERLAPWITH_H3K4ME3 == 1) %>%
    dplyr::select(
      gene,
      c(logmean_expression_ws, logmean_expression_mg, logmean_expression_sdl),
      MG_H3K27ME3_OVERLAPWITH_H3K4ME3
    ) %>%
    pivot_longer(
      cols = c(logmean_expression_ws, logmean_expression_mg, logmean_expression_sdl),
      names_to = "condition",
      values_to = "expression"
    ) %>%
    mutate(
      category =  paste0(case_when(
        MG_H3K27ME3_OVERLAPWITH_H3K4ME3 == 1 ~ "WTIH H3K27ME3 and H3K4ME3 at MG",
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
  
  ggsave(filename = file.path(resDir,paste0(id, 
                                            ".expression",".allStages", ".png")),
         plot = g_heat, width = 12, height = 12, units = "cm")

}



# sdl
{
  id="SDL_H3K27ME3_OVERLAPWITH_H3K4ME3"
  heatmap_df <- tpm_matrix %>%
    filter(SDL_H3K27ME3_OVERLAPWITH_H3K4ME3 == 1) %>%
    dplyr::select(
      gene,
      c(logmean_expression_ws, logmean_expression_mg, logmean_expression_sdl),
      SDL_H3K27ME3_OVERLAPWITH_H3K4ME3
    ) %>%
    pivot_longer(
      cols = c(logmean_expression_ws, logmean_expression_mg, logmean_expression_sdl),
      names_to = "condition",
      values_to = "expression"
    ) %>%
    mutate(
      category =  paste0(case_when(
        SDL_H3K27ME3_OVERLAPWITH_H3K4ME3 == 1 ~ "WTIH H3K27ME3 and H3K4ME3 at SDL",
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
  
  ggsave(filename = file.path(resDir,paste0(id, 
                                            ".expression",".allStages", ".png")),
         plot = g_heat, width = 12, height = 12, units = "cm")
  
}














