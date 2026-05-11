# Plot the expression of H2Aub peaks overlapped with and without H3K27me3

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

resDir <- "../results/plotExpression_H2AubIntersectH3K27me3Peaks"
if (!file.exists(resDir)){dir.create(resDir)}

geneList_all <- keys(org.At.tair.db, keytype = "TAIR") # all genes in database

# read in the annotation from ChIPseeker

samplefiles <- list.files("../results/ChIPseeker_H2AubIntersectH3K27me3Peaks", pattern= ".tsv", full.names=T)
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
                       MG_H2AUB_OVERLAPWITH_H3K27ME3=ifelse(gene %in% geneList$MG_H2AUB_OVERLAPWITH_H3K27ME3, 1, 0),
                       MG_H2AUB_WITHOUT_H3K27ME3=ifelse(gene %in% geneList$MG_H2AUB_WITHOUT_H3K27ME3, 1, 0),
                       SDL_H2AUB_OVERLAPWITH_H3K27ME3=ifelse(gene %in% geneList$SDL_H2AUB_OVERLAPWITH_H3K27ME3, 1, 0),
                       SDL_H2AUB_WITHOUT_H3K27ME3=ifelse(gene %in% geneList$SDL_H2AUB_WITHOUT_H3K27ME3, 1, 0),
                       WS_H2AUB_OVERLAPWITH_H3K27ME3=ifelse(gene %in% geneList$WS_H2AUB_OVERLAPWITH_H3K27ME3, 1, 0),
                       WS_H2AUB_WITHOUT_H3K27ME3=ifelse(gene %in% geneList$WS_H2AUB_WITHOUT_H3K27ME3, 1, 0))
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
peakset_to_plot <- "MG_H2AUB_OVERLAPWITH_H3K27ME3"

templist <- c(rep("logmean_expression_mg", 3), rep("logmean_expression_sdl", 4), 
              rep("logmean_expression_ws", 3))

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

# mg
{
  EXP_WITH_H3K27ME3 = tpm_matrix$logmean_expression_mg[tpm_matrix$MG_H2AUB_OVERLAPWITH_H3K27ME3 == 1]
  EXP_WITHOUT_H3K27ME3 = tpm_matrix$logmean_expression_mg[tpm_matrix$MG_H2AUB_WITHOUT_H3K27ME3 == 1]
  EXP_ALL_GENES = tpm_matrix$logmean_expression_mg
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
    ggtitle("H2Aub Marked Genes' Expression in MG, With or Without H3K27me3") +
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
  ggsave(filename = file.path(boxplotDir,paste0("MG", "_H2AubGenes_withAndWithoutH3K27me3", ".png")),
         plot = g, width = 12, height = 12, units = "cm")
}

#sdl
{
  EXP_WITH_H3K27ME3 = tpm_matrix$logmean_expression_sdl[tpm_matrix$SDL_H2AUB_OVERLAPWITH_H3K27ME3 == 1]
  EXP_WITHOUT_H3K27ME3 = tpm_matrix$logmean_expression_sdl[tpm_matrix$SDL_H2AUB_WITHOUT_H3K27ME3 == 1]
  EXP_ALL_GENES = tpm_matrix$logmean_expression_sdl
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
    ggtitle("H2Aub Marked Genes' Expression in SDL, With or Without H3K27me3") +
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
  ggsave(filename = file.path(boxplotDir,paste0("SDL", "_H2AubGenes_withAndWithoutH3K27me3", ".png")),
         plot = g, width = 12, height = 12, units = "cm")
}


# plot heatmap for all three stages
# mg
{
  heatmap_df <- tpm_matrix %>%
    filter(MG_H2AUB_OVERLAPWITH_H3K27ME3 == 1 | MG_H2AUB_WITHOUT_H3K27ME3 == 1) %>%
    dplyr::select(
      gene,
      c(logmean_expression_ws, logmean_expression_mg, logmean_expression_sdl),
      MG_H2AUB_OVERLAPWITH_H3K27ME3,
      MG_H2AUB_WITHOUT_H3K27ME3
    ) %>%
    pivot_longer(
      cols = c(logmean_expression_ws, logmean_expression_mg, logmean_expression_sdl),
      names_to = "condition",
      values_to = "expression"
    ) %>%
    mutate(
      category =  paste0(case_when(
        MG_H2AUB_OVERLAPWITH_H3K27ME3 == 1 ~ "With H3K27me3",
        MG_H2AUB_WITHOUT_H3K27ME3 == 1 ~ "Without H3K27me3"
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
    )
  
  g_heat
  
}

boxplot_df <- heatmap_df %>%
  filter(category == "With H3K27me3")

boxplot_df <- boxplot_df %>%
  mutate(
    condition = factor(
      condition,
      levels = c("logmean_expression_ws",
                 "logmean_expression_mg",
                 "logmean_expression_sdl")
    )
  )

g_box <- ggplot(
  boxplot_df,
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

g_box

{
  # mg
  mg_all_H2Aub_genes <- c(geneList$MG_H2AUB_OVERLAPWITH_H3K27ME3, geneList$MG_H2AUB_WITHOUT_H3K27ME3) %>% unique
  ws_all_H2Aub_genes <- c(geneList$WS_H2AUB_OVERLAPWITH_H3K27ME3, geneList$WS_H2AUB_WITHOUT_H3K27ME3) %>% unique
  
  
  # peaks with increased H3K27me3 from em to sdl
  deseq_res_2 <- file.path(paste0("../results/diffBind_ChIPseeker_", "MGvsSDL_H3K27ME3_newDataset"), 
                           "DESeq2Res_withAnnotations.tsv") %>%
    read.table(header = TRUE, fill = TRUE, sep = "\t", quote = "")
  stage <- c()
  stage[1] <- names(deseq_res_2)[7] %>% str_remove("Conc_")
  stage[2] <- names(deseq_res_2)[8] %>% str_remove("Conc_")
  FDR_th = 0.05
  peakList <- deseq_res_2 %>% filter(FDR < FDR_th) %>%
    filter(.data[[paste0("Conc_", stage[2])]] > .data[[paste0("Conc_", stage[1])]])
  
  geneList_fromDe <- peakList$ASSOCIATED_GENES %>% str_split(", ") %>% unlist %>% unique
  geneList_fromDe <- geneList_fromDe[!is.na(geneList_fromDe)]
  
  resDir
  color_platte <- c("#CC79A7", "#56B4E9", "#F0E442")
  
  venn.diagram(list(mg_all_H2Aub_genes, geneList_fromDe), 
               filename = file.path(resDir, "testPlot.png"),
               category = c("MG H2AUb Genes", "SDL H3K27me3 Up Genes"),
               disable.logging = T,
               fill = c(color_platte[1], color_platte[2]),
               fontfamily = "Arial",
               cat.fontfamily = "Arial",
               hyper.test = T, total.population = length(geneList_all), cat.dist = c(0.06, 0.06), margin = 0.1 )
  
  venn.diagram(list(ws_all_H2Aub_genes, geneList_fromDe), 
               filename = file.path(resDir, "wsH2AubGenes_intersect_sdlH2K27me3UpGenes.png"),
               category = c("WS H2AUb Genes", "SDL H3K27me3 Up Genes"),
               disable.logging = T,
               fill = c(color_platte[1], color_platte[2]),
               fontfamily = "Arial",
               cat.fontfamily = "Arial",
               hyper.test = T, total.population = length(geneList_all), cat.dist = c(0.06, 0.06), margin = 0.1 )
}

venn.diagram(list(ws_all_H2Aub_genes, mg_all_H2Aub_genes), 
             filename = file.path(resDir, "wsH2AubGenes_intersect_mgH2AubGenes.png"),
             category = c("WS H2AUb Genes", "MG H2AUb Genes"),
             disable.logging = T,
             fill = c(color_platte[1], color_platte[2]),
             fontfamily = "Arial",
             cat.fontfamily = "Arial",
             hyper.test = T, total.population = length(geneList_all), cat.dist = c(0.06, 0.06), margin = 0.1 )




venn.diagram(list(sample(geneList_all, length(geneList_fromDe)), sample(geneList_all, length(mg_all_H2Aub_genes))), 
             filename = file.path(resDir, "simulation.png"),
             category = c("MG H2AUb Genes", "SDL H3K27me3 Up Genes"),
             disable.logging = T,
             fill = c(color_platte[1], color_platte[2]),
             fontfamily = "Arial",
             cat.fontfamily = "Arial",
             hyper.test = T, total.population = length(geneList_all), cat.dist = c(0.06, 0.06), margin = 0.1 )


suppressMessages({
  library(clusterProfiler)
  library(org.At.tair.db)
})


peakList <- deseq_res_2 %>% filter(FDR < FDR_th) %>%
  filter(.data[[paste0("Conc_", "SDL")]] > .data[[paste0("Conc_", "MG")]])

geneList_sdlK27me3up <- peakList$ASSOCIATED_GENES %>% str_split(", ") %>% unlist %>% unique
geneList_sdlK27me3up <- geneList_sdlK27me3up[!is.na(geneList_sdlK27me3up)]

all_H3K27me3_genes <- deseq_res_2$ASSOCIATED_GENES %>% str_split(", ") %>% unlist %>% unique

# against all sdl H3K27me3 up genes
egos <- enrichGO(gene = geneList_sdlK27me3up[(geneList_sdlK27me3up %in% mg_all_H2Aub_genes)], 
                 keyType = "TAIR",
                 OrgDb = org.At.tair.db, 
                 ont = "BP", 
                 pAdjustMethod = "BH",
                 qvalueCutoff = 0.05, 
                 readable = TRUE,
                 pool = TRUE,
                 universe = geneList_sdlK27me3up)
dotplot(egos)
data.frame(egos) %>% view

# against all H2Aub genes
egos <- enrichGO(gene = mg_all_H2Aub_genes[(mg_all_H2Aub_genes %in% geneList_sdlK27me3up)], 
                 keyType = "TAIR",
                 OrgDb = org.At.tair.db, 
                 ont = "BP", 
                 pAdjustMethod = "BH",
                 qvalueCutoff = 0.05, 
                 readable = TRUE,
                 pool = TRUE,
                 universe = mg_all_H2Aub_genes)
dotplot(egos)
data.frame(egos) %>% view


# 
deseq_res_ws_mg_K27me3 <- file.path(paste0("../results/diffBind_ChIPseeker_", "WSvsMG_H3K27ME3"), 
                         "DESeq2Res_withAnnotations.tsv") %>%
  read.table(header = TRUE, fill = TRUE, sep = "\t", quote = "")

FDR_th = 0.05
peakList <- deseq_res_ws_mg_K27me3 %>% filter(FDR < FDR_th) %>%
  filter(.data[[paste0("Conc_", "MG")]] < .data[[paste0("Conc_", "WS")]])

geneList_wsK27me3up <- peakList$ASSOCIATED_GENES %>% str_split(", ") %>% unlist %>% unique
geneList_wsK27me3up <- geneList_wsK27me3up[!is.na(geneList_wsK27me3up)]

egos <- enrichGO(gene = geneList_wsK27me3up[(geneList_wsK27me3up %in% mg_all_H2Aub_genes)], 
                 keyType = "TAIR",
                 OrgDb = org.At.tair.db, 
                 ont = "BP", 
                 pAdjustMethod = "BH",
                 qvalueCutoff = 0.05, 
                 readable = TRUE,
                 pool = TRUE,
                 universe = geneList_wsK27me3up)
dotplot(egos)


# 
deseq_res_k27_mgH2aubGenes <- deseq_res_2 %>% separate_longer_delim(cols = ASSOCIATED_GENES, delim = ", ") %>%
  filter(ASSOCIATED_GENES %in% mg_all_H2Aub_genes)

df_longer <- deseq_res_k27_mgH2aubGenes %>%
  pivot_longer(
    cols = c(Conc_MG, Conc_SDL),
    names_to = "condition",
    values_to = "Conc_H3K27me3"
  )

ggplot(df_longer, aes(x = condition, y = GENE_NAMES, fill = Conc_H3K27me3)) +
  geom_tile() +
  scale_fill_viridis_c(option = "magma") +
  theme_minimal() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )
































