# plot the expression of all marked genes of a histone modification across different stages

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

resDir <- "../results/plotExpression_allPeaks"
if (!file.exists(resDir)){dir.create(resDir)}

geneList_all <- keys(org.At.tair.db, keytype = "TAIR") # all genes in database

# read in the annotation from ChIPseeker

samplefiles <- list.files("../results/ChIPseeker_consensusPeaks", pattern= ".tsv", full.names=T)
samplenames <- str_remove(basename(samplefiles),".consensusPeaks.withAnnotation.tsv")
names(samplefiles) <- samplenames

peak_annotation_df <- vector(mode = "list", length = length(samplefiles))
for(i in 1:length(samplefiles)){
  peak_annotation_df[[i]] <- read.table(samplefiles[[i]], header = TRUE, sep = "\t", row.names = 1) %>% as.data.frame()
}
names(peak_annotation_df) <- samplenames

###
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

# K27Ac
tpm_matrix <- mutate(tpm_matrix, 
                     MG_H2AUB=ifelse(gene %in% geneList$MG_H2AUB, 1, 0),
                     MG_K27ME3=ifelse(gene %in% geneList$MG_K27ME3, 1, 0),
                     MG_K4ME3=ifelse(gene %in% geneList$MG_K4ME3, 1, 0),
                     SDL_H2AUB_REPMERGED=ifelse(gene %in% geneList$SDL_H2AUB_REPMERGED, 1, 0),
                     SDL_H2AUB=ifelse(gene %in% geneList$SDL_H2AUB, 1, 0),
                     SDL_K4ME3=ifelse(gene %in% geneList$SDL_K4ME3, 1, 0),
                     SDL14D_K27ME3=ifelse(gene %in% geneList$SDL_H3K27ME3_newDataset, 1, 0),
                     WS_H2AUB=ifelse(gene %in% geneList$WS_H2AUB, 1, 0),
                     WS_K27ME3=ifelse(gene %in% geneList$WS_K27ME3, 1, 0),
                     WS_K4ME3=ifelse(gene %in% geneList$WS_K4ME3, 1, 0)) # 1=marked, 0=not marked

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
peakset_to_plot <- "WS_H2AUB"





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

# Plot to compare three histone modification at each stage
# ws
EXP_WS_H2AUB = tpm_matrix$logmean_expression_ws[tpm_matrix$WS_H2AUB == 1]
EXP_WS_K4ME3 = tpm_matrix$logmean_expression_ws[tpm_matrix$WS_K4ME3 == 1]
EXP_WS_K27ME3 = tpm_matrix$logmean_expression_ws[tpm_matrix$WS_K27ME3 == 1]
EXP_ALL_GENES = tpm_matrix$logmean_expression_ws


df_exp <- bind_rows(
  data.frame(group = "WS_K4ME3",  expression = EXP_WS_K4ME3),
  data.frame(group = "WS_H2AUB",  expression = EXP_WS_H2AUB),
  data.frame(group = "WS_K27ME3", expression = EXP_WS_K27ME3),
  data.frame(group = "ALL_GENES", expression = EXP_ALL_GENES)
)

df_exp$group <- factor(
  df_exp$group,
  levels = c("WS_K4ME3", "WS_H2AUB", "WS_K27ME3", "ALL_GENES")
)

g <- ggplot(df_exp, aes(x = group, y = expression, fill = group)) +
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
  xlab("hist")+
  ggtitle("Expression in WS")
  #  ylim(0,2)+
  ylab("log(TPM+1)")

my_comparisons <- list(
  c("WS_K27ME3", "ALL_GENES"),
  c("WS_H2AUB", "ALL_GENES"),
  c("WS_K4ME3", "ALL_GENES"),
  c("WS_H2AUB", "WS_K4ME3"),
  c("WS_H2AUB", "WS_K27ME3"),
  c("WS_K4ME3", "WS_K27ME3")
)

# to get p value manually # example
wilcox.test(EXP_WS_K4ME3, EXP_ALL_GENES, exact = TRUE)$p.value %>% log10()

g <- g +
  stat_compare_means(comparisons = my_comparisons,label = "p", method = "wilcox.test")

ggsave(filename = file.path(boxplotDir,paste0("WS", "_acrossHist", ".png")),
       plot = g, width = 12, height = 12, units = "cm")


#
# mg
EXP_MG_H2AUB = tpm_matrix$logmean_expression_mg[tpm_matrix$MG_H2AUB == 1]
EXP_MG_K4ME3 = tpm_matrix$logmean_expression_mg[tpm_matrix$MG_K4ME3 == 1]
EXP_MG_K27ME3 = tpm_matrix$logmean_expression_mg[tpm_matrix$MG_K27ME3 == 1]
EXP_ALL_GENES = tpm_matrix$logmean_expression_mg


df_exp <- bind_rows(
  data.frame(group = "MG_K4ME3",  expression = EXP_MG_K4ME3),
  data.frame(group = "MG_H2AUB",  expression = EXP_MG_H2AUB),
  data.frame(group = "MG_K27ME3", expression = EXP_MG_K27ME3),
  data.frame(group = "ALL_GENES", expression = EXP_ALL_GENES)
)

df_exp$group <- factor(
  df_exp$group,
  levels = c("MG_K4ME3", "MG_H2AUB", "MG_K27ME3", "ALL_GENES")
)

g <- ggplot(df_exp, aes(x = group, y = expression, fill = group)) +
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
  xlab("hist")+
  ggtitle("Expression in MG") +
#  ylim(0,2)+
ylab("log(TPM+1)")

my_comparisons <- list(
  c("MG_K27ME3", "ALL_GENES"),
  c("MG_H2AUB", "ALL_GENES"),
  c("MG_K4ME3", "ALL_GENES"),
  c("MG_H2AUB", "MG_K4ME3"),
  c("MG_H2AUB", "MG_K27ME3"),
  c("MG_K4ME3", "MG_K27ME3")
)

# to get p value manually # example
wilcox.test(EXP_MG_K4ME3, EXP_ALL_GENES, exact = TRUE)$p.value %>% log10()

g <- g +
  stat_compare_means(comparisons = my_comparisons,label = "p", method = "wilcox.test")

ggsave(filename = file.path(boxplotDir,paste0("MG", "_acrossHist", ".png")),
       plot = g, width = 12, height = 12, units = "cm")

#

# SDL
EXP_SDL_H2AUB = tpm_matrix$logmean_expression_sdl[tpm_matrix$SDL_H2AUB_REPMERGED == 1]
EXP_SDL_K4ME3 = tpm_matrix$logmean_expression_sdl[tpm_matrix$SDL_K4ME3 == 1]
EXP_SDL_K27ME3 = tpm_matrix$logmean_expression_sdl[tpm_matrix$SDL14D_K27ME3 == 1]
EXP_ALL_GENES = tpm_matrix$logmean_expression_sdl


df_exp <- bind_rows(
  data.frame(group = "SDL_K4ME3",  expression = EXP_SDL_K4ME3),
  data.frame(group = "SDL_H2AUB",  expression = EXP_SDL_H2AUB),
  data.frame(group = "SDL_K27ME3", expression = EXP_SDL_K27ME3),
  data.frame(group = "ALL_GENES", expression = EXP_ALL_GENES)
)

df_exp$group <- factor(
  df_exp$group,
  levels = c("SDL_K4ME3", "SDL_H2AUB", "SDL_K27ME3", "ALL_GENES")
)

g <- ggplot(df_exp, aes(x = group, y = expression, fill = group)) +
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
  xlab("hist")+
  ggtitle("Expression in SDL") +
#  ylim(0,2)+
ylab("log(TPM+1)")

my_comparisons <- list(
  c("SDL_K27ME3", "ALL_GENES"),
  c("SDL_H2AUB", "ALL_GENES"),
  c("SDL_K4ME3", "ALL_GENES"),
  c("SDL_H2AUB", "SDL_K4ME3"),
  c("SDL_H2AUB", "SDL_K27ME3"),
  c("SDL_K4ME3", "SDL_K27ME3")
)

# to get p value manually # example
wilcox.test(EXP_SDL_K4ME3, EXP_ALL_GENES, exact = TRUE)$p.value %>% log10()

g <- g +
  stat_compare_means(comparisons = my_comparisons,label = "p", method = "wilcox.test")

ggsave(filename = file.path(boxplotDir,paste0("SDL", "_acrossHist", ".png")),
       plot = g, width = 12, height = 12, units = "cm")




