#
# This script is to use DEseq2 result table to make scatter plots and MAplots to show show the comparison between two stages
# scatter plot: X and Y axis are the intensity of histone modifcation at two stages, each scatter represent  a peak
# MA plot: log Fold Change vs Mean Intensity
# with peaks that are significantly different highlighted, and genes of interest annotated
# similiar to 30979870:fig1e

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
  # library(ggpubr)
  # library(ggupset)
})

# SCATTER PLOT #

# make a result dir
resDir_scatter <- "../results/XYscatterPlot"
if(!file.exists(resDir_scatter)){dir.create(resDir_scatter)}

# define color platte
# color_platte <- c("#999999", "#D55E00", "#0072B2") # NS, Up in Stage1, Up in Stage2
color_platte <- c("grey", "blue", "red")

# scatter plot function
scatter_plot <- function(inputFilePath, hist, resDir_scatter, # input: annotated DEseq2 resultes
                         color_platte,
                         flipStages = F){
  FDR_th=0.05
  Fold_th=0
  DEseq2_res_annotated <- inputFilePath %>%
    read.table(header = TRUE, fill = TRUE, sep = "\t", quote = "")
  # extract stages compared
  stage1 <- names(DEseq2_res_annotated)[7] %>% str_remove("Conc_")
  stage2 <- names(DEseq2_res_annotated)[8] %>% str_remove("Conc_")
  # define significant peaks
  DEseq2_res_annotated <- DEseq2_res_annotated %>% mutate(significance=if_else((FDR<FDR_th) & (abs(Fold) > Fold_th), 
                                                                               if_else(.data[[paste0("Conc_", stage1)]] > 
                                                                                       .data[[paste0("Conc_", stage2)]],
                                                                                       paste0("Up in ", stage1),
                                                                                       paste0("Up in ", stage2)), "NS"))
  
  upPeaksStage1 <- sum(DEseq2_res_annotated$significance == paste0("Up in ", stage1), na.rm = TRUE)
  upPeaksStage2 <- sum(DEseq2_res_annotated$significance == paste0("Up in ", stage2), na.rm = TRUE)
  
  if(flipStages){
  label_txt <- paste0(
    "Up in ", stage1, ": ", upPeaksStage1, "\n",
    "Up in ", stage2, ": ", upPeaksStage2, "\n",
    "FDR < ", FDR_th, "\n") }else{ #    "Abs(log2FC) > ", Fold_th
    label_txt <- paste0(
      "Up in ", stage2, ": ", upPeaksStage2, "\n",
      "Up in ", stage1, ": ", upPeaksStage1, "\n",
      "FDR < ", FDR_th, "\n") }

  
  print(label_txt)
  # plot
  p <- ggplot() +
    geom_point(data = DEseq2_res_annotated %>% filter(significance == "NS"), aes(x=.data[[paste0("Conc_", stage1)]], y=.data[[paste0("Conc_", stage2)]]),
               alpha = 0.4, color = color_platte[1]) + # plot NS points
    geom_point(data= DEseq2_res_annotated %>% filter(significance == paste0("Up in ", stage1)), aes(x=.data[[paste0("Conc_", stage1)]], y=.data[[paste0("Conc_", stage2)]]),
               alpha = 0.7, color = color_platte[if_else(!flipStages, 2, 3)]) + # plot Stage1 Up points # SWITCH the color accordingly if X and Y axis are flipped
    geom_point(data = DEseq2_res_annotated %>% filter(significance == paste0("Up in ", stage2)), aes(x=.data[[paste0("Conc_", stage1)]], y=.data[[paste0("Conc_", stage2)]]),
               alpha = 0.7, color = color_platte[if_else(!flipStages, 3, 2)]) #+ # plot Stage2 Up points
  #guides(color="none")
 
  
  p <- p + theme_classic(base_family = "Arial", base_size = 8) +
    scale_x_continuous(labels = scales::number_format(accuracy = 0.1, decimal.mark = '.')) +
    scale_y_continuous(labels = scales::number_format(accuracy = 0.1, decimal.mark = '.'))
  if(flipStages){
    p <- p + coord_flip()
  }
  
  p <- p +
    annotate(
      "text",
      x = -Inf, y = Inf,
      label = label_txt,
      hjust = -0.05, vjust = 1.05,
      size = 4
    ) + coord_cartesian(clip = "off")
  
  w = 18
  h = 18
  ggsave(filename = file.path(resDir_scatter, paste0("XYscatterPlot_", stage1, "vs", stage2, "_", hist,
                                             "_", w, "x", h, ".png")),
         units = "cm",
         width = w,
         height = h)
  
  return(p)
}

# Plot
hist = "H2AUB"
scatter_plot(file.path("../results/diffBind_ChIPseeker_WSvsSDL_H2AUB_newDataset", 
                       "DESeq2Res_withAnnotations.tsv"), hist, resDir_scatter, color_platte, flipStages = F)
#scatter_plot(file.path("../results/diffBind_ChIPseeker_WSvsSDL_K4ME3",
#                       "DESeq2Res_withAnnotations.tsv"), hist, resDir_scatter, color_platte)

# plot for H2Aub
hist = "H2AUB"
scatter_plot(file.path("../results/diffBind_ChIPseeker_WSvsSDL_H2AUB",
                       "DESeq2Res_withAnnotations.tsv"), hist, resDir_scatter, color_platte, flipStages = T)
scatter_plot(file.path("../results/diffBind_ChIPseeker_WSvsMG_H2AUB",
                       "DESeq2Res_withAnnotations.tsv"), hist, resDir_scatter, color_platte, flipStages = T)

# Plot for H3K27me3 # Use unannotated DEseq2 results for now. I will add annotation later


# MA PLOT #

resDir_MA <- "../results/MAPlot"
if(!file.exists(resDir_MA)){dir.create(resDir_MA)}

# MA plot function
MA_plot <- function(inputFilePath, hist, resDir_MA, # input: annotated DEseq2 resultes
                         color_platte){
  FDR_th = 0.05
  Fold_th = 1
  DEseq2_res_annotated <- inputFilePath %>%
    read.table(header = TRUE, fill = TRUE, sep = "\t", quote = "")
  # extract stages compared
  stage1 <- names(DEseq2_res_annotated)[7] %>% str_remove("Conc_")
  stage2 <- names(DEseq2_res_annotated)[8] %>% str_remove("Conc_")
  # define significant peaks
  DEseq2_res_annotated <- DEseq2_res_annotated %>% mutate(significance=if_else( ((FDR<FDR_th) & (abs(Fold) > Fold_th ) ), 
                                                                               if_else(.data[[paste0("Conc_", stage1)]] > 
                                                                                         .data[[paste0("Conc_", stage2)]],
                                                                                       paste0("Up in ", stage1),
                                                                                       paste0("Up in ", stage2)), "NS"))
  
  # plot
  p <- ggplot() +
    geom_point(data = DEseq2_res_annotated %>% filter(significance == "NS"), aes(x=.data[["Conc"]], y=.data[["Fold"]]),
               alpha = 0.4, color = color_platte[1]) + # plot NS points
    geom_point(data= DEseq2_res_annotated %>% filter(significance == "Up in WS"), aes(x=.data[["Conc"]], y=.data[["Fold"]]),
               alpha = 0.7, color = color_platte[2]) + # plot Stage1 Up points
    geom_point(data = DEseq2_res_annotated %>% filter(significance == "Up in MG"), aes(x=.data[["Conc"]], y=.data[["Fold"]]),
               alpha = 0.7, color = color_platte[3]) # plot Stage2 Up points

  p <- p + theme_classic(base_family = "Arial", base_size = 8) +
    scale_x_continuous(labels = scales::number_format(accuracy = 0.1, decimal.mark = '.')) +
    scale_y_continuous(labels = scales::number_format(accuracy = 0.1, decimal.mark = '.'))
  w = 18
  h = 18
  ggsave(filename = file.path(resDir_MA, paste0("MAPlot_", stage1, "vs", stage2, "_", hist,
                                                     "_", w, "x", h, ".png")),
         units = "cm",
         width = w,
         height = h)
  
  return(p)
}

# Plot
hist = "H2AUB"
MA_plot(file.path("../results/diffBind_ChIPseeker_MGvsSDL_H2AUB_newDataset",
                       "DESeq2Res_withAnnotations.tsv"), hist, resDir_MA, color_platte)

# plot for H2Aub
hist = "H2AUB"
MA_plot(file.path("../results/diffBind_ChIPseeker_WSvsSDL_H2AUB",
                       "DESeq2Res_withAnnotations.tsv"), hist, resDir_MA, color_platte)
MA_plot(file.path("../results/diffBind_ChIPseeker_WSvsMG_H2AUB",
                       "DESeq2Res_withAnnotations.tsv"), hist, resDir_MA, color_platte)











