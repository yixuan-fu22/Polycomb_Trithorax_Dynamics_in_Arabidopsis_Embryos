# main program for making box or heat plots with expression values

### Load renv
library(renv)
renv::activate()

### load packages and functions
source("./plot_expression/load_libraries.R")
#renv::snapshot()
source("./plot_expression/make_boxplot.R")
source("./plot_expression/make_heatplot.R")
source("./plot_expression/plot_expression.R")

# load normalized expression matrix
normalized_expression <- read.table(file.path("..","data_fromOtherProjects","tpm_fromRNAseqOfEmbryosSeedlingsAndFlowers","embryos_sdl_flw.deseq2normalized.tpm.csv"),sep = ",",
                                    header = TRUE)

# input <- read.table(file.path("input","gene_list.txt"))$V1 # read from file
# plot_expression(input = input, expression = normalized_expression, 
#                plottype = "boxplot",
#  savefolder = "../results/Rplots/Normalized",
#  ylab = "TPM, DEseq2 Normalization")


# Plotting log transformation
# log_normalized_expression <- normalized_expression %>%
#   dplyr::select(-c(tissue,tissue_short,sample)) %>% +1 %>% log2() %>%
#   cbind(., normalized_expression %>% dplyr::select(c(tissue,tissue_short,sample))) %>%
#   select(tissue, tissue_short, sample, everything())
log_normalized_expression <- normalized_expression %>%
  dplyr::select(tissue, tissue_short, sample, everything()) %>%
  mutate(across(where(is.numeric), ~log2(. + 1)))


#input <- c("AT1G21970", "AT5G47670", "AT3G50870", "AT3G15030", "AT4G27160", "AT3G24650")



#input <- c("AT2G26760", "AT5G67260", "AT3G12280", "AT1G63100")
# PIN6/ATGH9B1/REM35/ATMCM2/CDC2B/ARR10/ACT2/AtCER6/ATRBR1/ATEXP6/SKU5/RGI5/ACT8/JMJ16/FAS1/BPP1/AtLNG1/ATEXPB3/GlcAT14A

saveFolder <- "../results/expressionPlots_individualGenes"

input <- c("AT1G77110")

plot_expression(input = input, expression = log_normalized_expression, 
                plottype = "heatmap",
                savefolder = saveFolder,
                ylab = "log2(DEseq2-normalized TPM + 1)")

#input = c("LEC1", "LEC2", "ABI3", "STM", "FUS3")
input = "CYCB1;4"

plot_expression(input = input, expression = log_normalized_expression %>% filter(tissue_short %in% c("LT", "MG")), 
                plottype = "boxplot",
                savefolder = saveFolder,
                ylab = "log2(DEseq2-normalized TPM + 1)")

#plot_expression(input = input, expression = log_normalized_expression, 
#                plottype = "boxplot",
#                savefolder = saveFolder,
#                ylab = "log2(DEseq2-normalized TPM + 1)")

 

