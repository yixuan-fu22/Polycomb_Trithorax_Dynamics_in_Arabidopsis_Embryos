# main program for making box or heat plots with expression values

### Load renv
library(renv)
#renv::init(bioconductor = 3.22)
renv::activate()

### load packages and functions
source("./load_libraries.R")
#renv::snapshot()
source("./make_boxplot.R")
source("./make_heatplot.R")
source("./plot_expression.R")

p <- plot_expression # alias for function
# pre-loading
expression <- load_matrix()

input <- read.table(file.path("input","gene_list.txt"))$V1 # read from file
# p(input = input, expression = expression, plottype = "boxplot",
#  savefolder = "../results/Rplots/motifEnrichmentGenes_290125")

# plot

p(input = c("TRB3"), expression = expression, plottype = "boxplot")

# read input # for the purpose of presentation

(readline() %>% strsplit(" "))[[1]] %>% p()

(readline() %>% strsplit(" "))[[1]] %>% p(plottype = "heatmap")

