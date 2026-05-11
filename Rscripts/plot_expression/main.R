# main program for making box or heat plots with expression values

### Load renv
library(renv)
renv::activate()
renv::snapshot()

### load packages and functions
source("./load_libraries.R")
#renv::snapshot()
source("./make_boxplot.R")
source("./make_heatplot.R")
source("./plot_expression.R")

p <- plot_expression # alias for function
# pre-loading
expression <- load_matrix()

input <- c("AT5G45830") # read from file
p(input = input, plottype = "boxplot",
  savefolder = "../results/Rplots/selectedGenes")

# read input # for the purpose of presentation

(readline() %>% strsplit(" "))[[1]] %>% p()

(readline() %>% strsplit(" "))[[1]] %>% p(plottype = "heatmap")

