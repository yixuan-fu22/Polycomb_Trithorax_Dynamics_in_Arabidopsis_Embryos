# This Script contains functions taht read expression (tpm) matrix made by Sleuth.Rmd, 
# and make box/heat plots to visualize expression in different tissues
# Rverison 4.2




load_matrix <- function(){
  # read in expression matrixs
  # embryo's
  em_expression <- read.table(file.path("..","results","tpm","embryos.tpm.csv"),sep = ",", 
                              header = TRUE)
  # sdl's
  sdl_expression <- read.table(file.path("..","results","tpm","sdl.tpm.csv"),sep = ",",
                               header = TRUE)
  # flower's
  flw_expression <- read.table(file.path("..","results","tpm","flw.tpm.csv"),sep = ",",
                               header = TRUE)
  # more samples
  #
  # combine
  expression <- rbind(em_expression, sdl_expression, flw_expression)
  return(expression)
}

get_symbol_tair_table <- function(gene_of_interests){
  tairSYMBOL <- org.At.tairSYMBOL %>% as.data.frame()
  return(tairSYMBOL %>% filter(symbol %in% gene_of_interests | gene_id %in% gene_of_interests) %>%  group_by(gene_id))
}

get_tair_and_symbol <- function(input){
  symbol_tair_table <- get_symbol_tair_table(input) %>% unique()
  tair.number <- character(length = length(input))
  gene.symbol <- character(length = length(input))
  for(i in 1:length(input)){
    if(grepl("^AT\\dG\\d+$", input[i]) == TRUE){
      tair.number[i] <- input[i]
      # find gene symbol
      gene.symbol[i] <- (symbol_tair_table %>% filter(gene_id == input[i]))$symbol %>% paste(collapse = ";")
    }else{
      gene.symbol[i] <- input[i]
      # find tair number
      temp <- (symbol_tair_table %>% filter(symbol == input[i]))$gene_id
      tair.number[i] <- temp[1]
      if(length(temp) > 1){
        print(paste0("more than 1 gene found for ", input[i],"; the first one was plotted"))
      }
      if(is.na(temp[1] %>% as.character)){
        print(paste0("no gene found for ", input[i]))
      }
    }
  }
  #remove NAs from the list
  gene.symbol <- gene.symbol[!is.na(tair.number)]
  tair.number <- tair.number[!is.na(tair.number)]
  return(list(tair.number=tair.number,gene.symbol=gene.symbol))
}


plot_expression <- function(input, expression, plottype = "boxplot", 
                            savefolder="../results/Rplots/boxplots", ylab="TPM"){
  tair_and_symbol <- get_tair_and_symbol(input)

  if(plottype == "boxplot"){
    boxplot <- make_boxplot(tair.number = tair_and_symbol$tair.number,
                            gene.symbol = tair_and_symbol$gene.symbol,
                            expression = expression,
                            color_by_group = TRUE,
                            path = savefolder,
                            ylab = ylab)
    boxplot
    return(boxplot)
  }else if(plottype == "heatmap"){
    heatplot <- make_heatplot(tair.number = tair_and_symbol$tair.number,
                              gene.symbol = tair_and_symbol$gene.symbol,
                              expression = expression,
                              path = savefolder)
    heatplot
    return(heatplot)
  }else{
    print("unknown plot type")
  }
}

