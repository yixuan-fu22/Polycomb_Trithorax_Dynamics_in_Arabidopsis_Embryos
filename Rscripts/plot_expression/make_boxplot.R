suppressMessages({
  library(dplyr)
  library(ggplot2)
  #library(ggpubr)
})

# get_tissue_type <- function(t){
#   if(t %in% c("pre-G","G","EH","LH","ET","LT","BC","MG")){
#     return("embryos")
#   }else if(t %in% c("SDL")){
#     return("seedlings")
#   }else if(t %in% c("FLW")){
#     return("flowers")
#   }else{
#     return("unknown")
#   }
# }

get_tissue_type <- function(t_list){
  type <- character(length = length(t_list))
  for(i in 1:length(t_list)){
    if(t_list[i] %in% c("pre-G","G","EH","LH","ET","LT","BC","MG")){
      type[i] <- "embryos"
    }else if(t_list[i] %in% c("SDL","SDL2")){
      type[i] <-"seedlings"
    }else if(t_list[i] %in% c("FLW")){
      type[i] <-"flowers"
    }else{
      type[i] <-"unknown"
    }
  }
  return(type)
}

make_boxplot <- function(tair.number,
                         gene.symbol,
                         expression,
                         level=c("pre-G","G","EH","LH","ET","LT","BC","MG","SDL","SDL2","FLW"),
                         #level=c("LT","MG", "SDL"),
                         path="../results/Rplots/boxplots",
                         yrange=c(0,NA),
                         color_by_group=FALSE,
                         save=TRUE,
                         ylab="TPM"){
  
  for(i in 1:length(tair.number)){
    selectedGene <- expression %>%
      dplyr::select(sample,tissue_short ,matches(tair.number[i])) %>%
      mutate(tissue_short=factor(tissue_short, levels=level)) %>%
      mutate(type=get_tissue_type(tissue_short)) %>%
      mutate(type=factor(type, levels = c("embryos","seedlings","flowers")))
    
    if(color_by_group){
      boxplot <- selectedGene %>%
        ggplot(aes(x=tissue_short, y=as.numeric(expression[[tair.number[i]]]),
                   color=type, fill=type)) + 
        geom_boxplot(alpha=0.2) +
        scale_fill_manual(values=c("orange", "green","red")) +
        scale_color_manual(values=c("red", "blue","purple")) +
        theme_bw(base_size = 14) +
        xlab("")+
        ylab(ylab)+
        ylim(yrange)+
        ggtitle(paste0(gene.symbol[i]," (",tair.number[i], ")"))+
        guides(color = FALSE, fill = FALSE)
      if(save){
        ggsave(file.path(path,paste0(gene.symbol[i], "_", tair.number[i], ".box.png")), boxplot, width = 6, 
               height = 4, units = "in") 
      }
    }else if(!color_by_group){
      boxplot <- selectedGene %>%
        ggplot(aes(x=tissue_short, y=as.numeric(expression[[tair.number[i]]]))) + 
        geom_boxplot(color="red", fill="orange", alpha=0.2) +
        theme_bw(base_size = 14) +
        xlab("")+
        ylab(ylab)+
        ylim(yrange)+
        ggtitle(paste0(gene.symbol[i]," (",tair.number[i], ")"))
      if(save){
        ggsave(file.path(path,paste0(gene.symbol[i], "_", tair.number[i], ".box.png")), boxplot, width = 6, 
               height = 4, units = "in") 
      }
    }
  }
  return(boxplot)
}