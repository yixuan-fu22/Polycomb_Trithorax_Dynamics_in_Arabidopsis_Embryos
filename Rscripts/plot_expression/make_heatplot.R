suppressMessages({
  library(dplyr)
  library(ggplot2)
})

make_heatplot <- function(tair.number,
                         gene.symbol,
                         expression,
                         level=c("pre-G","G","EH","LH","ET","LT","BC","MG", "SDL", "SDL2", "FLW"),
                         path="../results/Rplots/heatplots",
                         yrange=c(0,NA),
                         save=TRUE){
  for (i in 1:length(tair.number)) {
    heatplot <-
      expression %>%
      mutate(tissue_short=factor(tissue_short, levels=level)) %>%
      ggplot(aes(x=tissue_short, y=gene.symbol[i],fill=as.numeric(expression[[tair.number[i]]]))) +
      geom_tile()+
      xlab("")+
      ylab("")+
#      theme_ipsum_es()+
      labs(fill="TPM")+
      theme(plot.background = element_rect(fill = 'white'))+
#      gradient_color(c("blue", "white", "red"))+
      scale_fill_gradient(
        low = "blue",
        high = "red"
      )
    # ylim(yrange)+
    # ggtitle(paste0(gene.symbol," (",tair.number, ")"))
    
    #plots <- plots %>% append(heatplot)
    if(save){
      ggsave(file.path(path,paste0(gene.symbol[i], "_", tair.number[i], ".heat.png")), heatplot, 
             width = 10, height = 2*length(tair.number[i]), units = "in")
    }
  }
  return(heatplot) # return the last plot
}










