plt_confusion_matrix <- function(confMatrix){
  
  
  plt <- as.data.frame(confMatrix$table)
  nClass <- unique(plt$Reference)
  
  if (length(nClass)==4){
    x_legend <- c("C1","C2","C3","C4")
    y_legend <- c("C4","C3","C2","C1")
  } else if (length(nClass)==3){
    x_legend <- c("C1","C2","C3")
    y_legend <- c("C3","C2","C1")
  }else if (length(nClass)==2){
    x_legend <- c("C4","C1/C2/C3")
    y_legend <- c("C1/C2/C3","C4")
  }
  plt$Prediction <- factor(plt$Prediction, levels=rev(levels(plt$Prediction)))
  
  conf_plot <- ggplot(plt, aes(Reference,Prediction, fill= Freq)) + 
    geom_tile() +
    geom_text(aes(label=Freq), size = 20) +
    scale_fill_gradient(low="white", high="green4") +
    labs(x = "Reference",y = "Prediction") +
    scale_x_discrete(labels=x_legend) + 
    scale_y_discrete(labels=y_legend) + 
    theme(axis.text = element_text(size = 60),
          axis.title = element_text(size = 60),
          legend.key.size =  unit(2,'cm'),
          legend.text = element_text(size = 40),
          legend.title = element_text(size = 60),
          axis.text.y = element_text(angle = 90, hjust = 0.5))
    # # ggtitle(label = "Confusion Matrix",
    #         subtitle = paste(length(nClass), "cluster"))+
    # theme(
    #   plot.title = element_text(size = 18, face = "bold", color = "darkblue"),
    #   plot.subtitle = element_text(size = 16, face = "bold.italic", color = "deepskyblue"),
    #   axis.title = element_text(size = 18, family = "Helvetica"),
    #   axis.text = element_text(size = 18, family = "Helvetica"),
    #   legend.title = element_text(size = 16, family = "Helvetica"),
    #   legend.text = element_text(size = 16, family = "Helvetica"))
  
  return(conf_plot)
}