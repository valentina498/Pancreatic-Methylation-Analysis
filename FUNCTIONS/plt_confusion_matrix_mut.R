plt_confusion_matrix_mut <- function(confMatrix, project, mut){
  
  
  plt <- as.data.frame(confMatrix$table)
  nClass <- unique(plt$Reference)
  
  x_legend <- c("Mutated","WT")
  y_legend <- c("WT","Mutated")

  plt$Prediction <- factor(plt$Prediction, levels=rev(levels(plt$Prediction)))
  
  conf_plot <- ggplot(plt, aes(Reference,Prediction, fill= Freq)) +
  geom_tile() +
  geom_text(aes(label=Freq)) +
  scale_fill_gradient(low="white", high="green4") +
  labs(x = "Reference",y = "Prediction") +
  scale_x_discrete(labels=x_legend) +
  scale_y_discrete(labels=y_legend) +
  geom_point() +
  theme(text = element_text(size = 20), axis.text.x = element_text(angle = 90, hjust = 1))
  ggtitle(
    label = "Confusion Matrix",
    subtitle = paste(mut, "prediction on", project) +
      theme(
        plot.title = element_text(size = 16, face = "bold", color = "darkblue"),
        plot.subtitle = element_text(size = 11, face = "bold.italic", color = "deepskyblue")
        )
    ) +

  return(conf_plot)
}

plt_confusion_matrix_mut_REV <- function(confMatrix, project, mut){
  
  
  plt <- as.data.frame(confMatrix$table)
  nClass <- unique(plt$Reference)
  
  x_legend <- c("WT","Mutated")
  y_legend <- c("Mutated","WT")
  
  plt$Prediction <- factor(plt$Prediction, levels=rev(levels(plt$Prediction)))
  
  # conf_plot <- ggplot(plt, aes(Reference,Prediction, fill= Freq)) +
  #   geom_tile() +
  #   geom_text(aes(label=Freq)) +
  #   scale_fill_gradient(low="white", high="green4") +
  #   labs(x = "Reference",y = "Prediction") +
  #   scale_x_discrete(labels=x_legend) +
  #   scale_y_discrete(labels=y_legend) +
  #   #geom_point() +
  #   theme(text = element_text(size = 20), axis.text.x = element_text(angle = 90, hjust = 1)) +
  #   ggtitle(label = "Confusion Matrix",
  #           subtitle = paste(mut, "prediction on", project)
  # )
  
  # conf_plot <- ggplot(plt, aes(Reference, Prediction, fill = Freq)) +
  #   geom_tile() +
  #   geom_text(aes(label = Freq)) +
  #   scale_fill_gradient(low = "white", high = "green4") +
  #   labs(x = "Reference", y = "Prediction") +
  #   scale_x_discrete(labels = x_legend) +
  #   scale_y_discrete(labels = y_legend) +
  #   theme(text = element_text(size = 20),
  #         axis.text.x = element_text(angle = 90, hjust = 1)) +
  #   ggtitle(label = "Confusion Matrix",
  #           subtitle = paste(mut, "prediction on", project))
    
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
          plot.title = element_text(size = 40))+
      ggtitle(label = paste(mut, "prediction on", project)
    )
  
  
    return(conf_plot)
}