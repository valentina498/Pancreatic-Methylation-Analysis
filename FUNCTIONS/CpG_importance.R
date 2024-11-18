CpG_importance <- function(feature_importance, nCluster){
  
  if (nCluster == "cc4"){
    k = 4
  } else if (nCluster == "cc3"){
    k = 3
  } else if (nCluster == "cc2"){
    k = 2
  }
  #feature_importance <- importance(rf_model)
  feature_importance_df <- as.data.frame(feature_importance)
  feature_importance_df$Feature <- rownames(feature_importance_df)
  
  sorted_gini <- feature_importance_df[order(-feature_importance_df$MeanDecreaseGini), ]
  sorted_acc <- feature_importance_df[order(-feature_importance_df$MeanDecreaseAccuracy), ]
  
  plt_gini <- ggplot(sorted_gini[c(1:50),], aes(x=reorder(Feature, MeanDecreaseGini), y=MeanDecreaseGini),fill=MeanDecreaseGini) +
    geom_bar(stat='identity') +
    coord_flip() +
    #scale_fill_gradient(low = "white", high = "#8B0000") +
    xlab('Feature') +
    ylab('Importance (MeanDecreaseGini)') +
    ggtitle(label = "Feature importance: Mean Decrease Gini",
            subtitle = paste(k, "cluster")) +
    theme(
      plot.title = element_text(size = 16, face = "bold", color = "darkblue"),
      plot.subtitle = element_text(size = 11, face = "bold.italic", color = "deepskyblue"))
  
  plt_acc <- ggplot(sorted_acc[c(1:50),], aes(x=reorder(Feature, MeanDecreaseAccuracy), y=MeanDecreaseAccuracy),fill=MeanDecreaseAccuracy) +
    geom_bar(stat='identity') +
    coord_flip() +
    #scale_fill_gradient(low = "white", high = "#8B0000") +
    xlab('Feature') +
    ylab('Importance (MeanDecreaseAccuracy)') +
    ggtitle(label = "Feature importance: Mean Decrease Accuracy",
            subtitle = paste(k, "cluster"))+
    theme(
      plot.title = element_text(size = 16, face = "bold", color = "darkblue"),
      plot.subtitle = element_text(size = 11, face = "bold.italic", color = "deepskyblue"))
  
  plot <- grid.arrange(plt_gini, plt_acc, ncol = 2)
  
  res = list(plot = plot, sorted_gini = sorted_gini, sorted_acc = sorted_acc)
  return(res)
  
}