######################### MUTATION PREDICTION  #################################

library(TCGAbiolinks)
library(SummarizedExperiment)
library(sesame)
library(MASS)
library(DT)
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(grid)
library(reshape2)
library("survminer")
library(maftools)
library(cluster)
library(RColorBrewer)

#RF
library(randomForest)
library(caret)
library(PRROC)
library(pROC)

library(extrafont)
loadfonts(device = "win")

## FUNCTIONS
source("FUNCTIONS\\plt_confusion_matrix_mut.R")

## UTILS
load(file = "RESULTS\\DATAFRAME\\clinical_tcga_geneStrat.Rda")
row.names(clinical_tcga_geneStrat) <- NULL

load(file = "RESULTS\\DATAFRAME\\clinical_cptac_geneStrat.Rda")
row.names(clinical_cptac_geneStrat) <- NULL

load(file = "RESULTS\\DATAFRAME\\methyl_tcga_scale.Rda")
load(file = "RESULTS\\DATAFRAME\\methyl_cptac_scale.Rda")

##
mut_toPredict <- c("KRAS", "TP53", "KRAS_and_TP53")
name_mut <- c("KRAS", "TP53", "KRAS and TP53")
methyl_tcga <- as.data.frame(t(methyl_tcga_scale))
methyl_cptac <- as.data.frame(t(methyl_cptac_scale))
rownames(methyl_cptac) <- colnames(methyl_cptac_scale)

plt_roc <- list()
conf_plot_test <- list()
pr_plot <- list()

for (l in (1:length(mut_toPredict))){
  
  for (i in (1:nrow(methyl_tcga))){
    methyl_tcga$mut[i] <- clinical_tcga_geneStrat[clinical_tcga_geneStrat$submitter_id %in% rownames(methyl_tcga)[i],mut_toPredict[l]]
  }
  methyl_tcga$mut <- as.factor(methyl_tcga$mut)
  
  # train_data <- methyl_tcga[idx[[mut_toPredict[m]]],]
  # val_data <- methyl_tcga[-idx[[mut_toPredict[m]]],]
  
  
  for (i in (1:nrow(clinical_cptac_geneStrat))){
    methyl_cptac$mut[i] <- clinical_cptac_geneStrat[clinical_cptac_geneStrat$submitter_id %in% rownames(methyl_cptac)[i],mut_toPredict[l]]
  }
  
  set.seed(50)
  
  model_Rforest <- randomForest(mut ~ ., data = methyl_tcga, 
                                ntree = 100, importance = TRUE,verbose = TRUE)
  
  train_prob_forest <- predict(model_Rforest, methyl_tcga, type = "prob")[,2]
  test_prob_forest <- predict(model_Rforest, methyl_cptac, type = "prob")[,2]
  
  
  ############# CURVE ROC ##############
  roc_train <- roc(methyl_tcga$mut, train_prob_forest)
  roc_test <- roc(methyl_cptac$mut, test_prob_forest)
  
  
  # Calcolo AUC
  auc_train <- as.numeric(auc(roc_train))
  df_roc_train <- data.frame(
    specificity = rev(roc_train$specificities),
    sensitivity = rev(roc_train$sensitivities)
  )
  
  auc_test <- as.numeric(auc(roc_test))
  df_roc_test <- data.frame(
    specificity = rev(roc_test$specificities),
    sensitivity = rev(roc_test$sensitivities)
  )
  
  
  plt_roc[[l]] <- ggplot() +
    #geom_line(data = df_roc_train, aes(x = 1 - specificity, y = sensitivity), color = "dodgerblue4", size = 1.5) +
    #geom_line(data = df_roc_val, aes(x = 1 - specificity, y = sensitivity), color = "darkorange", linewidth = 1.5) +
    geom_line(data = df_roc_test, aes(x = 1 - specificity, y = sensitivity), color = "springgreen4", linewidth = 1.5) +
    #geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red", linewidth = 1.5) +
    geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1), linetype = "dashed", color = "red", linewidth = 1.5) +
    
    #scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    #scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    
    coord_equal() +
    # xlim(0, 1) +
    # ylim(0, 1) +
    
    labs(title = paste("AUROC", name_mut[l], "pred"),
         x = "1 - Specificity",
         y = "Sensitivity") +
    theme_minimal() +
    theme(axis.title = element_text(size = 20),
          axis.text = element_text(size = 20),
          plot.title = element_text( size = 22, hjust = 0.5, family = "Helvetica"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(color = "black"),
          axis.ticks = element_line(color = "black"),
          panel.border = element_rect(color = "black", fill = NA, size = 0.5) # Aggiunge il riquadro
          )+
    
    #annotate("text", x = 0.8, y = 0.2, label = paste("Train AUC =", round(auc_train, 2)), color = "dodgerblue4",size  = 6)+
    #annotate("text", x = 0.7, y = 0.15, label = paste("Validation AUC =", round(auc_val, 2)), color = "darkorange",size  = 6)+
    annotate("text", x = 0.7, y = 0.1, label = paste("Test AUROC =", round(auc_test, 2)), color = "springgreen4",size  = 7)
  
  ################# MATRICE DI CONFUSIONE ################
  # prediction on test set
  threshold <- 0.35
  #predictions <- ifelse(test_prob_forest > optimal_threshold$threshold, 1, 0)
  predictions <- ifelse(test_prob_forest > threshold, 1, 0)
  
  # Confusion matrix
  #confusion_matrix <- table(Predicted = predictions, Actual = methyl_cptac$mut)
  confusion_matrix <- confusionMatrix(factor(predictions), factor(methyl_cptac$mut),positive = "1")
  conf_plot_test[[l]] <- plt_confusion_matrix_mut_REV(confusion_matrix,"CPTAC-PDA",mut_toPredict[l])
  
  print(confusion_matrix)
  
  
  ################# PR CURVE ####################
  
  pr_train <- pr.curve(scores.class0 = train_prob_forest[methyl_tcga$mut == 1],
                       scores.class1 = train_prob_forest[methyl_tcga$mut == 0],
                       curve = TRUE)
  pr_test <- pr.curve(scores.class0 = test_prob_forest[methyl_cptac$mut == 1],
                      scores.class1 = test_prob_forest[methyl_cptac$mut == 0],
                      curve = TRUE)
  
  # PR curve
  pr_data <- data.frame(Recall = pr_test$curve[, 1], Precision = pr_test$curve[, 2])
  
  # Chance level
  chance_level <- mean(methyl_cptac$mut == 1)
  
  # Graph
  pr_plot[[l]] <- ggplot(pr_data, aes(x = Recall, y = Precision)) +
    geom_line(color = "blue", linewidth = 1.2) +
    geom_hline(yintercept = chance_level, linetype = "dashed", color = "black", linewidth = 1) +
    labs(title = paste("AUPRC", name_mut[l], "prediction"), x = "Recall", y = "Precision") +
    theme_minimal()+
    theme(axis.title = element_text(size = 20),
          axis.text = element_text(size = 20),
          plot.title = element_text( size = 22, hjust = 0.5, family = "Helvetica"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(color = "black"),
          axis.ticks = element_line(color = "black"),
          panel.border = element_rect(color = "black", fill = NA, size = 0.5) # Aggiunge il riquadro
    )+
    #scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) + # Limite asse x senza spazi
    #scale_y_continuous(limits = c(0, 1), expand = c(0, 0))+   # Limite asse y senza spazi
    annotate("text", x = 0.7, y = 0.1, label = paste("Test AUPRC =", round(pr_test$auc.integral, 2)), color = "blue",size  = 7)+
    annotate("text", x = 0.7, y = 0.05, label = paste("Chance Level =", round(chance_level, 2)), color = "black",size  = 7) 
  
  methyl_tcga <- as.data.frame(t(methyl_tcga_scale))
  
}

subplt_roc <- grid.arrange(grobs = plt_roc, ncol = 3)
subplt_CM <- grid.arrange(grobs = conf_plot_test, ncol = 3)
subplt_PR <- grid.arrange(grobs = pr_plot, ncol = 3)

ggsave(paste0("RESULTS\\IMMAGINI\\RF_mutClass\\RF_AUC_mutClassifier_NEW.png"), subplt_roc,
       width = 25,
       height = 9)

ggsave(paste0("RESULTS\\IMMAGINI\\RF_mutClass\\RF_CM_mutClassifier_NEW.png"), subplt_CM,
       width = 45,
       height = 10)

ggsave(paste0("RESULTS\\IMMAGINI\\RF_mutClass\\RF_PR_mutClassifier_NEW.png"), subplt_PR,
       width = 25,
       height = 9)

