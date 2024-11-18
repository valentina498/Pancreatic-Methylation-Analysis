####################### C4 vs OTHER ########################
library(TCGAbiolinks)
library(SummarizedExperiment)
library(sesame)
library(MASS)
library(DT)
library(ConsensusClusterPlus)
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(reshape2)
library("survminer")
library(survival)
library(maftools)
library(cluster)
library(RColorBrewer)


##UTILS
load(file = paste0("RESULTS\\DATAFRAME\\clinical_tcga_CC_KRAS_and_TP53.Rda"))
load(file = "RESULTS\\DATAFRAME\\eval_CC.Rda")


load(file = "RESULTS\\DATAFRAME\\clinical_tcga_geneStrat.Rda")
row.names(clinical_tcga_geneStrat) <- NULL

load(file = "RESULTS\\DATAFRAME\\clinical_cptac_geneStrat.Rda")
row.names(clinical_tcga_geneStrat) <- NULL

load(file = "RESULTS\\DATAFRAME\\methyl_tcga_scale.Rda")
load(file = "RESULTS\\DATAFRAME\\methyl_cptac_scale.Rda")
methyl_cptac_scale <- methyl_cptac_scale[, !duplicated(colnames(methyl_cptac_scale))]

##FUNCTION
source("FUNCTIONS\\cluster_annotation.R")
source("./FUNCTIONS/cluster_annotation_class.R")
source("FUNCTIONS\\plt_confusion_matrix.R")
source("FUNCTIONS\\CpG_importance.R")

## Choosing of best mutation comb from clustering parameter
best_CC <- "cc4_pam_euclidean"
clinical_tcga_CC[,best_CC] <- ifelse(clinical_tcga_CC[,best_CC] == "4", "C4","Co")

## Kaplan-Meier con C4 vs other
test_comb <- "KRAS_and_TP53"
name_folder <- paste0(test_comb,"\\cons_pam_euclidean\\")

fit <- survfit(Surv(as.numeric(clinical_tcga_CC$overall_survival), clinical_tcga_CC$deceased) ~ clinical_tcga_CC[,best_CC])
plot <- ggsurvplot(fit,
           data = clinical_tcga_CC,
           pval = T,
           
           risk.table = T,
           risk.table.y.text = FALSE,
           xlab ="Time in Months",
           break.time.by = 12,
           legend.labs =
             c("C4","Co"),
           xlim = c(0,84),
           
           title = "TCGA-PAAD", 
           subtitle = paste(best_CC),
           font.title = c(16, "bold", "darkblue"),
           font.subtitle = c(11, "bold.italic", "deepskyblue"))

sub_ann <- cluster_annotation(clinical_tcga_CC, "KRAS_and_TP53","scale", "pam", "euclidean", name_folder)
ggsave(paste0("RESULTS\\IMMAGINI\\ConsensusClustering\\",name_folder,"c2_annotation.png"), plot = sub_ann,
       width = 19,
       height = 21)


########################## RANDOM FOREST: CLUSTER prediction ###############################
library(randomForest)
library(caret)
library(Hmisc) ##cindex

library(extrafont)
loadfonts(device = "win")

##FLAGS
calc_RF = FALSE


#### CLINICAL DATA PREPARATION ####
# TCGA
clinical_tcga_RF <- clinical_tcga_CC[,c(1:21)]
clinical_tcga_RF$KRAS_and_TP53 <- clinical_tcga_CC$KRAS_and_TP53
clinical_tcga_RF$cluster <- clinical_tcga_CC[,best_CC]

# CPTAC
clinical_cptac_RF <- clinical_cptac_geneStrat[,c(1:21)]
clinical_cptac_RF$KRAS_and_TP53 <- clinical_cptac_geneStrat$KRAS_and_TP53
clinical_cptac_RF$KRAS <- as.factor(ifelse(clinical_cptac_RF$KRAS == 1, "Mutant","WT"))
clinical_cptac_RF$TP53 <- as.factor(ifelse(clinical_cptac_RF$TP53 == 1, "Mutant","WT"))
clinical_cptac_RF$KRAS_and_TP53 <- as.factor(ifelse(clinical_cptac_RF$KRAS_and_TP53 == 1, "Mutant","WT"))

#### METH DATA PREPARATION ####

## Add LABEL column in methyl_tcga_scale
methyl_tcga_RF <- as.data.frame(t(as.matrix(methyl_tcga_scale)))
for (i in (1:nrow(clinical_tcga_CC))){
  methyl_tcga_RF$cluster[i] <- clinical_tcga_CC[clinical_tcga_CC$submitter_id == rownames(methyl_tcga_RF)[i],best_CC]
}
methyl_tcga_RF$cluster <- as.factor(methyl_tcga_RF$cluster)
#methyl_tcga_RF <- as.matrix(methyl_tcga_RF)

## Transpose methyl_cptac_scale
methyl_cptac_RF <- as.data.frame(t(as.matrix(methyl_cptac_scale)))
rownames(methyl_cptac_RF) <- colnames(methyl_cptac_scale)

if (calc_RF == TRUE){
  
  ## TRAINING OF RANDOM FOREST MODEL
  rf_model <- randomForest(cluster ~ ., data = methyl_tcga_RF, 
                           ntree = 100, mtry = floor(sqrt(ncol(methyl_tcga_RF))), 
                           importance = TRUE, verbose = TRUE)
  
  saveRDS(rf_model, file = paste0("./RESULTS/DATAFRAME/Classifier/RF/RF_model_KRAS_and_TP53.rds"))
  
  print("Training done!")
  
} else {
  rf_model <- readRDS(file = paste0("./RESULTS/DATAFRAME/Classifier/RF/RF_model_KRAS_and_TP53.rds"))
  print("Loading done!")
}




### PREDICTION ON TRAIN
train_pred_forest <- as.data.frame(rf_model$predicted)
train_pred_forest$pred <- as.character(train_pred_forest[,"rf_model$predicted"])

for (p in (1:nrow(clinical_tcga_RF))){
  clinical_tcga_RF$cluster_pred[p] <- train_pred_forest[rownames(train_pred_forest) == clinical_tcga_RF$submitter_id[p],2]
}
clinical_tcga_RF$cluster <- as.factor(clinical_tcga_RF$cluster)
clinical_tcga_RF$cluster_pred <- as.factor(clinical_tcga_RF$cluster_pred)

#### FEATURE IMPORTANCE ####
feature_importance <- randomForest::importance(rf_model)
res <- CpG_importance(feature_importance, "cc2")
plot_importance <- res$plot

feature_sort <- data.frame(matrix(nrow=200, ncol=0))
sorted_gini <- res$sorted_gini
sorted_acc <- res$sorted_acc

col_gini <- paste0("cc2","_gini")
feature_sort[,col_gini] <- rownames(sorted_gini)[1:200]

col_acc <- paste0("cc2","_acc")
feature_sort[,col_acc] <- rownames(sorted_acc)[1:200]
save(feature_sort, file = "RESULTS/DATAFRAME/Classifier/RF/RF_feature_importance.Rda")

#### CONFUSION MATRIX OF RF MODEL ####
confMatrix <- confusionMatrix(clinical_tcga_RF$cluster_pred, clinical_tcga_RF$cluster)
print(confMatrix)

conf_plot <- plt_confusion_matrix(confMatrix)
ggsave(paste0("./RESULTS/IMMAGINI/ConsensusClustering/KRAS_and_TP53/cons_pam_euclidean/RF/RF_confusionM_cc2.png"), conf_plot,
       width = 15,
       height = 15)

### PREDICTION ON TEST
test_pred_forest <- predict(rf_model, methyl_cptac_RF)
test_pred_forest <- as.data.frame(test_pred_forest)
test_pred_forest$test_pred_forest <- as.character(test_pred_forest$test_pred_forest)
for (p in (1:nrow(clinical_cptac_RF))){
  clinical_cptac_RF$cluster_pred[p] <- test_pred_forest[rownames(test_pred_forest) == clinical_cptac_RF$submitter_id[p],1]
}
clinical_cptac_RF$cluster_pred <- as.factor(clinical_cptac_RF$cluster_pred)




############################### KAPLAN-MEIER ###################################
####### comparison of gt and pred on training dataset and prediction on test set

##TCGA-PAAD
fit_cluster <- survfit(Surv(as.numeric(clinical_tcga_RF$overall_survival), 
                            clinical_tcga_RF$deceased) ~ clinical_tcga_RF[,"cluster"])
fit_pred <- survfit(Surv(as.numeric(clinical_tcga_RF$overall_survival), 
                         clinical_tcga_RF$deceased) ~ clinical_tcga_RF[,"cluster_pred"])
fit_summary <- summary(fit_pred)

cindex_cluster <- rcorr.cens(as.numeric(clinical_tcga_RF[,"cluster"]), Surv(as.numeric(clinical_tcga_RF$overall_survival), 
                                                                            clinical_tcga_RF$deceased))
cindex_pred <- rcorr.cens(as.numeric(clinical_tcga_RF[,"cluster_pred"]), Surv(as.numeric(clinical_tcga_RF$overall_survival), 
                                                                              clinical_tcga_RF$deceased))
print(cindex_cluster)
print(cindex_pred)
plot_cluster <- ggsurvplot(fit_cluster,
                         data = clinical_tcga_RF,
                         pval = T,
                         
                         risk.table = T,
                         risk.table.y.text = FALSE,
                         xlab ="Time in Months",
                         break.time.by = 12,
                         legend.labs =
                           c("C4","C1/C2/C3"),
                         xlim = c(0,72),
                         
                         title = "TCGA-PAAD",
                         legend.title = "",
                         palette = c("#C77CFF", "black")
                         #font.title = c(16, "bold", "darkblue"),
                         #font.subtitle = c(11, "bold.italic", "deepskyblue")
                         )
plot_cluster$plot <- plot_cluster$plot + theme(
  plot.title = element_text(family = "Helvetica",size = 24, hjust = 0.5),
  legend.title = element_text(family = "Helvetica",size = 24),  # Dimensione del font del titolo della legenda
  legend.text = element_text(family = "Helvetica",size = 24),   # Dimensione del font del testo della legenda
  axis.title = element_text(family = "Helvetica",size = 24),    # Dimensione del font dei nomi degli assi
  axis.text = element_text(family = "Helvetica",size = 24)      # Dimensione del font dei valori degli assi
)

plot_cluster$table <- plot_cluster$table + theme(
  axis.title = element_text(family = "Helvetica",size = 24),
  axis.text = element_text(family = "Helvetica",size = 24)
)


plot_train <- ggsurvplot(fit_pred,
                        data = clinical_tcga_RF,
                        pval = T,
                        
                        risk.table = T,
                        risk.table.y.text = FALSE,
                        xlab ="Time in Months",
                        break.time.by = 12,
                        legend.labs =
                          c("C4","C1/C2/C3"),
                        xlim = c(0,72),
                        
                        title = "TCGA-PAAD prediction",
                        legend.title = "",
                        palette = c("#C77CFF", "black")
                        #font.title = c(16, "bold", "darkblue"),
                        #font.subtitle = c(11, "bold.italic", "deepskyblue")
                        )
plot_train$plot <- plot_train$plot + theme(
  plot.title = element_text(family = "Helvetica",size = 24, hjust = 0.5),
  legend.title = element_text(family = "Helvetica",size = 24),  # Dimensione del font del titolo della legenda
  legend.text = element_text(family = "Helvetica",size = 24),   # Dimensione del font del testo della legenda
  axis.title = element_text(family = "Helvetica",size = 24),    # Dimensione del font dei nomi degli assi
  axis.text = element_text(family = "Helvetica",size = 24)      # Dimensione del font dei valori degli assi
)

plot_train$table <- plot_train$table + theme(
  axis.title = element_text(family = "Helvetica",size = 24),
  axis.text = element_text(family = "Helvetica",size = 24)
)


##CPTAC-3
fit_pred_test <- survfit(Surv(as.numeric(clinical_cptac_RF$overall_survival), 
                               clinical_cptac_RF$deceased) ~ clinical_cptac_RF[,"cluster_pred"],
                          data = clinical_cptac_RF)
fit_summary <- summary(fit_pred_test)
cindex_test <- rcorr.cens(as.numeric(clinical_cptac_RF[,"cluster_pred"]), Surv(as.numeric(clinical_cptac_RF$overall_survival), 
                                                                                clinical_cptac_RF$deceased))
print(cindex_test)

plot_test <- ggsurvplot(fit_pred_test,
                   data = clinical_cptac_RF,
                   pval = T,
                   
                   risk.table = T,
                   risk.table.y.text = FALSE,
                   xlab ="Time in Months",
                   break.time.by = 12,
                   legend.labs =
                     c("C4","C1/C2/C3"),
                   xlim = c(0,72),
                   
                   title = "CPTAC-3 prediction",
                   legend.title = "",
                   palette = c("#C77CFF", "black")
                   #font.title = c(16, "bold", "darkblue"),
                   #font.subtitle = c(11, "bold.italic", "deepskyblue")
                   )
plot_test$plot <- plot_test$plot + theme(
  plot.title = element_text(family = "Helvetica",size = 24, hjust = 0.5),
  legend.title = element_text(family = "Helvetica",size = 24),  # Dimensione del font del titolo della legenda
  legend.text = element_text(family = "Helvetica",size = 24),   # Dimensione del font del testo della legenda
  axis.title = element_text(family = "Helvetica",size = 24),    # Dimensione del font dei nomi degli assi
  axis.text = element_text(family = "Helvetica",size = 24)      # Dimensione del font dei valori degli assi
)

plot_test$table <- plot_test$table + theme(
  axis.title = element_text(family = "Helvetica",size = 24),
  axis.text = element_text(family = "Helvetica",size = 24)
)

plot_KM_tot <- arrange_ggsurvplots(list(plot_cluster,plot_train,plot_test), print = TRUE,
                            ncol = 3, risk.table.height = 0.3)
ggsave(paste0("RESULTS\\IMMAGINI\\ConsensusClustering\\",name_folder,"/RF/KM_cluster_train_test.png"), plot_KM_tot,
       width = 15,
       height = 7)

##################### CLUSTER ANNOTATION ######################################
sub_ann_tcga <- cluster_annotation_class("cc2",clinical_tcga_RF,"TCGA-PAAD","KRAS_and_TP53")
sub_ann_cptac <- cluster_annotation_class("cc2",clinical_cptac_RF,"CPTAC","KRAS_and_TP53")


###### PLOT ######
ggsave(paste0("RESULTS\\IMMAGINI\\ConsensusClustering\\",name_folder,"/RF/annotation_trainRF_tcga.png"), sub_ann_tcga$subplot_annotation,
       width = 19,
       height = 21)
ggsave(paste0("RESULTS\\IMMAGINI\\ConsensusClustering\\",name_folder,"/RF/annotation_testRF_cptac.png"), sub_ann_cptac$subplot_annotation,
       width = 19,
       height = 21)


