############################# CTREE C_4 vs C_other #############################
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

## CTREE
library(riskRegression)
library(survival)
library(partykit)
library(coin)
library(Hmisc)  # C-index
library(cowplot)
library(randomForestSRC)
library(extrafont)
loadfonts(device = "win")


## UTILS
load(file = "RESULTS\\DATAFRAME\\clinical_tcga_geneStrat.Rda")
row.names(clinical_tcga_geneStrat) <- NULL

load(file = "RESULTS\\DATAFRAME\\clinical_cptac_geneStrat.Rda")
row.names(clinical_tcga_geneStrat) <- NULL

load(file = "RESULTS\\DATAFRAME\\methyl_tcga_scale.Rda")
load(file = "RESULTS\\DATAFRAME\\methyl_cptac_scale.Rda")
methyl_cptac_scale <- methyl_cptac_scale[, !duplicated(colnames(methyl_cptac_scale))]

load(file = "RESULTS/DATAFRAME/Classifier/RF/RF_feature_importance.Rda")

load(file = "RESULTS\\DATAFRAME\\train_mean_zScore.Rda")
load(file = "RESULTS\\DATAFRAME\\train_sd_zScore.Rda")

## Clinical data preparation
# TCGA
clinical_tcga_CT <- clinical_tcga_geneStrat[,c(1:21)]
clinical_tcga_CT$KRAS_and_TP53 <- clinical_tcga_geneStrat$KRAS_and_TP53

# CPTAC
clinical_cptac_CT <- clinical_cptac_geneStrat[,c(1:21)]
clinical_cptac_CT$KRAS_and_TP53 <- clinical_cptac_geneStrat$KRAS_and_TP53
clinical_cptac_CT$KRAS <- as.factor(ifelse(clinical_cptac_CT$KRAS == 1, "Mutant","WT"))
clinical_cptac_CT$TP53 <- as.factor(ifelse(clinical_cptac_CT$TP53 == 1, "Mutant","WT"))
clinical_cptac_CT$KRAS_and_TP53 <- as.factor(ifelse(clinical_cptac_CT$KRAS_and_TP53 == 1, "Mutant","WT"))

## Meth data preparation (features selection)
methyl_tcga_CT <- as.data.frame(t(as.matrix(methyl_tcga_scale)))
methyl_cptac_CT <- as.data.frame(t(as.matrix(methyl_cptac_scale)))

rank_type <- 1 #choose 1 = gini e 2 = accuracy
n_features <- 10 #number from 1 to 200 (numbers of features to training ctree model)

methyl_tcga_CT <- methyl_tcga_CT[,feature_sort[1:n_features,rank_type]]
methyl_cptac_CT <- methyl_cptac_CT[,feature_sort[1:n_features,rank_type]]

for (i in (1:nrow(methyl_tcga_CT))){
  methyl_tcga_CT[i,"overall_survival"] <- clinical_tcga_CT[clinical_tcga_CT$submitter_id == rownames(methyl_tcga_CT)[i],"overall_survival"]
  methyl_tcga_CT[i,"deceased"] <- clinical_tcga_CT[clinical_tcga_CT$submitter_id == rownames(methyl_tcga_CT)[i],"deceased"]
}

for (i in (1:nrow(methyl_cptac_CT))){
  methyl_cptac_CT[i,"overall_survival"] <- clinical_cptac_CT[clinical_cptac_CT$submitter_id == rownames(methyl_cptac_CT)[i],"overall_survival"]
  methyl_cptac_CT[i,"deceased"] <- clinical_cptac_CT[clinical_cptac_CT$submitter_id == rownames(methyl_cptac_CT)[i],"deceased"]
}
################# TRAINING ON TCGA #######################
ctree_surv <- ctree(Surv(overall_survival, deceased) ~. , 
                    data = methyl_tcga_CT)

print(ctree_surv)
plot(ctree_surv)
png("./RESULTS/IMMAGINI/CTREE/ctree_plot.png", width = 800, height = 600)  # Specifica larghezza e altezza desiderate
plot(ctree_surv)
dev.off()

thresh_zscore <- -0.4869
cpg_node <- "cg16854533"
thresh_normal <- thresh_zscore * train_sd[cpg_node] + train_mean[cpg_node]

################# PREDICTION ON TCGA AND CPTAC #######################
pred_tcga <- predict(ctree_surv, newdata = methyl_tcga_CT, type = "response")
pred_cptac <- predict(ctree_surv, newdata = methyl_cptac_CT, type = "response")

pred_median <- sort(unique(pred_tcga))
pred_label <- c("G1", "G2")

for (i in (1:nrow(methyl_tcga_CT))){
  methyl_tcga_CT[i,"pred_CT"] <- pred_tcga[names(pred_tcga) == rownames(methyl_tcga_CT)[i]]
}
methyl_tcga_CT$pred_CT <- ifelse(methyl_tcga_CT$pred_CT == "16", pred_label[1],pred_label[2])

for (i in (1:nrow(methyl_cptac_CT))){
  methyl_cptac_CT[i,"pred_CT"] <- pred_cptac[names(pred_cptac) == rownames(methyl_cptac_CT)[i]]
}
methyl_cptac_CT$pred_CT <- ifelse(methyl_cptac_CT$pred_CT == "16", pred_label[1],pred_label[2])


plt_ggsurv <- list()

# Calc and plotting survival curve of TRAIN and TEST
fit <- survfit(Surv(overall_survival, deceased) ~ pred_CT, data = methyl_tcga_CT)
print(fit)
plt_ggsurv[[1]] <- ggsurvplot(fit, 
                              data = methyl_tcga_CT, 
                              pval = T,
                              pval.coord = c(0, 0.15),
                              break.time.by = 12,
                              risk.table = T,
                              risk.table.y.text = FALSE,
                              xlim = c(0,74),
                              #legend.labs=c("Group 1", "Group 2", "Group 3"),
                              legend.labs=c("Hypomethylated", "Hypermethylated"),
                              palette = c("#C84D4C","#44A043"),
                              title = paste("TCGA-PAAD"),
                              legend.title=""
                              # legend.labs =c(pred_label)
)
plt_ggsurv[[1]]$plot <- plt_ggsurv[[1]]$plot + theme(
  plot.title = element_text(hjust = 0.5, family = "Helvetica"),
  legend.title = element_text(family = "Helvetica",size = 14),  # Dimensione del font del titolo della legenda
  legend.text = element_text(family = "Helvetica",size = 14),   # Dimensione del font del testo della legenda
  axis.title = element_text(family = "Helvetica",size = 18),    # Dimensione del font dei nomi degli assi
  axis.text = element_text(family = "Helvetica",size = 18)      # Dimensione del font dei valori degli assi
)

plt_ggsurv[[1]]$table <- plt_ggsurv[[1]]$table + theme(
  axis.title = element_text(family = "Helvetica",size = 18),
  axis.text = element_text(family = "Helvetica",size = 18)
)

fit <- survfit(Surv(overall_survival, deceased) ~ pred_CT, data = methyl_cptac_CT)
print(fit)
plt_ggsurv[[2]] <- ggsurvplot(fit, 
                              data = methyl_cptac_CT, 
                              pval = T,
                              pval.coord = c(0, 0.15),
                              break.time.by = 12,
                              risk.table = T,
                              risk.table.y.text = FALSE,
                              xlim = c(0,62),
                              #legend.labs=c("Group 1", "Group 2", "Group 3"),
                              legend.labs=c("Hypomethylated", "Hypermethylated"),  
                              palette = c("#C84D4C","#44A043"),
                              title = paste("CPTAC-PDA"),
                              legend.title=""
                              # legend.labs = c(pred_label)
)
plt_ggsurv[[2]]$plot <- plt_ggsurv[[2]]$plot + theme(
  plot.title = element_text(hjust = 0.5, family = "Helvetica"),
  legend.title = element_text(family = "Helvetica",size = 14),  # Dimensione del font del titolo della legenda
  legend.text = element_text(family = "Helvetica",size = 14),   # Dimensione del font del testo della legenda
  axis.title = element_text(family = "Helvetica",size = 18),    # Dimensione del font dei nomi degli assi
  axis.text = element_text(family = "Helvetica",size = 18)      # Dimensione del font dei valori degli assi
)

plt_ggsurv[[2]]$table <- plt_ggsurv[[2]]$table + theme(
  axis.title = element_text(family = "Helvetica",size = 18),
  axis.text = element_text(family = "Helvetica",size = 18)
)
plot <- arrange_ggsurvplots(list(plt_ggsurv[[1]],plt_ggsurv[[2]]), print = TRUE,
                            ncol = 2, risk.table.height = 0.3)
ggsave(paste0("./RESULTS/IMMAGINI/CTREE/KM_ctree_prediction_2Group.png"), plot,
       width = 10,
       height = 6)

