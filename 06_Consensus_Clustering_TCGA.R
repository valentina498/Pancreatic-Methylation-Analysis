####################### CONSENSUS CLUSTERING TCGA-PAAD #########################
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
load(file = "RESULTS\\DATAFRAME\\clinical_tcga_geneStrat.Rda")
row.names(clinical_tcga_geneStrat) <- NULL

load(file = "RESULTS\\DATAFRAME\\methyl_tcga_scale.Rda")

##FUNCTION
source("FUNCTIONS\\calc_silh_pval_KM.R")
source("FUNCTIONS\\CC_kaplan_meier.R")
source("FUNCTIONS\\cluster_annotation.R")

## Clustering params
test_comb <- "KRAS_and_TP53"
v_clusterAlg <- c("hc","pam","km")
v_distance <- c("pearson","spearman","euclidean","maximum","canberra")

## Empty df

col <- c("mut","clusterAlg","dataType","distance",
         "silCC3","pvalCC3","numItem_CC3",
         "silCC4","pvalCC4","numItem_CC4")
eval_param <- data.frame(matrix(nrow=0, ncol=length(col)))
colnames(eval_param) <- col

## Clinical df preparation
for (i in (20:ncol(clinical_tcga_geneStrat))){
  clinical_tcga_geneStrat[,i] <- ifelse(clinical_tcga_geneStrat[,i]=="1","Mutant","WT")
}
  
methyl_consensus <- as.matrix(methyl_tcga_scale)
clinical_tcga_CC <- clinical_tcga_geneStrat

clinical_tcga_CC$age_strat <- ifelse(clinical_tcga_CC$age_strat == "Over 65", ">65", "≤ 65")
clinical_tcga_CC$gender <- ifelse(clinical_tcga_CC$gender == "female", "Female", "Male")

for (i in (1: length(v_clusterAlg))){
  if(i == 3){
    v_distance <- c("euclidean")
  } else {
    v_distance <- c("pearson","spearman","euclidean","maximum","canberra")
  }
  
  for (j in (1:length(v_distance))){
    
    res <- calc_silh_pval_KM_REV(methyl_consensus,test_comb,v_clusterAlg[i],v_distance[j],
                                 FALSE, clinical_tcga_geneStrat, 20, TRUE)
    
    clinical_tcga_CC[[paste0("cc3_",v_clusterAlg[i],"_",v_distance[j])]] <- res$clinical_CC[,paste0("cc3_",v_clusterAlg[i],"_",v_distance[j])]
    clinical_tcga_CC[[paste0("cc4_",v_clusterAlg[i],"_",v_distance[j])]] <- res$clinical_CC[,paste0("cc4_",v_clusterAlg[i],"_",v_distance[j])]
    save(clinical_tcga_CC, file = paste0("RESULTS\\DATAFRAME\\clinical_tcga_CC_",test_comb,".Rda"))
    
    ### Results acquisition from res
    cc3 <- res$results[[3]]
    cc4 <- res$results[[4]]
    temp <- data.frame(mut = test_comb,clusterAlg=v_clusterAlg[i], dataType = "scale", distance = v_distance[j], 
                       silCC3 = mean(res$cc3Sil[,3]), pvalCC3 = res$cc3Pval, numItem_CC3 = res$numCC3,
                       silCC4 = mean(res$cc4Sil[,3]), pvalCC4= res$cc4Pval, numItem_CC4 = res$numCC4)
    
    eval_param <- rbind(eval_param, temp)
    
    
    ################### Kaplan-Meier ################
    ###
    name_folder <- paste0(test_comb,"\\cons_",v_clusterAlg[i],"_",v_distance[j],"\\")
    ###
    
    
    clinical <- clinical_tcga_CC
    
    string <- paste0("cc3_",v_clusterAlg[i],"_",v_distance[j])
    fit_cc3 <- survfit(Surv(as.numeric(clinical$overall_survival), clinical$deceased) ~ clinical[,string])
    
    string <- paste0("cc4_",v_clusterAlg[i],"_",v_distance[j])
    fit_cc4 <- survfit(Surv(as.numeric(clinical$overall_survival), clinical$deceased) ~ clinical[,string])
    
    
    CC_kaplan_meier_REV(test_comb, v_clusterAlg[i], v_distance[j],
                        clinical, fit_cc3, fit_cc4, 
                        name_folder,TRUE)
    
    ################### Other plot #################
    
    
    ######### Annotation ########
    sub_ann_normal <- cluster_annotation(clinical_tcga_CC, test_comb,"scale", v_clusterAlg[i], v_distance[j], name_folder)
    
    ggsave(paste0("RESULTS\\IMMAGINI\\ConsensusClustering\\",name_folder,"cc4_annotation.png"), plot = sub_ann_normal,
           width = 18,
           height = 21)
    
    
    
    while (dev.cur() != 1) dev.off()
  }
  
}
save(eval_param, file = "RESULTS\\DATAFRAME\\eval_CC.Rda")


