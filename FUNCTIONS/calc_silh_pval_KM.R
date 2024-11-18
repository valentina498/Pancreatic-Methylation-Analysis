########################## FUNZIONI ########################
calc_silh_pval_KM <- function(methyl,methyl_scale,test_comb,algCC,distCC,display, clinical, thresh_num, flag_save){
  
  library("survminer")
  library(survival)
  ## methyl = matrice di metilazione con pazienti sulle colonne e siti CpG sulle righe
  ## dataType = "scale" o "normal" se i dati sono rispettivamente standardizzati o no
  ## algCC = algoritmo da utilizzare per fare clustering ("hc", "km", "pam")
  ## distCC = parametro per definire la distanza tra gli item ("pearson","spearman","euclidean","binary","maximum","canberra","minkowski")
  ##                                                           N.B. l'algoritmo "km" supporta solo la distanza "euclidean"
  ## display = TRUE o FALSE per mostrare o meno l'avanzamento di ConsensusClusterPlus
  ## clinical = matrice dei dati clinici per OS 
  
  ######################### DATI NORMALI #######################################
  results_N = ConsensusClusterPlus(methyl, maxK=6,
                                 reps=500,pItem=0.8,pFeature=1,
                                 clusterAlg = algCC, distance = distCC,
                                 seed=50, 
                                 plot="png",
                                 title = paste0("RESULTS\\IMMAGINI\\ConsensusClustering\\",test_comb,"\\cons_normal_",algCC,"_",distCC),
                                 verbose = display)
  
  ######################### DATI STANDARDIZZATI ################################
  results_S = ConsensusClusterPlus(methyl_scale, maxK=6,
                                   reps=500,pItem=0.8,pFeature=1,
                                   clusterAlg = algCC, distance = distCC,
                                   seed=50, 
                                   plot="png",
                                   title = paste0("RESULTS\\IMMAGINI\\ConsensusClustering\\",test_comb,"\\cons_scale_",algCC,"_",distCC),
                                   verbose = display)
  
  
  #### INSERISCI I RISULTATI DEI CLUSTER IN CLINICAL PER COSTRUIRE LE KM
  result_N_cc3 <- results_N[[3]]$consensusClass
  result_N_cc4 <- results_N[[4]]$consensusClass
  
  result_S_cc3 <- results_S[[3]]$consensusClass
  result_S_cc4 <- results_S[[4]]$consensusClass
  
  
  for (k in (1:nrow(clinical))){
    clinical[clinical$submitter_id == names(result_N_cc3)[k],paste0("cc3_normal_",algCC,"_",distCC)] <- result_N_cc3[names(result_N_cc3)[k]]
    clinical[clinical$submitter_id == names(result_N_cc4)[k],paste0("cc4_normal_",algCC,"_",distCC)] <- result_N_cc4[names(result_N_cc4)[k]]
    
    clinical[clinical$submitter_id == names(result_S_cc3)[k],paste0("cc3_scale_",algCC,"_",distCC)] <- result_S_cc3[names(result_S_cc3)[k]]
    clinical[clinical$submitter_id == names(result_S_cc4)[k],paste0("cc4_scale_",algCC,"_",distCC)] <- result_S_cc4[names(result_S_cc4)[k]]
  }
  
  ######################### DATI NORMALI #######################################
  #### CALCOLA LE METRICHE DA INSERISE POI IN eval_param ####
  
  # Calcola la silhouette
  cc3_N <- results_N[[3]]
  cc4_N <- results_N[[4]]
  cc3Sil_N = silhouette(x = cc3_N[[3]], 
                      dist = as.matrix(1- cc3_N[[4]])) #per ottenere la matrice di distanza e non di similarità si sottrae la matrice di similarità ad 1
  cc4Sil_N = silhouette(x = cc4_N[[3]], 
                      dist = as.matrix(1- cc4_N[[4]])) 
  
  # Calcola il p-value
  string3 <- paste0("cc3_normal_",algCC,"_",distCC)
  string4 <- paste0("cc4_normal_",algCC,"_",distCC)
  
  log_rank3 <- survdiff(formula = Surv(overall_survival, deceased) ~ clinical[,string3], data = clinical)
  log_rank4 <- survdiff(formula = Surv(overall_survival, deceased) ~ clinical[,string4], data = clinical)
  
  cc3Pval_N <- log_rank3$pvalue
  cc4Pval_N <- log_rank4$pvalue
  
  # Calcola la numerosità dei cluster
  numItem_CC3_N <- "good"
  numItem_CC4_N <- "good"
  if(any(table(cc3_N$consensusClass)<thresh_num) == TRUE){
    numItem_CC3_N <- min(table(cc3_N$consensusClass))
  }
  
  if(any(table(cc4_N$consensusClass)<thresh_num) == TRUE){
    numItem_CC4_N <-  min(table(cc4_N$consensusClass))
  }
  
  ######################### DATI STANDARDIZZATI ###############################
  #### CALCOLA LE METRICHE DA INSERISE POI IN eval_param ####
  
  # Calcola la silhouette
  cc3_S <- results_S[[3]]
  cc4_S <- results_S[[4]]
  cc3Sil_S = silhouette(x = cc3_S[[3]], 
                        dist = as.matrix(1- cc3_S[[4]])) #per ottenere la matrice di distanza e non di similarità si sottrae la matrice di similarità ad 1
  cc4Sil_S = silhouette(x = cc4_S[[3]], 
                        dist = as.matrix(1- cc4_S[[4]])) 
  
  # Calcola il p-value
  string3 <- paste0("cc3_scale_",algCC,"_",distCC)
  string4 <- paste0("cc4_scale_",algCC,"_",distCC)
  
  log_rank3 <- survdiff(formula = Surv(overall_survival, deceased) ~ clinical[,string3], data = clinical)
  log_rank4 <- survdiff(formula = Surv(overall_survival, deceased) ~ clinical[,string4], data = clinical)
  
  cc3Pval_S <- log_rank3$pvalue
  cc4Pval_S <- log_rank4$pvalue
  
  # Calcola la numerosità dei cluster
  numItem_CC3_S <- "good"
  numItem_CC4_S <- "good"
  if(any(table(cc3_S$consensusClass)<thresh_num) == TRUE){
    numItem_CC3_S <- min(table(cc3_S$consensusClass))
  }
  
  if(any(table(cc4_S$consensusClass)<thresh_num) == TRUE){
    numItem_CC4_S <-  min(table(cc4_S$consensusClass))
  }
  
  risultati_calc_sil_p <- list(results_N = results_N,
                               results_S = results_S,
                               
                               cc3Sil_N = cc3Sil_N, cc4Sil_N = cc4Sil_N,
                               cc3Sil_S = cc3Sil_S, cc4Sil_S = cc4Sil_S,
                               
                               cc3Pval_N = cc3Pval_N, cc4Pval_N = cc4Pval_N,
                               cc3Pval_S = cc3Pval_S, cc4Pval_S = cc4Pval_S,
                               
                               numCC3_N = numItem_CC3_N, numCC4_N = numItem_CC4_N,
                               numCC3_S = numItem_CC3_S, numCC4_S = numItem_CC4_S,
                               
                               clinical_CC = clinical)
  return(risultati_calc_sil_p)
}

calc_silh_pval_KM_REV <- function(methyl,test_comb,algCC,distCC,display, clinical, thresh_num, flag_save){
  
  library("survminer")
  library(survival)
  ## methyl = matrice di metilazione con pazienti sulle colonne e siti CpG sulle righe
  ## dataType = "scale" o "normal" se i dati sono rispettivamente standardizzati o no
  ## algCC = algoritmo da utilizzare per fare clustering ("hc", "km", "pam")
  ## distCC = parametro per definire la distanza tra gli item ("pearson","spearman","euclidean","binary","maximum","canberra","minkowski")
  ##                                                           N.B. l'algoritmo "km" supporta solo la distanza "euclidean"
  ## display = TRUE o FALSE per mostrare o meno l'avanzamento di ConsensusClusterPlus
  ## clinical = matrice dei dati clinici per OS 
  
  ######################### DATI NORMALI #######################################
  results = ConsensusClusterPlus(methyl, maxK=6,
                                   reps=500,pItem=0.8,pFeature=1,
                                   clusterAlg = algCC, distance = distCC,
                                   seed=50, 
                                   plot="png",
                                   title = paste0("RESULTS\\IMMAGINI\\ConsensusClustering\\",test_comb,"\\cons_",algCC,"_",distCC),
                                   verbose = display)
  
  
  #### INSERISCI I RISULTATI DEI CLUSTER IN CLINICAL PER COSTRUIRE LE KM
  result_cc3 <- results[[3]]$consensusClass
  result_cc4 <- results[[4]]$consensusClass
  
  
  
  for (k in (1:nrow(clinical))){
    clinical[clinical$submitter_id == names(result_cc3)[k],paste0("cc3_",algCC,"_",distCC)] <- result_cc3[names(result_cc3)[k]]
    clinical[clinical$submitter_id == names(result_cc4)[k],paste0("cc4_",algCC,"_",distCC)] <- result_cc4[names(result_cc4)[k]]
    
  }
  
  #### CALCOLA LE METRICHE DA INSERISE POI IN eval_param ####
  
  # Calcola la silhouette
  cc3 <- results[[3]]
  cc4 <- results[[4]]
  cc3Sil = silhouette(x = cc3[[3]], 
                        dist = as.matrix(1- cc3[[4]])) #per ottenere la matrice di distanza e non di similarità si sottrae la matrice di similarità ad 1
  cc4Sil = silhouette(x = cc4[[3]], 
                        dist = as.matrix(1- cc4[[4]])) 
  
  # Calcola il p-value
  string3 <- paste0("cc3_",algCC,"_",distCC)
  string4 <- paste0("cc4_",algCC,"_",distCC)
  
  log_rank3 <- survdiff(formula = Surv(overall_survival, deceased) ~ clinical[,string3], data = clinical)
  log_rank4 <- survdiff(formula = Surv(overall_survival, deceased) ~ clinical[,string4], data = clinical)
  
  cc3Pval <- log_rank3$pvalue
  cc4Pval <- log_rank4$pvalue
  
  # Calcola la numerosità dei cluster
  numItem_CC3 <- "good"
  numItem_CC4 <- "good"
  if(any(table(cc3$consensusClass)<thresh_num) == TRUE){
    numItem_CC3 <- min(table(cc3$consensusClass))
  }
  
  if(any(table(cc4$consensusClass)<thresh_num) == TRUE){
    numItem_CC4 <-  min(table(cc4$consensusClass))
  }
  
  
  risultati_calc_sil_p <- list(results = results,
                               
                               cc3Sil = cc3Sil, cc4Sil = cc4Sil,
                               
                               cc3Pval = cc3Pval, cc4Pval = cc4Pval,
                               
                               numCC3 = numItem_CC3, numCC4 = numItem_CC4,
                               
                               clinical_CC = clinical)
  return(risultati_calc_sil_p)
}