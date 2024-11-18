##################### ANALISI DEI DATI DI MUTAZIONE ############################
library(TCGAbiolinks)
library(SummarizedExperiment)
library(sesame)
library(MASS)
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(reshape2)
library("survminer")
library(survival)
library(maftools)

## UTILS
load(file = "RESULTS\\DATAFRAME\\clinical_tcga_strat.Rda")
load(file = "RESULTS\\DATAFRAME\\clinical_cptac_strat.Rda")


load(file = "RESULTS\\DATAFRAME\\SNV_tcga_maf.Rda")
load(file = "RESULTS\\DATAFRAME\\SNV_cptac_maf.Rda")

## FUNCTION
source("FUNCTIONS\\matrix_mut.R")
source("FUNCTIONS\\calc_gene_sig.R")
source("FUNCTIONS\\comb_mutation.R")
source("FUNCTIONS\\plot_pval_geneComb.R")

## FLAG
flag_save = TRUE
flag_gene_sig = FALSE
flag_create_mutMatrix = FALSE
flag_prepare_comb = FALSE
flag_comb_sig = TRUE

###### ESPLORAZIONE MUTAZIONI ######

# TCGA-PAAD
maf_tcga <- read.maf(maf = snv_tcga)
plotmafSummary(maf = maf_tcga,
               addStat = 'median',
               rmOutlier = TRUE,
               dashboard = TRUE)
oncoplot(maf = maf_tcga, top = 20, removeNonMutated = TRUE)


# CPTAC-3
maf_cptac <- read.maf(maf = snv_cptac)
plotmafSummary(maf = maf_cptac,
               addStat = 'median',
               rmOutlier = TRUE,
               dashboard = TRUE)
oncoplot(maf = maf_cptac, top = 20, removeNonMutated = TRUE)

###### CREAZIONE DELLE MUTATION MATRIX ######
## sono matrici con geni sulle colonne, pazienti sulle righe e si trova 1 se il gene è mutato
## e 0 se il gene non è mutato

## La funzione matrix_mut calcola anche gli N geni più mutati nel progetto e riempie
## "clinical" con le colonne relative a quei geni

if (flag_create_mutMatrix == TRUE){
  
  ## TCGA-PAAD
  res_tcga <- matrix_mut(snv_tcga, clinical_tcga_strat, 20)
  mutMatrix_tcga <- res_tcga$mutMatrix
  clinical_tcga_geneStrat <- res_tcga$clinical_geneStrata
  
  save(mutMatrix_tcga, file = "RESULTS\\DATAFRAME\\mutMatrix_tcga.Rda")
  save(clinical_tcga_geneStrat, file = "RESULTS\\DATAFRAME\\clinical_tcga_geneStrat.Rda")
  
  ## CPTAC-3
  res_cptac <- matrix_mut(snv_cptac, clinical_cptac_strat, 20)
  mutMatrix_cptac <- res_cptac$mutMatrix
  clinical_cptac_geneStrat <- res_cptac$clinical_geneStrata
  
  save(mutMatrix_cptac, file = "RESULTS\\DATAFRAME\\mutMatrix_cptac.Rda")
  save(clinical_cptac_geneStrat, file = "RESULTS\\DATAFRAME\\clinical_cptac_geneStrat.Rda")
  
  
} else {
  
  load(file = "RESULTS\\DATAFRAME\\mutMatrix_tcga.Rda")
  load(file = "RESULTS\\DATAFRAME\\clinical_tcga_geneStrat.Rda")
  
  load(file = "RESULTS\\DATAFRAME\\mutMatrix_cptac.Rda")
  load(file = "RESULTS\\DATAFRAME\\clinical_cptac_geneStrat.Rda")
}

## Prendo solo i primi 20 geni top
top_genes_cptac <- colnames(clinical_cptac_geneStrat[20:39])
top_genes_tcga <- colnames(clinical_tcga_geneStrat[20:39])
top_genes_tot <- intersect(top_genes_cptac, top_genes_tcga)


###### CALCOLA I GENI SIGNIFICATIVI IN TCGA tra quelli in comune tra tcga e cptac ######
if (flag_gene_sig == TRUE){
  ##### NON FUNZIONA
  # fit_cptac <- survfit(Surv(as.numeric(clinical_cptac_geneStrat$overall_survival), clinical_cptac_geneStrat$deceased) ~ clinical_cptac_geneStrat[,gene])
  # fit_tcga <- survfit(Surv(as.numeric(clinical_tcga_geneStrat$overall_survival), clinical_tcga_geneStrat$deceased) ~ clinical_tcga_geneStrat[,gene])
  res <- calc_gene_sig(top_genes_tot, 
                                fit_tcga, clinical_tcga_geneStrat, 
                                TRUE, 0.1)
  gene_sig_tcga <- res$gene_sig
  gene_pval <- res$gene_pval
  save(gene_sig_tcga, file = "RESULTS\\DATAFRAME\\gene_sig_tcga.Rda")
} else {
  load(file = "RESULTS\\DATAFRAME\\gene_sig_tcga.Rda")
}

###### PREPARA I DF CLINICAL AGGIUNGENDO LE COLONNE DELLE COMBINAZIONI DI MUTAZIONI ######
if (flag_prepare_comb == TRUE){
  
  ## TCGA-PAAD
  clinical_tcga_geneStrat <- comb_mutation(clinical_tcga_geneStrat)
  
  ## CPTAC-3
  clinical_cptac_geneStrat <- comb_mutation(clinical_cptac_geneStrat)
  
  save(clinical_tcga_geneStrat, file = "RESULTS\\DATAFRAME\\clinical_tcga_geneStrat.Rda")
  save(clinical_cptac_geneStrat, file = "RESULTS\\DATAFRAME\\clinical_cptac_geneStrat.Rda")
} else{
  
  load(file = "RESULTS\\DATAFRAME\\clinical_tcga_geneStrat.Rda")
  load(file = "RESULTS\\DATAFRAME\\clinical_cptac_geneStrat.Rda")
}

###### CALCOLA LE COMBINAZIONI DI MUT SIGNIFICATIVE IN TCGA ######
mut_combination <- c("KRAS", "TP53","CDKN2A",
                     "KRAS_and_TP53","KRAS_and_CDKN2A","TP53_and_CDKN2A",
                     "KRAS_or_TP53","KRAS_or_CDKN2A","TP53_or_CDKN2A",
                     "KRAS_and_TP53_and_CDKN2A", "KRAS_or_TP53_or_CDKN2A")

mut_title <- c("KRAS", "TP53","CDKN2A",
               "KRAS and TP53","KRAS and CDKN2A","TP53 and CDKN2A",
               "KRAS or TP53","KRAS or CDKN2A","TP53 or CDKN2A",
               "KRAS and TP53 and CDKN2A", "KRAS or TP53 or CDKN2A")
if (flag_comb_sig == TRUE){
  
  clinical_tcga <- clinical_tcga_geneStrat
  clinical_cptac <- clinical_cptac_geneStrat
  
  mut_comb_data <- data.frame(mut = mut_combination)
  
  
  for (i in (1:length(mut_combination))){
    gene_comb <- mut_combination[i]
    gene_title <- mut_title[i]
    
    fit_cptac <- survfit(Surv(as.numeric(clinical_cptac$overall_survival), clinical_cptac$deceased) ~ clinical_cptac[,gene_comb])
    fit_tcga <- survfit(Surv(as.numeric(clinical_tcga$overall_survival), clinical_tcga$deceased) ~ clinical_tcga[,gene_comb])
    
    comb_pval <- plot_pval_geneComb(clinical_tcga, fit_tcga, clinical_cptac, fit_cptac,
                                             gene_comb, gene_title,TRUE, i)
    
    mut_comb_data$pval[i] <- comb_pval
    mut_comb_data$N_wt[i] <- table(clinical_tcga[,gene_comb])["0"]
    mut_comb_data$N_mut[i] <- table(clinical_tcga[,gene_comb])["1"]

    
  }
  save(mut_comb_data, file = "RESULTS\\DATAFRAME\\mut_comb_data.Rda")
} else {
  load(file = "RESULTS\\DATAFRAME\\mut_comb_data.Rda")
}


