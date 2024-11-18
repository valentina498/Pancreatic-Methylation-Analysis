############################################
########   CpG Probes Filtering   ########## 

library(ChAMPdata)

#################
###   STEP 1  ###
#######   Merging of probe in common between Illumina EPIC (CPTAC) and Illumina 450K (TCGA)   #########
do_merge = FALSE
save_probe = FALSE

if(do_merge == TRUE){
  data<-read.csv("C:\\Users\\valen\\Downloads\\humanmethylation450_15017482_v1-2.csv",header=T,sep=",",na.strings="",row.names=NULL,quote="",skip=7)
  probe_tcga <- read.csv("UTILS\\HM450.hg38.manifest.gencode.v36.tsv", sep='\t')
  probe_cptac <- read.csv("UTILS\\EPIC.hg38.manifest.gencode.v36.tsv", sep='\t')
  
  # Merging of probes EPIC and probes 450K info
  probe_merge <- merge(probe_tcga, probe_cptac)
  probe <- data[data$IlmnID %in% probe_merge$probeID,colnames(data) %in% c("CHR","IlmnID","UCSC_RefGene_Name","UCSC_RefGene_Accession","UCSC_RefGene_Group")]
  colnames(probe)<-c("Probe","Chrm","Gene","Accession","Region")
  
  # Remove ChrX and ChrY probes
  probe <- probe[!(probe$Chrm %in% c("X","Y")),]
  
  probe<-probe[!(is.na(probe$Gene)),]
  
  if (save_probe == TRUE){
    save(probe, file="RESULTS\\DATAFRAME\\probe_merge_EPIC_450K.Rda")
  }
} else {
  load(file="RESULTS\\DATAFRAME\\probe_merge_EPIC_450K.Rda")
}
## STEP 1 results : probes in common between 450K and EPIC, 
## Excluding the CpG on the X and Y chromosomes, 
## and those that do not have the gene information


###   STEP 2  ###
####### Exclution of NaN Values #######

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

## UTILS
load(file = "RESULTS\\DATAFRAME\\clinical_tcga_geneStrat.Rda")
row.names(clinical_tcga_geneStrat) <- NULL

load(file = "RESULTS\\DATAFRAME\\clinical_cptac_geneStrat.Rda")
row.names(clinical_cptac_geneStrat) <- NULL

tcga_data <- readRDS(file = "RESULTS\\DATAFRAME\\Meth_tcga_data.rds")
load(file = "RESULTS\\DATAFRAME\\MethDF_tcga.Rda")
cptac_data <- readRDS(file = "RESULTS\\DATAFRAME\\Meth_cptac_data.rds")
load(file = "RESULTS\\DATAFRAME\\MethDF_cptac.Rda")

load(file = "RESULTS\\DATAFRAME\\mut_comb_data.Rda")

## FUNCTION
source("FUNCTIONS\\filter_CpG_met.R")

#### Best mutations' combination
#### with lowest p value
best_comb <- mut_comb_data$mut[mut_comb_data$pval==min(mut_comb_data$pval)]
mut_title <- c("mutated KRAS and TP53")

### Filtering of methylation dataset
probe_merge_id <- probe$Probe
tcga_data_filtered <- filter_CpG_met(tcga_data, probe_merge_id, best_comb, clinical_tcga_geneStrat)
methyl_tcga <- as.data.frame(assay(tcga_data_filtered))


## STEP 2 results: CpG probes with NaN<30%


### STEP 3 
#### Differential methylated CpG site analysis

calc_dmc = FALSE

if (calc_dmc == TRUE){
  
  dmc <- TCGAanalyze_DMC(
    data = tcga_data_filtered,
    groupCol = best_comb,
    p.cut = 0.05,
    diffmean.cut = 0.15,
    adj.method = "BH",
    
    save = FALSE,
    
    #Parametri di figura
    plot.filename = paste0("RESULTS\\IMMAGINI\\VolcanoPlot_DMC\\",i,"_",best_comb,"_DCM.png"),
    title = paste0("Differential Methylated CpG (TCGA-PAAD)\n",mut_title),
    legend = "Legend",
    #label = c("Not Significant","Hypermethylated", "Hypomethylated"),
    
    cores = 1 # if set to 1 there will be a progress bar
  )
  save(dmc, file = paste0("RESULTS\\DATAFRAME\\DMC\\",best_comb,"_DCM.Rda"))
  
} else{
  
  load(file = paste0("RESULTS\\DATAFRAME\\DMC\\",best_comb,"_DCM.Rda"))
  dmc_hypo_hyper <- dmc[dmc$status != "Not Significant",]
  dmc_hypo <- dmc[dmc$status == "Hypomethylated in Mutant",]
  dmc_hyper <- dmc[dmc$status == "Hypermethylated in Mutant",]
  
}
## STEP 3 results: DMC in KRAS/TP53 co-mutation


### STEP 4
#### Standardization of DMC

Zscore =TRUE
if(Zscore == TRUE){
  methyl_tcga_normal <- as.matrix(methyl_tcga[rownames(methyl_tcga) %in% rownames(dmc_hypo_hyper), ])
  
  methyl_tcga_scaleT <- scale(t(methyl_tcga_normal), center = TRUE)
  train_mean <- attr(methyl_tcga_scaleT, "scaled:center")
  train_sd <- attr(methyl_tcga_scaleT, "scaled:scale")
  methyl_tcga_scale <- as.data.frame(t(methyl_tcga_scaleT))
  
  methyl_cptac_normal <- methyl_cptac[rownames(methyl_cptac) %in% rownames(dmc_hypo_hyper),]
  methyl_cptac_normal[is.na(methyl_cptac_normal)] <- 0.00000001

  methyl_cptac_scale <- as.data.frame(t(scale(t(methyl_cptac_normal), center = train_mean,scale = train_sd)))
  methyl_cptac_scale <- methyl_cptac_scale[, !duplicated(colnames(methyl_cptac_scale))]
  
  save(methyl_tcga_scale, file = "RESULTS\\DATAFRAME\\methyl_tcga_scale.Rda")
  save(methyl_cptac_scale, file = "RESULTS\\DATAFRAME\\methyl_cptac_scale.Rda")
  save(train_mean, file = "RESULTS\\DATAFRAME\\train_mean_zScore.Rda")
  save(train_sd, file = "RESULTS\\DATAFRAME\\train_sd_zScore.Rda")
} else{
  load(file = "RESULTS\\DATAFRAME\\methyl_tcga_scale.Rda")
  load(file = "RESULTS\\DATAFRAME\\methyl_cptac_scale.Rda")
  load(file = "RESULTS\\DATAFRAME\\train_mean_zScore.Rda")
  load(file = "RESULTS\\DATAFRAME\\train_sd_zScore.Rda")
}

