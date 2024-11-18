#################### PREPARAZIONE DEI DATASET ##################################
#################          TCGA-PAAD           #################################

## LIBRERIE
library(TCGAbiolinks)
library(SummarizedExperiment)
library(sesame)
library(MASS)
library(tidyverse)

setwd("C:\\Users\\valen\\OneDrive - Politecnico di Bari\\Mongelli\\Risultati\\Progetto_revisione")

## FUNZIONI
source("FUNCTIONS\\solution_prepare_SNV.R")

## FLAG
flag_query_SNV = FALSE
flag_query_METH = FALSE
flag_download = FALSE
 

## pazienti con dati clinici
## pazienti con dati di metilazione
## pazienti con dati di mutazione

##### DATI CLINICI #####
clinical_tcga <- GDCquery_clinic("TCGA-PAAD")
clinical_tcga <- as.data.frame(clinical_tcga)

# PREPARAZIONE DEL DATASET PER ANALISI OVERALL SURVIVAL

## (a) Check sulle colonne utili per valutare "Overall Survival"
## (b) Conversione "age" da giorni ad anni
## (c) Conversione overall survival da giorni a mesi
## (d) Filtraggio dei pazienti: si considerano solo i pazienti con "Infiltrating duct carcinoma"


# (a) Check sulle colonne utili per valutare "Overall Survival"
clinical_tcga <- clinical_tcga[clinical_tcga$vital_status != "Not Reported",c("project","submitter_id", # ID PROGETTO E PAZIENTE
                                                                              "vital_status","days_to_last_follow_up","days_to_death", # colonne x OVERALL SURVIVAL
                                                                              "gender","age_at_diagnosis",
                                                                              "ajcc_pathologic_stage","ajcc_pathologic_t","ajcc_pathologic_n","ajcc_pathologic_m",
                                                                              "primary_diagnosis")]

# (b) Conversione "age" da giorni ad anni
clinical_tcga$age_at_diagnosis <- round(as.numeric(clinical_tcga$age_at_diagnosis)/365)

# (c) Conversione overall survival da giorni a mesi
clinical_tcga$days_to_last_follow_up <- round(as.numeric(clinical_tcga$days_to_last_follow_up)/30.42)
clinical_tcga$days_to_death <- round(as.numeric(clinical_tcga$days_to_death)/30.42)
clinical_tcga$overall_survival <- ifelse(clinical_tcga$vital_status == "Alive",
                                         clinical_tcga$days_to_last_follow_up,
                                         clinical_tcga$days_to_death)

# (d) Filtraggio dei pazienti: si considerano solo i pazienti con "Infiltrating duct carcinoma"
clinical_tcga <- clinical_tcga[clinical_tcga$primary_diagnosis == "Infiltrating duct carcinoma, NOS",]


case_id_clinical <- clinical_tcga$submitter_id


#####

#### DATI MUTAZIONI ####

##### SCARICAMENTO E PREPARAZIONE DATI SNP #####

# Creazione della QUERY
if (flag_query_SNV == TRUE){
  SNV_query_tcga <- GDCquery(project = "TCGA-PAAD",
                             data.category = "Simple Nucleotide Variation",
                             data.format = "maf",
                             access = "open")
  
  output_snv_tcga <- SNV_query_tcga[[1]][[1]]
  # Download file SNV
  if (flag_download == TRUE){
    GDCdownload(query_SNV_tcga, files.per.chunk = 1)
  }
  
  
  # Soluzione bug GDCprepare
  # (ci sono alcuni file maf vuoti o con le colonne "Tumor_Seq_Allele2" NaN)
  snv_tcga <- solution_prepare_SNV(project_gdc = "TCGA-PAAD", 
                                   flag_solution_prepare = "TRUE")
  
  #snv_tcga$Tumor_Sample_Barcode <- substr(snv_tcga$Tumor_Sample_Barcode, 1, 12)
  save(snv_tcga, file="RESULTS\\DATAFRAME\\SNV_tcga_maf.Rda")
  
} else {
  
  load(file="RESULTS\\DATAFRAME\\SNV_tcga_maf.Rda")
  
}

case_id_snv <- unique(substr(snv_tcga$Tumor_Sample_Barcode, 1, 12))


#####

#### DATI METILAZIONE ####

##### SCARICAMENTO E PREPARAZIONE DATI METILAZIONE #####

if (flag_query_METH==TRUE){
  METH_query_tcga <- GDCquery(
    project = c("TCGA-PAAD"),
    data.category = "DNA Methylation",
    platform = "Illumina Human Methylation 450",
    data.type = "Methylation Beta Value",
    sample.type = "Primary Tumor")
  
  
  output_meth_tcga <- getResults(METH_query_tcga)
  
  # Download FILE
  if (flag_download==TRUE){
    GDCdownload(METH_query_tcga)
  }
  
  # Data PREAPARE
  ### TUMOR
  tcga_data <- GDCprepare(METH_query_tcga, summarizedExperiment = TRUE)
  
  methyl_mat <- assay(tcga_data)
  methyl_tcga <- as.data.frame(methyl_mat)
  
  saveRDS(tcga_data, file="RESULTS/DATAFRAME/Meth_tcga_data.rds")
  save(methyl_tcga, file = "RESULTS\\DATAFRAME\\MethDf_tcga.Rda")
} else{
  tcga_data <- readRDS(file="RESULTS/DATAFRAME/Meth_tcga_data.rds")
  load(file = "RESULTS\\DATAFRAME\\MethDf_tcga.Rda")
}

colnames(methyl_tcga) <- substr(colnames(methyl_tcga),1,12)
case_id_meth <- unique(colnames(methyl_tcga))



###### TOTALE CASES ######

case_id_tcga <- intersect(intersect(case_id_clinical, case_id_snv), case_id_meth)
case_id_tcga <- case_id_tcga[case_id_tcga != "TCGA-IB-7651"] #N.B. tolgo questo paziente perchè ha più di 8000 SNV

## Ricalcolo tutti i df considerando solo i case_id_tcga

clinical_tcga <- clinical_tcga[clinical_tcga$submitter_id %in% case_id_tcga,]

snv_tcga <- snv_tcga[substr(snv_tcga$Tumor_Sample_Barcode, 1, 12) %in% case_id_tcga,]
snv_tcga$Tumor_Sample_Barcode <- substr(snv_tcga$Tumor_Sample_Barcode, 1, 12)

tcga_data <- tcga_data[, substr(rownames(colData(tcga_data)), 1, 12) %in% case_id_tcga]
rownames(colData(tcga_data)) <- substr(rownames(colData(tcga_data)), 1, 12)

methyl_mat <- assay(tcga_data)
methyl_tcga <- as.data.frame(methyl_mat)

#### SALVATAGGIO DI TUTTI I DATASET UTILI
save(case_id_tcga, file = "RESULTS\\DATAFRAME\\Case_id_tcga.Rda")

save(clinical_tcga, file = "RESULTS\\DATAFRAME\\Clinical_tcga.Rda")

save(snv_tcga, file="RESULTS\\DATAFRAME\\SNV_tcga_maf.Rda")

saveRDS(tcga_data, file="RESULTS/DATAFRAME/Meth_tcga_data.rds")
save(methyl_tcga, file = "RESULTS\\DATAFRAME\\MethDf_tcga.Rda")





