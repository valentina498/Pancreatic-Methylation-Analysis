#################### PREPARAZIONE DEI DATASET ##################################
#################          CPTAC-3             #################################

## LIBRERIE
library(TCGAbiolinks)
library(SummarizedExperiment)
library(sesame)
library(MASS)
library(tidyverse)

setwd("C:\\Users\\valen\\OneDrive - Politecnico di Bari\\Mongelli\\Risultati\\Progetto_revisione")
## UTILS
aliquot <- read.csv(file = "UTILS\\aliquot_cptac_last.tsv", sep="\t")
sample <- read.csv(file = "UTILS\\sample_cptac_last.tsv", sep="\t")

df_case_aliquot <- aliquot[,c("case_submitter_id","sample_submitter_id","aliquot_submitter_id")]
df_case_sample <- sample[sample$biospecimen_anatomic_site == "Pancreas",c("case_submitter_id","sample_submitter_id","sample_type")]

df_aliquot_sample <- merge(df_case_aliquot, df_case_sample, by.df_case_aliquot = sample_submitter_id)
df_aliquot_sample <- df_aliquot_sample[df_aliquot_sample$sample_type == "Primary Tumor", ]
df_aliquot_sample$aliquot_submitter_id <- substr(df_aliquot_sample$aliquot_submitter_id, 1, 13)

## FLAG
flag_query_clinic <- FALSE
flag_query_SNV = FALSE
flag_query_METH = FALSE
flag_download = FALSE

## pazienti con dati clinici
## pazienti con dati di metilazione
## pazienti con dati di mutazione

#### DATI CLINICI ####
if(flag_query_clinic == TRUE){
  clinical_cptac <- GDCquery_clinic("CPTAC-3")
  clinical_cptac <- clinical_cptac[clinical_cptac$primary_site == "Pancreas",]
  save(clinical_cptac, file = "RESULTS\\DATAFRAME\\clinical_cptac_fromQuery.Rda")
} else{
  load(file = "RESULTS\\DATAFRAME\\clinical_cptac_fromQuery.Rda")
}

# PREPARAZIONE DEL DATASET PER ANALISI OVERALL SURVIVAL

## (a) Check sulle colonne utili per valutare "Overall Survival"
## (b) Conversione "age" da giorni ad anni
## (c) Conversione overall survival da giorni a mesi

# (a) Check sulle colonne utili per valutare "Overall Survival"
clinical_cptac <- clinical_cptac[clinical_cptac$vital_status != "Not Reported",c("project","submitter_id", # ID PROGETTO E PAZIENTE
                                                                                 "vital_status","days_to_last_follow_up","days_to_death", # colonne x OVERALL SURVIVAL
                                                                                 "gender","age_at_diagnosis",
                                                                                 "ajcc_pathologic_stage","ajcc_pathologic_t","ajcc_pathologic_n","ajcc_pathologic_m",
                                                                                 "primary_diagnosis")]

# (b) Conversione "age" da giorni ad anni
clinical_cptac$age_at_diagnosis <- round(as.numeric(clinical_cptac$age_at_diagnosis)/365)

# (c) Conversione overall survival da giorni a mesi
clinical_cptac$days_to_last_follow_up <- round(as.numeric(clinical_cptac$days_to_last_follow_up)/30.42)
clinical_cptac$days_to_death <- round(as.numeric(clinical_cptac$days_to_death)/30.42)
clinical_cptac$overall_survival <- ifelse(clinical_cptac$vital_status == "Alive",
                                          clinical_cptac$days_to_last_follow_up,
                                          clinical_cptac$days_to_death)

case_id_clinical <- clinical_cptac$submitter_id


#### DATI MUTAZIONI ####

ids <- unique(df_aliquot_sample$aliquot_submitter_id)
###

if (flag_query_SNV == TRUE){
  SNV_query_cptac <- GDCquery(
    project = "CPTAC-3", 
    data.category = "Simple Nucleotide Variation", 
    data.format = "maf",
    access = "open",
    barcode = ids)
  output_snv_cptac <- SNV_query_cptac[[1]][[1]]
  
  # Download dei file
  if (flag_download == TRUE){
    GDCdownload(SNV_query_cptac, files.per.chunk = 1)
  }
  
  snv_cptac <- GDCprepare(SNV_query_cptac, summarizedExperiment = TRUE)
  
  snv_cptac$Tumor_Sample_Barcode <- substr(snv_cptac$Tumor_Sample_Barcode, 1, 13)
  
  for (i in (1:nrow(df_aliquot_sample))){
    idx <- which(snv_cptac$Tumor_Sample_Barcode == df_aliquot_sample$aliquot_submitter_id[i])
    snv_cptac$Tumor_Sample_Barcode[idx] <- df_aliquot_sample$case_submitter_id[i]
  }
  
  save(snv_cptac, file = "RESULTS\\DATAFRAME\\SNV_cptac_maf.Rda")
} else{
  load(file = "RESULTS\\DATAFRAME\\SNV_cptac_maf.Rda")
}

case_id_snv <- unique(snv_cptac$Tumor_Sample_Barcode)


#### DATI METILAZIONE ####

if(flag_query_METH == TRUE){
  METH_query_cptac <- GDCquery(project = "CPTAC-3",
                               data.category = "DNA Methylation",
                               data.type = "Methylation Beta Value",
                               platform = "Illumina Methylation Epic",
                               access = "open",
                               barcode = ids)
  output_meth_cptac <- getResults(METH_query_cptac)
  output_meth_cptac <- subset(output_meth_cptac, 
                              !(output_meth_cptac$sample_type %in% c("Solid Tissue Normal",
                                                                     "Solid Tissue Normal;Solid Tissue Normal")))
  METH_query_cptac[[1]][[1]] <- output_meth_cptac
  # Download FILE
  if (flag_download == TRUE){
    GDCdownload(METH_query_cptac)
  }
  
  cptac_data <- GDCprepare(METH_query_cptac, summarizedExperiment = TRUE)
  
  methyl_mat <- assay(cptac_data)
  methyl_cptac <- as.data.frame(methyl_mat)
  colnames(methyl_cptac) <- substr(colnames(methyl_cptac), 1, 13)
  rownames(colData(cptac_data)) <- substr(rownames(colData(cptac_data)), 1, 13)
  
  for (i in (1:nrow(df_aliquot_sample))){
    idx <- which(colnames(methyl_cptac) == df_aliquot_sample$aliquot_submitter_id[i])
    colnames(methyl_cptac)[idx] <- df_aliquot_sample$case_submitter_id[i]
    
    idx_data <- which(rownames(colData(cptac_data)) == df_aliquot_sample$aliquot_submitter_id[i])
    rownames(colData(cptac_data))[idx_data] <- df_aliquot_sample$case_submitter_id[i]
  }
  
  saveRDS(cptac_data, file="RESULTS/DATAFRAME/Meth_cptac_data.rds")
  save(methyl_cptac, file = "RESULTS\\DATAFRAME\\MethDf_cptac.Rda")
} else {
  cptac_data <- readRDS(file="RESULTS/DATAFRAME/Meth_cptac_data.rds")
  load(file = "RESULTS\\DATAFRAME\\MethDf_cptac.Rda")
}

case_id_meth <- unique(colnames(methyl_cptac))


###### TOTALE CASES ######

case_id_cptac <- intersect(intersect(case_id_clinical, case_id_snv), case_id_meth)


## Ricalcolo tutti i df considerando solo i case_id_cptac

clinical_cptac <- clinical_cptac[clinical_cptac$submitter_id %in% case_id_cptac,]

snv_cptac <- snv_cptac[snv_cptac$Tumor_Sample_Barcode %in% case_id_cptac,]


cptac_data <- cptac_data[, rownames(colData(cptac_data), 1, 12) %in% case_id_cptac]
methyl_mat <- assay(cptac_data)
methyl_cptac <- as.data.frame(methyl_mat)

#### SALVATAGGIO DI TUTTI I DATASET UTILI
save(case_id_cptac, file = "RESULTS\\DATAFRAME\\Case_id_cptac.Rda")

save(clinical_cptac, file = "RESULTS\\DATAFRAME\\Clinical_cptac.Rda")

save(snv_cptac, file="RESULTS\\DATAFRAME\\SNV_cptac_maf.Rda")

saveRDS(cptac_data, file="RESULTS/DATAFRAME/Meth_cptac_data.rds")
save(methyl_cptac, file = "RESULTS\\DATAFRAME\\MethDf_cptac.Rda")
