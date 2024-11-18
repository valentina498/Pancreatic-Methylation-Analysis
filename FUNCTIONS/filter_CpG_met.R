filter_CpG_met <- function(data, probe_merge, new_col, clinical_geneStrat){
  
  # a) Filtro sulle sonde: si prendono solo le sonde selezionate secondo i criteri:
  #        -- sonde in comune tra Illumina EPIC e 450K
  #        -- esclusione delle sonde sui cromosomi sessuali X e Y
  # b) Filtro sulle sonde che assumono valori nulli in più del 30% dei campioni 
  #    (sono mantenute solo le sonde con un numero di valori mancanti NA minore del 30%)
  # c) Sostituzione dei NA con 0
  
  feature_criterion1 <- rownames(assay(data)) %in% probe_merge
  data_filtered <- data[feature_criterion1, ]
  
  feature_criterion2 <- rowSums(is.na(assay(data_filtered))) <= 40 #30% dei pazienti
  data_filtered <- data_filtered[feature_criterion2, ]
  
  
  assay(data_filtered)[is.na(assay(data_filtered))] <- 0.00000001
  
  for (j in (1:length(new_col))){
    for (i in (1:nrow(clinical_geneStrat))){
      data_filtered@colData[i,new_col[j]] <- clinical_geneStrat[clinical_geneStrat$submitter_id == rownames(data_filtered@colData)[i], new_col[j]]
    }
  }
  return(data_filtered)
}