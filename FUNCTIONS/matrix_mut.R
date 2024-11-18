matrix_mut <- function(snv_total, clinical, N){
  
  gene_unique <- unique(snv_total$Hugo_Symbol)      # Nomi COLONNE
  cases <- unique(snv_total$Tumor_Sample_Barcode) # Nomi RIGHE
  
  mutation_matrix <- matrix(0, nrow=length(cases), ncol=length(gene_unique))
  rownames(mutation_matrix) <- cases
  colnames(mutation_matrix) <- gene_unique
  mutation_df <- as.data.frame(mutation_matrix)
  
  for (i in (1:length(cases))){
    maf_case <- snv_total[snv_total$Tumor_Sample_Barcode==cases[i],]
    
    gene_case <- unique(maf_case$Hugo_Symbol)
    mutation_df[cases[i],gene_case] <- 1

  }
  
  mafr <- read.maf(maf = snv_total)
  freq_mut <- getGeneSummary(mafr)
  top_mut <- freq_mut$Hugo_Symbol[1:N]
  clinical <- clinical[clinical$submitter_id %in% cases,]
  
  for (i in (1:20)){
    gene <- top_mut[i]
    clinical[gene] <- NA
    
    for (j in (1:length(cases))){
      row <- which(clinical$submitter_id %in% cases[j])
      clinical[row, gene] <- mutation_df[cases[j], gene]
    }
  }
  
  risultati <- list(mutMatrix = mutation_df, clinical_geneStrata = clinical)
  return(risultati)
}