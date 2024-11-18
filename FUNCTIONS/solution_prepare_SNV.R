solution_prepare_SNV <- function(project_gdc, flag_solution_prepare){
  
  #FUNZIONE PER RISOLVERE IL PROBLEMA NELLA PREPARAZIONE DEL FILE MAF IN TCGA
  
  if (flag_solution_prepare==TRUE){
    query_SNV_prova <- GDCquery(project = project_gdc,
                                data.category = "Simple Nucleotide Variation",
                                data.format = "maf",
                                access = "open")
    output_query_prova <- query_SNV_prova[[1]][[1]]
    
    output_query_temp <- output_query_prova[1,]
    query_SNV_prova[[1]][[1]] <- output_query_temp
    
    
    
    maf_total <- GDCprepare(query_SNV_prova, summarizedExperiment = FALSE)
    
    
    
    case_null <- c()
    for (i in (2:180)){
      output_query_temp <- output_query_prova[i,]
      query_SNV_prova[[1]][[1]] <- output_query_temp
      
      maf <- GDCprepare(query_SNV_prova, summarizedExperiment = TRUE)
      
      if (ncol(maf)==140){
        maf_total <- merge(maf_total, maf, all = TRUE)
      } else {
        case_null <- c(case_null, output_query_temp$cases[i])
      }
      print(i)
    }
    
    #save(maf_total, file = paste0(directory,"\\SNV_tcga_maf.Rda"))
  }
  
  #load(file = paste0(directory,"\\SNV_tcga_maf.Rda"))
  
  return(maf_total)
}