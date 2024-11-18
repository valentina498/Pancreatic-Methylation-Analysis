comb_mutation <- function(clinical){
  
  ## AND a coppie
  # mutation KRAS AND TP53
  clinical[,"KRAS_and_TP53"] <- ifelse(clinical$KRAS == 1 & 
                                 clinical$TP53 == 1, 1, 0)
  # mutation KRAS AND CDKN2A
  clinical[,"KRAS_and_CDKN2A"] <- ifelse(clinical$KRAS == 1 &
                                 clinical$CDKN2A == 1, 1, 0)
  # mutation TP53 AND CDKN2A, wild type KRAS
  clinical[,"TP53_and_CDKN2A"] <- ifelse(clinical$TP53 == 1 &
                                 clinical$CDKN2A == 1, 1, 0)
  
  ## OR a coppie
  # mutation KRAS AND TP53
  clinical[,"KRAS_or_TP53"] <- ifelse(clinical$KRAS == 1 | 
                                         clinical$TP53 == 1, 1, 0)
  # mutation KRAS AND CDKN2A
  clinical[,"KRAS_or_CDKN2A"] <- ifelse(clinical$KRAS == 1 |
                                           clinical$CDKN2A == 1, 1, 0)
  # mutation TP53 AND CDKN2A, wild type KRAS
  clinical[,"TP53_or_CDKN2A"] <- ifelse(clinical$TP53 == 1 |
                                           clinical$CDKN2A == 1, 1, 0)
  

  #### AND a tre
  clinical[,"KRAS_and_TP53_and_CDKN2A"] <- ifelse(clinical$KRAS == 1 & 
                                   clinical$TP53 == 1 &
                                   clinical$CDKN2A == 1, 1, 0)
  
  #### OR a tre
  clinical[,"KRAS_or_TP53_or_CDKN2A"] <- ifelse(clinical$KRAS == 1 | 
                                  clinical$TP53 == 1 |
                                  clinical$CDKN2A == 1, 1, 0)
  
  return(clinical)
  
}