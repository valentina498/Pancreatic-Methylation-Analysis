calc_gene_sig <- function(top_genes_tot, fit_tcga, clinical_tcga, flag_save_img, p_val_thresh){
  
  ## proc_choice = "A" -> calcola i geni mutati significativi
  ## proc_choice = "B" -> calcola le combinazioni di mutazioni significative
  
  gene_sig_tcga <- c() #vettore che conterrà i geni che sono significativi per OS in tcga
  gene_pval <- c()
  xlim = c(0,84)
  for (i in (1:length(top_genes_tot))){
    gene <- top_genes_tot[i]
    # plot1 <- ggsurvplot(fit_cptac,
    #                     data = clinical_cptac,
    #                     pval = T,
    #                     
    #                     risk.table = T,
    #                     risk.table.y.text = FALSE,
    #                     xlab ="Time in Months",
    #                     break.time.by = 12,
    #                     legend.labs =
    #                       c("WildType","Mutant"),
    #                     xlim = xlim,
    #                     
    #                     title = "CPTAC-3", subtitle = paste0("Mutated: ",gene),
    #                     font.title = c(16, "bold", "darkblue"),
    #                     font.subtitle = c(11, "bold.italic", "deepskyblue"))
    # 
    # 
    # plot2 <- ggsurvplot(fit_tcga,
    #                     data = clinical_tcga,
    #                     pval = T,
    #                     
    #                     risk.table = T,
    #                     risk.table.y.text = FALSE,
    #                     xlab ="Time in Months",
    #                     break.time.by = 12,
    #                     legend.labs =
    #                       c("WildType","Mutant"),
    #                     xlim = xlim,
    #                     
    #                     title = "TCGA-PAAD", subtitle = paste0("Mutated: ",gene),
    #                     font.title = c(16, "bold", "darkblue"),
    #                     font.subtitle = c(11, "bold.italic", "deepskyblue"))
    # 
    # plot <- arrange_ggsurvplots(list(plot2,plot1), print = TRUE,
    #                             ncol = 2, nrow = 1, risk.table.height = 0.3)
    # 
    # if (flag_save_img == TRUE){
    #   ggsave(paste0("RESULTS\\IMMAGINI\\Mut_KM\\",gene,"_KM.png"), plot,
    #          width = 12,
    #          height = 7)
    # }
    log_rank <- survdiff(Surv(overall_survival, deceased) ~ clinical_tcga[,gene], clinical_tcga)
    
    if (log_rank[["pvalue"]] < p_val_thresh){
      gene_sig_tcga <- c(gene_sig_tcga, gene)
      #gene_pval <- c(gene_pval, log_rank[["pvalue"]])
    }
    
    gene_pval <- c(gene_pval, paste(gene,"=",log_rank[["pvalue"]]))
    
  } 
  res <- list(gene_sig = gene_sig_tcga, gene_pval = gene_pval)

  return(res)
  
}