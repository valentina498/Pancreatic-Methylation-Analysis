plot_pval_geneComb <- function(clinical_tcga, fit_tcga, clinical_cptac, fit_cptac, 
                               gene_comb, gene_title, flag_save_img, i){
  
  
  plot1 <- ggsurvplot(fit_cptac,
                      data = clinical_cptac,
                      pval = T,
                      
                      risk.table = T,
                      risk.table.y.text = FALSE,
                      xlab ="Time in Months",
                      break.time.by = 12,
                      legend.labs =
                        c("Any WT","Both Mutant"),
                      xlim = c(0,72),
                      
                      title = "CPTAC-3", subtitle = gene_title,
                      font.title = c(20),
                      font.subtitle = c(20))
    
  plot1$plot <- plot1$plot + theme(
    legend.title = element_text(family = "Helvetica",size = 18),  # Dimensione del font del titolo della legenda
    legend.text = element_text(family = "Helvetica",size = 18),   # Dimensione del font del testo della legenda
    axis.title = element_text(family = "Helvetica",size = 18),    # Dimensione del font dei nomi degli assi
    axis.text = element_text(family = "Helvetica",size = 18)      # Dimensione del font dei valori degli assi
  )
  
  plot1$table <- plot1$table + theme(
    axis.title = element_text(family = "Helvetica",size = 18),
    axis.text = element_text(family = "Helvetica",size = 18)
  )
  
  
  plot2 <- ggsurvplot(fit_tcga,
                      data = clinical_tcga,
                      pval = T,
                      
                      risk.table = T,
                      risk.table.y.text = FALSE,
                      xlab ="Time in Months",
                      break.time.by = 12,
                      legend.labs =
                        c("Any WT","Both Mutant"),
                      xlim = c(0,72),
                      title = "TCGA-PAAD", 
                      subtitle = gene_title,
                      font.title = c(20),
                      font.subtitle = c(20))
  
  plot2$plot <- plot2$plot + theme(
    legend.title = element_text(family = "Helvetica",size = 18),  # Dimensione del font del titolo della legenda
    legend.text = element_text(family = "Helvetica",size = 18),   # Dimensione del font del testo della legenda
    axis.title = element_text(family = "Helvetica",size = 18),    # Dimensione del font dei nomi degli assi
    axis.text = element_text(family = "Helvetica",size = 18)      # Dimensione del font dei valori degli assi
  )
  
  plot2$table <- plot2$table + theme(
    axis.title = element_text(family = "Helvetica",size = 18),
    axis.text = element_text(family = "Helvetica",size = 18)
  )
  
  plot <- arrange_ggsurvplots(list(plot2,plot1), print = TRUE,
                              ncol = 2, nrow = 1, risk.table.height = 0.3)
  
  if (flag_save_img == TRUE){
    ggsave(paste0("RESULTS\\IMMAGINI\\Mut_KM\\PLOT KM tutte le combinazioni\\",i,"_",gene_comb,"_KM_modifica.png"), plot,
           width = 12,
           height = 9)
  }
  log_rank <- survdiff(Surv(overall_survival, deceased) ~ clinical_tcga[,gene_comb], data = clinical_tcga)
  
  comb_pval <- log_rank[["pvalue"]]
  #if (log_rank[["pvalue"]] < p_val_thresh){
  #  gene_comb_sig_tcga <- c(gene_comb_sig_tcga, gene_comb)
  #}
  
  return(comb_pval)
}