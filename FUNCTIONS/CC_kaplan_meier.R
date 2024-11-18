CC_kaplan_meier <- function(test_comb, v_clusterAlg, v_distance, clinical, 
                            fit_cc3_normal, fit_cc4_normal, fit_cc3_scale, fit_cc4_scale, 
                            name_folder_normal, name_folder_scale,flag_save){
  
  plot1 <- ggsurvplot(fit_cc3_normal,
                      data = clinical,
                      #pval = T,
                      
                      risk.table = T,
                      risk.table.y.text = FALSE,
                      xlab ="Time in Months",
                      break.time.by = 12,
                      legend.labs =
                        c("Class 1","Class 2","Class 3"),
                      xlim = c(0,84),
                      
                      title = "TCGA-PAAD", 
                      subtitle = paste("CC3 normal",v_clusterAlg,"(",v_distance,")"),
                      font.title = c(16, "bold", "darkblue"),
                      font.subtitle = c(11, "bold.italic", "deepskyblue"))
  
  plot2 <- ggsurvplot(fit_cc4_normal,
                      data = clinical,
                      #pval = T,
                      
                      risk.table = T,
                      risk.table.y.text = FALSE,
                      xlab ="Time in Months",
                      break.time.by = 12,
                      legend.labs =
                        c("Class 1","Class 2","Class 3","Class 4"),
                      xlim = c(0,84),
                      
                      title = "TCGA-PAAD", 
                      subtitle = paste("CC4 normal",v_clusterAlg,"(",v_distance,")"),
                      font.title = c(16, "bold", "darkblue"),
                      font.subtitle = c(11, "bold.italic", "deepskyblue"))
  
  plot <- arrange_ggsurvplots(list(plot1,plot2), print = TRUE,
                              ncol = 2, risk.table.height = 0.3)
  
  name_folder <- paste0(test_comb[k],"\\cons_normal_",v_clusterAlg,"_",v_distance,"\\")
  if (flag_save == TRUE){
    ggsave(paste0("RESULTS\\IMMAGINI\\ConsensusClustering\\",name_folder_normal,"CC_",v_clusterAlg,"_",v_distance,".png"), plot,
           width = 12,
           height = 8)
  }
  
  
  plot3 <- ggsurvplot(fit_cc3_scale,
                      data = clinical,
                      #pval = T,
                      
                      risk.table = T,
                      risk.table.y.text = FALSE,
                      xlab ="Time in Months",
                      break.time.by = 12,
                      legend.labs =
                        c("Class 1","Class 2","Class 3"),
                      xlim = c(0,84),
                      
                      title = "TCGA-PAAD", 
                      subtitle = paste("CC3 scale",v_clusterAlg,"(",v_distance,")"),
                      font.title = c(16, "bold", "darkblue"),
                      font.subtitle = c(11, "bold.italic", "deepskyblue"))
  
  plot4 <- ggsurvplot(fit_cc4_scale,
                      data = clinical,
                      #pval = T,
                      
                      risk.table = T,
                      risk.table.y.text = FALSE,
                      xlab ="Time in Months",
                      break.time.by = 12,
                      legend.labs =
                        c("Class 1","Class 2","Class 3","Class 4"),
                      xlim = c(0,84),
                      
                      title = "TCGA-PAAD", 
                      subtitle = paste("CC4 scale",v_clusterAlg,"(",v_distance,")"),
                      font.title = c(16, "bold", "darkblue"),
                      font.subtitle = c(11, "bold.italic", "deepskyblue"))
  
  plot <- arrange_ggsurvplots(list(plot3,plot4), print = TRUE,
                              ncol = 2, risk.table.height = 0.3)
  
  name_folder <- paste0(test_comb[k],"\\cons_scale_",v_clusterAlg,"_",v_distance,"\\")
  if (flag_save == TRUE){
    ggsave(paste0("RESULTS\\IMMAGINI\\ConsensusClustering\\",name_folder_scale,"CC_",v_clusterAlg,"_",v_distance,".png"), plot,
           width = 12,
           height = 8)
  }
  
  
}

CC_kaplan_meier_REV <- function(test_comb, v_clusterAlg, v_distance, clinical, 
                            fit_cc3_normal, fit_cc4_normal, 
                            name_folder_normal,flag_save){
  
  plot1 <- ggsurvplot(fit_cc3_normal,
                      data = clinical,
                      #pval = T,
                      
                      risk.table = T,
                      risk.table.y.text = FALSE,
                      xlab ="Time in Months",
                      break.time.by = 12,
                      legend.labs =
                        c("Class 1","Class 2","Class 3"),
                      xlim = c(0,84),
                      
                      title = "TCGA-PAAD", 
                      subtitle = paste("CC3",v_clusterAlg,"(",v_distance,")"),
                      font.title = c(16, "bold", "darkblue"),
                      font.subtitle = c(11, "bold.italic", "deepskyblue"))
  
  plot2 <- ggsurvplot(fit_cc4_normal,
                      data = clinical,
                      #pval = T,
                      
                      risk.table = T,
                      risk.table.y.text = FALSE,
                      xlab ="Time in Months",
                      break.time.by = 12,
                      legend.labs =
                        c("Class 1","Class 2","Class 3","Class 4"),
                      xlim = c(0,84),
                      
                      title = "TCGA-PAAD", 
                      subtitle = paste("CC4 ",v_clusterAlg,"(",v_distance,")"),
                      font.title = c(16, "bold", "darkblue"),
                      font.subtitle = c(11, "bold.italic", "deepskyblue"))
  
  plot <- arrange_ggsurvplots(list(plot1,plot2), print = TRUE,
                              ncol = 2, risk.table.height = 0.3)
  
  name_folder <- paste0(test_comb,"\\cons_",v_clusterAlg,"_",v_distance,"\\")
  if (flag_save == TRUE){
    ggsave(paste0("RESULTS\\IMMAGINI\\ConsensusClustering\\",name_folder,"CC_",v_clusterAlg,"_",v_distance,".png"), plot,
           width = 12,
           height = 8)
  }
  
}