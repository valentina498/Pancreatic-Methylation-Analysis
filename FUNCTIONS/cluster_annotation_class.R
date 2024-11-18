cluster_annotation_class <- function(nCluster, clinical_res, project, test_comb){
  
  annotation_CC <- c("gender","age_strat","stage_strat","T_strat","N_strat","M_strat","KRAS","TP53",test_comb)
  ann_title <- c("Gender","Age","Tumor stage","T stage","N stage","M stage","KRAS","TP53","KRAS&TP53")
  
  if(project == "TCGA-PAAD"){
    palettes <- list(
      "gender" = c("#2E9FDF", "red"), #gender
      "age_strat" = c("darkorange", "green4"), #age
      "stage_strat" = c("gold", "orange","red1","red4"), #stage
      "T_strat" = brewer.pal(4, "Set1"), #T
      "N_strat" = c("slategray1", "skyblue2","steelblue4"), #N
      "M_strat" = c("plum1", "plum3","plum4"), #M
      "KRAS" = c("skyblue4","skyblue2"),
      "TP53" = c("skyblue4","skyblue2"),
      "KRAS_and_TP53" = c("skyblue4","skyblue2")
    )
    
    
  } else {
    palettes <- list(
      "gender" = c("#2E9FDF", "red"), #gender
      "age_strat" = c("darkorange", "green4"), #age
      "stage_strat" = c("gold", "orange","red1","red4"), #stage
      "T_strat" = brewer.pal(5, "Set1"), #T
      "N_strat" = c("slategray1", "skyblue2","skyblue3","steelblue4"), #N
      "M_strat" = c("plum1", "plum3","plum4"), #M
      "KRAS" = c("skyblue4","skyblue2"),
      "TP53" = c("skyblue4","skyblue2"),
      "KRAS_and_TP53" = c("skyblue4","skyblue2")
    )
  }
  
  
  
  
  ##### PLOT DI ANNOTAZIONE CON I DATI CLINICI
  plots_ann <- list()
  
  for (p in (1:length(annotation_CC))){
    string_ann <- annotation_CC[p]
    col_name = c("cluster_pred")
    
    clinical_res_p = clinical_res[,c(col_name,string_ann)]
    if(p == 9){
      plots_ann[[p]] <- ggplot(clinical_res_p, aes(x = .data[[col_name]], fill = .data[[string_ann]])) +
        geom_bar(position = "fill") +
        scale_fill_manual(values = palettes[[string_ann]],
                          labels = c("Mutant" = "Both Mut", "WT" = "Any WT")
        ) +
        labs(title = ann_title[p],
             x = "Cluster",
             fill = ann_title[p]) +
        theme(
          legend.position = "top",
          plot.title = element_text(size = 24 ,hjust = 0.5),
          legend.title = element_text(size = 16),
          legend.text = element_text(size = 16),
          axis.title = element_text(size = 16),
          axis.text = element_text(size = 16),
          axis.title.y = element_blank()
        ) +
        scale_y_continuous(labels = scales::percent)
    } else{
      plots_ann[[p]] <- ggplot(clinical_res_p, aes(x = .data[[col_name]], fill = .data[[string_ann]])) +
        geom_bar(position = "fill") +
        scale_fill_manual(values = palettes[[string_ann]]
        ) +
        labs(title = ann_title[p],
             x = "Cluster",
             fill = ann_title[p]) +
        theme(
          legend.position = "top",
          plot.title = element_text(size = 24 ,hjust = 0.5),
          legend.title = element_text(size = 16),
          legend.text = element_text(size = 16),
          axis.title = element_text(size = 16),
          axis.text = element_text(size = 16),
          axis.title.y = element_blank()
        ) +
        scale_y_continuous(labels = scales::percent)
    }
    
    
    #print(p)
    
  }
  
  subplot_annotation <- grid.arrange(grobs = plots_ann, ncol = 3, top=paste("Prediction on",project))
  
  
  
  #### PLOT ANNOTAZIONE feature mutation profile 
  
  string_ann <- test_comb
  clinical_res_p = clinical_res[,c(col_name,string_ann)]
  
  plot_comb <- ggplot(clinical_res_p, aes(x = .data[[col_name]], fill = .data[[string_ann]])) +
    geom_bar(position = "fill") +
    scale_fill_manual(values = c("skyblue4","skyblue2")) +
    labs(title = ann_title[p],
         x = "Cluster",
         fill = ann_title[p]) +
    theme(
      legend.position = "top",
      plot.title = element_text(size = 24 ,hjust = 0.5),
      legend.title = element_text(size = 16),
      legend.text = element_text(size = 16),
      axis.title = element_text(size = 16),
      axis.text = element_text(size = 16),
      axis.title.y = element_blank()
    ) +
    scale_y_continuous(labels = scales::percent)
  
  res <- list(subplot_annotation = subplot_annotation, plot_comb = plot_comb)
  
  return(res)
  
}