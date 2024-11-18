cluster_annotation <- function(clinical_res, test_comb,dataType,algCC,distCC, name_plot){
  # annotation_CC <- c("gender","age_strat","stage_strat","T_strat","N_strat","M_strat","KRAS","TP53","CDKN2A",test_comb)
  # ann_title <- c("Gender","Age","Tumor stage","T stage","N stage","M stage","KRAS","TP53","CDKN2A", test_comb)
  # 
  # palettes <- list(
  #   "gender" = c("#2E9FDF", "red"), #gender
  #   "age_strat" = c("darkorange", "green4"), #age
  #   "stage_strat" = c("gold", "orange","red1","red4"), #stage
  #   "T_strat" = brewer.pal(4, "Set1"), #T
  #   "N_strat" = c("slategray1", "skyblue2","steelblue4"), #N
  #   "M_strat" = c("plum1", "plum3","plum4"), #M
  #   "KRAS" = c("skyblue4","skyblue2"),
  #   "TP53" = c("skyblue4","skyblue2"),
  #   #"CDKN2A" = c("skyblue4","skyblue2")
  #   "KRAS_and_TP53"= c("skyblue4","skyblue2")
  # )
  # 
  # n_cc <- c("cc3", "cc4")
  # n_cc_title <- c("K cluster = 3", "K cluster = 4")
  # for (l in (1:length(n_cc))){
  #   plots_ann <- list()
  #   col_name = paste0(n_cc[l],"_",dataType,"_",algCC,"_",distCC)
  #   
  #   ##### PLOT DI ANNOTAZIONE CON I DATI CLINICI
  #   for (p in (1:(length(annotation_CC)-1))){
  #     string_ann <- annotation_CC[p]
  # 
  #     clinical_res_p = clinical_res[,c(col_name,string_ann)]
  #     plots_ann[[p]] <- ggplot(clinical_res_p, aes(x = .data[[col_name]], fill = .data[[string_ann]])) +
  #       geom_bar(position = "fill") +
  #       scale_fill_manual(values = palettes[[string_ann]]) +
  #       labs(title = ann_title[p],
  #            x = "Cluster",
  #            y = "Propotion",
  #            fill = ann_title[p]) +
  #       theme(
  #         plot.title = element_text(face = "bold.italic", size = 14, colour = "deepskyblue")
  #       ) + 
  #       scale_y_continuous(labels = scales::percent)
  #     
  #   }
  # 
  #   subplot_annotation <- grid.arrange(grobs = plots_ann, ncol = 3)
  #   
  #   ggsave(paste0("RESULTS\\IMMAGINI\\ConsensusClustering\\",name_plot,n_cc[l],"_annotation.png"), plot = subplot_annotation,
  #          width = 14,
  #          height = 14)
  #   
  #   ## Annotazione della combinazione di mutazioni scelta con "test_comb"
  #   p = length(annotation_CC)
  #   string_ann <- annotation_CC[p]
  #   clinical_res_p = clinical_res[,c(col_name,string_ann)]
  #   
  #   plot_comb <- ggplot(clinical_res_p, aes(x = .data[[col_name]], fill = .data[[string_ann]])) +
  #     geom_bar(position = "fill") +
  #     scale_fill_manual(values = c("skyblue4","skyblue2")) +
  #     labs(title = ann_title[p],
  #          x = "Cluster",
  #          y = "Propotion",
  #          fill = ann_title[p]) +
  #     theme(
  #       plot.title = element_text(face = "bold.italic", size = 14, colour = "deepskyblue")
  #     ) + 
  #     scale_y_continuous(labels = scales::percent)
  #   
  #   ggsave(paste0("RESULTS\\IMMAGINI\\ConsensusClustering\\",name_plot,n_cc[l],"_mutProfile.png"), plot = plot_comb,
  #          width = 14,
  #          height = 14)
  #   
    
    
  # }
  
  #bar_width = 0.5
  annotation_CC <- c("gender","age_strat","stage_strat","T_strat","N_strat","M_strat","KRAS","TP53",test_comb)
  ann_title <- c("\nGender\n","\nAge\n","\nTumor stage\n","\nT stage\n","\nN stage\n","\nM stage\n","\nKRAS\n","\nTP53\n", "\nKRAS/TP53\n")

  palettes <- list(
    "gender" = c("#2E9FDF", "red"), #gender
    "age_strat" = c("darkorange", "green4"), #age
    "stage_strat" = c("gold", "orange","red1","red4"), #stage
    "T_strat" = brewer.pal(4, "Set1"), #T
    "N_strat" = c("slategray1", "skyblue2","steelblue4"), #N
    "M_strat" = c("plum1", "plum3","plum4"), #M
    "KRAS" = c("skyblue4","skyblue2"),
    "TP53" = c("skyblue4","skyblue2"),
    "KRAS_and_TP53"= c("skyblue4","skyblue2")
  )

  #n_cc <- c("cc3", "cc4")
  n_cc <- c("cc4")
  n_cc_title <- c("K cluster = 4")
  
  
  for (l in (1:length(n_cc))){
    plots_ann <- list()
    col_name = paste0(n_cc[l],"_",algCC,"_",distCC)

    ##### PLOT DI ANNOTAZIONE CON I DATI CLINICI
    for (p in (1:(length(annotation_CC)))){
      string_ann <- annotation_CC[p]

      clinical_res_p = clinical_res[,c(col_name,string_ann)]
      
      if(p == 9){
        plots_ann[[p]] <- ggplot(clinical_res_p, aes(x = .data[[col_name]], fill = .data[[string_ann]])) +
          geom_bar(position = "fill") +
          scale_fill_manual(values = palettes[[string_ann]],
                            labels = c("Mutant" = "Both Mut", "WT" = "Any WT"))+
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
          scale_y_continuous(labels = scales::percent)+
          guides(fill = guide_legend(title = NULL))
      } else{
        plots_ann[[p]] <- ggplot(clinical_res_p, aes(x = .data[[col_name]], fill = .data[[string_ann]])) +
          geom_bar(position = "fill") +
          scale_fill_manual(values = palettes[[string_ann]]) +
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
          scale_y_continuous(labels = scales::percent)+
          guides(fill = guide_legend(title = NULL))
      }
      

    }

    subplot_annotation <- grid.arrange(grobs = plots_ann, ncol = 3)

    
  }
  return(subplot_annotation)
  
}


