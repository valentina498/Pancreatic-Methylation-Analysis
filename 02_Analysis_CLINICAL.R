################################################################################
######################### ANALISI DEI DATI CLINICI #############################
########################### TCGA-PAAD vs CPTAC-3 ###############################

library(ggplot2)
library(gridExtra)
library(reshape2)
library("survminer")
library(survival)
library(RColorBrewer)

setwd("C:\\Users\\valen\\OneDrive - Politecnico di Bari\\Mongelli\\Risultati\\Progetto_revisione")

## UTILS
load(file="./RESULTS/DATAFRAME/Clinical_cptac.Rda")
load(file="./RESULTS/DATAFRAME/Clinical_tcga.Rda")

clinical_cptac$deceased <- ifelse(clinical_cptac$vital_status == "Dead", TRUE, FALSE)
clinical_tcga$deceased <- ifelse(clinical_tcga$vital_status == "Dead", TRUE, FALSE)

## FLAG
flag_save = TRUE        #salvataggio DATAFRAME CON STRATIFICAZIONE
flag_save_image = TRUE #salvataggio immagini FREQUENZA
flag_save_km = TRUE    #salvataggio immagini KAPLAN-MEIER

## (A) STATISTICA DESCRITTIVA ##
################################
# 1) GENDER
# 2) AGE at Diagnosis
# 3) TUMOR STAGE AJCC
# 4) TUMOR STAGE AJCC-T
# 5) TUMOR STAGE AJCC-N
# 6) TUMOR STAGE AJCC-M
# 7) Primary Diagnosis
theme_set(theme_gray())
### (1) GENDER ###
gender <- data.frame(gender = c("Male", "Female"),
                     TCGA = c(sum(clinical_tcga$gender=="male"),sum(clinical_tcga$gender=="female")),
                     CPTAC = c(sum(clinical_cptac$gender=="male"),sum(clinical_cptac$gender=="female"))) %>%
  melt(id.vars="gender",variable.name="Project")

plot <- ggplot(gender, aes(x=Project,y=value,fill=gender)) +
  geom_col(position=position_dodge())+
  labs(title = "Gender",
       fill = "Gender")+
  scale_fill_manual(values = c("red","#2E9FDF"))+
  theme(
    plot.title = element_text(face = "bold.italic", size = 18, colour = "deepskyblue"),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14)
  )
plot


if (flag_save_image == TRUE){
  ggsave("./RESULTS/IMMAGINI/Clinical_freq/Gender_freq1.png", plot,
         width = 7,
         height = 7)
}

### (2) AGE at Diagnosis ###

age_cptac <- ggplot(clinical_cptac, aes(y = age_at_diagnosis)) +
  geom_boxplot() +
  ylim(30, 90) +
  labs(title = "Age at diagnosis (CPTAC-3)",
       y = "Age")+
  theme(
    plot.title = element_text(face = "bold.italic", size = 18, colour = "deepskyblue"),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14)
  )


age_tcga <- ggplot(clinical_tcga, aes(y = age_at_diagnosis)) +
  geom_boxplot() +
  ylim(30, 90) +
  labs(title = "Age at diagnosis (TCGA)",
       y = "Age")+
  theme(
    plot.title = element_text(face = "bold.italic", size = 18, colour = "deepskyblue"),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14)
  )

subplot <- grid.arrange(age_tcga, age_cptac, ncol = 2)

boxplot_stats_cptac <- boxplot.stats(clinical_cptac$age_at_diagnosis)
boxplot_stats_tcga <- boxplot.stats(clinical_tcga$age_at_diagnosis)

if (flag_save_image == TRUE){
  ggsave("RESULTS\\IMMAGINI\\Clinical_freq\\Age_freq.png", subplot,
         width = 9,
         height = 7)
  
}

### (3) TUMOR STAGE AJCC ###
stage_ajcc <- data.frame(stage_ajcc = c("Stage I", "Stage II","Stage III", "Stage IV"),
                         TCGA = c(sum(clinical_tcga$ajcc_pathologic_stage %in% c("Stage I","Stage IA","Stage IB")),
                                  sum(clinical_tcga$ajcc_pathologic_stage %in% c("Stage IIA","Stage IIB")),
                                  sum(clinical_tcga$ajcc_pathologic_stage %in% c("Stage III")),
                                  sum(clinical_tcga$ajcc_pathologic_stage %in% c("Stage IV"))),
                         CPTAC = c(sum(clinical_cptac$ajcc_pathologic_stage %in% c("Stage IA","Stage IB")),
                                   sum(clinical_cptac$ajcc_pathologic_stage %in% c("Stage IIA","Stage IIB")),
                                   sum(clinical_cptac$ajcc_pathologic_stage %in% c("Stage III")),
                                   sum(clinical_cptac$ajcc_pathologic_stage %in% c("Stage IV")))) %>%
  melt(id.vars="stage_ajcc",variable.name="Project")

plot <- ggplot(stage_ajcc, aes(x=Project,y=value,fill=stage_ajcc)) +
  geom_col(position=position_dodge())+
  scale_fill_manual(values = c("gold", "orange","red1","red4"))+
  theme(
    plot.title = element_text(face = "bold.italic", size = 18, colour = "deepskyblue"),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14)
  )+
  labs(title = "Tumor Stage (ajcc)",
       fill = "Tumor Stage (ajcc)")

if (flag_save_image == TRUE){
  ggsave("RESULTS\\IMMAGINI\\Clinical_freq\\Tumor_stage_ajcc.png", plot,
         width = 7,
         height = 7)
}

### (4) TUMOR STAGE AJCC-T ###
stage_ajcc_T <- data.frame(stage_ajcc_T = c("T1", "T2","T3", "T4","TX"),
                           TCGA = c(sum(clinical_tcga$ajcc_pathologic_t %in% c("T1")),
                                    sum(clinical_tcga$ajcc_pathologic_t %in% c("T2")),
                                    sum(clinical_tcga$ajcc_pathologic_t %in% c("T3")),
                                    sum(clinical_tcga$ajcc_pathologic_t %in% c("T4")),
                                    sum(clinical_tcga$ajcc_pathologic_t %in% c("TX"))),
                           CPTAC = c(sum(clinical_cptac$ajcc_pathologic_t %in% c("T1","T1a","T1c")),
                                     sum(clinical_cptac$ajcc_pathologic_t %in% c("T2")),
                                     sum(clinical_cptac$ajcc_pathologic_t %in% c("T3")),
                                     sum(clinical_cptac$ajcc_pathologic_t %in% c("T4")),
                                     sum(clinical_cptac$ajcc_pathologic_t %in% c("TX")))) %>%
  melt(id.vars="stage_ajcc_T",variable.name="Project")

plot <- ggplot(stage_ajcc_T, aes(x=Project,y=value,fill=stage_ajcc_T)) +
  geom_col(position=position_dodge())+
  scale_fill_manual(values = brewer.pal(5, "Set1"))+
  theme(
    plot.title = element_text(face = "bold.italic", size = 18, colour = "deepskyblue"),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14)
  )+
  labs(title = "Tumor Stage (ajcc T)",
       fill = "Tumor Stage (ajcc T)")

if (flag_save_image == TRUE){
  ggsave("RESULTS\\IMMAGINI\\Clinical_freq\\Tumor_stage_T.png", plot,
         width = 7,
         height = 7)
}

### (5) TUMOR STAGE AJCC-N ###
stage_ajcc_N <- data.frame(stage_ajcc_N = c("N0", "N1","N2", "NX"),
                           TCGA = c(sum(clinical_tcga$ajcc_pathologic_n %in% c("N0")),
                                    sum(clinical_tcga$ajcc_pathologic_n %in% c("N1","N1b")),
                                    sum(clinical_tcga$ajcc_pathologic_n %in% c("N2")),
                                    sum(clinical_tcga$ajcc_pathologic_n %in% c("NX"))),
                           CPTAC = c(sum(clinical_cptac$ajcc_pathologic_n %in% c("N0")),
                                     sum(clinical_cptac$ajcc_pathologic_n %in% c("N1")),
                                     sum(clinical_cptac$ajcc_pathologic_n %in% c("N2")),
                                     sum(clinical_cptac$ajcc_pathologic_n %in% c("NX")))) %>%
  melt(id.vars="stage_ajcc_N",variable.name="Project")

plot <- ggplot(stage_ajcc_N, aes(x=Project,y=value,fill=stage_ajcc_N)) +
  geom_col(position=position_dodge())+
  scale_fill_manual(values = c("slategray1", "skyblue2","skyblue3","steelblue4"))+
  theme(
    plot.title = element_text(face = "bold.italic", size = 18, colour = "deepskyblue"),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14)
  )+
  labs(title = "Tumor Stage (ajcc N)",
       fill = "Tumor Stage (ajcc N)")

if (flag_save_image == TRUE){
  ggsave("RESULTS\\IMMAGINI\\Clinical_freq\\Tumor_stage_N.png", plot,
         width = 7,
         height = 7)
}

### (6) TUMOR STAGE AJCC-M ###
stage_ajcc_M <- data.frame(stage_ajcc_M = c("M0", "M1","MX"),
                           TCGA = c(sum(clinical_tcga$ajcc_pathologic_m %in% c("M0")),
                                    sum(clinical_tcga$ajcc_pathologic_m %in% c("M1")),
                                    sum(clinical_tcga$ajcc_pathologic_m %in% c("MX"))),
                           CPTAC = c(sum(clinical_cptac$ajcc_pathologic_m %in% c("M0")),
                                     sum(clinical_cptac$ajcc_pathologic_m %in% c("M1","M1b")),
                                     sum(clinical_cptac$ajcc_pathologic_m %in% c("MX")))) %>%
  melt(id.vars="stage_ajcc_M",variable.name="Project")

plot <- ggplot(stage_ajcc_M, aes(x=Project,y=value,fill=stage_ajcc_M)) +
  geom_col(position=position_dodge())+
  scale_fill_manual(values = c("plum1", "plum3","plum4"))+
  theme(
    plot.title = element_text(face = "bold.italic", size = 18, colour = "deepskyblue"),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14)
  )+
  labs(title = "Tumor Stage (ajcc M)",
       fill = "Tumor Stage (ajcc M)")

if (flag_save_image == TRUE){
  ggsave("RESULTS\\IMMAGINI\\Clinical_freq\\Tumor_stage_M.png", plot,
         width = 7,
         height = 7)
}

################################################################################################################################


## (B) KAPLAN MEIER (overall survival)##
########################################

#### Stratification ####
## AGE ##
clinical_cptac$age_strat <- ifelse(clinical_cptac$age_at_diagnosis < 65, "<65","≥65")
clinical_tcga$age_strat <- ifelse(clinical_tcga$age_at_diagnosis < 65, "<65","≥65")

## ajcc Stage ##

# CPTAC-3#
for (i in (1:nrow(clinical_cptac))){
  if (clinical_cptac[i,"ajcc_pathologic_stage"] %in% c("Stage IA", "Stage IB")){
    clinical_cptac[i,"stage_strat"] <- "Stage I"
  } else if (clinical_cptac[i,"ajcc_pathologic_stage"] %in% c("Stage IIA", "Stage IIB")){
    clinical_cptac[i,"stage_strat"] <- "Stage II"
  } else if (clinical_cptac[i,"ajcc_pathologic_stage"] %in% c("Stage III")){
    clinical_cptac[i,"stage_strat"] <- "Stage III"
  } else if (clinical_cptac[i,"ajcc_pathologic_stage"] %in% c("Stage IV")){
    clinical_cptac[i,"stage_strat"] <- "Stage IV"
  } 
}

#TCGA-PAAD#
for (i in (1:nrow(clinical_tcga))){
  if (clinical_tcga[i,"ajcc_pathologic_stage"] %in% c("Stage I","Stage IA", "Stage IB")){
    clinical_tcga[i,"stage_strat"] <- "Stage I"
  } else if (clinical_tcga[i,"ajcc_pathologic_stage"] %in% c("Stage IIA", "Stage IIB")){
    clinical_tcga[i,"stage_strat"] <- "Stage II"
  } else if (clinical_tcga[i,"ajcc_pathologic_stage"] %in% c("Stage III")){
    clinical_tcga[i,"stage_strat"] <- "Stage III"
  } else if (clinical_tcga[i,"ajcc_pathologic_stage"] %in% c("Stage IV")){
    clinical_tcga[i,"stage_strat"] <- "Stage IV"
  } 
}

## T Stage

#CPTAC-3#
for (i in (1:nrow(clinical_cptac))){
  if (clinical_cptac[i,"ajcc_pathologic_t"] %in% c("T1","T1a","T1c")){
    clinical_cptac[i,"T_strat"] <- "T1-a-c"
  } else if (clinical_cptac[i,"ajcc_pathologic_t"] %in% c("T2")){
    clinical_cptac[i,"T_strat"] <- "T2"
  } else if (clinical_cptac[i,"ajcc_pathologic_t"] %in% c("T3")){
    clinical_cptac[i,"T_strat"] <- "T3"
  } else if (clinical_cptac[i,"ajcc_pathologic_t"] %in% c("T4")){
    clinical_cptac[i,"T_strat"] <- "T4"
  } else if (clinical_cptac[i,"ajcc_pathologic_t"] %in% c("TX","Unknown")){
    clinical_cptac[i,"T_strat"] <- "TX"
  } 
}

#TCGA-PAAD#
for (i in (1:nrow(clinical_tcga))){
  if (clinical_tcga[i,"ajcc_pathologic_t"] %in% c("T1")){
    clinical_tcga[i,"T_strat"] <- "T1-a-c"
  } else if (clinical_tcga[i,"ajcc_pathologic_t"] %in% c("T2")){
    clinical_tcga[i,"T_strat"] <- "T2"
  } else if (clinical_tcga[i,"ajcc_pathologic_t"] %in% c("T3")){
    clinical_tcga[i,"T_strat"] <- "T3"
  } else if (clinical_tcga[i,"ajcc_pathologic_t"] %in% c("T4")){
    clinical_tcga[i,"T_strat"] <- "T4"
  } else if (clinical_tcga[i,"ajcc_pathologic_t"] %in% c("TX")){
    clinical_tcga[i,"T_strat"] <- "TX"
  } 
}

## N Stage

#CPTAC-3#
for (i in (1:nrow(clinical_cptac))){
  if (clinical_cptac[i,"ajcc_pathologic_n"] %in% c("N0")){
    clinical_cptac[i,"N_strat"] <- "N0"
  } else if (clinical_cptac[i,"ajcc_pathologic_n"] %in% c("N1")){
    clinical_cptac[i,"N_strat"] <- "N1"
  } else if (clinical_cptac[i,"ajcc_pathologic_n"] %in% c("N2")){
    clinical_cptac[i,"N_strat"] <- "N2"
  } else if (clinical_cptac[i,"ajcc_pathologic_n"] %in% c("NX","Unknown")){
    clinical_cptac[i,"N_strat"] <- "NX"
  } 
}

#TCGA-PAAD#
for (i in (1:nrow(clinical_tcga))){
  if (clinical_tcga[i,"ajcc_pathologic_n"] %in% c("N0")){
    clinical_tcga[i,"N_strat"] <- "N0"
  } else if (clinical_tcga[i,"ajcc_pathologic_n"] %in% c("N1")){
    clinical_tcga[i,"N_strat"] <- "N1"
  } else if (clinical_tcga[i,"ajcc_pathologic_n"] %in% c("N2")){
    clinical_tcga[i,"N_strat"] <- "N2"
  } else if (clinical_tcga[i,"ajcc_pathologic_n"] %in% c("NX")){
    clinical_tcga[i,"N_strat"] <- "NX"
  } 
}

## M stage

#CPTAC-3#
for (i in (1:nrow(clinical_cptac))){
  if (clinical_cptac[i,"ajcc_pathologic_m"] %in% c("M0")){
    clinical_cptac[i,"M_strat"] <- "M0"
  } else if (clinical_cptac[i,"ajcc_pathologic_m"] %in% c("M1","M1b")){
    clinical_cptac[i,"M_strat"] <- "M1"
  } else if (clinical_cptac[i,"ajcc_pathologic_m"] %in% c("MX")){
    clinical_cptac[i,"M_strat"] <- "MX"
  } 
}

#TCGA-PAAD#
for (i in (1:nrow(clinical_tcga))){
  if (clinical_tcga[i,"ajcc_pathologic_m"] %in% c("M0")){
    clinical_tcga[i,"M_strat"] <- "M0"
  } else if (clinical_tcga[i,"ajcc_pathologic_m"] %in% c("M1","M1b")){
    clinical_tcga[i,"M_strat"] <- "M1"
  } else if (clinical_tcga[i,"ajcc_pathologic_m"] %in% c("MX")){
    clinical_tcga[i,"M_strat"] <- "MX"
  } 
}



clinical_cptac_strat <- clinical_cptac
clinical_tcga_strat <- clinical_tcga
if (flag_save==TRUE){
  save(clinical_cptac_strat, file = "RESULTS\\DATAFRAME\\clinical_cptac_strat.Rda")
  save(clinical_tcga_strat, file = "RESULTS\\DATAFRAME\\clinical_tcga_strat.Rda")
}


#### PLOT ####

# xlim = c(0,84)
### K-M Gender ###

centr = theme_grey() + theme(plot.title = element_text(hjust = 0.5, size = 14))

fit_cptac <- survfit(Surv(as.numeric(clinical_cptac$overall_survival), clinical_cptac$deceased) ~ clinical_cptac$gender)
plot1 <- ggsurvplot(fit_cptac,
                    data = clinical_cptac,
                    pval = T,
                    risk.table = F,
                    risk.table.y.text = FALSE,
                    palette = c("#2E9FDF", "red"),
                    xlab ="Time in Months",
                    break.time.by = 12,
                    legend.labs =
                      c("Male", "Female"),
                    xlim = c(0,62),
                    size = 1.5,
                    legend.title = "",
                    ggtheme = centr)+
  ggtitle("CPTAC, Gender")

fit_tcga <- survfit(Surv(as.numeric(clinical_tcga$overall_survival), clinical_tcga$deceased) ~ clinical_tcga$gender)
plot2 <- ggsurvplot(fit_tcga,
                    data = clinical_tcga,
                    pval = T,
                    xlim = c(0,74),
                    risk.table = F,
                    risk.table.y.text = FALSE,
                    palette = c("#2E9FDF", "red"),
                    xlab ="Time in Months",
                    break.time.by = 12,
                    legend.labs =
                      c("Male", "Female"),
                    size = 1.5,
                    legend.title = "",
                    ggtheme = centr)+
  ggtitle("TCGA, Gender")

plot <- arrange_ggsurvplots(list(plot2,plot1), print = TRUE,
                            ncol = 2, nrow = 1, risk.table.height = 0.3)
if (flag_save_km == TRUE){
  ggsave("./RESULTS/IMMAGINI/Clinical_KM/Gender_KM.png", plot,
         width = 12,
         height = 7)
}

### K-M Age ###
fit_cptac <- survfit(Surv(as.numeric(clinical_cptac$overall_survival), clinical_cptac$deceased) ~ clinical_cptac$age_strat)
plot1 <- ggsurvplot(fit_cptac,
                    data = clinical_cptac,
                    pval = T,
                    risk.table = F,
                    risk.table.y.text = FALSE,
                    palette = c("darkorange", "green4"),
                    xlab ="Time in Months",
                    break.time.by = 12,
                    legend.labs = c("≥65","<65"),
                    xlim = c(0,62),
                    size = 1.5,
                    legend.title = "",
                    ggtheme = centr)+
  ggtitle("CPTAC, Age ad diagnosis")

fit_tcga <- survfit(Surv(as.numeric(clinical_tcga$overall_survival), clinical_tcga$deceased) ~ clinical_tcga$age_strat)
plot2 <- ggsurvplot(fit_tcga,
                    data = clinical_tcga,
                    pval = T,
                    risk.table = F,
                    risk.table.y.text = FALSE,
                    palette = c("darkorange", "green4"),
                    xlab ="Time in Months",
                    break.time.by = 12,
                    legend.labs = c("≥65","<65"),
                    xlim = c(0,74),
                    size = 1.5,
                    legend.title = "",
                    ggtheme = centr) +
  ggtitle("TCGA, Age ad diagnosis")

plot <- arrange_ggsurvplots(list(plot2,plot1), print = TRUE,
                            ncol = 2, nrow = 1, risk.table.height = 0.3)
if (flag_save_km == TRUE){
  ggsave("./RESULTS/IMMAGINI/Clinical_KM/Age_KM.png", plot,
         width = 12,
         height = 7)
}

### K-M ajcc STAGE ###
fit_cptac <- survfit(Surv(as.numeric(clinical_cptac$overall_survival), clinical_cptac$deceased) ~ clinical_cptac$stage_strat)

plot1 <- ggsurvplot(fit_cptac,
                    data = clinical_cptac,
                    pval = T,
                    risk.table = F,
                    risk.table.y.text = FALSE,
                    palette = c("gold", "orange","red1","red4"),
                    xlab ="Time in Months",
                    break.time.by = 12,
                    legend.labs = c("Stage I","Stage II", "Stage III", "Stage IV"),
                    xlim = c(0,62),
                    size = 1.5,
                    legend.title = "",
                    ggtheme = centr)+
ggtitle("CPTAC, AJCC STAGE")


fit_tcga <- survfit(Surv(as.numeric(clinical_tcga$overall_survival), clinical_tcga$deceased) ~ clinical_tcga$stage_strat)
plot2 <- ggsurvplot(fit_tcga,
                    data = clinical_tcga,
                    pval = T,
                    risk.table = F,
                    risk.table.y.text = FALSE,
                    palette = c("gold", "orange","red1","red4"),
                    xlab ="Time in Months",
                    break.time.by = 12,
                    legend.labs = c("Stage I","Stage II", "Stage III", "Stage IV"),
                    xlim = c(0,74),
                    size = 1.5,
                    legend.title = "",
                    ggtheme = centr)+
  ggtitle("TCGA, AJCC STAGE")

plot <- arrange_ggsurvplots(list(plot2,plot1), print = TRUE,
                            ncol = 2, nrow = 1, risk.table.height = 0.3)

if (flag_save_km == TRUE){
  ggsave("./RESULTS/IMMAGINI/Clinical_KM/Stage_KM.png", plot,
         width = 12,
         height = 7)
}

### K-M ajcc STAGE T ###
fit_cptac <- survfit(Surv(as.numeric(clinical_cptac$overall_survival), clinical_cptac$deceased) ~ clinical_cptac$T_strat)
plot1 <- ggsurvplot(fit_cptac,
                    data = clinical_cptac,
                    pval = T,
                    risk.table = F,
                    risk.table.y.text = FALSE,
                    xlab ="Time in Months",
                    break.time.by = 12,
                    legend.labs = c("T1","T2", "T3", "T4", "TX"),
                    xlim = c(0,62),
                    size = 1.5,
                    legend.title = "",
                    ggtheme = centr)+
  ggtitle("CPTAC, AJCC STAGE T")

fit_tcga <- survfit(Surv(as.numeric(clinical_tcga$overall_survival), clinical_tcga$deceased) ~ clinical_tcga$T_strat)
plot2 <- ggsurvplot(fit_tcga,
                    data = clinical_tcga,
                    pval = T,
                    risk.table = F,
                    risk.table.y.text = FALSE,
                    xlab ="Time in Months",
                    break.time.by = 12,
                    legend.labs = c("T1","T2", "T3", "T4"),
                    xlim = c(0,74),
                    size = 1.5,
                    legend.title = "",
                    ggtheme = centr)+
  ggtitle("TCGA, AJCC STAGE T")

plot <- arrange_ggsurvplots(list(plot2,plot1), print = TRUE,
                            ncol = 2, nrow = 1, risk.table.height = 0.3)

if (flag_save_km == TRUE){
  ggsave("./RESULTS/IMMAGINI/Clinical_KM/stageT_KM.png", plot,
         width = 12,
         height = 7)
}


### K-M ajcc STAGE N ###
fit_cptac <- survfit(Surv(as.numeric(clinical_cptac$overall_survival), clinical_cptac$deceased) ~ clinical_cptac$N_strat)
plot1 <- ggsurvplot(fit_cptac,
                    data = clinical_cptac,
                    pval = T,
                    risk.table = F,
                    risk.table.y.text = FALSE,
                    palette = c("slategray1", "skyblue2","steelblue3","steelblue4"),
                    xlab ="Time in Months",
                    break.time.by = 12,
                    legend.labs = c("N0","N1", "N2", "NX"),
                    xlim = c(0,62),
                    size = 1.5,
                    legend.title = "",
                    ggtheme = centr)+
  ggtitle("CPTAC, AJCC STAGE N")

fit_tcga <- survfit(Surv(as.numeric(clinical_tcga$overall_survival), clinical_tcga$deceased) ~ clinical_tcga$N_strat)
plot2 <- ggsurvplot(fit_tcga,
                    data = clinical_tcga,
                    pval = T,
                    risk.table = T,
                    risk.table.y.text = FALSE,
                    palette = c("slategray1", "skyblue2"),
                    xlab ="Time in Months",
                    break.time.by = 12,
                    legend.labs = c("N0","N1"),
                    xlim = c(0,74),
                    size = 1.5,
                    legend.title = "",
                    ggtheme = centr) +
  ggtitle("TCGA (AJCC STAGE N)")

plot <- arrange_ggsurvplots(list(plot2,plot1), print = TRUE,
                            ncol = 2, nrow = 1, risk.table.height = 0.3)
if (flag_save_km == TRUE){
  ggsave("./RESULTS/IMMAGINI/Clinical_KM/stageN_KM.png", plot,
         width = 12,
         height = 7)
}


# ### K-M ajcc STAGE M ###
# fit_cptac <- survfit(Surv(as.numeric(clinical_cptac$overall_survival), clinical_cptac$deceased) ~ clinical_cptac$M_strat)
# plot1 <- ggsurvplot(fit_cptac,
#                     data = clinical_cptac,
#                     pval = T,
#                     risk.table = F,
#                     risk.table.y.text = FALSE,
#                     palette = c("plum1", "plum3","plum4"),
#                     xlab ="Time in Months",
#                     break.time.by = 12,
#                     legend.labs = c("M0","M1", "MX"),
#                     xlim = c(0,62),
#                     size = 1.5,
#                     legend.title = "",
#                     ggtheme = centr)+
#   ggtitle("CPTAC AJCC STAGE M")
# 
# fit_tcga <- survfit(Surv(as.numeric(clinical_tcga$overall_survival), clinical_tcga$deceased) ~ clinical_tcga$M_strat)
# plot2 <- ggsurvplot(fit_tcga,
#                     data = clinical_tcga,
#                     pval = T,
#                     risk.table = F,
#                     risk.table.y.text = FALSE,
#                     palette = c("plum1", "plum3","plum4"),
#                     xlab ="Time in Months",
#                     break.time.by = 12,
#                     legend.labs =c("M0","M1", "MX"),
#                     xlim = c(0,74),
#                     size = 1.5,
#                     legend.title = "",
#                     ggtheme = centr)+
#   ggtitle("TCGA (AJCC STAGE M)")
# 
# plot <- arrange_ggsurvplots(list(plot2,plot1), print = TRUE,
#                             ncol = 2, nrow = 1, risk.table.height = 0.3)
# 
# if (flag_save_km == TRUE){
#   ggsave("./RESULTS/IMMAGINI/Clinical_KM/stageM_KM.png", plot,
#          width = 12,
#          height = 7)
# }
# 
# ### K-M PRIMARY-DIAGNOSIS ###
# fit_cptac <- survfit(Surv(as.numeric(clinical_cptac$overall_survival), clinical_cptac$deceased) ~ clinical_cptac$primary_diagnosis)
# plot1 <- ggsurvplot(fit_cptac,
#                     data = clinical_cptac,
#                     pval = T,
# 
#                     risk.table = T,
#                     risk.table.y.text = FALSE,
# 
#                     palette = c("red3"),
#                     xlab ="Time in Months",
#                     break.time.by = 12,
#                     xlim = xlim,
# 
#                     legend.labs =
#                       c("Infiltrating duct carcinoma"))+
#   ggtitle("CPTAC (Primary diagnosis)")
# 
# fit_tcga <- survfit(Surv(as.numeric(clinical_tcga$overall_survival), clinical_tcga$deceased) ~ clinical_tcga$primary_diagnosis)
# plot2 <- ggsurvplot(fit_tcga,
#                     data = clinical_tcga,
#                     pval = T,
# 
#                     risk.table = T,
#                     risk.table.y.text = FALSE,
# 
#                     palette = c("red3"),
#                     xlab ="Time in Months",
#                     break.time.by = 12,
#                     xlim = xlim,
# 
#                     legend.labs =
#                       c("Infiltrating duct carcinoma"))+
#   ggtitle("TCGA (Primary diagnosis)")
# 
# plot <- arrange_ggsurvplots(list(plot2,plot1), print = TRUE,
#                             ncol = 2, nrow = 1, risk.table.height = 0.3)
# 
# if (flag_save_km == TRUE){
#   ggsave("RESULTS\\IMMAGINI\\Clinical_KM\\PrimaryD_KM.png", plot,
#          width = 12,
#          height = 7)
# }
################################################################################################################################

