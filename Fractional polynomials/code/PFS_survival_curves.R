# Darolutamide in HSPC NMA #
# Survival extrapolation FP models
# PFS outcome

# 1. Packages used --------------------------------------------------------
pkgs <- c("tidyverse", "readxl", "here", "haven", "survival")
lapply(pkgs, library, character.only = T)
source(here("01_scripts", "utils.R"))
set.seed(03082022)

# 2. Import data --------------------------------------------------------
load(here("02_data", "arasens_ipd.rda"))
arasens <- arasens_ipd$CROD
arasens_daro <- filter(arasens, treatment == "Darolutamide+docetaxel arm")

# Colors for the plot
color=c("#FF0000", "#00FF00", "#B200ED", "#FFA500", "#0000FF", "#fb9a99", "#e31a1c", 
                 "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a", "#581845", "#b15928", "#FFEA00",
                 "#000000")
                 
#3. Survival curve darolutamide ---------------------------
# Second-order model with lowest DIC (5 models)
fp_data_1 <- read.csv(here("02_data", "PFS_results_FP_sim14_108.csv")) 
fp_data_2 <- read.csv(here("02_data", "PFS_results_FP_sim11_108.csv")) 
fp_data_3 <- read.csv(here("02_data", "PFS_results_FP_sim19_108.csv")) 
fp_data_4 <- read.csv(here("02_data", "PFS_results_FP_sim10_108.csv")) 
fp_data_5 <- read.csv(here("02_data", "PFS_results_FP_sim15_108.csv")) 

# Export survival at 108 months
daro_surv_model_1 <- fp_data_1[which(fp_data_1$X == "S[1,108]"),c(6,4,8)]
daro_surv_model_2 <- fp_data_2[which(fp_data_2$X == "S[1,108]"),c(6,4,8)]
daro_surv_model_3 <- fp_data_3[which(fp_data_3$X == "S[1,108]"),c(6,4,8)]
daro_surv_model_4 <- fp_data_4[which(fp_data_4$X == "S[1,108]"),c(6,4,8)]
daro_surv_model_5 <- fp_data_5[which(fp_data_5$X == "S[1,108]"),c(6,4,8)]
daro_surv_models <- bind_rows(daro_surv_model_1,
                              daro_surv_model_2,
                              daro_surv_model_3,
                              daro_surv_model_4,
                              daro_surv_model_5)

daro_surv_models <- apply(daro_surv_models, 1, format_results)
write.csv(daro_surv_models, here("05_tables", "daro_surv_models_PFS.csv"))

# Survival data to plot
graph_daro <- data.frame(time=c(0:108),
                         model_1=c(1, fp_data_1$mean[1:108]),
                         model_2=c(1, fp_data_2$mean[1:108]),
                         model_3=c(1, fp_data_3$mean[1:108]),
                         model_4=c(1, fp_data_4$mean[1:108]),
                         model_5=c(1, fp_data_5$mean[1:108]))

# KM curve
KM.est_daro <- survfit(Surv(time, event) ~ 1, data= arasens_daro, 
                       type="kaplan-meier", conf.int=FALSE)

jpeg(file = here("04_figures", "extrapolation_darolutamide_PFS_v3.jpg"), width = 900, height = 800, res = 120)

plot(KM.est_daro, xlab="Time (months)", ylab="PFS", xaxt="n", yaxt="n", main=" ",
     xlim = c(0,108), ylim=c(0,1),
     mark.time=FALSE, col=color[15], conf.int=F)

title("Progression-free survival (PFS) \n2nd order FP models")
axis(2, at=c(0, 0.2, 0.4, 0.6, 0.8, 1))
axis(1, at= seq(0, 108, 20))

## Add prediction lines#
lines(graph_daro$time, graph_daro$model_1, col=color[1])
lines(graph_daro$time, graph_daro$model_2, col=color[2])
lines(graph_daro$time, graph_daro$model_3, col=color[3])
lines(graph_daro$time, graph_daro$model_4, col=color[4])
lines(graph_daro$time, graph_daro$model_5, col=color[5])

# legend
legend(x = 70, y = 1,
       c("P1 = -0.5, P2 = -0.5", "P1 = -1, P2 = 0.5", "P1 = 0, P2 = 0.5",  
         "P1 = -1, P2 = 0", "P1 = -0.5, P2 = 0", "ARASENS"),
       col=c(color[1], color[2], color[3], color[4], color[5], color[15]),
       lty=c(1,1,1,1), ncol=1, text.width=6, box.lty=0)

# save plot
dev.off()


#4. Survival curve enza+ADT ---------------------------
# Second-order model with lowest DIC (5 models) ENZAMET as comparator trial
fp_data_1_ENZA <- read.csv(here("02_data", "PFS_results_FP_sim14_ENZA.csv")) 
fp_data_2_ENZA <- read.csv(here("02_data", "PFS_results_FP_sim11_ENZA.csv")) 
fp_data_3_ENZA <- read.csv(here("02_data", "PFS_results_FP_sim19_ENZA.csv")) 
fp_data_4_ENZA <- read.csv(here("02_data", "PFS_results_FP_sim10_ENZA.csv")) 
fp_data_5_ENZA <- read.csv(here("02_data", "PFS_results_FP_sim15_ENZA.csv")) 

                         
# Export survival at 108 months
enza_surv_model_1 <- fp_data_1[which(fp_data_1_ENZA$X == "S[3,108]"),c(6,4,8)]
enza_surv_model_2 <- fp_data_2[which(fp_data_1_ENZA$X == "S[3,108]"),c(6,4,8)]
enza_surv_model_3 <- fp_data_3[which(fp_data_1_ENZA$X == "S[3,108]"),c(6,4,8)]
enza_surv_model_4 <- fp_data_4[which(fp_data_1_ENZA$X == "S[3,108]"),c(6,4,8)]
enza_surv_model_5 <- fp_data_5[which(fp_data_1_ENZA$X == "S[3,108]"),c(6,4,8)]
enza_surv_models <- bind_rows(enza_surv_model_1,
                              enza_surv_model_2,
                              enza_surv_model_3,
                              enza_surv_model_4,
                              enza_surv_model_5)

enza_surv_models <- apply(enza_surv_models, 1, format_results)
write.csv(enza_surv_models, here("05_tables", "enza_surv_models_PFS.csv"))

# IPDs
ENZAMET_enza <- read_csv(here("02_data", "IPD reconstructed", "PFS_ENZAMET_enzalutamide_ipd.csv"))
ARCHES_enza <- read_csv(here("02_data", "IPD reconstructed", "PFS_ARCHES_ENZA+ADT_ipd.csv"))
VAI_enza <- read_csv(here("02_data", "IPD reconstructed", "PFS_Vaishampayan 2021_enzalutamide_ipd.csv"))

graph_enza <- data.frame(time=c(0:108),
                         model_1=c(1, fp_data_1_ENZA$mean[325:432]),
                         model_2=c(1, fp_data_2_ENZA$mean[325:432]),
                         model_3=c(1, fp_data_3_ENZA$mean[325:432]),
                         model_4=c(1, fp_data_4_ENZA$mean[325:432]),
                         model_5=c(1, fp_data_5_ENZA$mean[325:432]))

# KM curve
KM.est_enza <- survfit(Surv(time, event) ~ 1, data= ARCHES_enza, 
                       type="kaplan-meier", conf.int=FALSE)
KM.est_enzamet <- survfit(Surv(time, event) ~ 1, data= ENZAMET_enza, 
                       type="kaplan-meier", conf.int=FALSE)
KM.est_vai <- survfit(Surv(time, event) ~ 1, data= VAI_enza, 
                       type="kaplan-meier", conf.int=FALSE)

jpeg(file = here("04_figures", "extrapolation_enza+ADT_PFS_v4.jpg"), width = 900, height = 800, res = 120)

plot(KM.est_enza, xlab="Time (months)", ylab="PFS", xaxt="n", yaxt="n", main=" ",
     xlim = c(0,108), ylim=c(0,1),
     mark.time=FALSE, col=color[15], conf.int=F, lwd = 2)

title("Progression-free survival (PFS) \n2nd order FP models")
axis(2, at=c(0, 0.2, 0.4, 0.6, 0.8, 1))
axis(1, at= seq(0, 108, 20))

## Add prediction lines#
lines(graph_enza$time, graph_enza$model_1, col=color[1])
lines(graph_enza$time, graph_enza$model_2, col=color[2])
lines(graph_enza$time, graph_enza$model_3, col=color[3])
lines(graph_enza$time, graph_enza$model_4, col=color[4])
lines(graph_enza$time, graph_enza$model_5, col=color[5])
lines(KM.est_enzamet, col = color[15], conf.int = F, lty = 2)
lines(KM.est_vai, col = color[15], conf.int = F, lty = 3)

# legend
legend(x = 70, y = 1,
       c("P1 = -0.5, P2 = -0.5", "P1 = -1, P2 = 0.5", "P1 = 0, P2 = 0.5",  
         "P1 = -1, P2 = 0", "P1 = -0.5, P2 = 0", "ARCHES", "ENZAMET", "Vaishampayan"),
       col=c(color[1], color[2], color[3], color[4], color[5], color[15], color[15], color[15]),
       lty=c(1,1,1,1,1,1,2,3), lwd=c(1,1,1,1,1,2,1,1),
       ncol=1, text.width=6, box.lty=0)

# save plot
dev.off()


#5. Survival curve docetaxel+ADT ---------------------------
# Second-order model with lowest DIC (5 models) CHAARTED as comparator trial
fp_data_1_CHAAR <- read.csv(here("02_data", "PFS_results_FP_sim14_CHAART.csv")) 
fp_data_2_CHAAR <- read.csv(here("02_data", "PFS_results_FP_sim11_CHAART.csv")) 
fp_data_3_CHAAR <- read.csv(here("02_data", "PFS_results_FP_sim19_CHAART.csv")) 
fp_data_4_CHAAR <- read.csv(here("02_data", "PFS_results_FP_sim10_CHAART.csv")) 
fp_data_5_CHAAR <- read.csv(here("02_data", "PFS_results_FP_sim15_CHAART.csv")) 

# Export survival at 108 months
doc_surv_model_1 <- fp_data_1[which(fp_data_1_CHAAR$X == "S[2,108]"),c(6,4,8)]
doc_surv_model_2 <- fp_data_2[which(fp_data_1_CHAAR$X == "S[2,108]"),c(6,4,8)]
doc_surv_model_3 <- fp_data_3[which(fp_data_1_CHAAR$X == "S[2,108]"),c(6,4,8)]
doc_surv_model_4 <- fp_data_4[which(fp_data_1_CHAAR$X == "S[2,108]"),c(6,4,8)]
doc_surv_model_5 <- fp_data_5[which(fp_data_1_CHAAR$X == "S[2,108]"),c(6,4,8)]
doc_surv_models <- bind_rows(doc_surv_model_1,
                              doc_surv_model_2,
                              doc_surv_model_3,
                              doc_surv_model_4,
                              doc_surv_model_5)

doc_surv_models <- apply(doc_surv_models, 1, format_results)
write.csv(doc_surv_models, here("05_tables", "doc_surv_models_PFS.csv"))

# IPDs
graph_doc <- data.frame(time=c(0:108),
                        model_1=c(1, fp_data_1_CHAAR$mean[109:216]),
                        model_2=c(1, fp_data_2_CHAAR$mean[109:216]),
                        model_3=c(1, fp_data_3_CHAAR$mean[109:216]),
                        model_4=c(1, fp_data_4_CHAAR$mean[109:216]),
                        model_5=c(1, fp_data_5_CHAAR$mean[109:216]))

# ARASENS
arasens_pcb <- filter(arasens, treatment == "Placebo+docetaxel arm")

# CHAARTED
CHAARTED_doc <- read_csv(here("02_data", "IPD reconstructed", "PFS_CHAARTED_ADT+Doc_ipd.csv"))

# PEACE-1
PEACE_doc <- read_csv(here("02_data", "IPD reconstructed", "PFS_PEACE_ADT+Doc_ipd.csv"))

# STAMPEDE-3
STAMP3_doc <- read_csv(here("02_data", "IPD reconstructed", "PFS_STAMPEDE-3_ADT+Doc_ipd.csv"))

# STAMPEDE-4
STAMP4_doc <- read_csv(here("02_data", "IPD reconstructed", "PFS_STAMPEDE-4_ADT+DocP_ipd.csv"))

# GETUG
GETUG_doc <- read_csv(here("02_data", "IPD reconstructed", "PFS_GETUG_ADT+Doc_ipd.csv"))

# KM curves
KM.est_doc <- survfit(Surv(time, event) ~ 1, data= CHAARTED_doc, 
                      type="kaplan-meier", conf.int=FALSE)

KM.est_doc2 <- survfit(Surv(time, event) ~ 1, data= arasens_pcb, 
                      type="kaplan-meier", conf.int=FALSE)

KM.est_doc3 <- survfit(Surv(time, event) ~ 1, data= PEACE_doc, 
                      type="kaplan-meier", conf.int=FALSE)

KM.est_doc4 <- survfit(Surv(time, event) ~ 1, data= STAMP3_doc, 
                      type="kaplan-meier", conf.int=FALSE)

KM.est_doc5 <- survfit(Surv(time, event) ~ 1, data= STAMP4_doc, 
                      type="kaplan-meier", conf.int=FALSE)

KM.est_doc6 <- survfit(Surv(time, event) ~ 1, data= GETUG_doc, 
                      type="kaplan-meier", conf.int=FALSE)

jpeg(file = here("04_figures", "extrapolation_doc+ADT_PFS_v4.jpg"), width = 900, height = 800, res = 120)

plot(KM.est_doc, xlab="Time (months)", ylab="PFS", xaxt="n", yaxt="n", main=" ",
     xlim = c(0,108), ylim=c(0,1),
     mark.time=FALSE, col=color[15], conf.int=F, lwd = 2)

title("Progression-free survival (PFS) \n2nd order FP models")
axis(2, at=c(0, 0.2, 0.4, 0.6, 0.8, 1))
axis(1, at= seq(0, 108, 20))

## Add prediction lines#
lines(graph_doc$time, graph_doc$model_1, col=color[1])
lines(graph_doc$time, graph_doc$model_2, col=color[2])
lines(graph_doc$time, graph_doc$model_3, col=color[3])
lines(graph_doc$time, graph_doc$model_4, col=color[4])
lines(graph_doc$time, graph_doc$model_5, col=color[5])
lines(KM.est_doc2, conf.int = F, col = color[15], lty = 2)
lines(KM.est_doc3, conf.int = F, col = color[15], lty = 3)
lines(KM.est_doc4, conf.int = F, col = color[15], lty = 4)
lines(KM.est_doc5, conf.int = F, col = color[15], lty = 5)
lines(KM.est_doc6, conf.int = F, col = color[15], lty = 6)


# legend
legend(x = 60, y = 1,
       c("P1 = -0.5, P2 = -0.5", "P1 = -1, P2 = 0.5", "P1 = 0, P2 = 0.5",  
         "P1 = -1, P2 = 0", "P1 = -0.5, P2 = 0", "CHAARTED", "ARASENS",
         "PEACE-1", "STAMPEDE arms C vs A", "STAMPEDE arms C vs G", "GETUG"),
       col=c(color[1], color[2], color[3], color[4], color[5], color[15], color[15], color[15],
             color[15], color[15], color[15]),
       lty=c(1,1,1,1,1,1,2,3,4,5,6), lwd=c(1,1,1,1,1,2,1,1,1,1,1),
       ncol=1, text.width=6, box.lty=0)

# save plot
dev.off()

#6. Survival curve abi ---------------------------
# Second-order model with lowest DIC (5 models) LATITUDE as comparator trial
fp_data_1_LATITUDE <- read.csv(here("02_data", "PFS_results_FP_sim14_LATITUDE.csv")) 
fp_data_2_LATITUDE <- read.csv(here("02_data", "PFS_results_FP_sim11_LATITUDE.csv")) 
fp_data_3_LATITUDE <- read.csv(here("02_data", "PFS_results_FP_sim19_LATITUDE.csv")) 
fp_data_4_LATITUDE <- read.csv(here("02_data", "PFS_results_FP_sim10_LATITUDE.csv")) 
fp_data_5_LATITUDE <- read.csv(here("02_data", "PFS_results_FP_sim15_LATITUDE.csv")) 


# Export survival at 108 months
abi_surv_model_1 <- fp_data_1[which(fp_data_1_LATITUDE$X == "S[5,108]"),c(6,4,8)]
abi_surv_model_2 <- fp_data_2[which(fp_data_1_LATITUDE$X == "S[5,108]"),c(6,4,8)]
abi_surv_model_3 <- fp_data_3[which(fp_data_1_LATITUDE$X == "S[5,108]"),c(6,4,8)]
abi_surv_model_4 <- fp_data_4[which(fp_data_1_LATITUDE$X == "S[5,108]"),c(6,4,8)]
abi_surv_model_5 <- fp_data_5[which(fp_data_1_LATITUDE$X == "S[5,108]"),c(6,4,8)]
abi_surv_models <- bind_rows(abi_surv_model_1,
                             abi_surv_model_2,
                             abi_surv_model_3,
                             abi_surv_model_4,
                             abi_surv_model_5)

abi_surv_models <- apply(abi_surv_models, 1, format_results)
write.csv(abi_surv_models, here("05_tables", "abi_surv_models_PFS.csv"))


# IPDs
graph_abi <- data.frame(time=c(0:108),
                        model_1=c(1, fp_data_1$mean[433:540]),
                        model_2=c(1, fp_data_2$mean[433:540]),
                        model_3=c(1, fp_data_3$mean[433:540]),
                        model_4=c(1, fp_data_4$mean[433:540]),
                        model_5=c(1, fp_data_5$mean[433:540]))

LATITUDE_abi <- read_csv(here("02_data", "IPD reconstructed", "PFS_LATITUDE_AA+P+ADT_ipd.csv"))
STAMP_abi <- read_csv(here("02_data", "IPD reconstructed", "PFS_STAMPEDE-4_ADT+AAP_ipd.csv"))

# KM curve
KM.est_abi <- survfit(Surv(time, event) ~ 1, data= LATITUDE_abi, 
                       type="kaplan-meier", conf.int=FALSE)
KM.est_abi2 <- survfit(Surv(time, event) ~ 1, data= STAMP_abi, 
                       type="kaplan-meier", conf.int=FALSE)


jpeg(file = here("04_figures", "extrapolation_abiraterone_PFS_v4.jpg"), width = 900, height = 800, res = 120)

plot(KM.est_abi, xlab="Time (months)", ylab="PFS", xaxt="n", yaxt="n", main=" ",
     xlim = c(0,108), ylim=c(0,1),
     mark.time=FALSE, col=color[15], conf.int=F, lwd = 2)

title("Progression-free survival (PFS) \n2nd order FP models")
axis(2, at=c(0, 0.2, 0.4, 0.6, 0.8, 1))
axis(1, at= seq(0, 108, 20))

## Add prediction lines#
lines(graph_abi$time, graph_abi$model_1, col=color[1])
lines(graph_abi$time, graph_abi$model_2, col=color[2])
lines(graph_abi$time, graph_abi$model_3, col=color[3])
lines(graph_abi$time, graph_abi$model_4, col=color[4])
lines(graph_abi$time, graph_abi$model_5, col=color[5])
lines(KM.est_abi2, conf.int = F, col = color[15], lty = 2)

# legend
legend(x = 60, y = 1,
       c("P1 = -0.5, P2 = -0.5", "P1 = -1, P2 = 0.5", "P1 = 0, P2 = 0.5",  
         "P1 = -1, P2 = 0", "P1 = -0.5, P2 = 0", "LATITUDE", "STAMPEDE arms C vs G"),
       col=c(color[1], color[2], color[3], color[4], color[5], color[15], color[15]),
       lty=c(1,1,1,1, 1,1,2), lwd=c(1,1,1,1, 1,2,1), ncol=1, text.width=6, box.lty=0)

# save plot
dev.off()


# Survival at 108 months for each treatment-----------------------
# It uses the sim14 model (P1: -0.5, P2: -0.5)
surv_daro <- fp_data_1[which(fp_data_1$X == "S[1,108]"),c(6,4,8)]
surv_doc <- fp_data_1[which(fp_data_1$X == "S[2,108]"),c(6,4,8)]
surv_enza <- fp_data_1[which(fp_data_1$X == "S[3,108]"),c(6,4,8)]
surv_adt <- fp_data_1[which(fp_data_1$X == "S[4,108]"),c(6,4,8)]
surv_abi <- fp_data_1[which(fp_data_1$X == "S[5,108]"),c(6,4,8)]
surv_apa <- fp_data_1[which(fp_data_1$X == "S[6,108]"),c(6,4,8)]
surv_enzadoc <- fp_data_1[which(fp_data_1$X == "S[7,108]"),c(6,4,8)]
surv_abidoc  <- fp_data_1[which(fp_data_1$X == "S[8,108]"),c(6,4,8)]

surv_108 <- bind_rows(surv_daro,
                      surv_doc,
                      surv_enza,
                      surv_adt,
                      surv_abi,
                      surv_apa,
                      surv_enzadoc,
                      surv_abidoc)

surv_108 <- apply(surv_108, 1, summary.stats)
write.csv(surv_108, here("05_tables", "median_survival_PFS.csv"))
