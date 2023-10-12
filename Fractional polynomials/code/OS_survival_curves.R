# Darolutamide in HSPC NMA #
# Survival extrapolation FP models
# OS outcome

# 1. Packages used --------------------------------------------------------
pkgs <- c("tidyverse", "readxl", "here", "haven", "survival")
lapply(pkgs, library, character.only = T)
source(here("Fractional polynomials/code", "utils.R"))
set.seed(03082022)

# Colors for the plots
color=c("#FF0000", "#00FF00", "#B200ED", "#FFA500", "#0000FF", "#fb9a99", "#e31a1c", 
                 "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a", "#581845", "#b15928", "#FFEA00",
                 "#000000")

# 2. Import data --------------------------------------------------------
# 2.1 IPDs
load(here("Fractional polynomials/data", "arasens_ipd.rda"))
arasens <- arasens_ipd$OS
arasens_pcb <- filter(arasens, treatment == "Placebo+docetaxel arm")
arasens_daro <- filter(arasens, treatment == "Darolutamide+docetaxel arm")
ENZAMET_enza_doc <- read_csv(here("Fractional polynomials/data", "IPD reconstructed", "OS_ENZAMET_enzalutamide+Doc_ipd.csv"))
doc <- read_csv(here("Fractional polynomials/data", "IPD reconstructed", "OS_PEACE_ADT+Doc_ipd.csv"))
abi <- read_csv(here("Fractional polynomials/data", "IPD reconstructed", "OS_STAMPEDE-2_AA+ADT_ipd.csv"))
abi2 <- read_csv(here("Fractional polynomials/data", "IPD reconstructed", "OS_LATITUDE_AA+P+ADT_ipd.csv"))


#2.2 Second-order model with lowest DIC (5 models)
fp_data_1 <- read.csv(here("Fractional polynomials/data", "OS_results_FP_sim14_108.csv"))  # model with P1: -0.5, P2: -0.5
fp_data_2 <- read.csv(here("Fractional polynomials/data", "OS_results_FP_sim15_108.csv"))  # model with P1: -0.5, P2:  0
fp_data_3 <- read.csv(here("Fractional polynomials/data", "OS_results_FP_sim11_108.csv"))  # model with P1: -1,   P2:  0.5
fp_data_4 <- read.csv(here("Fractional polynomials/data", "OS_results_FP_sim10_108.csv"))  # model with P1: -1,   P2:  0
fp_data_5 <- read.csv(here("Fractional polynomials/data", "OS_results_FP_sim09_108.csv"))  # model with P1: -1,   P2: -0.5

#3. Export survival at 108 months-------------------------------------
#3.1 Darolutamide
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
write.csv(daro_surv_models, here("Fractional polynomials/data", "daro_surv_models_OS.csv"))

#3.2 Enzalutamide
enza_surv_model_1 <- fp_data_1[which(fp_data_1$X == "S[4,108]"),c(6,4,8)]
enza_surv_model_2 <- fp_data_2[which(fp_data_2$X == "S[4,108]"),c(6,4,8)]
enza_surv_model_3 <- fp_data_3[which(fp_data_3$X == "S[4,108]"),c(6,4,8)]
enza_surv_model_4 <- fp_data_4[which(fp_data_4$X == "S[4,108]"),c(6,4,8)]
enza_surv_model_5 <- fp_data_5[which(fp_data_5$X == "S[4,108]"),c(6,4,8)]
enza_surv_models <- bind_rows(enza_surv_model_1,
                              enza_surv_model_2,
                              enza_surv_model_3,
                              enza_surv_model_4,
                              enza_surv_model_5)

enza_surv_models <- apply(enza_surv_models, 1, format_results)
write.csv(enza_surv_models, here("Fractional polynomials/data", "enza_surv_models_OS.csv"))

# 3.3 Docetaxel
doc_surv_model_1 <- fp_data_1[which(fp_data_1$X == "S[2,108]"),c(6,4,8)]
doc_surv_model_2 <- fp_data_2[which(fp_data_2$X == "S[2,108]"),c(6,4,8)]
doc_surv_model_3 <- fp_data_3[which(fp_data_3$X == "S[2,108]"),c(6,4,8)]
doc_surv_model_4 <- fp_data_4[which(fp_data_4$X == "S[2,108]"),c(6,4,8)]
doc_surv_model_5 <- fp_data_5[which(fp_data_5$X == "S[2,108]"),c(6,4,8)]
doc_surv_models <- bind_rows( doc_surv_model_1,
                              doc_surv_model_2,
                              doc_surv_model_3,
                              doc_surv_model_4,
                              doc_surv_model_5)

doc_surv_models <- apply(doc_surv_models, 1, format_results)
write.csv(doc_surv_models, here("Fractional polynomials/data", "doc_surv_models_OS.csv"))

# 3.4 Abiratrone
abi_surv_model_1 <- fp_data_1[which(fp_data_1$X == "S[5,108]"),c(6,4,8)]
abi_surv_model_2 <- fp_data_2[which(fp_data_1$X == "S[5,108]"),c(6,4,8)]
abi_surv_model_3 <- fp_data_3[which(fp_data_1$X == "S[5,108]"),c(6,4,8)]
abi_surv_model_4 <- fp_data_4[which(fp_data_1$X == "S[5,108]"),c(6,4,8)]
abi_surv_model_5 <- fp_data_5[which(fp_data_1$X == "S[5,108]"),c(6,4,8)]
abi_surv_models <- bind_rows( abi_surv_model_1,
                              abi_surv_model_2,
                              abi_surv_model_3,
                              abi_surv_model_4,
                              abi_surv_model_5)

abi_surv_models <- apply(abi_surv_models, 1, format_results)
write.csv(abi_surv_models, here("Fractional polynomials/data", "abi_surv_models_OS.csv"))

# 3.5 All the treatments
# From sim14 model (P1: -0.5, P2: -0.5)
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
write.csv(surv_108, here("Fractional polynomials/data", "survival_108months_OS.csv"))

#4. Survival curve darolutamide ---------------------------
# FP extrapolation
graph_daro <- data.frame(time=c(0:108),
                         model_1=c(1, fp_data_1$mean[1:108]),
                         model_2=c(1, fp_data_2$mean[1:108]),
                         model_3=c(1, fp_data_3$mean[1:108]),
                         model_4=c(1, fp_data_4$mean[1:108]),
                         model_5=c(1, fp_data_5$mean[1:108]))
# KM curve
KM.est_daro <- survfit(Surv(time, event) ~ 1, data= arasens_daro, 
                       type="kaplan-meier", conf.int=FALSE)

#jpeg(file = here("04_figures", "extrapolation_darolutamide_OS.jpg"), width = 900, height = 800, res = 120)

plot(KM.est_daro, xlab="Time (months)", ylab="OS", xaxt="n", yaxt="n", main=" ",
     xlim = c(0,108), ylim=c(0,1),
     mark.time=FALSE, col=color[15], conf.int=F)

title("Overall survival (OS) \n2nd order FP models")
axis(2, at=c(0, 0.2, 0.4, 0.6, 0.8, 1))
axis(1, at= seq(0, 108, 20))

## Add prediction lines#
lines(graph_daro$time, graph_daro$model_1, col=color[1])
lines(graph_daro$time, graph_daro$model_2, col=color[2])
lines(graph_daro$time, graph_daro$model_3, col=color[3])
lines(graph_daro$time, graph_daro$model_4, col=color[4])
lines(graph_daro$time, graph_daro$model_5, col=color[5])

# legend
legend(x = 60, y = 1,
       c("P1: -0.5, P2: -0.5",
         "P1: -0.5, P2:  0",
         "P1: -1,   P2:  0.5",
         "P1: -1,   P2:  0",
         "P1: -1,   P2: -0.5 ",
         "ARASENS"),
       col=c(color[1], color[2], color[3], color[4], color[5], color[15]),
       lty=c(1,1,1,1), ncol=1, text.width=6, box.lty=0)

# save plot
#dev.off()


#5. Survival curve enza+ADT ---------------------------
#5.1 Second-order model with lowest DIC (5 models) ARCHES trial of interest
fp_data_1_ARCHES <- read.csv(here("Fractional polynomials/data", "OS_results_FP_sim14_ARCHES.csv"))  # model with P1: -0.5, P2: -0.5
fp_data_2_ARCHES <- read.csv(here("Fractional polynomials/data", "OS_results_FP_sim15_ARCHES.csv"))  # model with P1: -0.5, P2:  0
fp_data_3_ARCHES <- read.csv(here("Fractional polynomials/data", "OS_results_FP_sim11_ARCHES.csv"))  # model with P1: -1,   P2:  0.5
fp_data_4_ARCHES <- read.csv(here("Fractional polynomials/data", "OS_results_FP_sim10_ARCHES.csv"))  # model with P1: -1,   P2:  0
fp_data_5_ARCHES <- read.csv(here("Fractional polynomials/data", "OS_results_FP_sim09_ARCHES.csv"))  # model with P1: -1,   P2: -0.5


# FP extrapolation
graph_enza <- data.frame(time=c(0:108),
                         model_1=c(1, fp_data_1_ARCHES$mean[325:432]),
                         model_2=c(1, fp_data_2_ARCHES$mean[325:432]),
                         model_3=c(1, fp_data_3_ARCHES$mean[325:432]),
                         model_4=c(1, fp_data_4_ARCHES$mean[325:432]),
                         model_5=c(1, fp_data_5_ARCHES$mean[325:432]))

# IPDs

# ENZAMET
ENZAMET_enza <- read_csv(here("Fractional polynomials/data", "IPD reconstructed", "OS_ENZAMET_enzalutamide_ipd.csv"))

# Vaishampayan
enza_vai <- read_csv(here("Fractional polynomials/data", "IPD reconstructed", "OS_Vaishampayan 2021_enzalutamide_ipd.csv"))

# ARCHES
enza_arches <- read_csv(here("Fractional polynomials/data", "IPD reconstructed", "OS_ARCHES_ENZA+ADT_ipd.csv"))

# KM curves
KM.est_enza <- survfit(Surv(time, event) ~ 1, data= ENZAMET_enza, 
                       type="kaplan-meier", conf.int=FALSE)
KM.est_enza2 <- survfit(Surv(time, event) ~ 1, data= enza_vai, 
                       type="kaplan-meier", conf.int=FALSE)
KM.est_enza3 <- survfit(Surv(time, event) ~ 1, data= enza_arches, 
                       type="kaplan-meier", conf.int=FALSE)

#jpeg(file = here("04_figures", "extrapolation_enza+ADT_OS_v4.jpg"), width = 900, height = 800, res = 120)

plot(KM.est_enza, xlab="Time (months)", ylab="OS", xaxt="n", yaxt="n", main=" ",
     xlim = c(0,108), ylim=c(0,1),
     mark.time=FALSE, col=color[15], conf.int=F, lty = 2)

title("Overall survival (OS) \n2nd order FP models")
axis(2, at=c(0, 0.2, 0.4, 0.6, 0.8, 1))
axis(1, at= seq(0, 108, 20))

## Add prediction lines#
lines(KM.est_enza2, conf.int = F, col = color[15], lty = 3)
lines(KM.est_enza3, conf.int = F, col = color[15], lty = 1, lwd = 2)
lines(graph_enza$time, graph_enza$model_1, col=color[1])
lines(graph_enza$time, graph_enza$model_2, col=color[2])
lines(graph_enza$time, graph_enza$model_3, col=color[3])
lines(graph_enza$time, graph_enza$model_4, col=color[4])
lines(graph_enza$time, graph_enza$model_5, col=color[5])

# legend
legend(x = 60, y = 1,
       c("P1: -0.5, P2: -0.5",
         "P1: -0.5, P2:  0",
         "P1: -1,   P2:  0.5",
         "P1: -1,   P2:  0",
         "P1: -1,   P2: -0.5 ",
         "ENZAMET", "Vaishampayan 2021", "ARCHES"),
       col=c(color[1], color[2], color[3], color[4], color[5], color[15], color[15], color[15]),
       lty=c(1,1,1,1,1,2,3,1), lwd=c(1,1,1,1,1,1,1,2), ncol=1, text.width=6, box.lty=0)

# save plot
#dev.off()


#6. Survival curve docetaxel+ADT ---------------------------
# FP extrapolation
graph_doc <- data.frame(time=c(0:108),
                        model_1=c(1, fp_data_1$mean[109:216]),
                        model_2=c(1, fp_data_2$mean[109:216]),
                        model_3=c(1, fp_data_3$mean[109:216]),
                        model_4=c(1, fp_data_4$mean[109:216]),
                        model_5=c(1, fp_data_5$mean[109:216]))

# IPD data
# ARASENS
arasens_pcb <- filter(arasens, treatment == "Placebo+docetaxel arm")

# PEACE-1
peace <- read_csv(here("Fractional polynomials/data", "IPD reconstructed", "OS_PEACE_ADT+Doc_ipd.csv"))

# STAMPEDE-3
stamp3 <- read_csv(here("Fractional polynomials/data", "IPD reconstructed", "OS_STAMPEDE-3_ADT+Doc_ipd.csv"))

# STAMPEDE-4
stamp4 <- read_csv(here("Fractional polynomials/data", "IPD reconstructed", "OS_STAMPEDE-4_SOC+DocP_ipd.csv"))

# GETUG
getug <- read_csv(here("Fractional polynomials/data", "IPD reconstructed", "OS_GETUG_ADT+Doc_ipd.csv"))

# CHAARTED
chaart <- read_csv(here("Fractional polynomials/data", "IPD reconstructed", "OS_CHAARTED_ADT+Doc_ipd.csv"))

# KM curves
KM.est_doc <- survfit(Surv(time, event) ~ 1, data= peace, 
                       type="kaplan-meier", conf.int=FALSE)

KM.est_doc_arasens <- survfit(Surv(time, event) ~ 1, data= arasens_pcb, 
                       type="kaplan-meier", conf.int=FALSE)

KM.est_doc_stamp3 <- survfit(Surv(time, event) ~ 1, data= stamp3, 
                       type="kaplan-meier", conf.int=FALSE)

KM.est_doc_stamp4 <- survfit(Surv(time, event) ~ 1, data= stamp4, 
                       type="kaplan-meier", conf.int=FALSE)

KM.est_doc_getug <- survfit(Surv(time, event) ~ 1, data= getug, 
                       type="kaplan-meier", conf.int=FALSE)

KM.est_doc_peace <- survfit(Surv(time, event) ~ 1, data= peace, 
                       type="kaplan-meier", conf.int=FALSE)

KM.est_doc_chaart <- survfit(Surv(time, event) ~ 1, data= chaart, 
                       type="kaplan-meier", conf.int=FALSE)

#jpeg(file = here("04_figures", "extrapolation_doc+ADT_OS_v4.jpg"), width = 900, height = 800, res = 120)

plot(KM.est_doc, xlab="Time (months)", ylab="OS", xaxt="n", yaxt="n", main=" ",
     xlim = c(0,108), ylim=c(0,1),
     mark.time=FALSE, col=color[15], conf.int=F)

title("Overall survival (OS) \n2nd order FP models")
axis(2, at=c(0, 0.2, 0.4, 0.6, 0.8, 1))
axis(1, at= seq(0, 108, 20))

## Add prediction lines#
lines(graph_doc$time, graph_doc$model_1, col=color[1])
lines(graph_doc$time, graph_doc$model_2, col=color[2])
lines(graph_doc$time, graph_doc$model_3, col=color[3])
lines(graph_doc$time, graph_doc$model_4, col=color[4])
lines(graph_doc$time, graph_doc$model_5, col=color[5])
lines(KM.est_doc_arasens, conf.int = F, col = color[15], lty = 2, lwd = 2)
lines(KM.est_doc_getug, conf.int = F,   col = color[15], lty = 3)
lines(KM.est_doc_chaart, conf.int = F,  col = color[15], lty = 4)
lines(KM.est_doc_stamp3, conf.int = F,  col = color[15], lty = 5)
lines(KM.est_doc_stamp4, conf.int = F,  col = color[15], lty = 6)
lines(KM.est_doc_peace, conf.int = F,   col = color[15], lty = 7)


# legend
legend(x = 0, y = 0.5,
       c("P1: -0.5, P2: -0.5",
         "P1: -0.5, P2:  0",
         "P1: -1,   P2:  0.5",
         "P1: -1,   P2:  0",
         "P1: -1,   P2: -0.5 ", "ARASENS", "GETUG", "CHAARTED",
         "STAMPEDE arms C vs A", "STAMPEDE arms C vs G", "PEACE-1"),
       col=c(color[1], color[2], color[3], color[4], color[5], color[15], color[15], color[15],
             color[15], color[15], color[15]),
       lty=c(1,1,1,1,1,2,3,4,5,6,7), lwd=c(1,1,1,1,1,2,1,1,1,1,1),
       ncol=1, text.width=6, box.lty=0)

# save plot
#dev.off()

#7. Survival curve abi ---------------------------
# abi data
graph_abi <- data.frame(time=c(0:108),
                        model_1=c(1, fp_data_1$mean[433:540]),
                        model_2=c(1, fp_data_2$mean[433:540]),
                        model_3=c(1, fp_data_3$mean[433:540]),
                        model_4=c(1, fp_data_4$mean[433:540]),
                        model_5=c(1, fp_data_5$mean[433:540]))

# IPDs

# STAMPEDE-2
abi_stamp2 <- read_csv(here("Fractional polynomials/data", "IPD reconstructed", "OS_STAMPEDE-2_AA+ADT_ipd.csv"))

# STAMPEDE-4
abi_stamp4 <- read_csv(here("Fractional polynomials/data", "IPD reconstructed", "OS_STAMPEDE-4_SOC+AAP_ipd.csv"))

# LATITUDE
abi_lat <- read_csv(here("Fractional polynomials/data", "IPD reconstructed", "OS_LATITUDE_AA+P+ADT_ipd.csv"))


# KM curves
KM.est_abi <- survfit(Surv(time, event) ~ 1, data= abi_lat, 
                       type="kaplan-meier", conf.int=FALSE)
KM.est_abi2 <- survfit(Surv(time, event) ~ 1, data= abi_stamp2, 
                      type="kaplan-meier", conf.int=FALSE)
KM.est_abi3 <- survfit(Surv(time, event) ~ 1, data= abi_stamp4, 
                      type="kaplan-meier", conf.int=FALSE)

#jpeg(file = here("04_figures", "extrapolation_abi+ADT_OS_v4.jpg"), width = 900, height = 800, res = 120)

plot(KM.est_abi, xlab="Time (months)", ylab="OS", xaxt="n", yaxt="n", main=" ",
     xlim = c(0,108), ylim=c(0,1),
     mark.time=FALSE, col=color[15], conf.int=F, lwd = 2)

title("Overall survival (OS) \n2nd order FP models")
axis(2, at=c(0, 0.2, 0.4, 0.6, 0.8, 1))
axis(1, at= seq(0, 108, 20))

## Add prediction lines#
lines(KM.est_abi2, col = color[15], conf.int = F, lty = 2)
lines(KM.est_abi3, col = color[15], conf.int = F, lty = 3)
lines(graph_abi$time, graph_abi$model_1, col=color[1])
lines(graph_abi$time, graph_abi$model_2, col=color[2])
lines(graph_abi$time, graph_abi$model_3, col=color[3])
lines(graph_abi$time, graph_abi$model_4, col=color[4])
lines(graph_abi$time, graph_abi$model_5, col=color[5])

# legend
legend(x = 60, y = 1,
       c("P1: -0.5, P2: -0.5",
         "P1: -0.5, P2:  0",
         "P1: -1,   P2:  0.5",
         "P1: -1,   P2:  0",
         "P1: -1,   P2: -0.5 ",
         "LATITUDE", "STAMPEDE arms G vs A", "STAMPEDE arms C vs G"),
       col=c(color[1], color[2], color[3], color[4], color[5], color[15], color[15], color[15]),
       lty=c(1,1,1,1, 1, 1,2,3), lwd=c(1,1,1,1, 1, 2,1,1), 
       ncol=1, text.width=6, box.lty=0)

# save plot
#dev.off()


#8. Survival curve apa ---------------------------
# apa data
graph_apa <- data.frame(time=c(0:108),
                        model_1=c(1, fp_data_1$mean[541:648]),
                        model_2=c(1, fp_data_2$mean[541:648]),
                        model_3=c(1, fp_data_3$mean[541:648]),
                        model_4=c(1, fp_data_4$mean[541:648]),
                        model_5=c(1, fp_data_5$mean[541:648]))

#IPDs

# TITAN
apa <- read_csv(here("Fractional polynomials/data", "IPD reconstructed", "OS_TITAN_apalutamide_ipd.csv"))

# KM curve
KM.est_apa <- survfit(Surv(time, event) ~ 1, data= apa, 
                       type="kaplan-meier", conf.int=FALSE)

#jpeg(file = here("04_figures", "extrapolation_apalutamide.jpg"), width = 900, height = 800, res = 120)

plot(KM.est_apa, xlab="Time (months)", ylab="OS", xaxt="n", yaxt="n", main=" ",
     xlim = c(0,108), ylim=c(0,1),
     mark.time=FALSE, col=color[15], conf.int=F)

title("Overall survival (OS) \n2nd order FP models")
axis(2, at=c(0, 0.2, 0.4, 0.6, 0.8, 1))
axis(1, at= seq(0, 108, 20))

## Add prediction lines#
lines(graph_apa$time, graph_apa$model_1, col=color[1])
lines(graph_apa$time, graph_apa$model_2, col=color[2])
lines(graph_apa$time, graph_apa$model_3, col=color[3])
lines(graph_apa$time, graph_apa$model_4, col=color[4])
lines(graph_apa$time, graph_apa$model_5, col=color[5])

# legend
legend(x = 60, y = 1,
       c("P1: -0.5, P2: -0.5",
         "P1: -0.5, P2:  0",
         "P1: -1,   P2:  0.5",
         "P1: -1,   P2:  0",
         "P1: -1,   P2: -0.5",
         "TITAN"),
       col=c(color[1], color[2], color[3], color[4], color[5], color[15]),
       lty=c(1,1,1,1), ncol=1, text.width=6, box.lty=0)

# save plot
#dev.off()

#9. Model averaging----------------------------------
#9.1 Import model arrays  
fp_data_1_array <- read.csv(here("Fractional polynomials/data", "OS_array_FP_sim14.csv"))  # model with P1: -0.5, P2: -0.5
fp_data_2_array <- read.csv(here("Fractional polynomials/data", "OS_array_FP_sim15.csv"))  # model with P1: -0.5, P2:  0
fp_data_3_array <- read.csv(here("Fractional polynomials/data", "OS_array_FP_sim11.csv"))  # model with P1: -1,   P2:  0.5
fp_data_4_array <- read.csv(here("Fractional polynomials/data", "OS_array_FP_sim10.csv"))  # model with P1: -1,   P2:  0
fp_data_5_array <- read.csv(here("Fractional polynomials/data", "OS_array_FP_sim09.csv"))  # model with P1: -1,   P2: -0.5

#9.2 Darolutamide
summary.stats((fp_data_1_array[,2]+ fp_data_2_array[,2]+ fp_data_3_array[,2]+
                 fp_data_4_array[,2]+ fp_data_5_array[,2])/5)

#9.3 Docetaxel
summary.stats((fp_data_1_array[,3]+ fp_data_2_array[,3]+ fp_data_3_array[,3]+
                 fp_data_4_array[,3]+ fp_data_5_array[,3])/5)

#9.4 ADT
summary.stats((fp_data_1_array[,4]+ fp_data_2_array[,4]+ fp_data_3_array[,4]+
                 fp_data_4_array[,4]+ fp_data_5_array[,4])/5)

#9.5 Enzalutamide
summary.stats((fp_data_1_array[,5]+ fp_data_2_array[,5]+ fp_data_3_array[,5]+
                 fp_data_4_array[,5]+ fp_data_5_array[,5])/5)

#9.6 Abiraterone
summary.stats((fp_data_1_array[,6]+ fp_data_2_array[,6]+ fp_data_3_array[,6]+
                 fp_data_4_array[,6]+ fp_data_5_array[,6])/5)

#9.7 Apalutamide
summary.stats((fp_data_1_array[,7]+ fp_data_2_array[,7]+ fp_data_3_array[,7]+
                 fp_data_4_array[,7]+ fp_data_5_array[,7])/5)

#9.8 Enza+docetaxel
summary.stats((fp_data_1_array[,8]+ fp_data_2_array[,8]+ fp_data_3_array[,8]+
                 fp_data_4_array[,8]+ fp_data_5_array[,8])/5)

#9.9 Abiraterone+docetaxel
summary.stats((fp_data_1_array[,9]+ fp_data_2_array[,9]+ fp_data_3_array[,9]+
                 fp_data_4_array[,9]+ fp_data_5_array[,9])/5)

for (i in arrays) {
  arrays[[i]] <- 
}
arrays <- list(
fp_data_1_array,
fp_data_2_array,
fp_data_3_array,
fp_data_4_array,
fp_data_5_array)

arrays
