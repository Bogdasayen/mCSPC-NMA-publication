# Script to create data ready for FP analyses 
# PFS outcome

# 1. Packages used --------------------------------------------------------
pkgs <- c("tidyverse", "readxl", "here", "survival", "survminer", "Hmisc")
lapply(pkgs, library, character.only = T)
source(here("code", "set_FP_data.R")) 

# 2. Import data ---------------------------------------------------------
km_data_directory <- here("02_data", "KM digitized")
ipd_directory <- here("02_data", "IPD reconstructed")
FP_data_directory <- here("02_data", "FP data")
outcome_names <- "PFS"
study_names <- list()
treatment_names <- list()
FP_data <- list()

#2.1 ARCHES   -------------------------------------------------------------  
study_names[["PFS"]] <- c("ARCHES")
treatment_names[["ARCHES"]] <- c("ADT", "ENZA+ADT")
set_FP_data(outcome_names   = outcome_names, 
            study_names     = study_names, 
            treatment_names = treatment_names)

#2.2 CHAARTED   -------------------------------------------------------------  
study_names[["PFS"]] <- c("CHAARTED")
treatment_names[["CHAARTED"]] <- c("ADT", "ADT+Doc")
set_FP_data(outcome_names   = outcome_names, 
            study_names     = study_names, 
            treatment_names = treatment_names)

#2.3 ENZAMET  Doc -------------------------------------------------------------  
study_names[["PFS"]] <- c("ENZAMET")
treatment_names[["ENZAMET"]] <- c("ADT+Doc", "enzalutamide+Doc")
set_FP_data(outcome_names   = outcome_names, 
            study_names     = study_names, 
            treatment_names = treatment_names)

#2.4 ENZAMET  noDoc -------------------------------------------------------------  
study_names[["PFS"]] <- c("ENZAMET")
treatment_names[["ENZAMET"]] <- c("ADT", "enzalutamide")
set_FP_data(outcome_names   = outcome_names, 
            study_names     = study_names, 
            treatment_names = treatment_names)

#2.5 GETUG   -------------------------------------------------------------  
study_names[["PFS"]] <- c("GETUG")
treatment_names[["GETUG"]] <- c("ADT", "ADT+Doc")
set_FP_data(outcome_names   = outcome_names, 
            study_names     = study_names, 
            treatment_names = treatment_names)

#2.6 LATITUDE   -------------------------------------------------------------  
study_names[["PFS"]] <- c("LATITUDE")
treatment_names[["LATITUDE"]] <- c("ADT", "AA+P+ADT")
set_FP_data(outcome_names   = outcome_names, 
            study_names     = study_names, 
            treatment_names = treatment_names)

#2.7 PEACE-1   -------------------------------------------------------------  
study_names[["PFS"]] <- c("PEACE")
treatment_names[["PEACE"]] <- c("ADT+Doc", "AA+ADT+Doc")
set_FP_data(outcome_names   = outcome_names, 
            study_names     = study_names, 
            treatment_names = treatment_names)

#2.8 STAMPEDE-3   -------------------------------------------------------------  
study_names[["PFS"]] <- c("STAMPEDE-3")
treatment_names[["STAMPEDE-3"]] <- c("ADT", "ADT+Doc")
set_FP_data(outcome_names   = outcome_names, 
            study_names     = study_names, 
            treatment_names = treatment_names)

#2.9 STAMPEDE-4   -------------------------------------------------------------  
study_names[["PFS"]] <- c("STAMPEDE-4")
treatment_names[["STAMPEDE-4"]] <- c("ADT+DocP", "ADT+AAP")
set_FP_data(outcome_names   = outcome_names, 
            study_names     = study_names, 
            treatment_names = treatment_names)

#2.10 TITAN   -------------------------------------------------------------  
study_names[["PFS"]] <- c("TITAN")
treatment_names[["TITAN"]] <- c("ADT", "apalutamide")
set_FP_data(outcome_names   = outcome_names, 
            study_names     = study_names, 
            treatment_names = treatment_names)

#2.11 Vaishampayan 2021   -------------------------------------------------------------  
study_names[["PFS"]] <- c("Vaishampayan 2021")
treatment_names[["Vaishampayan 2021"]] <- c("bicalutamide", "enzalutamide")
set_FP_data(outcome_names   = outcome_names, 
            study_names     = study_names, 
            treatment_names = treatment_names)

#2.12 ARASENS ----------------------------------------------
arasens <- read_csv(here("02_data", "arasens_PFS.csv"))
arasens_dara <- filter(arasens, treatment == "Darolutamide+docetaxel arm")
arasens_pcb <- filter(arasens, treatment == "Placebo+docetaxel arm")

# Daralutamide arm
# Based on CROD_ITT_km.jpeg file on SharePoint (therefore timepoints every 10 months)
natrisk_dara <- read_excel(here("02_data", "KM digitized", "ARASENS daralutamide PFS.xlsx"), 
                           sheet = 2)
s_dara <- rep(1, nrow(natrisk_dara)-1)
a_dara <- rep(1, nrow(natrisk_dara)-1)
dt_dara <- rep(10, nrow(natrisk_dara)-1)
cutpoints <- seq(from = 0, to = 60, by = 10)
r_dara <- arasens_dara %>% mutate(ranges = cut2(time, cutpoints)) %>% # number of events by time interval 
  group_by(ranges) %>% 
  summarise(r = sum(event))# %>% select(r)
r_dara <- r_dara$r
z_dara <- natrisk_dara$Natrisk[-length(natrisk_dara$Natrisk)] # number of patients at risk at each interval
time_j_dara <- natrisk_dara$Time[-1] # start of every time interval

FP_data_dara <- qpcR:::cbind.na(s_dara, a_dara, dt_dara, r_dara, z_dara, time_j_dara)
colnames(FP_data_dara) <- c("s", "a", "dt", "r", "z", "time_j")
write.csv(FP_data_dara[-nrow(FP_data_dara),],
          here("02_data", "FP data", "PFS_ARASENS_daralutamide_FPdata.csv"))

# Placebo arm
natrisk_pcb <- read_excel(here("02_data", "KM digitized", "ARASENS ADT+Doc PFS.xlsx"), 
                          sheet = 2)
s_pcb <- rep(1, nrow(natrisk_pcb)-1)
a_pcb <- rep(1, nrow(natrisk_pcb)-1)
dt_pcb <- rep(10, nrow(natrisk_pcb)-1)
cutpoints_pcb <- seq(from = 0, to = 60, by = 10)
r_pcb <- arasens_pcb %>% mutate(ranges = cut2(time, cutpoints_pcb)) %>% # number of events by time interval 
  group_by(ranges) %>% 
  summarise(r = sum(event))# %>% select(r)
r_pcb <- r_pcb$r
z_pcb <- natrisk_pcb$Natrisk[-length(natrisk_pcb$Natrisk)] # number of patients at risk at each interval
time_j_pcb <- natrisk_pcb$Time[-1] # start of every time interval

FP_data_pcb <- qpcR:::cbind.na(s_pcb, a_pcb, dt_pcb, r_pcb, z_pcb, time_j_pcb)
FP_data_pcb <- FP_data_pcb[-7,]
colnames(FP_data_pcb) <- c("s", "a", "dt", "r", "z", "time_j")
write.csv(FP_data_pcb, here("02_data", "FP data", "PFS_ARASENS_ADT+Doc_FPdata.csv"))

#3. Merge all datasets-------------------------------------
# ARASENS; arm1: Daro; arm2: ADT+Doc
study1_a1 <- read_csv(here("02_data", "FP data", 
                           "PFS_ARASENS_daralutamide_FPdata.csv"), col_select = -1)
study1_a1

study1_a2 <- read_csv(here("02_data", "FP data", 
                           "PFS_ARASENS_ADT+Doc_FPdata.csv"), col_select = -1)
study1_a2$a <- 2

# ARCHES; arm1: Enza; arm2: ADT
study2_a1 <- read_csv(here("02_data", "FP data", 
                           "PFS_ARCHES_ENZA+ADT_FPdata.csv"), col_select = -1)
study2_a1$s <- 2

study2_a2 <- read_csv(here("02_data", "FP data", 
                           "PFS_ARCHES_ADT_FPdata.csv"), col_select = -1)
study2_a2$s <- 2
study2_a2$a <- 2

# CHAARTED; arm1: Doc+ADT; arm2: ADT
study3_a1 <- read_csv(here("02_data", "FP data", 
                           "PFS_CHAARTED_ADT+Doc_FPdata.csv"), col_select = -1)
study3_a1$s <- 3

study3_a2 <- read_csv(here("02_data", "FP data", 
                           "PFS_CHAARTED_ADT_FPdata.csv"), col_select = -1)
study3_a2$s <- 3
study3_a2$a <- 2

# ENZAMET no doc; arm1: enza+ADT, arm2: ADT
study4_a1 <- read_csv(here("02_data", "FP data", 
                           "PFS_ENZAMET_enzalutamide_FPdata.csv"), col_select = -1)
study4_a1$s <- 4

study4_a2 <- read_csv(here("02_data", "FP data", 
                           "PFS_ENZAMET_ADT_FPdata.csv"), col_select = -1)
study4_a2$s <- 4
study4_a2$a <- 2

# ENZAMET doc; arm1: doc+ADT vs arm2: Enzalutamide + Doc + ADT
study5_a1 <- read_csv(here("02_data", "FP data", 
                           "PFS_ENZAMET_ADT+Doc_FPdata.csv"), col_select = -1)
study5_a1$s <- 5

study5_a2 <- read_csv(here("02_data", "FP data", 
                           "PFS_ENZAMET_enzalutamide+Doc_FPdata.csv"), col_select = -1)
study5_a2$s <- 5
study5_a2$a <- 2

# GETUG; arm1: ADT+Doc; arm2: ADT
study6_a1 <- read_csv(here("02_data", "FP data", 
                           "PFS_GETUG_ADT+Doc_FPdata.csv"), col_select = -1)
study6_a1$s <- 6

study6_a2 <- read_csv(here("02_data", "FP data", 
                           "PFS_GETUG_ADT_FPdata.csv"), col_select = -1)
study6_a2$s <- 6
study6_a2$a <- 2

# LATITUDE; arm1: abi; arm2: ADT
study7_a1 <- read_csv(here("02_data", "FP data", 
                           "PFS_LATITUDE_AA+P+ADT_FPdata.csv"), col_select = -1)
study7_a1$s <- 7

study7_a2 <- read_csv(here("02_data", "FP data", 
                           "PFS_LATITUDE_ADT_FPdata.csv"), col_select = -1)
study7_a2$s <- 7
study7_a2$a <- 2

# PEACE, arm1: Doc+ADT vs arm2: AA+Doc+ADT
study8_a1 <- read_csv(here("02_data", "FP data", 
                           "PFS_PEACE_ADT+Doc_FPdata.csv"), col_select = -1)
study8_a1$s <- 8

study8_a2 <- read_csv(here("02_data", "FP data", 
                           "PFS_PEACE_AA+ADT+Doc_FPdata.csv"), col_select = -1)
study8_a2$s <- 8
study8_a2$a <- 2

# STAMPEDE-3; arm1: Doc+ADT; arm2: ADT
study9_a1 <- read_csv(here("02_data", "FP data", 
                            "PFS_STAMPEDE-3_ADT+Doc_FPdata.csv"), col_select = -1)
study9_a1$s <- 9

study9_a2 <- read_csv(here("02_data", "FP data", 
                            "PFS_STAMPEDE-3_ADT_FPdata.csv"), col_select = -1)
study9_a2$s <- 9
study9_a2$a <- 2

# STAMPEDE-4; arm1: abi+ADT; arm2: Doc+ADT
study10_a1 <- read_csv(here("02_data", "FP data", 
                           "PFS_STAMPEDE-4_ADT+AAP_FPdata.csv"), col_select = -1)
study10_a1$s <- 10

study10_a2 <- read_csv(here("02_data", "FP data", 
                           "PFS_STAMPEDE-4_ADT+DocP_FPdata.csv"), col_select = -1)
study10_a2$s <- 10
study10_a2$a <- 2

# TITAN; arm1: apalutamide; arm2: ADT
study11_a1 <- read_csv(here("02_data", "FP data", 
                            "PFS_TITAN_apalutamide_FPdata.csv"), col_select = -1)
study11_a1$s <- 11

study11_a2 <- read_csv(here("02_data", "FP data", 
                            "PFS_TITAN_ADT_FPdata.csv"), col_select = -1)
study11_a2$s <- 11
study11_a2$a <- 2

# Vaishampayan 2021; arm1: enza+ADT; arm2: ADT [double-check this with Howard/Philip]
study12_a1 <- read_csv(here("02_data", "FP data", 
                            "PFS_Vaishampayan 2021_enzalutamide_FPdata.csv"), col_select = -1)
study12_a1$s <- 12

study12_a2 <- read_csv(here("02_data", "FP data", 
                            "PFS_Vaishampayan 2021_bicalutamide_FPdata.csv"), col_select = -1)
study12_a2$s <- 12
study12_a2$a <- 2

#4 Export dataset to be analysed-------------------------
studies <- rbind(study1_a1, study1_a2, study2_a1, study2_a2, study3_a1, study3_a2, 
                 study4_a1, study4_a2, study5_a1, study5_a2, study6_a1, study6_a2,
                 study7_a1, study7_a2, study8_a1, study8_a2, study9_a1, study9_a2,
                 study10_a1, study10_a2, study11_a1, study11_a2, study12_a1, study12_a2)
write_csv(studies, here("02_data", "FP data", "FP_data_PFS.csv"))
