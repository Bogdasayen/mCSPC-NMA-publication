# Script to create data ready for FP analyses 
# OS outcome

# 1. Packages used --------------------------------------------------------
pkgs <- c("tidyverse", "readxl", "here", "survival", "survminer", "Hmisc")
lapply(pkgs, library, character.only = T)
source(here("code", "set_FP_data.R")) 

# 2. Import data ---------------------------------------------------------
km_data_directory <- here("02_data", "KM digitized")
ipd_directory <- here("02_data", "IPD reconstructed")
FP_data_directory <- here("02_data", "FP data")
outcome_names <- "OS"
study_names <- list()
treatment_names <- list()
FP_data <- list()

#1 ENZAMET  Doc -------------------------------------------------------------  
study_names[["OS"]] <- c("ENZAMET")
treatment_names[["ENZAMET"]] <- c("ADT+Doc", "enzalutamide+Doc")
set_FP_data(outcome_names   = outcome_names, 
            study_names     = study_names, 
            treatment_names = treatment_names)

#2 ENZAMET  noDoc -------------------------------------------------------------  
study_names[["OS"]] <- c("ENZAMET")
treatment_names[["ENZAMET"]] <- c("ADT", "enzalutamide")
set_FP_data(outcome_names   = outcome_names, 
            study_names     = study_names, 
            treatment_names = treatment_names)

#3 CHAARTED   -------------------------------------------------------------  
study_names[["OS"]] <- c("CHAARTED")
treatment_names[["CHAARTED"]] <- c("ADT", "ADT+Doc")
set_FP_data(outcome_names   = outcome_names, 
            study_names     = study_names, 
            treatment_names = treatment_names)

#4 GETUG   -------------------------------------------------------------  
study_names[["OS"]] <- c("GETUG")
treatment_names[["GETUG"]] <- c("ADT", "ADT+Doc")
set_FP_data(outcome_names   = outcome_names, 
            study_names     = study_names, 
            treatment_names = treatment_names)

#5 LATITUDE   -------------------------------------------------------------  
study_names[["OS"]] <- c("LATITUDE")
treatment_names[["LATITUDE"]] <- c("ADT", "AA+P+ADT")
set_FP_data(outcome_names   = outcome_names, 
            study_names     = study_names, 
            treatment_names = treatment_names)

#6 ARCHES   -------------------------------------------------------------  
study_names[["OS"]] <- c("ARCHES")
treatment_names[["ARCHES"]] <- c("ADT", "ENZA+ADT")
set_FP_data(outcome_names   = outcome_names, 
            study_names     = study_names, 
            treatment_names = treatment_names)

#7 PEACE-1   -------------------------------------------------------------  
study_names[["OS"]] <- c("PEACE")
treatment_names[["PEACE"]] <- c("ADT+Doc", "AA+ADT+Doc")
set_FP_data(outcome_names   = outcome_names, 
            study_names     = study_names, 
            treatment_names = treatment_names)

#8 STAMPEDE-2   -------------------------------------------------------------  
study_names[["OS"]] <- c("STAMPEDE-2")
treatment_names[["STAMPEDE-2"]] <- c("AA+ADT", "ADT")
set_FP_data(outcome_names   = outcome_names, 
            study_names     = study_names, 
            treatment_names = treatment_names)

#9 STAMPEDE-3   -------------------------------------------------------------  
study_names[["OS"]] <- c("STAMPEDE-3")
treatment_names[["STAMPEDE-3"]] <- c("ADT", "ADT+Doc")
set_FP_data(outcome_names   = outcome_names, 
            study_names     = study_names, 
            treatment_names = treatment_names)

#10 STAMPEDE-4   -------------------------------------------------------------  
study_names[["OS"]] <- c("STAMPEDE-4")
treatment_names[["STAMPEDE-4"]] <- c("SOC+DocP", "SOC+AAP")
set_FP_data(outcome_names   = outcome_names, 
            study_names     = study_names, 
            treatment_names = treatment_names)


#11 TITAN   -------------------------------------------------------------  
study_names[["OS"]] <- c("TITAN")
treatment_names[["TITAN"]] <- c("ADT", "apalutamide")
set_FP_data(outcome_names   = outcome_names, 
            study_names     = study_names, 
            treatment_names = treatment_names)

#12 Vaishampayan 2021   -------------------------------------------------------------  
study_names[["OS"]] <- c("Vaishampayan 2021")
treatment_names[["Vaishampayan 2021"]] <- c("bicalutamide", "enzalutamide")
set_FP_data(outcome_names   = outcome_names, 
            study_names     = study_names, 
            treatment_names = treatment_names)

#13 ARASENS ----------------------------------------------
arasens <- read_csv(here("02_data", "arasens_OS.csv"))
arasens_daro <- filter(arasens, treatment == "Darolutamide+docetaxel arm")
arasens_pcb <- filter(arasens, treatment == "Placebo+docetaxel arm")

# Darolutamide arm
natrisk_dara <- read_excel(here("02_data", "KM digitized", "ARASENS daralutamide OS.xlsx"), 
                           sheet = 2)
s_dara <- rep(1, nrow(natrisk_dara)-1)
a_dara <- rep(1, nrow(natrisk_dara)-1)
dt_dara <- rep(3, nrow(natrisk_dara)-1)
cutpoints <- seq(from = 0, to = 54, by = 3)
r_dara <- arasens_daro %>% mutate(ranges = cut2(time, cutpoints)) %>% # number of events by time interval 
  group_by(ranges) %>% 
  summarise(r = sum(event))# %>% select(r)
r_dara <- r_dara$r
z_dara <- natrisk_dara$Natrisk[-length(natrisk_dara$Natrisk)] # number of patients at risk at each interval
time_j_dara <- natrisk_dara$Time[-1] # start of every time interval

FP_data_dara <- qpcR:::cbind.na(s_dara, a_dara, dt_dara, r_dara, z_dara, time_j_dara)
colnames(FP_data_dara) <- c("s", "a", "dt", "r", "z", "time_j")
write.csv(FP_data_dara[-nrow(FP_data_dara),],
          here("02_data", "FP data", "OS_ARASENS_daralutamide_FPdata.csv"))

# Placebo arm
natrisk_pcb <- read_excel(here("02_data", "KM digitized", "ARASENS ADT+Doc OS.xlsx"), 
                           sheet = 2)
s_pcb <- rep(1, nrow(natrisk_pcb)-1)
a_pcb <- rep(1, nrow(natrisk_pcb)-1) # this should be 2
dt_pcb <- rep(3, nrow(natrisk_pcb)-1)
cutpoints_pcb <- seq(from = 0, to = 57, by = 3)
r_pcb <- arasens_pcb %>% mutate(ranges = cut2(time, cutpoints_pcb)) %>% # number of events by time interval 
  group_by(ranges) %>% 
  summarise(r = sum(event))# %>% select(r)
r_pcb <- r_pcb$r
z_pcb <- natrisk_pcb$Natrisk[-length(natrisk_pcb$Natrisk)] # number of patients at risk at each interval
time_j_pcb <- natrisk_pcb$Time[-1] # start of every time interval

FP_data_pcb <- qpcR:::cbind.na(s_pcb, a_pcb, dt_pcb, r_pcb, z_pcb, time_j_pcb)
FP_data_pcb <- FP_data_pcb[-c(20,21),]
colnames(FP_data_pcb) <- c("s", "a", "dt", "r", "z", "time_j")
write.csv(FP_data_pcb, here("02_data", "FP data", "OS_ARASENS_ADT+Doc_FPdata.csv"))

#11. Merge all datasets-------------------------------------
# ARASENS; arm1: Daro; arm2: ADT+Doc
study1_a1 <- read_csv(here("02_data", "FP data", 
                           "OS_ARASENS_daralutamide_FPdata.csv"), col_select = -1)
study1_a1

study1_a2 <- read_csv(here("02_data", "FP data", 
                           "OS_ARASENS_ADT+Doc_FPdata.csv"), col_select = -1)
study1_a2$a <- 2

# ARCHES; arm1: enza+ADT; arm2: ADT
study2_a1 <- read_csv(here("02_data", "FP data", 
                           "OS_ARCHES_ENZA+ADT_FPdata.csv"), col_select = -1)
study2_a1$s <- 2

study2_a2 <- read_csv(here("02_data", "FP data", 
                           "OS_ARCHES_ADT_FPdata.csv"), col_select = -1)
study2_a2$s <- 2
study2_a2$a <- 2

# CHAARTED; arm1: Doc+ADT; arm2: ADT
study3_a1 <- read_csv(here("02_data", "FP data", 
                           "OS_CHAARTED_ADT+Doc_FPdata.csv"), col_select = -1)
study3_a1$s <- 3

study3_a2 <- read_csv(here("02_data", "FP data", 
                           "OS_CHAARTED_ADT_FPdata.csv"), col_select = -1)
study3_a2$s <- 3
study3_a2$a <- 2

# ENZAMET; arm1:enza+ADT; arm2:ADT
study4_a1 <- read_csv(here("02_data", "FP data", 
                           "OS_ENZAMET_enzalutamide_FPdata.csv"), col_select = -1)
study4_a1$s <- 4

study4_a2 <- read_csv(here("02_data", "FP data", 
                           "OS_ENZAMET_ADT_FPdata.csv"), col_select = -1)
study4_a2$s <- 4
study4_a2$a <- 2

# ENZAMET Doc; arm1:ADT+Doc; arm2: enza+Doc
study5_a1 <- read_csv(here("02_data", "FP data", 
                           "OS_ENZAMET_ADT+Doc_FPdata.csv"), col_select = -1)
study5_a1$s <- 5

study5_a2 <- read_csv(here("02_data", "FP data", 
                           "OS_ENZAMET_enzalutamide+Doc_FPdata.csv"), col_select = -1)
study5_a2$s <- 5
study5_a2$a <- 2

# GETUG; arm1: ADT+Doc; arm2: ADT
study6_a1 <- read_csv(here("02_data", "FP data", 
                           "OS_GETUG_ADT+Doc_FPdata.csv"), col_select = -1)
study6_a1$s <- 6

study6_a2 <- read_csv(here("02_data", "FP data", 
                           "OS_GETUG_ADT_FPdata.csv"), col_select = -1)
study6_a2$s <- 6
study6_a2$a <- 2

# LATITUDE; arm1: abi; arm2: ADT
study7_a1 <- read_csv(here("02_data", "FP data", 
                           "OS_LATITUDE_AA+P+ADT_FPdata.csv"), col_select = -1)
study7_a1$s <- 7

study7_a2 <- read_csv(here("02_data", "FP data", 
                           "OS_LATITUDE_ADT_FPdata.csv"), col_select = -1)
study7_a2$s <- 7
study7_a2$a <- 2

# PEACE; arm1: ADT+Doc; arm2: abi
study8_a1 <- read_csv(here("02_data", "FP data", 
                           "OS_PEACE_ADT+Doc_FPdata.csv"), col_select = -1)
study8_a1$s <- 8

study8_a2 <- read_csv(here("02_data", "FP data", 
                           "OS_PEACE_AA+ADT+Doc_FPdata.csv"), col_select = -1)
study8_a2$s <- 8
study8_a2$a <- 2

# STAMPEDE-2; arm1: abi+ADT; arm2: ADT
study9_a1 <- read_csv(here("02_data", "FP data", 
                           "OS_STAMPEDE-2_AA+ADT_FPdata.csv"), col_select = -1)
study9_a1$s <- 9

study9_a2 <- read_csv(here("02_data", "FP data", 
                           "OS_STAMPEDE-2_ADT_FPdata.csv"), col_select = -1)
study9_a2$s <- 9
study9_a2$a <- 2

# STAMPEDE-3; arm1: Doc+ADT; arm2: ADT
study10_a1 <- read_csv(here("02_data", "FP data", 
                           "OS_STAMPEDE-3_ADT+Doc_FPdata.csv"), col_select = -1)
study10_a1$s <- 10

study10_a2 <- read_csv(here("02_data", "FP data", 
                           "OS_STAMPEDE-3_ADT_FPdata.csv"), col_select = -1)
study10_a2$s <- 10
study10_a2$a <- 2

# STAMPEDE-4; arm1: abi+ADT; arm2: Doc+ADT
study11_a1 <- read_csv(here("02_data", "FP data", 
                            "OS_STAMPEDE-4_SOC+AAP_FPdata.csv"), col_select = -1)
study11_a1$s <- 11

study11_a2 <- read_csv(here("02_data", "FP data", 
                            "OS_STAMPEDE-4_SOC+DocP_FPdata.csv"), col_select = -1)
study11_a2$s <- 11
study11_a2$a <- 2

# TITAN; arm1: apalutamide; arm2: ADT
study12_a1 <- read_csv(here("02_data", "FP data", 
                           "OS_TITAN_apalutamide_FPdata.csv"), col_select = -1)
study12_a1$s <- 12

study12_a2 <- read_csv(here("02_data", "FP data", 
                           "OS_TITAN_ADT_FPdata.csv"), col_select = -1)
study12_a2$s <- 12
study12_a2$a <- 2

# Vaishampayan 2021; arm1: enza+ADT; arm2: ADT [double-check this with Howard/Philip]
study13_a1 <- read_csv(here("02_data", "FP data", 
                           "OS_Vaishampayan 2021_enzalutamide_FPdata.csv"), col_select = -1)
study13_a1$s <- 13

study13_a2 <- read_csv(here("02_data", "FP data", 
                           "OS_Vaishampayan 2021_bicalutamide_FPdata.csv"), col_select = -1)
study13_a2$s <- 13
study13_a2$a <- 2

# Create 1 dataset to be analysed
studies <- rbind(study1_a1, study1_a2, study2_a1, study2_a2, study3_a1, study3_a2,
                 study4_a1, study4_a2, study5_a1, study5_a2, study6_a1, study6_a2,
                 study7_a1, study7_a2, study8_a1, study8_a2, study9_a1, study9_a2,
                 study10_a1, study10_a2, study11_a1, study11_a2, study12_a1, study12_a2,
                 study13_a1, study13_a2)
write_csv(studies, here("02_data", "FP data", "FP_data_OS.csv"))
