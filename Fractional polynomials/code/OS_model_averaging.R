# Code for model averaging
# Outcome: OS

# 1. Packages used --------------------------------------------------------
pkgs <- c("tidyverse", "readxl", "here", "haven", "R2OpenBUGS")
lapply(pkgs, library, character.only = T)
set.seed(03082022)
source(here("code", "utils.R")) 

# 2. Import data --------------------------------------------------------
# These data is only to match treatment codes, not actual HRs
OS_data <- read_excel(here("02_data", "mHSPC OS trt codes.xlsx"))

# Fractional polynomial data created in Os_data_preparation.R
FP_data <- read_csv(here("02_data", "FP data", "FP_data_OS.csv"))
which(FP_data$z - FP_data$r < 0) # check if any r is higher than natrisk at each interval

# Transform data into a list for 1st order FP
first_FP <- list(
  P1   = 0,                          # P1 for 1st order FP
  N    = nrow(FP_data),             # intervals derived from the reconstructed IPD
  nt   = max(OS_data$t1, OS_data$t2),                         # number of treatments
  ns   = nrow(OS_data),                         # number of studies
  maxt = max(FP_data$time_j),       # maximum follow-up (months)
  mean = c(0, 0),                    # priors for the mean effect
  prec = structure(.Data= c(1.00000E-04, 0, 0, 1.00000E-04), 
                   .Dim=c(2, 2)),   # priors for the precision
  t = matrix(c(OS_data$t1, OS_data$t2), # treatments IDs
             nrow = nrow(OS_data), ncol = 2),
  na = OS_data$na,               # number of arms
  s = FP_data$s,                # study ID
  r = FP_data$r,                # number of events
  z = FP_data$z,                # patients at risk
  a = FP_data$a,                # treatment arm
  time = FP_data$time_j,        # time (months)
  dt = FP_data$dt               # interval length (months)
)

first_FP

# Transform data for 2nd order FP model
second_FP <- first_FP
second_FP$P2 <- 0
second_FP$maxt <- 60 # 60 months requested by Bayer
second_FP$mean <- c(0, 0, 0)
second_FP$prec <- structure(.Data= c(1.00000E-04, 0, 0,
                                     0, 1.00000E-02, 0,
                                     0, 0, 1.00000E-02), .Dim=c(3, 3))

# call the model
second_FP_model_FE <- here("03_models", "second_FP_model_FE.txt")

#Initial Values 
second_inits <- function(){
  #chain 1
  list(d= structure(.Data= rep(c(NA, rep(-0.1, first_FP$nt - 1)), 3),
                    .Dim=c(first_FP$nt, 3)), 
       mu= structure(.Data= c(rep(0.1, 3 * first_FP$ns)), 
                     .Dim=c(first_FP$ns, 3)))
  #chain 2
  list(d= structure(.Data= rep(c(NA, rep(1, first_FP$nt - 1)), 3),
                    .Dim=c(first_FP$nt, 3)),  
       mu = structure(.Data= c(rep(1, 3 * first_FP$ns)),
                      .Dim=c(first_FP$ns, 3)))
}



#1.Hazard ratios at timepoints of interest----------------
HR_params <- c("HR_1year", "HR_2year", "HR_3year", "HR_4year", "HR_5year", "RMST")

#1.1 Second-order FP, FE model, P1 = -0.5, P2 = -0.5 ----------------
second_FP$P1 <- -0.5
second_FP$P2 <- -0.5
second_FP_sim1_HR <- bugs(data = second_FP, inits = second_inits, model.file = second_FP_model_FE, 
                         n.iter = 40000, parameters.to.save = HR_params, n.burnin = 20000, 
                         n.chains = 2, 
                         debug = FALSE)

# Extract HR at timepoints of interest for the base case model
HR_export_36 <- as_tibble(round(1/second_FP_sim1_HR$summary[grep("HR_3year", rownames(second_FP_sim1_HR$summary)), c(1, 7, 3)], digits = 2))
HR_export_36 <- mutate(HR_export_36, string = paste0(mean, " (", `97.5%` ,", ", `2.5%`, ")"))
HR_export_48 <- as_tibble(round(1/second_FP_sim1_HR$summary[grep("HR_4year", rownames(second_FP_sim1_HR$summary)), c(1, 7, 3)], digits = 2))
HR_export_48 <- mutate(HR_export_48, string = paste0(mean, " (", `97.5%` ,", ", `2.5%`, ")"))
HR_export_60 <- as_tibble(round(1/second_FP_sim1_HR$summary[grep("HR_5year", rownames(second_FP_sim1_HR$summary)), c(1, 7, 3)], digits = 2))
HR_export_60 <- mutate(HR_export_60, string = paste0(mean, " (", `97.5%` ,", ", `2.5%`, ")"))

HR_export <- bind_cols(HR_export_36, HR_export_48, HR_export_60)
write_csv(HR_export, here("05_tables", "OS_HR_basecase.csv"))

# Select a sample of MCMC chains to produce model averaging results
second_FP_sim1_HR_array <- second_FP_sim1_HR$sims.array[,1,]
second_FP_sim1_HR_array <- apply(second_FP_sim1_HR_array, 2, sample, size = 500)

#1.2 Second-order FP, FE model, P1 = -0.5, P2 = 0 ----------------
second_FP$P1 <- -0.5
second_FP$P2 <- 0
second_FP_sim2_HR <- bugs(data = second_FP, inits = second_inits, model.file = second_FP_model_FE, 
                         n.iter = 40000, parameters.to.save = HR_params, n.burnin = 20000, 
                         n.chains = 2, 
                         debug = FALSE)
second_FP_sim2_HR_array <- second_FP_sim2_HR$sims.array[,1,]
second_FP_sim2_HR_array <- apply(second_FP_sim2_HR_array, 2, sample, size = 500)

#1.3 Second-order FP, FE model, P1 = -1, P2 = 0.5 ----------------
second_FP$P1 <- -1
second_FP$P2 <- 0.5
second_FP_sim3_HR <- bugs(data = second_FP, inits = second_inits, model.file = second_FP_model_FE, 
                         n.iter = 40000, parameters.to.save = HR_params, n.burnin = 20000, 
                         n.chains = 2, 
                         debug = FALSE)
second_FP_sim3_HR_array <- second_FP_sim3_HR$sims.array[,1,]
second_FP_sim3_HR_array <- apply(second_FP_sim3_HR_array, 2, sample, size = 500)

#1.4 Second-order FP, FE model, P1 = -1, P2 = 0 ----------------
second_FP$P1 <- -1
second_FP$P2 <- 0
second_FP_sim4_HR <- bugs(data = second_FP, inits = second_inits, model.file = second_FP_model_FE, 
                         n.iter = 40000, parameters.to.save = HR_params, n.burnin = 20000, 
                         n.chains = 2, 
                         debug = FALSE)
second_FP_sim4_HR_array <- second_FP_sim4_HR$sims.array[,1,]
second_FP_sim4_HR_array <- apply(second_FP_sim4_HR_array, 2, sample, size = 500)

#1.5 Second-order FP, FE model, P1 = -1, P2 = -0.5 ----------------
second_FP$P1 <- -1
second_FP$P2 <- -0.5
second_FP_sim5_HR <- bugs(data = second_FP, inits = second_inits, model.file = second_FP_model_FE, 
                         n.iter = 40000, parameters.to.save = HR_params, n.burnin = 20000, 
                         n.chains = 2, 
                         debug = FALSE)
second_FP_sim5_HR_array <- second_FP_sim5_HR$sims.array[,1,]
second_FP_sim5_HR_array <- apply(second_FP_sim5_HR_array, 2, sample, size = 500)

# Functions to call several times (it calculates averaged HRs)
model_average <- function(bugs_arrays, time_point, treatment){
  summary.stats(1/((bugs_arrays[[1]][, time_point][,treatment] + 
                    bugs_arrays[[2]][, time_point][,treatment] +
                    bugs_arrays[[3]][, time_point][,treatment] + 
                    bugs_arrays[[4]][, time_point][,treatment] +
                    bugs_arrays[[5]][, time_point][,treatment])/length(bugs_arrays)))
}

# list of MCMC arrays 
arrays <- list(second_FP_sim1_HR_array, second_FP_sim2_HR_array, second_FP_sim3_HR_array,
               second_FP_sim4_HR_array, second_FP_sim5_HR_array)

# HR at 3 years
three_years <- c("HR_3year[2]", "HR_3year[3]", "HR_3year[4]", 
               "HR_3year[5]", "HR_3year[6]", "HR_3year[7]", "HR_3year[8]")

doc_HR_36 <-     model_average(bugs_arrays = arrays, time_point = three_years, treatment = 1)
enza_HR_36 <-    model_average(bugs_arrays = arrays, time_point = three_years, treatment = 2)
adt_HR_36 <-     model_average(bugs_arrays = arrays, time_point = three_years, treatment = 3)
abi_HR_36 <-     model_average(bugs_arrays = arrays, time_point = three_years, treatment = 4)
apa_HR_36 <-     model_average(bugs_arrays = arrays, time_point = three_years, treatment = 5)
enzadoc_HR_36 <- model_average(bugs_arrays = arrays, time_point = three_years, treatment = 6)
abidoc_HR_36 <-  model_average(bugs_arrays = arrays, time_point = three_years, treatment = 7)

HR_36 <- c(doc_HR_36,enza_HR_36,adt_HR_36,abi_HR_36,apa_HR_36,enzadoc_HR_36,abidoc_HR_36)

# HR at 4 years
four_years <- c("HR_4year[2]", "HR_4year[3]", "HR_4year[4]", 
                 "HR_4year[5]", "HR_4year[6]", "HR_4year[7]", "HR_4year[8]")

doc_HR_48 <-     model_average(bugs_arrays = arrays, time_point = four_years, treatment = 1)
enza_HR_48 <-    model_average(bugs_arrays = arrays, time_point = four_years, treatment = 2)
adt_HR_48 <-     model_average(bugs_arrays = arrays, time_point = four_years, treatment = 3)
abi_HR_48 <-     model_average(bugs_arrays = arrays, time_point = four_years, treatment = 4)
apa_HR_48 <-     model_average(bugs_arrays = arrays, time_point = four_years, treatment = 5)
enzadoc_HR_48 <- model_average(bugs_arrays = arrays, time_point = four_years, treatment = 6)
abidoc_HR_48 <-  model_average(bugs_arrays = arrays, time_point = four_years, treatment = 7)

HR_48 <- c(doc_HR_48,enza_HR_48,adt_HR_48,abi_HR_48,apa_HR_48,enzadoc_HR_48,abidoc_HR_48)

# HR at 5 years
five_years <- c("HR_5year[2]", "HR_5year[3]", "HR_5year[4]", 
                "HR_5year[5]", "HR_5year[6]", "HR_5year[7]", "HR_5year[8]")

doc_HR_60     <- model_average(bugs_arrays = arrays, time_point = five_years, treatment = 1)
enza_HR_60    <- model_average(bugs_arrays = arrays, time_point = five_years, treatment = 2)
adt_HR_60     <- model_average(bugs_arrays = arrays, time_point = five_years, treatment = 3)
abi_HR_60     <- model_average(bugs_arrays = arrays, time_point = five_years, treatment = 4)
apa_HR_60     <- model_average(bugs_arrays = arrays, time_point = five_years, treatment = 5)
enzadoc_HR_60 <- model_average(bugs_arrays = arrays, time_point = five_years, treatment = 6)
abidoc_HR_60  <- model_average(bugs_arrays = arrays, time_point = five_years, treatment = 7)

HR_60 <- c(doc_HR_60, enza_HR_60, adt_HR_60, abi_HR_60, apa_HR_60, enzadoc_HR_60, abidoc_HR_60)

# HR to be exported
HRs <- cbind(HR_36, HR_48, HR_60)
write_csv(as_tibble(HRs), here("05_tables", "model_average_HRs_OS_v2.csv"))

#2. RMST at 60 months----------------------------------
# Base case with model with the lowest DIC
# Absolute RSMT
rmst_daro    <- summary.stats(second_FP_sim1_HR$sims.array[,1,"RMST[1,60]"])
rmst_doc     <- summary.stats(second_FP_sim1_HR$sims.array[,1,"RMST[2,60]"])
rmst_enza    <- summary.stats(second_FP_sim1_HR$sims.array[,1,"RMST[3,60]"])
rmst_adt     <- summary.stats(second_FP_sim1_HR$sims.array[,1,"RMST[4,60]"])
rmst_abi     <- summary.stats(second_FP_sim1_HR$sims.array[,1,"RMST[5,60]"])
rmst_apa     <- summary.stats(second_FP_sim1_HR$sims.array[,1,"RMST[6,60]"])
rmst_enzadoc <- summary.stats(second_FP_sim1_HR$sims.array[,1,"RMST[7,60]"])
rmst_abidoc  <- summary.stats(second_FP_sim1_HR$sims.array[,1,"RMST[8,60]"])

rmst <- c(rmst_daro, rmst_doc, rmst_enza, rmst_adt, rmst_abi, rmst_apa,
          rmst_enzadoc, rmst_abidoc)

# Relative RSMT to reference treatment
diff_doc <- summary.stats(second_FP_sim1_HR$sims.array[,1,"RMST[2,60]"] - 
                          second_FP_sim1_HR$sims.array[,1,"RMST[1,60]"])
diff_enza <- summary.stats(second_FP_sim1_HR$sims.array[,1,"RMST[3,60]"] - 
                           second_FP_sim1_HR$sims.array[,1,"RMST[1,60]"])
diff_adt <- summary.stats(second_FP_sim1_HR$sims.array[,1,"RMST[4,60]"] - 
                          second_FP_sim1_HR$sims.array[,1,"RMST[1,60]"])
diff_abi <- summary.stats(second_FP_sim1_HR$sims.array[,1,"RMST[5,60]"] - 
                          second_FP_sim1_HR$sims.array[,1,"RMST[1,60]"])
diff_apa <- summary.stats(second_FP_sim1_HR$sims.array[,1,"RMST[6,60]"] - 
                          second_FP_sim1_HR$sims.array[,1,"RMST[1,60]"])
diff_enzadoc <- summary.stats(second_FP_sim1_HR$sims.array[,1,"RMST[7,60]"] - 
                              second_FP_sim1_HR$sims.array[,1,"RMST[1,60]"])
diff_abidoc <- summary.stats(second_FP_sim1_HR$sims.array[,1,"RMST[8,60]"] - 
                             second_FP_sim1_HR$sims.array[,1,"RMST[1,60]"])

diff <- c(NA, diff_doc, diff_enza, diff_adt, diff_abi, diff_apa,
          diff_enzadoc, diff_abidoc)

# Export results
write.csv(bind_cols(rmst, diff), here("05_tables", "rmst_os_v3.csv"))

# Functions to call several times (slightly different from HRs)
model_average_rmst <- function(bugs_arrays, params, treatment){
  summary.stats((bugs_arrays[[1]][, params][,treatment] + 
                 bugs_arrays[[2]][, params][,treatment] +
                 bugs_arrays[[3]][, params][,treatment] + 
                 bugs_arrays[[4]][, params][,treatment] +
                 bugs_arrays[[5]][, params][,treatment])/length(bugs_arrays))
}

# parameters of interest
params <- c("RMST[1,60]","RMST[2,60]","RMST[3,60]","RMST[4,60]",
            "RMST[5,60]","RMST[6,60]","RMST[7,60]","RMST[8,60]")


# RMST per treatment
rmst_daro    <- model_average_rmst(arrays, params = params, treatment = 1)
rmst_doc     <- model_average_rmst(arrays, params = params, treatment = 2)
rmst_enza    <- model_average_rmst(arrays, params = params, treatment = 3)
rmst_adt     <- model_average_rmst(arrays, params = params, treatment = 4)
rmst_abi     <- model_average_rmst(arrays, params = params, treatment = 5)
rmst_apa     <- model_average_rmst(arrays, params = params, treatment = 6)
rmst_enzadoc <- model_average_rmst(arrays, params = params, treatment = 7)
rmst_abidoc  <- model_average_rmst(arrays, params = params, treatment = 8)

rmst <- c(rmst_daro, rmst_doc, rmst_enza, rmst_adt, rmst_abi, rmst_apa,
          rmst_enzadoc, rmst_abidoc)

# Relative RSMT to reference treatment
# Function to estimate the relative RMST (slightly different from HRs)
model_average_rmst_rel <- function(bugs_arrays, params, treatment, reference){
  summary.stats((bugs_arrays[[1]][, params][,treatment] + 
                   bugs_arrays[[2]][, params][,treatment] +
                   bugs_arrays[[3]][, params][,treatment] + 
                   bugs_arrays[[4]][, params][,treatment] +
                   bugs_arrays[[5]][, params][,treatment])/length(bugs_arrays) -
                  (bugs_arrays[[1]][, params][,reference] + 
                   bugs_arrays[[2]][, params][,reference] +
                   bugs_arrays[[3]][, params][,reference] + 
                   bugs_arrays[[4]][, params][,reference] +
                   bugs_arrays[[5]][, params][,reference])/length(bugs_arrays))
}

diff_doc <-     model_average_rmst_rel(arrays, params = params, treatment = 2, reference = 1)
diff_enza <-    model_average_rmst_rel(arrays, params = params, treatment = 3, reference = 1)
diff_adt <-     model_average_rmst_rel(arrays, params = params, treatment = 4, reference = 1)
diff_abi <-     model_average_rmst_rel(arrays, params = params, treatment = 5, reference = 1)
diff_apa <-     model_average_rmst_rel(arrays, params = params, treatment = 6, reference = 1)
diff_enzadoc <- model_average_rmst_rel(arrays, params = params, treatment = 7, reference = 1)
diff_abidoc <-  model_average_rmst_rel(arrays, params = params, treatment = 8, reference = 1)

diff <- c(NA, diff_doc, diff_enza, diff_adt, diff_abi, diff_apa,
          diff_enzadoc, diff_abidoc)

write.csv(bind_cols(rmst, diff), here("05_tables", "model_average_rmst_os_v2.csv"))


#3. Survival (%OS) at 60 months-----------------------------------
# Import arrays from the "S" parameter
fp_data_1_array <- read_csv(here("02_data", "OS_array_FP_sim14.csv"), col_select = -1)  # model with P1: -0.5, P2: -0.5
fp_data_2_array <- read_csv(here("02_data", "OS_array_FP_sim15.csv"), col_select = -1)  # model with P1: -0.5, P2:  0
fp_data_3_array <- read_csv(here("02_data", "OS_array_FP_sim11.csv"), col_select = -1)  # model with P1: -1,   P2:  0.5
fp_data_4_array <- read_csv(here("02_data", "OS_array_FP_sim10.csv"), col_select = -1)  # model with P1: -1,   P2:  0
fp_data_5_array <- read_csv(here("02_data", "OS_array_FP_sim09.csv"), col_select = -1)  # model with P1: -1,   P2: -0.5

# Base case model
fp_data_1 <- read.csv(here("02_data", "OS_results_FP_sim14.csv"))  # model with P1: -0.5, P2: -0.5
surv_daro    <- filter(fp_data_1, X == "S[1,60]")
surv_doc     <- filter(fp_data_1, X == "S[2,60]")
surv_enza    <- filter(fp_data_1, X == "S[3,60]")
surv_adt     <- filter(fp_data_1, X == "S[4,60]")
surv_abi     <- filter(fp_data_1, X == "S[5,60]")
surv_apa     <- filter(fp_data_1, X == "S[6,60]")
surv_enzadoc <- filter(fp_data_1, X == "S[7,60]")
surv_abidoc  <- filter(fp_data_1, X == "S[8,60]")

surv_base <- bind_rows(surv_daro, surv_doc, surv_enza, surv_adt, surv_abi, surv_apa, surv_enzadoc, surv_abidoc)
surv_base <- select(surv_base, c(mean, X2.5., X97.5.)) %>% round(digits = 2)
surv_base <- mutate(surv_base, string = paste0(mean, " (", X2.5. ,", ", X97.5., ")"))
write_csv(surv_base, here("05_tables", "basecase_surv_60months.csv"))

# function to calculate average survival at 60 months
model_average_surv <- function(bugs_arrays, treatment){
  average <- (bugs_arrays[[1]][, treatment] + 
                    bugs_arrays[[2]][, treatment] +
                    bugs_arrays[[3]][, treatment] + 
                    bugs_arrays[[4]][, treatment] +
                    bugs_arrays[[5]][, treatment])/length(bugs_arrays)
  summary.stats(average[,1])
}

# list with MCMC arrays
arrays_surv <- list(fp_data_1_array,fp_data_2_array,fp_data_3_array,fp_data_4_array,fp_data_5_array)

# average survival at 108 months per treatment
surv_daro    <- model_average_surv(bugs_arrays = arrays_surv, treatment = 1)
surv_doc     <- model_average_surv(bugs_arrays = arrays_surv, treatment = 2)
surv_enza    <- model_average_surv(bugs_arrays = arrays_surv, treatment = 3)
surv_adt     <- model_average_surv(bugs_arrays = arrays_surv, treatment = 4)
surv_abi     <- model_average_surv(bugs_arrays = arrays_surv, treatment = 5)
surv_apa     <- model_average_surv(bugs_arrays = arrays_surv, treatment = 6)
surv_enzadoc <- model_average_surv(bugs_arrays = arrays_surv, treatment = 7)
surv_abidoc  <- model_average_surv(bugs_arrays = arrays_surv, treatment = 8)

# export results
surv_average <- c(surv_daro, surv_doc, surv_enza, surv_adt, surv_abi, surv_apa, surv_enzadoc, surv_abidoc)
write_csv(as_tibble(surv_average), here("05_tables", "model_average_surv_60months.csv"))
