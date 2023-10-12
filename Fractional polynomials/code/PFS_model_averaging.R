# Code for model averaging
# Outcome: PFS

# 1. Packages used --------------------------------------------------------
pkgs <- c("tidyverse", "readxl", "here", "haven", "R2OpenBUGS")
lapply(pkgs, library, character.only = T)
set.seed(03082022)
source(here("Fractional polynomials/code", "utils.R")) 

# 2. Import data --------------------------------------------------------
# These data is only to match treatment codes, not actual HRs
PFS_data <- read_excel(here("Fractional polynomials/data", "mHSPC PFS trt codes.xlsx"))

# Fractional polynomial data created in PFS_data_preparation.R
FP_data <- read_csv(here("Fractional polynomials/data", "FP data", "FP_data_PFS.csv"))
which(FP_data$z - FP_data$r < 0) # check if any r is higher than natrisk at each interval

# Transform data into a list for 1st order FP
first_FP <- list(
  P1   = 0,                          # P1 for 1st order FP
  N    = nrow(FP_data),             # intervals derived from the reconstructed IPD
  nt   = max(PFS_data$t1, PFS_data$t2),                         # number of treatments
  ns   = nrow(PFS_data),                         # number of studies
  maxt = max(FP_data$time_j),       # maximum follow-up (months)
  mean = c(0, 0),                    # priors for the mean effect
  prec = structure(.Data= c(1.00000E-04, 0, 0, 1.00000E-04), 
                   .Dim=c(2, 2)),   # priors for the precision
  t = matrix(c(PFS_data$t1, PFS_data$t2), # treatments IDs
             nrow = nrow(PFS_data), ncol = 2),
  na = PFS_data$na,               # number of arms
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
second_FP$maxt <- 36 # 36 months requested by Bayer
second_FP$mean <- c(0, 0, 0)
second_FP$prec <- structure(.Data= c(1.00000E-04, 0, 0,
                                     0, 1.00000E-02, 0,
                                     0, 0, 1.00000E-02), .Dim=c(3, 3))

# call the model
second_FP_model_FE <- here("Fractional polynomials/code", "second_FP_model_FE.txt")

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
HR_params <- c("HR_1year", "HR_2year", "HR_3year", "RMST")

#1.1 Second-order FP, FE model, P1 = -0.5, P2 = -0.5 ----------------
second_FP$P1 <- -0.5
second_FP$P2 <- -0.5
second_FP_sim1_HR <- bugs(data = second_FP, inits = second_inits, model.file = second_FP_model_FE, 
                          n.iter = 40000, parameters.to.save = HR_params, n.burnin = 20000, 
                          n.chains = 2, 
                          debug = FALSE)

# Extract HR at timepoints of interest for the base case model
HR_export_12 <- as_tibble(round(1/second_FP_sim1_HR$summary[grep("HR_1year", rownames(second_FP_sim1_HR$summary)), c(1, 7, 3)], digits = 2))
HR_export_12 <- mutate(HR_export_12, string = paste0(mean, " (", `97.5%` ,", ", `2.5%`, ")"))
HR_export_24 <- as_tibble(round(1/second_FP_sim1_HR$summary[grep("HR_2year", rownames(second_FP_sim1_HR$summary)), c(1, 7, 3)], digits = 2))
HR_export_24 <- mutate(HR_export_24, string = paste0(mean, " (", `97.5%` ,", ", `2.5%`, ")"))
HR_export_36 <- as_tibble(round(1/second_FP_sim1_HR$summary[grep("HR_3year", rownames(second_FP_sim1_HR$summary)), c(1, 7, 3)], digits = 2))
HR_export_36 <- mutate(HR_export_36, string = paste0(mean, " (", `97.5%` ,", ", `2.5%`, ")"))

HR_export <- bind_cols(HR_export_12, HR_export_24, HR_export_36)
write_csv(HR_export, here("Fractional polynomials/data", "PFS_HR_basecase.csv"))

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

# HR at 1 year
one_year <- c("HR_1year[2]", "HR_1year[3]", "HR_1year[4]", 
              "HR_1year[5]", "HR_1year[6]", "HR_1year[7]", "HR_1year[8]")

doc_HR_12 <-     model_average(bugs_arrays = arrays, time_point = one_year, treatment = 1)
enza_HR_12 <-    model_average(bugs_arrays = arrays, time_point = one_year, treatment = 2)
adt_HR_12 <-     model_average(bugs_arrays = arrays, time_point = one_year, treatment = 3)
abi_HR_12 <-     model_average(bugs_arrays = arrays, time_point = one_year, treatment = 4)
apa_HR_12 <-     model_average(bugs_arrays = arrays, time_point = one_year, treatment = 5)
enzadoc_HR_12 <- model_average(bugs_arrays = arrays, time_point = one_year, treatment = 6)
abidoc_HR_12 <-  model_average(bugs_arrays = arrays, time_point = one_year, treatment = 7)

HR_12 <- c(doc_HR_12,enza_HR_12,adt_HR_12,abi_HR_12,apa_HR_12,enzadoc_HR_12,abidoc_HR_12)

# HR at 2 years
two_years <- c("HR_2year[2]", "HR_2year[3]", "HR_2year[4]", 
               "HR_2year[5]", "HR_2year[6]", "HR_2year[7]", "HR_2year[8]")

doc_HR_24 <- model_average(bugs_arrays = arrays, time_point = two_years, treatment = 1)
enza_HR_24 <- model_average(bugs_arrays = arrays, time_point = two_years, treatment = 2)
adt_HR_24 <- model_average(bugs_arrays = arrays, time_point = two_years, treatment = 3)
abi_HR_24 <- model_average(bugs_arrays = arrays, time_point = two_years, treatment = 4)
apa_HR_24 <- model_average(bugs_arrays = arrays, time_point = two_years, treatment = 5)
enzadoc_HR_24 <- model_average(bugs_arrays = arrays, time_point = two_years, treatment = 6)
abidoc_HR_24 <- model_average(bugs_arrays = arrays, time_point = two_years, treatment = 7)

HR_24 <- c(doc_HR_24, enza_HR_24, adt_HR_24, abi_HR_24, apa_HR_24, enzadoc_HR_24, abidoc_HR_24)

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

# HR to be exported
HRs <- cbind(HR_12, HR_24, HR_36)
write_csv(as_tibble(HRs), here("Fractional polynomials/data", "model_average_HRs_PFS.csv"))

#2. RMST at 36 months----------------------------------
#Base case with model with the lowest DIC
# Absolute RSMT
rmst_daro    <- summary.stats(second_FP_sim1_HR$sims.array[,1,"RMST[1,36]"])
rmst_doc     <- summary.stats(second_FP_sim1_HR$sims.array[,1,"RMST[2,36]"])
rmst_enza    <- summary.stats(second_FP_sim1_HR$sims.array[,1,"RMST[3,36]"])
rmst_adt     <- summary.stats(second_FP_sim1_HR$sims.array[,1,"RMST[4,36]"])
rmst_abi     <- summary.stats(second_FP_sim1_HR$sims.array[,1,"RMST[5,36]"])
rmst_apa     <- summary.stats(second_FP_sim1_HR$sims.array[,1,"RMST[6,36]"])
rmst_enzadoc <- summary.stats(second_FP_sim1_HR$sims.array[,1,"RMST[7,36]"])
rmst_abidoc  <- summary.stats(second_FP_sim1_HR$sims.array[,1,"RMST[8,36]"])

rmst <- c(rmst_daro, rmst_doc, rmst_enza, rmst_adt, rmst_abi, rmst_apa,
          rmst_enzadoc, rmst_abidoc)

# Relative RSMT to reference treatment
diff_doc <- summary.stats(second_FP_sim1_HR$sims.array[,1,"RMST[2,36]"] - 
                          second_FP_sim1_HR$sims.array[,1,"RMST[1,36]"])
diff_enza <- summary.stats(second_FP_sim1_HR$sims.array[,1,"RMST[3,36]"] - 
                           second_FP_sim1_HR$sims.array[,1,"RMST[1,36]"])
diff_adt <- summary.stats(second_FP_sim1_HR$sims.array[,1,"RMST[4,36]"] - 
                          second_FP_sim1_HR$sims.array[,1,"RMST[1,36]"])
diff_abi <- summary.stats(second_FP_sim1_HR$sims.array[,1,"RMST[5,36]"] - 
                          second_FP_sim1_HR$sims.array[,1,"RMST[1,36]"])
diff_apa <- summary.stats(second_FP_sim1_HR$sims.array[,1,"RMST[6,36]"] - 
                          second_FP_sim1_HR$sims.array[,1,"RMST[1,36]"])
diff_enzadoc <- summary.stats(second_FP_sim1_HR$sims.array[,1,"RMST[7,36]"] - 
                              second_FP_sim1_HR$sims.array[,1,"RMST[1,36]"])
diff_abidoc <- summary.stats(second_FP_sim1_HR$sims.array[,1,"RMST[8,36]"] - 
                             second_FP_sim1_HR$sims.array[,1,"RMST[1,36]"])

diff <- c(NA, diff_doc, diff_enza, diff_adt, diff_abi, diff_apa,
          diff_enzadoc, diff_abidoc)

# Export results
write.csv(bind_cols(rmst, diff), here("Fractional polynomials/data", "rmst_PFS.csv"))

# Functions to call several times (slightly different from HRs)
model_average_rmst <- function(bugs_arrays, params, treatment){
  summary.stats((bugs_arrays[[1]][, params][,treatment] + 
                   bugs_arrays[[2]][, params][,treatment] +
                   bugs_arrays[[3]][, params][,treatment] + 
                   bugs_arrays[[4]][, params][,treatment] +
                   bugs_arrays[[5]][, params][,treatment])/length(bugs_arrays))
}

# parameters of interest
params <- c("RMST[1,36]","RMST[2,36]","RMST[3,36]","RMST[4,36]",
            "RMST[5,36]","RMST[6,36]","RMST[7,36]","RMST[8,36]")


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

write.csv(bind_cols(rmst, diff), here("Fractional polynomials/data", "model_average_rmst_PFS.csv"))

#3. Survival (%PFS) at 36 months-----------------------------------
# Import arrays from the "S" parameter
fp_data_1_array <- read_csv(here("Fractional polynomials/data", "PFS_array_FP_sim14.csv"), col_select = -1)  # model with P1: -0.5, P2: -0.5
fp_data_2_array <- read_csv(here("Fractional polynomials/data", "PFS_array_FP_sim15.csv"), col_select = -1)  # model with P1: -0.5, P2:  0
fp_data_3_array <- read_csv(here("Fractional polynomials/data", "PFS_array_FP_sim11.csv"), col_select = -1)  # model with P1: -1,   P2:  0.5
fp_data_4_array <- read_csv(here("Fractional polynomials/data", "PFS_array_FP_sim10.csv"), col_select = -1)  # model with P1: -1,   P2:  0
fp_data_5_array <- read_csv(here("Fractional polynomials/data", "PFS_array_FP_sim19.csv"), col_select = -1)  # model with P1: -1,   P2: -0.5

# Base case model
fp_data_1 <- read.csv(here("Fractional polynomials/data", "PFS_results_FP_sim14.csv"))  # model with P1: -0.5, P2: -0.5
surv_daro    <- filter(fp_data_1, X == "S[1,36]")
surv_doc     <- filter(fp_data_1, X == "S[2,36]")
surv_enza    <- filter(fp_data_1, X == "S[3,36]")
surv_adt     <- filter(fp_data_1, X == "S[4,36]")
surv_abi     <- filter(fp_data_1, X == "S[5,36]")
surv_apa     <- filter(fp_data_1, X == "S[6,36]")
surv_enzadoc <- filter(fp_data_1, X == "S[7,36]")
surv_abidoc  <- filter(fp_data_1, X == "S[8,36]")

surv_base <- bind_rows(surv_daro, surv_doc, surv_enza, surv_adt, surv_abi, surv_apa, surv_enzadoc, surv_abidoc)
surv_base <- select(surv_base, c(mean, X2.5., X97.5.)) %>% round(digits = 2)
surv_base <- mutate(surv_base, string = paste0(mean, " (", X2.5. ,", ", X97.5., ")"))
write_csv(surv_base, here("Fractional polynomials/data", "basecase_surv_36months.csv"))

# function to calculate average survival at 108 months
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
write_csv(as_tibble(surv_average), here("Fractional polynomials/data", "model_average_surv_36months.csv"))
