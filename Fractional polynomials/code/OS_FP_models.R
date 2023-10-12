# Darulatamide in mHSPC NMA #
# Outcome: OS (KM curves data)

# 1. Packages used --------------------------------------------------------
pkgs <- c("tidyverse", "readxl", "here", "haven", "multinma", "R2OpenBUGS",
          "igraph")
lapply(pkgs, library, character.only = T)
set.seed(03082022)
#source(here ("01_scripts", "cross_or.R")) 
source(here("Fractional polynomials/code", "utils.R")) 

# 2. Import data --------------------------------------------------------
# These data is only to match treatment codes, not actual HRs
OS_data <- read_excel(here("Fractional polynomials/data", "mHSPC OS trt codes.xlsx"))

# Set up the evidence network
OS_long <- pivot_longer(OS_data, cols = c(t1_name, t2_name), names_to = "trt", values_to = "treatment")
OS_long <- pivot_longer(OS_data, cols = c(t1, t2), names_to = "trt", values_to = "treatment")
OS_long$sample_size <- c(651, 654, 574, 576, 397, 393, 563, 562, 563, 562, 192, 193, 597, 602, 355, 355, 960, 957, 362, 724, 377, 189, 525, 527, 36, 35)
OS_net <- set_agd_arm(OS_long, study = study, trt = treatment, y = y, se = se, sample_size = sample_size)

#jpeg(here("04_figures", "network-plot_base.jpg"), width = 800, height = 700, res = 120)
plot(OS_net, weight_nodes = TRUE, nudge = 0.1) +
  #ggraph::geom_edge_fan(aes(label = unique(OS_net$agd_arm$.study)), 
  #               angle_calc = 'along',
  #               label_dodge = unit(2.5, 'mm')) +
  ggplot2::theme(legend.position = "bottom",
                 legend.box = "vertical")
#dev.off()

# Import data from KM curves (prepared in OS_data-preparation-FP.R)
FP_data <- read_csv(here("Fractional polynomials/data", "FP_data_OS.csv"))
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

# 3.2 Fit Fixed Effects 1st order fractional polynomial model----------------
# Call the model
first_FP_model_FE <- here("Fractional polynomials/code", "first_FP_model_FE.txt")
  
#Initial Values 
inits <- function(){
  #chain 1
  list(d= structure(.Data= rep(c(NA, rep(-0.1, first_FP$nt - 1)), 2),
                    .Dim=c(first_FP$nt, 2)), 
       mu= structure(.Data= c(rep(0.1, 2 * first_FP$ns)), 
                     .Dim=c(first_FP$ns, 2)))
  #chain 2
  list(d= structure(.Data= rep(c(NA, rep(1, first_FP$nt - 1)), 2),
                    .Dim=c(first_FP$nt, 2)),  
       mu = structure(.Data= c(rep(1, 2 * first_FP$ns)),
                      .Dim=c(first_FP$ns, 2)))
}

# Parameters to save
params_FE <- c("d", "totresdev")

# First order FP P1 = -2---------------------------------------------
first_FP$P1 <- -2
first_FP_sim1 <- bugs(data = first_FP, inits = inits,
                             model.file = first_FP_model_FE, n.iter = 20000,
                             parameters.to.save = params_FE, n.burnin = 10000, n.chains = 2,
                             debug = FALSE)

first_FP_sim1

# First order FP P1 = -1---------------------------------------------
first_FP$P1 <- -1
first_FP_sim2 <- bugs(data = first_FP, inits = inits,
                      model.file = first_FP_model_FE, n.iter = 20000,
                      parameters.to.save = params_FE, n.burnin = 10000, n.chains = 2,
                      debug = FALSE)

first_FP_sim2

# First order FP P1 = -0.5---------------------------------------------
first_FP$P1 <- -0.5
first_FP_sim3 <- bugs(data = first_FP, inits = inits,
                      model.file = first_FP_model_FE, n.iter = 20000,
                      parameters.to.save = params_FE, n.burnin = 10000, n.chains = 2,
                      debug = FALSE)

first_FP_sim3

# First order FP P1 = 0---------------------------------------------
first_FP$P1 <- 0
first_FP_sim4 <- bugs(data = first_FP, inits = inits,
                      model.file = first_FP_model_FE, n.iter = 20000,
                      parameters.to.save = params_FE, n.burnin = 10000, n.chains = 2,
                      debug = FALSE)

first_FP_sim4

# First order FP P1 = 0.5---------------------------------------------
first_FP$P1 <- 0.5
first_FP_sim5 <- bugs(data = first_FP, inits = inits,
                      model.file = first_FP_model_FE, n.iter = 20000,
                      parameters.to.save = params_FE, n.burnin = 10000, n.chains = 2,
                      debug = FALSE)

first_FP_sim5

# First order FP P1 = 1---------------------------------------------
first_FP$P1 <- 1
first_FP_sim6 <- bugs(data = first_FP, inits = inits,
                      model.file = first_FP_model_FE, n.iter = 20000,
                      parameters.to.save = params_FE, n.burnin = 10000, n.chains = 2,
                      debug = FALSE)

first_FP_sim6

# First order FP P1 = 1---------------------------------------------
first_FP$P1 <- 2
first_FP_sim7 <- bugs(data = first_FP, inits = inits,
                      model.file = first_FP_model_FE, n.iter = 20000,
                      parameters.to.save = params_FE, n.burnin = 10000, n.chains = 2,
                      debug = FALSE)

first_FP_sim7

# 4. Random effects models  -------------------------------
# 4.1 Set up the evidence network

# 4.2 Fit random effects NMA model

# 4.3 Assess convergence and model fit

# 5. Fixed Effects 2nd order fractional polynomial model----------------
# call the model
second_FP_model_FE <- here("Fractional polynomials/code", "second_FP_model_FE.txt")

# data for 2nd order FP model
second_FP <- first_FP
second_FP$P2 <- 0
second_FP$mean <- c(0, 0, 0)
second_FP$prec <- structure(.Data= c(1.00000E-04, 0, 0,
                                     0, 1.00000E-02, 0,
                                     0, 0, 1.00000E-02), .Dim=c(3, 3))

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

# Monitor parameters
second_params_FE <- c("d", "totresdev")

# Second-order FP, FE model, P1 = -2, P2 = -2 ----------------
second_FP$P1 <- -2
second_FP$P2 <- -2
second_FP_sim_1 <- bugs(data = second_FP, inits = second_inits, model.file = second_FP_model_FE, 
                        n.iter = 20000, parameters.to.save = second_params_FE, n.burnin = 10000, 
                        n.chains = 2, 
                        debug = FALSE)

second_FP_sim_1

# Second-order FP, FE model, P1 = -2, P2 = -1 ----------------
second_FP$P1 <- -2
second_FP$P2 <- -1
second_FP_sim_2 <- bugs(data = second_FP, inits = second_inits, model.file = second_FP_model_FE, 
                        n.iter = 20000, parameters.to.save = second_params_FE, n.burnin = 10000, 
                        n.chains = 2, 
                        debug = FALSE)

second_FP_sim_2

# Second-order FP, FE model, P1 = -2, P2 = -0.5 ----------------
second_FP$P1 <- -2
second_FP$P2 <- -0.5
second_FP_sim_3 <- bugs(data = second_FP, inits = second_inits, model.file = second_FP_model_FE, 
                        n.iter = 20000, parameters.to.save = second_params_FE, n.burnin = 10000, 
                        n.chains = 2, 
                        debug = FALSE)

second_FP_sim_3

# Second-order FP, FE model, P1 = -2, P2 = 0 ----------------
second_FP$P1 <- -2
second_FP$P2 <- 0
second_FP_sim_4 <- bugs(data = second_FP, inits = second_inits, model.file = second_FP_model_FE, 
                        n.iter = 20000, parameters.to.save = second_params_FE, n.burnin = 10000, 
                        n.chains = 2, 
                        debug = FALSE)

second_FP_sim_4

# Second-order FP, FE model, P1 = -2, P2 = 0.5 ----------------
second_FP$P1 <- -2
second_FP$P2 <- 0.5
second_FP_sim_5 <- bugs(data = second_FP, inits = second_inits, model.file = second_FP_model_FE, 
                        n.iter = 20000, parameters.to.save = second_params_FE, n.burnin = 10000, 
                        n.chains = 2, 
                        debug = FALSE)

second_FP_sim_5

# Second-order FP, FE model, P1 = -2, P2 = 1 ----------------
second_FP$P1 <- -2
second_FP$P2 <- 1
second_FP_sim_6 <- bugs(data = second_FP, inits = second_inits, model.file = second_FP_model_FE, 
                        n.iter = 20000, parameters.to.save = second_params_FE, n.burnin = 10000, 
                        n.chains = 2, 
                        debug = FALSE)

second_FP_sim_6


# Second-order FP, FE model, P1 = -2, P2 = 2 ----------------
second_FP$P1 <- -2
second_FP$P2 <- 2
second_FP_sim_7 <- bugs(data = second_FP, inits = second_inits, model.file = second_FP_model_FE, 
                        n.iter = 20000, parameters.to.save = second_params_FE, n.burnin = 10000, 
                        n.chains = 2, 
                        debug = FALSE)

second_FP_sim_7

# Second-order FP, FE model, P1 = -1, P2 = -1 ----------------
second_FP$P1 <- -1
second_FP$P2 <- -1
second_FP_sim_8 <- bugs(data = second_FP, inits = second_inits, model.file = second_FP_model_FE, 
                        n.iter = 20000, parameters.to.save = second_params_FE, n.burnin = 10000, 
                        n.chains = 2, 
                        debug = FALSE)

second_FP_sim_8

# Second-order FP, FE model, P1 = -1, P2 = -0.5 ----------------
second_FP$P1 <- -1
second_FP$P2 <- -0.5
second_FP_sim_9 <- bugs(data = second_FP, inits = second_inits, model.file = second_FP_model_FE, 
                        n.iter = 20000, parameters.to.save = second_params_FE, n.burnin = 10000, 
                        n.chains = 2, 
                        debug = FALSE)

second_FP_sim_9

# Second-order FP, FE model, P1 = -1, P2 = 0 ----------------
second_FP$P1 <- -1
second_FP$P2 <- 0
second_FP_sim_10 <- bugs(data = second_FP, inits = second_inits, model.file = second_FP_model_FE, 
                         n.iter = 20000, parameters.to.save = second_params_FE, n.burnin = 10000, 
                         n.chains = 2, 
                         debug = FALSE)

second_FP_sim_10

# Second-order FP, FE model, P1 = -1, P2 = 0.5 ----------------
second_FP$P1 <- -1
second_FP$P2 <- 0.5
second_FP_sim_11 <- bugs(data = second_FP, inits = second_inits, model.file = second_FP_model_FE, 
                         n.iter = 20000, parameters.to.save = second_params_FE, n.burnin = 10000, 
                         n.chains = 2, 
                         debug = FALSE)

second_FP_sim_11

# Second-order FP, FE model, P1 = -1, P2 = 1 ----------------
second_FP$P1 <- -1
second_FP$P2 <- 1
second_FP_sim_12 <- bugs(data = second_FP, inits = second_inits, model.file = second_FP_model_FE, 
                         n.iter = 20000, parameters.to.save = second_params_FE, n.burnin = 10000, 
                         n.chains = 2, 
                         debug = FALSE)

second_FP_sim_12

# Second-order FP, FE model, P1 = -1, P2 = 2 ----------------
second_FP$P1 <- -1
second_FP$P2 <- 2
second_FP_sim_13 <- bugs(data = second_FP, inits = second_inits, model.file = second_FP_model_FE, 
                         n.iter = 20000, parameters.to.save = second_params_FE, n.burnin = 10000, 
                         n.chains = 2, 
                         debug = FALSE)

second_FP_sim_13

# Second-order FP, FE model, P1 = -0.5, P2 = -0.5 ----------------
second_FP$P1 <- -0.5
second_FP$P2 <- -0.5
second_FP_sim_14 <- bugs(data = second_FP, inits = second_inits, model.file = second_FP_model_FE, 
                         n.iter = 20000, parameters.to.save = second_params_FE, n.burnin = 10000, 
                         n.chains = 2, 
                         debug = FALSE)

second_FP_sim_14

# Second-order FP, FE model, P1 = -0.5, P2 = 0 ----------------
second_FP$P1 <- -0.5
second_FP$P2 <- 0
second_FP_sim_15 <- bugs(data = second_FP, inits = second_inits, model.file = second_FP_model_FE, 
                         n.iter = 20000, parameters.to.save = second_params_FE, n.burnin = 10000, 
                         n.chains = 2, 
                         debug = FALSE)

second_FP_sim_15

# Second-order FP, FE model, P1 = -0.5, P2 = 1 ----------------
second_FP$P1 <- -0.5
second_FP$P2 <- 1
second_FP_sim_16 <- bugs(data = second_FP, inits = second_inits, model.file = second_FP_model_FE, 
                         n.iter = 20000, parameters.to.save = second_params_FE, n.burnin = 10000, 
                         n.chains = 2, 
                         debug = FALSE)

second_FP_sim_16

# Second-order FP, FE model, P1 = -0.5, P2 = 2 ----------------
second_FP$P1 <- -0.5
second_FP$P2 <- 2
second_FP_sim_17 <- bugs(data = second_FP, inits = second_inits, model.file = second_FP_model_FE, 
                         n.iter = 20000, parameters.to.save = second_params_FE, n.burnin = 10000, 
                         n.chains = 2, 
                         debug = FALSE)

second_FP_sim_17

# Second-order FP, FE model, P1 = 0, P2 = 0 ----------------
second_FP$P1 <- 0
second_FP$P2 <- 0
second_FP_sim_18 <- bugs(data = second_FP, inits = second_inits, model.file = second_FP_model_FE, 
                         n.iter = 20000, parameters.to.save = second_params_FE, n.burnin = 10000, 
                         n.chains = 2, 
                         debug = FALSE)

second_FP_sim_18

# Second-order FP, FE model, P1 = 0, P2 = 0.5 ----------------
second_FP$P1 <- 0
second_FP$P2 <- 0.5
second_FP_sim_19 <- bugs(data = second_FP, inits = second_inits, model.file = second_FP_model_FE, 
                         n.iter = 20000, parameters.to.save = second_params_FE, n.burnin = 10000, 
                         n.chains = 2, 
                         debug = FALSE)

second_FP_sim_19

# Second-order FP, FE model, P1 = 0, P2 = 1 ----------------
second_FP$P1 <- 0
second_FP$P2 <- 1
second_FP_sim_20 <- bugs(data = second_FP, inits = second_inits, model.file = second_FP_model_FE, 
                         n.iter = 20000, parameters.to.save = second_params_FE, n.burnin = 10000, 
                         n.chains = 2, 
                         debug = FALSE)

second_FP_sim_20

# Second-order FP, FE model, P1 = 0, P2 = 2 ----------------
second_FP$P1 <- 0
second_FP$P2 <- 2
second_FP_sim_21 <- bugs(data = second_FP, inits = second_inits, model.file = second_FP_model_FE, 
                         n.iter = 20000, parameters.to.save = second_params_FE, n.burnin = 10000, 
                         n.chains = 2, 
                         debug = FALSE)

second_FP_sim_21

# Second-order FP, FE model, P1 = 1, P2 = 1 ----------------
second_FP$P1 <- 1
second_FP$P2 <- 1
second_FP_sim_22 <- bugs(data = second_FP, inits = second_inits, model.file = second_FP_model_FE, 
                         n.iter = 20000, parameters.to.save = second_params_FE, n.burnin = 10000, 
                         n.chains = 2, 
                         debug = FALSE)

second_FP_sim_22

# Second-order FP, FE model, P1 = 1, P2 = 2 ----------------
second_FP$P1 <- 1
second_FP$P2 <- 2
second_FP_sim_23 <- bugs(data = second_FP, inits = second_inits, model.file = second_FP_model_FE, 
                         n.iter = 20000, parameters.to.save = second_params_FE, n.burnin = 10000, 
                         n.chains = 2, 
                         debug = FALSE)

second_FP_sim_23

# Second-order FP, FE model, P1 = 2, P2 = 2 ----------------
second_FP$P1 <- 2
second_FP$P2 <- 2
second_FP_sim_24 <- bugs(data = second_FP, inits = second_inits, model.file = second_FP_model_FE, 
                         n.iter = 20000, parameters.to.save = second_params_FE, n.burnin = 10000, 
                         n.chains = 2, 
                         debug = FALSE)

second_FP_sim_24


# 7. Model selection---------------------------
# 7.1 First-order FP
first_FE_models <- list(first_FP_sim1,
                        first_FP_sim2,
                        first_FP_sim3,
                        first_FP_sim4,
                        first_FP_sim5,
                        first_FP_sim6,
                        first_FP_sim7)

first_FE_fit <- lapply(first_FE_models, extract_fit)
first_FE_fit <- unlist(first_FE_fit)
first_FE_fit <- round(matrix(first_FE_fit, nrow = 7, ncol = 4, byrow = TRUE))
p1 <- c(-2, -1, -0.5, 0, 0.5, 1, 2)
first_FE_fit <- cbind(first_FE_fit, p1)
colnames(first_FE_fit) <- c("Deviance", "pD", "DIC", "Residual deviance",  "P1")
first_FE_fit <- arrange(as_tibble(first_FE_fit), DIC)
write.csv(first_FE_fit, here("Fractional polynomials/data", "OS_first_FP_FE_fit.csv"))

# 7.2 Second-order FP
second_FE_models <- list(
  second_FP_sim_1,
  second_FP_sim_2,
  second_FP_sim_3,
  second_FP_sim_4,
  second_FP_sim_5,
  second_FP_sim_6,
  second_FP_sim_7,
  second_FP_sim_8,
  second_FP_sim_9,
  second_FP_sim_10,
  second_FP_sim_11,
  second_FP_sim_12,
  second_FP_sim_13,
  second_FP_sim_14,
  second_FP_sim_15,
  second_FP_sim_16,
  second_FP_sim_17,
  second_FP_sim_18,
  second_FP_sim_19,
  second_FP_sim_20,
  second_FP_sim_21,
  second_FP_sim_22,
  second_FP_sim_23,
  second_FP_sim_24)

second_FE_fit <- lapply(second_FE_models, extract_fit)
second_FE_fit <- unlist(second_FE_fit)
second_FE_fit <- round(matrix(second_FE_fit, nrow = 24, ncol = 4, byrow = TRUE))
p1 <- c(-2, -2, -2, -2, -2, -2, -2, -1, -1, -1, -1, -1, -1, -0.5, -0.5, -0.5, -0.5, 0, 0,   0, 0, 1, 1, 2)
p2 <- c(-2, -1, -0.5, 0, 0.5, 1, 2, -1, -0.5, 0, 0.5, 1, 2, -0.5,    0,    1,    2, 0, 0.5, 1, 2, 1, 2, 2)
second_FE_fit <- cbind(second_FE_fit, p1)
second_FE_fit <- cbind(second_FE_fit, p2)
colnames(second_FE_fit) <- c("Deviance", "pD", "DIC", "totresdev","P1", "P2")
second_FE_fit <- arrange(as_tibble(second_FE_fit), DIC)
write.csv(second_FE_fit, here("Fractional polynomials/data", "OS_second_FP_FE_fit.csv"))

#9. Export HRs for timepoints of interest--------------------------------

HR_params <- c("HR_2year", "HR_3year", "HR_4year")

# Second-order FP, FE model, P1 = -0.5, P2 = -0.5 ----------------
second_FP$P1 <- -0.5
second_FP$P2 <- -0.5
second_FP_sim_HR <- bugs(data = second_FP, inits = second_inits, model.file = second_FP_model_FE, 
                         n.iter = 40000, parameters.to.save = HR_params, n.burnin = 20000, 
                         n.chains = 2, 
                         debug = FALSE)

# HR at 2 years
HR_24 <- 1/second_FP_sim_HR$summary[c("HR_2year[2]", "HR_2year[3]", "HR_2year[4]", 
                             "HR_2year[5]", "HR_2year[6]", "HR_2year[7]", "HR_2year[8]"),]

HR_24 <- round(HR_24[,c(1,7,3)], digits = 2)
HR_24 <- paste0(HR_24[,1], " (", HR_24[,2], ",", HR_24[,3], ")")

# HR at 3 years
HR_36 <- 1/second_FP_sim_HR$summary[c("HR_3year[2]", "HR_3year[3]", "HR_3year[4]", 
                                      "HR_3year[5]", "HR_3year[6]", "HR_3year[7]", "HR_3year[8]"),]

HR_36 <- round(HR_36[,c(1,7,3)], digits = 2)
HR_36 <- paste0(HR_36[,1], " (", HR_36[,2], ",", HR_36[,3], ")")

# HR at 4 years
HR_48 <- 1/second_FP_sim_HR$summary[c("HR_4year[2]", "HR_4year[3]", "HR_4year[4]", 
                                      "HR_4year[5]", "HR_4year[6]", "HR_4year[7]", "HR_4year[8]"),]

HR_48 <- round(HR_48[,c(1,7,3)], digits = 2)
HR_48 <- paste0(HR_48[,1], " (", HR_48[,2], ",", HR_48[,3], ")")

# HR to be exported
HRs <- cbind(HR_24, HR_36, HR_48)
write_csv(as_tibble(HRs), here("Fractional polynomials/data", "time-varying_HRs_OS.csv"))

#10. Code to extract the RMST at 108 months-----------------------------
# Second-order FP, FE model, P1 = 0, P2 = 0.5 (model with the lowest DIC)
second_FP$P1 <- -0.5
second_FP$P2 <- -0.5
second_FP_sim_19 <- bugs(data = second_FP, inits = second_inits, model.file = second_FP_model_FE, 
                         n.iter = 20000, parameters.to.save = "RMST", n.burnin = 10000, 
                         n.chains = 2, 
                         debug = FALSE)

# Absolute RSMT
rmst_daro    <- summary.stats(second_FP_sim_19$sims.array[,1,"RMST[1,108]"])
rmst_doc     <- summary.stats(second_FP_sim_19$sims.array[,1,"RMST[2,108]"])
rmst_enza    <- summary.stats(second_FP_sim_19$sims.array[,1,"RMST[3,108]"])
rmst_adt     <- summary.stats(second_FP_sim_19$sims.array[,1,"RMST[4,108]"])
rmst_abi     <- summary.stats(second_FP_sim_19$sims.array[,1,"RMST[5,108]"])
rmst_apa     <- summary.stats(second_FP_sim_19$sims.array[,1,"RMST[6,108]"])
rmst_enzadoc <- summary.stats(second_FP_sim_19$sims.array[,1,"RMST[7,108]"])
rmst_abidoc  <- summary.stats(second_FP_sim_19$sims.array[,1,"RMST[8,108]"])

rmst <- c(rmst_daro, rmst_doc, rmst_enza, rmst_adt, rmst_abi, rmst_apa,
          rmst_enzadoc, rmst_abidoc)

# Relative RSMT to reference treatment
diff_doc <- summary.stats(second_FP_sim_19$sims.array[,1,"RMST[2,108]"] - 
                second_FP_sim_19$sims.array[,1,"RMST[1,108]"])
diff_enza <- summary.stats(second_FP_sim_19$sims.array[,1,"RMST[3,108]"] - 
                second_FP_sim_19$sims.array[,1,"RMST[1,108]"])
diff_adt <- summary.stats(second_FP_sim_19$sims.array[,1,"RMST[4,108]"] - 
                second_FP_sim_19$sims.array[,1,"RMST[1,108]"])
diff_abi <- summary.stats(second_FP_sim_19$sims.array[,1,"RMST[5,108]"] - 
                second_FP_sim_19$sims.array[,1,"RMST[1,108]"])
diff_apa <- summary.stats(second_FP_sim_19$sims.array[,1,"RMST[6,108]"] - 
                second_FP_sim_19$sims.array[,1,"RMST[1,108]"])
diff_enzadoc <- summary.stats(second_FP_sim_19$sims.array[,1,"RMST[7,108]"] - 
                second_FP_sim_19$sims.array[,1,"RMST[1,108]"])
diff_abidoc <- summary.stats(second_FP_sim_19$sims.array[,1,"RMST[8,108]"] - 
                second_FP_sim_19$sims.array[,1,"RMST[1,108]"])

diff <- c(NA, diff_doc, diff_enza, diff_adt, diff_abi, diff_apa,
          diff_enzadoc, diff_abidoc)

# Export results
write.csv(bind_cols(rmst, diff), here("Fractional polynomials/data", "rmst_os.csv"))


#11. Extrapolaion of models with the lowest DIC-------------------------------
# Second-order FP, FE model, P1 = -0.5, P2 = -0.5 ----------------
second_FP$P1 <- -0.5
second_FP$P2 <- -0.5
second_FP$maxt <- 108 # increasing maximum time to extrapolate survival
second_FP_sim_14 <- bugs(data = second_FP, inits = second_inits, model.file = second_FP_model_FE, 
                         n.iter = 40000, parameters.to.save = "S", n.burnin = 20000, 
                         n.chains = 2, 
                         debug = FALSE)

# export summary data for survival curves
results_second_FP_sim_14 <- second_FP_sim_14$summary
write.csv(results_second_FP_sim_14, here("Fractional polynomials/data", "OS_results_FP_sim14_108.csv"))

# export array for model averaging (only a sample to reduce computation time)
array_second_FP_sim_14 <- second_FP_sim_14$sims.array
time_point <- c("S[1,60]", "S[2,60]", "S[3,60]",
                "S[4,60]", "S[5,60]", "S[6,60]",
                "S[7,60]", "S[8,60]") # vector to extract timepoints of interest
array_second_FP_sim_14 <- array_second_FP_sim_14[, 1, time_point]
array_second_FP_sim_14 <- apply(array_second_FP_sim_14, 2, sample, size = 500) 
write.csv(array_second_FP_sim_14, here("Fractional polynomials/data", "OS_array_FP_sim14.csv"))

# Second-order FP, FE model, P1 = -0.5, P2 = 0------------------
second_FP$P1 <- -0.5
second_FP$P2 <- 0
second_FP$maxt <- 108 # increasing maximum time to extrapolate survival
second_FP_sim_15 <- bugs(data = second_FP, inits = second_inits, model.file = second_FP_model_FE, 
                         n.iter = 40000, parameters.to.save = "S", n.burnin = 20000, 
                         n.chains = 2, 
                         debug = FALSE)

results_second_FP_sim_15 <- second_FP_sim_15$summary
write.csv(results_second_FP_sim_15, here("Fractional polynomials/data", "OS_results_FP_sim15_108.csv"))

array_second_FP_sim_15 <- second_FP_sim_15$sims.array
array_second_FP_sim_15 <- array_second_FP_sim_15[, 1, time_point]
array_second_FP_sim_15 <- apply(array_second_FP_sim_15, 2, sample, size = 500)
write.csv(array_second_FP_sim_15, here("Fractional polynomials/data", "OS_array_FP_sim15.csv"))

# Second-order FP, FE model, P1 = -1, P2 = 0.5------------------
second_FP$P1 <- -1
second_FP$P2 <- 0.5
second_FP$maxt <- 108 # increasing maximum time to extrapolate survival
second_FP_sim_11 <- bugs(data = second_FP, inits = second_inits, model.file = second_FP_model_FE, 
                         n.iter = 40000, parameters.to.save = "S", n.burnin = 20000, 
                         n.chains = 2, 
                         debug = FALSE)

results_second_FP_sim_11 <- second_FP_sim_11$summary
write.csv(results_second_FP_sim_11, here("Fractional polynomials/data", "OS_results_FP_sim11_108.csv"))

array_second_FP_sim_11 <- second_FP_sim_11$sims.array
array_second_FP_sim_11 <- array_second_FP_sim_11[, 1, time_point]
array_second_FP_sim_11 <- apply(array_second_FP_sim_11, 2, sample, size = 500)
write.csv(array_second_FP_sim_11, here("Fractional polynomials/data", "OS_array_FP_sim11.csv"))

# Second-order FP, FE model, P1 = -1, P2 = 0------------------
second_FP$P1 <- -1
second_FP$P2 <- 0
second_FP$maxt <- 108 # increasing maximum time to extrapolate survival
second_FP_sim_10 <- bugs(data = second_FP, inits = second_inits, model.file = second_FP_model_FE, 
                         n.iter = 40000, parameters.to.save = "S", n.burnin = 20000, 
                         n.chains = 2, 
                         debug = FALSE)

results_second_FP_sim_10 <- second_FP_sim_10$summary
write.csv(results_second_FP_sim_10, here("Fractional polynomials/data", "OS_results_FP_sim10_108.csv"))

array_second_FP_sim_10 <- second_FP_sim_10$sims.array
array_second_FP_sim_10 <- array_second_FP_sim_10[, 1, time_point]
array_second_FP_sim_10 <- apply(array_second_FP_sim_10, 2, sample, size = 500)
write.csv(array_second_FP_sim_10, here("Fractional polynomials/data", "OS_array_FP_sim10.csv"))

# Second-order FP, FE model, P1 = -1, P2 = -0.5------------------
second_FP$P1 <- -1
second_FP$P2 <- -0.5
second_FP$maxt <- 108 # increasing maximum time to extrapolate survival
second_FP_sim_09 <- bugs(data = second_FP, inits = second_inits, model.file = second_FP_model_FE, 
                         n.iter = 40000, parameters.to.save = "S", n.burnin = 20000, 
                         n.chains = 2, 
                         debug = FALSE)

results_second_FP_sim_09 <- second_FP_sim_09$summary
write.csv(results_second_FP_sim_09, here("Fractional polynomials/data", "OS_results_FP_sim09_108.csv"))

array_second_FP_sim_09 <- second_FP_sim_09$sims.array
array_second_FP_sim_09 <- array_second_FP_sim_09[, 1, time_point]
array_second_FP_sim_09 <- apply(array_second_FP_sim_09, 2, sample, size = 500)
write.csv(array_second_FP_sim_09, here("Fractional polynomials/data", "OS_array_FP_sim09.csv"))

# ARCHES as the trial of interest-------------------------------------
# call the model
second_FP_model_FE <- here("Fractional polynomials/code", "second_FP_model_FE_ARCHES.txt")

# Second-order FP, FE model, P1 = -0.5, P2 = -0.5 ----------------
second_FP$P1 <- -0.5
second_FP$P2 <- -0.5
second_FP$maxt <- 108 # increasing maximum time to extrapolate survival
second_FP_sim_14 <- bugs(data = second_FP, inits = second_inits, model.file = second_FP_model_FE, 
                         n.iter = 40000, parameters.to.save = "S", n.burnin = 20000, 
                         n.chains = 2, 
                         debug = FALSE)

# export summary data for survival curves
results_second_FP_sim_14 <- second_FP_sim_14$summary
write.csv(results_second_FP_sim_14, here("Fractional polynomials/data", "OS_results_FP_sim14_ARCHES.csv"))

# Second-order FP, FE model, P1 = -0.5, P2 = 0------------------
second_FP$P1 <- -0.5
second_FP$P2 <- 0
second_FP$maxt <- 108 # increasing maximum time to extrapolate survival
second_FP_sim_15 <- bugs(data = second_FP, inits = second_inits, model.file = second_FP_model_FE, 
                         n.iter = 40000, parameters.to.save = "S", n.burnin = 20000, 
                         n.chains = 2, 
                         debug = FALSE)

results_second_FP_sim_15 <- second_FP_sim_15$summary
write.csv(results_second_FP_sim_15, here("Fractional polynomials/data", "OS_results_FP_sim15_ARCHES.csv"))

# Second-order FP, FE model, P1 = -1, P2 = 0.5------------------
second_FP$P1 <- -1
second_FP$P2 <- 0.5
second_FP$maxt <- 108 # increasing maximum time to extrapolate survival
second_FP_sim_11 <- bugs(data = second_FP, inits = second_inits, model.file = second_FP_model_FE, 
                         n.iter = 40000, parameters.to.save = "S", n.burnin = 20000, 
                         n.chains = 2, 
                         debug = FALSE)

results_second_FP_sim_11 <- second_FP_sim_11$summary
write.csv(results_second_FP_sim_11, here("Fractional polynomials/data", "OS_results_FP_sim11_ARCHES.csv"))

# Second-order FP, FE model, P1 = -1, P2 = 0------------------
second_FP$P1 <- -1
second_FP$P2 <- 0
second_FP$maxt <- 108 # increasing maximum time to extrapolate survival
second_FP_sim_10 <- bugs(data = second_FP, inits = second_inits, model.file = second_FP_model_FE, 
                         n.iter = 40000, parameters.to.save = "S", n.burnin = 20000, 
                         n.chains = 2, 
                         debug = FALSE)

results_second_FP_sim_10 <- second_FP_sim_10$summary
write.csv(results_second_FP_sim_10, here("Fractional polynomials/data", "OS_results_FP_sim10_ARCHES.csv"))

# Second-order FP, FE model, P1 = -1, P2 = -0.5------------------
second_FP$P1 <- -1
second_FP$P2 <- -0.5
second_FP$maxt <- 108 # increasing maximum time to extrapolate survival
second_FP_sim_09 <- bugs(data = second_FP, inits = second_inits, model.file = second_FP_model_FE, 
                         n.iter = 40000, parameters.to.save = "S", n.burnin = 20000, 
                         n.chains = 2, 
                         debug = FALSE)

results_second_FP_sim_09 <- second_FP_sim_09$summary
write.csv(results_second_FP_sim_09, here("Fractional polynomials/data", "OS_results_FP_sim09_ARCHES.csv"))

