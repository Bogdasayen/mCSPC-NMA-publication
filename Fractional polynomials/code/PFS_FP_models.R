# Darulotamide in mHSPC NMA #
# Outcome: PFS (KM curves data)

# 1. Packages used --------------------------------------------------------
pkgs <- c("tidyverse", "readxl", "here", "haven", "multinma", "R2OpenBUGS",
          "igraph")
lapply(pkgs, library, character.only = T)
set.seed(03082022)
#source(here ("01_scripts", "cross_or.R")) 
source(here("01_scripts", "utils.R")) 

# 2. Import data --------------------------------------------------------
PFS_data <- read_excel(here("02_data", "mHSPC PFS trt codes.xlsx"))

# Set up the evidence network
PFS_long <- pivot_longer(PFS_data, cols = c(t1, t2), names_to = "trt", values_to = "treatment")
PFS_net <- set_agd_arm(PFS_long, study = study, trt = treatment, y = y, se = se)
graph <- as.igraph(PFS_net)
plot(graph)

# Import data from KM curves
FP_data <- read_csv(here("02_data", "FP data", "FP_data_PFS.csv"))
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

# 3.2 Fit Fixed Effects 1st order fractional polynomial model----------------
# Call the model
first_FP_model_FE <- here("03_models", "first_FP_model_FE.txt")

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

# Parameters to monitor
first_params_FE <- c("d", "totresdev")

# First order FP P1 = -2---------------------------------------------
first_FP$P1 <- -2
first_FP_sim1 <- bugs(data = first_FP, inits = inits,
                      model.file = first_FP_model_FE, n.iter = 20000,
                      parameters.to.save = first_params_FE, n.burnin = 10000, n.chains = 2,
                      debug = FALSE)

first_FP_sim1

# First order FP P1 = -1---------------------------------------------
first_FP$P1 <- -1
first_FP_sim2 <- bugs(data = first_FP, inits = inits,
                      model.file = first_FP_model_FE, n.iter = 20000,
                      parameters.to.save = first_params_FE, n.burnin = 10000, n.chains = 2,
                      debug = FALSE)

first_FP_sim2

# First order FP P1 = -0.5---------------------------------------------
first_FP$P1 <- -0.5
first_FP_sim3 <- bugs(data = first_FP, inits = inits,
                      model.file = first_FP_model_FE, n.iter = 20000,
                      parameters.to.save = first_params_FE, n.burnin = 10000, n.chains = 2,
                      debug = FALSE)

first_FP_sim3

# First order FP P1 = 0---------------------------------------------
first_FP$P1 <- 0
first_FP_sim4 <- bugs(data = first_FP, inits = inits,
                      model.file = first_FP_model_FE, n.iter = 20000,
                      parameters.to.save = first_params_FE, n.burnin = 10000, n.chains = 2,
                      debug = FALSE)

first_FP_sim4

# First order FP P1 = 0.5---------------------------------------------
first_FP$P1 <- 0.5
first_FP_sim5 <- bugs(data = first_FP, inits = inits,
                      model.file = first_FP_model_FE, n.iter = 20000,
                      parameters.to.save = first_params_FE, n.burnin = 10000, n.chains = 2,
                      debug = FALSE)

first_FP_sim5

# First order FP P1 = 1---------------------------------------------
first_FP$P1 <- 1
first_FP_sim6 <- bugs(data = first_FP, inits = inits,
                      model.file = first_FP_model_FE, n.iter = 20000,
                      parameters.to.save = first_params_FE, n.burnin = 10000, n.chains = 2,
                      debug = FALSE)

first_FP_sim6

# First order FP P1 = 2---------------------------------------------
first_FP$P1 <- 2
first_FP_sim7 <- bugs(data = first_FP, inits = inits,
                      model.file = first_FP_model_FE, n.iter = 20000,
                      parameters.to.save = first_params_FE, n.burnin = 10000, n.chains = 2,
                      debug = FALSE)

first_FP_sim7

# 4. Random effects models  -------------------------------
# Call the model
first_FP_model_RE <- here("03_models", "first_FP_model_RE.txt")

#Initial Values 
inits_RE <- function(){
  #chain 1
  list(d= structure(.Data= rep(c(NA, rep(-0.1, first_FP$nt - 1)), 2),
                    .Dim=c(first_FP$nt, 2)), 
       mu= structure(.Data= c(rep(0.1, 2 * first_FP$ns)), 
                     .Dim=c(first_FP$ns, 2)),
       sd = 0.1)
  #chain 2
  list(d= structure(.Data= rep(c(NA, rep(1, first_FP$nt - 1)), 2),
                    .Dim=c(first_FP$nt, 2)),  
       mu = structure(.Data= c(rep(1, 2 * first_FP$ns)),
                      .Dim=c(first_FP$ns, 2)),
       sd = 1)
}

# Parameters to monitor
first_params_RE <- c("d", "totresdev", "sd")

# 4.2 Fit random effects NMA models
# First order FP P1 = -2---------------------------------------------
first_FP$P1 <- -2
first_FP_sim1_RE <- bugs(data = first_FP, inits = inits_RE,
                      model.file = first_FP_model_RE, n.iter = 20000,
                      parameters.to.save = first_params_RE, n.burnin = 10000, n.chains = 2,
                      debug = FALSE)

first_FP_sim1_RE

# First order FP P1 = -1---------------------------------------------
first_FP$P1 <- -1
first_FP_sim2_RE <- bugs(data = first_FP, inits = inits_RE,
                      model.file = first_FP_model_RE, n.iter = 20000,
                      parameters.to.save = first_params_RE, n.burnin = 10000, n.chains = 2,
                      debug = FALSE)

first_FP_sim2_RE

# First order FP P1 = -0.5---------------------------------------------
first_FP$P1 <- -0.5
first_FP_sim3_RE <- bugs(data = first_FP, inits = inits_RE,
                      model.file = first_FP_model_RE, n.iter = 20000,
                      parameters.to.save = first_params_RE, n.burnin = 10000, n.chains = 2,
                      debug = FALSE)

first_FP_sim3_RE

# First order FP P1 = 0---------------------------------------------
first_FP$P1 <- 0
first_FP_sim4_RE <- bugs(data = first_FP, inits = inits_RE,
                      model.file = first_FP_model_RE, n.iter = 20000,
                      parameters.to.save = first_params_RE, n.burnin = 10000, n.chains = 2,
                      debug = FALSE)

first_FP_sim4_RE

# First order FP P1 = 0.5---------------------------------------------
first_FP$P1 <- 0.5
first_FP_sim5_RE <- bugs(data = first_FP, inits = inits_RE,
                      model.file = first_FP_model_RE, n.iter = 20000,
                      parameters.to.save = first_params_RE, n.burnin = 10000, n.chains = 2,
                      debug = FALSE)

first_FP_sim5_RE

# First order FP P1 = 1---------------------------------------------
first_FP$P1 <- 1
first_FP_sim6_RE <- bugs(data = first_FP, inits = inits_RE,
                      model.file = first_FP_model_RE, n.iter = 20000,
                      parameters.to.save = first_params_RE, n.burnin = 10000, n.chains = 2,
                      debug = FALSE)

first_FP_sim6_RE

# First order FP P1 = 2---------------------------------------------
first_FP$P1 <- 2
first_FP_sim7_RE <- bugs(data = first_FP, inits = inits_RE,
                      model.file = first_FP_model_RE, n.iter = 20000,
                      parameters.to.save = first_params_RE, n.burnin = 10000, n.chains = 2,
                      debug = FALSE)

first_FP_sim7_RE


# 4.3 Assess convergence and model fit

# 5. Fixed Effects 2nd order fractional polynomial model----------------
# call the model
second_FP_model_FE <- here("03_models", "second_FP_model_FE.txt")

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

# Parameters to monitor
second_params_FE <- c("d", "totresdev")

# Second-order FP, FE model, P1 = -2, P2 = -2 ----------------
second_FP$P1 <- -2
second_FP$P2 <- -2
second_FP_sim_1 <- bugs(data = second_FP, inits = second_inits, model.file = second_FP_model_FE, 
                        n.iter = 20000, parameters.to.save = second_params_FE, n.burnin = 10000, 
                        n.chains = 2, 
                        debug = TRUE)

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


# 6. Second-order FP Random effects models  -------------------------------
# call the model
second_FP_model_RE <- here("03_models", "second_FP_model_RE.txt")

# data for 2nd order FP model
second_FP <- first_FP
second_FP$P2 <- 0
second_FP$mean <- c(0, 0, 0)
second_FP$prec <- structure(.Data= c(1.00000E-04, 0, 0, 
                                     0, 1.00000E-02, 0,
                                     0, 0, 1.00000E-02), .Dim=c(3, 3))

#Initial Values 
second_inits_RE <- function(){
  #chain 1
  list(d= structure(.Data= rep(c(NA, rep(-0.1, first_FP$nt - 1)), 3),
                    .Dim=c(first_FP$nt, 3)), 
       mu= structure(.Data= c(rep(0.1, 3 * first_FP$ns)), 
                     .Dim=c(first_FP$ns, 3)),
       sd = 0.1)
  #chain 2
  list(d= structure(.Data= rep(c(NA, rep(1, first_FP$nt - 1)), 3),
                    .Dim=c(first_FP$nt, 3)),  
       mu = structure(.Data= c(rep(1, 3 * first_FP$ns)),
                      .Dim=c(first_FP$ns, 3)),
       sd = 1)
}

# Parameters to monitor
second_params_RE <- c("d", "totresdev")

# Second-order FP, FE model, P1 = -2, P2 = -2 ----------------
second_FP$P1 <- -2
second_FP$P2 <- -2
second_FP_sim_1_RE <- bugs(data = second_FP, inits = second_inits_RE, model.file = second_FP_model_RE, 
                        n.iter = 20000, parameters.to.save = second_params_RE, n.burnin = 10000, 
                        n.chains = 2, 
                        debug = FALSE)

second_FP_sim_1_RE

# Second-order FP, FE model, P1 = -2, P2 = -1 ----------------
second_FP$P1 <- -2
second_FP$P2 <- -1
second_FP_sim_2_RE <- bugs(data = second_FP, inits = second_inits_RE, model.file = second_FP_model_RE, 
                        n.iter = 20000, parameters.to.save = second_params_RE, n.burnin = 10000, 
                        n.chains = 2, 
                        debug = FALSE)

second_FP_sim_2_RE

# Second-order FP, FE model, P1 = -2, P2 = -0.5 ----------------
second_FP$P1 <- -2
second_FP$P2 <- -0.5
second_FP_sim_3_RE <- bugs(data = second_FP, inits = second_inits_RE, model.file = second_FP_model_RE, 
                        n.iter = 20000, parameters.to.save = second_params_RE, n.burnin = 10000, 
                        n.chains = 2, 
                        debug = FALSE)

second_FP_sim_3_RE

# Second-order FP, FE model, P1 = -2, P2 = 0 ----------------
second_FP$P1 <- -2
second_FP$P2 <- 0
second_FP_sim_4_RE <- bugs(data = second_FP, inits = second_inits_RE, model.file = second_FP_model_RE, 
                        n.iter = 20000, parameters.to.save = second_params_RE, n.burnin = 10000, 
                        n.chains = 2, 
                        debug = FALSE)

second_FP_sim_4_RE

# Second-order FP, FE model, P1 = -2, P2 = 0.5 ----------------
second_FP$P1 <- -2
second_FP$P2 <- 0.5
second_FP_sim_5_RE <- bugs(data = second_FP, inits = second_inits_RE, model.file = second_FP_model_RE, 
                        n.iter = 20000, parameters.to.save = second_params_RE, n.burnin = 10000, 
                        n.chains = 2, 
                        debug = FALSE)

second_FP_sim_5_RE

# Second-order FP, FE model, P1 = -2, P2 = 1 ----------------
second_FP$P1 <- -2
second_FP$P2 <- 1
second_FP_sim_6_RE <- bugs(data = second_FP, inits = second_inits_RE, model.file = second_FP_model_RE, 
                        n.iter = 20000, parameters.to.save = second_params_RE, n.burnin = 10000, 
                        n.chains = 2, 
                        debug = FALSE)

second_FP_sim_6_RE


# Second-order FP, FE model, P1 = -2, P2 = 2 ----------------
second_FP$P1 <- -2
second_FP$P2 <- 2
second_FP_sim_7_RE <- bugs(data = second_FP, inits = second_inits_RE, model.file = second_FP_model_RE, 
                        n.iter = 20000, parameters.to.save = second_params_RE, n.burnin = 10000, 
                        n.chains = 2, 
                        debug = FALSE)

second_FP_sim_7_RE

# Second-order FP, FE model, P1 = -1, P2 = -1 ----------------
second_FP$P1 <- -1
second_FP$P2 <- -1
second_FP_sim_8_RE <- bugs(data = second_FP, inits = second_inits_RE, model.file = second_FP_model_RE, 
                        n.iter = 20000, parameters.to.save = second_params_RE, n.burnin = 10000, 
                        n.chains = 2, 
                        debug = FALSE)

second_FP_sim_8_RE

# Second-order FP, FE model, P1 = -1, P2 = -0.5 ----------------
second_FP$P1 <- -1
second_FP$P2 <- -0.5
second_FP_sim_9_RE <- bugs(data = second_FP, inits = second_inits_RE, model.file = second_FP_model_RE, 
                        n.iter = 20000, parameters.to.save = second_params_RE, n.burnin = 10000, 
                        n.chains = 2, 
                        debug = FALSE)

second_FP_sim_9_RE

# Second-order FP, FE model, P1 = -1, P2 = 0 ----------------
second_FP$P1 <- -1
second_FP$P2 <- 0
second_FP_sim_10_RE <- bugs(data = second_FP, inits = second_inits_RE, model.file = second_FP_model_RE, 
                         n.iter = 20000, parameters.to.save = second_params_RE, n.burnin = 10000, 
                         n.chains = 2, 
                         debug = FALSE)

second_FP_sim_10_RE

# Second-order FP, FE model, P1 = -1, P2 = 0.5 ----------------
second_FP$P1 <- -1
second_FP$P2 <- 0.5
second_FP_sim_11_RE <- bugs(data = second_FP, inits = second_inits_RE, model.file = second_FP_model_RE, 
                         n.iter = 20000, parameters.to.save = second_params_RE, n.burnin = 10000, 
                         n.chains = 2, 
                         debug = FALSE)

second_FP_sim_11_RE


# Second-order FP, FE model, P1 = -1, P2 = 1 ----------------
second_FP$P1 <- -1
second_FP$P2 <- 1
second_FP_sim_12_RE <- bugs(data = second_FP, inits = second_inits_RE, model.file = second_FP_model_RE, 
                         n.iter = 20000, parameters.to.save = second_params_RE, n.burnin = 10000, 
                         n.chains = 2, 
                         debug = FALSE)

second_FP_sim_12_RE

# Second-order FP, FE model, P1 = -1, P2 = 2 ----------------
second_FP$P1 <- -1
second_FP$P2 <- 2
second_FP_sim_13_RE <- bugs(data = second_FP, inits = second_inits_RE, model.file = second_FP_model_RE, 
                         n.iter = 20000, parameters.to.save = second_params_RE, n.burnin = 10000, 
                         n.chains = 2, 
                         debug = FALSE)

second_FP_sim_13_RE

# Second-order FP, FE model, P1 = -0.5, P2 = -0.5 ----------------
second_FP$P1 <- -0.5
second_FP$P2 <- -0.5
second_FP_sim_14_RE <- bugs(data = second_FP, inits = second_inits_RE, model.file = second_FP_model_RE, 
                         n.iter = 20000, parameters.to.save = second_params_RE, n.burnin = 10000, 
                         n.chains = 2, 
                         debug = FALSE)

second_FP_sim_14_RE

# Second-order FP, FE model, P1 = -0.5, P2 = 0 ----------------
second_FP$P1 <- -0.5
second_FP$P2 <- 0
second_FP_sim_15_RE <- bugs(data = second_FP, inits = second_inits_RE, model.file = second_FP_model_RE, 
                         n.iter = 20000, parameters.to.save = second_params_RE, n.burnin = 10000, 
                         n.chains = 2, 
                         debug = FALSE)

second_FP_sim_15_RE

# Second-order FP, FE model, P1 = -0.5, P2 = 1 ----------------
second_FP$P1 <- -0.5
second_FP$P2 <- 1
second_FP_sim_16_RE <- bugs(data = second_FP, inits = second_inits_RE, model.file = second_FP_model_RE, 
                         n.iter = 20000, parameters.to.save = second_params_RE, n.burnin = 10000, 
                         n.chains = 2, 
                         debug = FALSE)

second_FP_sim_16_RE

# Second-order FP, FE model, P1 = -0.5, P2 = 2 ----------------
second_FP$P1 <- -0.5
second_FP$P2 <- 2
second_FP_sim_17_RE <- bugs(data = second_FP, inits = second_inits_RE, model.file = second_FP_model_RE, 
                         n.iter = 20000, parameters.to.save = second_params_RE, n.burnin = 10000, 
                         n.chains = 2, 
                         debug = FALSE)

second_FP_sim_17_RE

# Second-order FP, FE model, P1 = 0, P2 = 0 ----------------
second_FP$P1 <- 0
second_FP$P2 <- 0
second_FP_sim_18_RE <- bugs(data = second_FP, inits = second_inits_RE, model.file = second_FP_model_RE, 
                         n.iter = 20000, parameters.to.save = second_params_RE, n.burnin = 10000, 
                         n.chains = 2, 
                         debug = FALSE)

second_FP_sim_18_RE

# Second-order FP, FE model, P1 = 0, P2 = 0.5 ----------------
second_FP$P1 <- 0
second_FP$P2 <- 0.5
second_FP_sim_19_RE <- bugs(data = second_FP, inits = second_inits_RE, model.file = second_FP_model_RE, 
                         n.iter = 20000, parameters.to.save = second_params_RE, n.burnin = 10000, 
                         n.chains = 2, 
                         debug = FALSE)

second_FP_sim_19_RE

# Second-order FP, FE model, P1 = 0, P2 = 1 ----------------
second_FP$P1 <- 0
second_FP$P2 <- 1
second_FP_sim_20_RE <- bugs(data = second_FP, inits = second_inits_RE, model.file = second_FP_model_RE, 
                         n.iter = 20000, parameters.to.save = second_params_RE, n.burnin = 10000, 
                         n.chains = 2, 
                         debug = FALSE)

second_FP_sim_20_RE

# Second-order FP, FE model, P1 = 0, P2 = 2 ----------------
second_FP$P1 <- 0
second_FP$P2 <- 2
second_FP_sim_21_RE <- bugs(data = second_FP, inits = second_inits_RE, model.file = second_FP_model_RE, 
                         n.iter = 20000, parameters.to.save = second_params_RE, n.burnin = 10000, 
                         n.chains = 2, 
                         debug = FALSE)

second_FP_sim_21_RE

# Second-order FP, FE model, P1 = 1, P2 = 1 ----------------
second_FP$P1 <- 1
second_FP$P2 <- 1
second_FP_sim_22_RE <- bugs(data = second_FP, inits = second_inits_RE, model.file = second_FP_model_RE, 
                         n.iter = 20000, parameters.to.save = second_params_RE, n.burnin = 10000, 
                         n.chains = 2, 
                         debug = FALSE)

second_FP_sim_22_RE

# Second-order FP, FE model, P1 = 1, P2 = 2 ----------------
second_FP$P1 <- 1
second_FP$P2 <- 2
second_FP_sim_23_RE <- bugs(data = second_FP, inits = second_inits_RE, model.file = second_FP_model_RE, 
                         n.iter = 20000, parameters.to.save = second_params_RE, n.burnin = 10000, 
                         n.chains = 2, 
                         debug = FALSE)

second_FP_sim_23_RE

# Second-order FP, FE model, P1 = 2, P2 = 2 ----------------
second_FP$P1 <- 2
second_FP$P2 <- 2
second_FP_sim_24_RE <- bugs(data = second_FP, inits = second_inits_RE, model.file = second_FP_model_RE, 
                         n.iter = 20000, parameters.to.save = second_params_RE, n.burnin = 10000, 
                         n.chains = 2, 
                         debug = FALSE)

second_FP_sim_24_RE


# 7. Model selection---------------------------
# 7.1 First-order FP FE models
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
colnames(first_FE_fit) <- c("Deviance", "pD", "DIC", "Residual deviance", "P1")
write.csv(first_FE_fit, here("05_tables", "PFS_first_FP_FE_fit.csv"))

# 7.2 First-order FP RE models
first_RE_models <- list(first_FP_sim1_RE,
                        first_FP_sim2_RE,
                        first_FP_sim3_RE,
                        first_FP_sim4_RE,
                        first_FP_sim5_RE,
                        first_FP_sim6_RE,
                        first_FP_sim7_RE)

first_RE_fit <- lapply(first_RE_models, extract_fit)
first_RE_fit <- unlist(first_RE_fit)
first_RE_fit <- round(matrix(first_RE_fit, nrow = 7, ncol = 4, byrow = TRUE))
p1 <- c(-2, -1, -0.5, 0, 0.5, 1, 2)
first_RE_fit <- cbind(first_RE_fit, p1)
colnames(first_RE_fit) <- c("Deviance", "pD", "DIC", "Residual deviance", "P1")
write.csv(first_RE_fit, here("05_tables", "PFS_first_FP_RE_fit.csv"))

# 7.3 Second-order FP FE models
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
colnames(second_FE_fit) <- c("Deviance", "pD", "DIC", "Residual deviance", "P1", "P2")
write.csv(second_FE_fit, here("05_tables", "PFS_second_FP_FE_fit.csv"))

# 7.4 Second-order FP RE models
second_RE_models <- list(
  second_FP_sim_1_RE,
  second_FP_sim_2_RE,
  second_FP_sim_3_RE,
  second_FP_sim_4_RE,
  second_FP_sim_5_RE,
  second_FP_sim_6_RE,
  second_FP_sim_7_RE,
  second_FP_sim_8_RE,
  second_FP_sim_9_RE,
  second_FP_sim_10_RE,
  second_FP_sim_11_RE,
  second_FP_sim_12_RE,
  second_FP_sim_13_RE,
  second_FP_sim_14_RE,
  second_FP_sim_15_RE,
  second_FP_sim_16_RE,
  second_FP_sim_17_RE,
  second_FP_sim_18_RE,
  second_FP_sim_19_RE,
  second_FP_sim_20_RE,
  second_FP_sim_21_RE,
  second_FP_sim_22_RE,
  second_FP_sim_23_RE,
  second_FP_sim_24_RE)

second_RE_fit <- lapply(second_RE_models, extract_fit)
second_RE_fit <- unlist(second_RE_fit)
second_RE_fit <- round(matrix(second_RE_fit, nrow = 24, ncol = 4, byrow = TRUE))
p1 <- c(-2, -2, -2, -2, -2, -2, -2, -1, -1, -1, -1, -1, -1, -0.5, -0.5, -0.5, -0.5, 0, 0,   0, 0, 1, 1, 2)
p2 <- c(-2, -1, -0.5, 0, 0.5, 1, 2, -1, -0.5, 0, 0.5, 1, 2, -0.5,    0,    1,    2, 0, 0.5, 1, 2, 1, 2, 2)
second_RE_fit <- cbind(second_RE_fit, p1)
second_RE_fit <- cbind(second_RE_fit, p2)
colnames(second_RE_fit) <- c("Deviance", "pD", "DIC", "Residual deviance", "P1", "P2")
write.csv(second_RE_fit, here("05_tables", "PFS_second_FP_RE_fit.csv"))

#8 Models with lowest DIC for extrapolation --------------------------
arrange(as_tibble(first_FE_fit), `Residual deviance`)
arrange(as_tibble(second_FE_fit), `Residual deviance`)

#8.1 ARASENS as the trial of interest----------------------------------
# Second-order FP, FE model, P1 = -0.5, P2 = -0.5 ----------------
second_FP$P1 <- -0.5
second_FP$P2 <- -0.5
second_FP$maxt <- 108 # increasing maximum time to extrapolate survival
second_FP_sim_14 <- bugs(data = second_FP, inits = second_inits, model.file = second_FP_model_FE, 
                         n.iter = 20000, parameters.to.save = "S", n.burnin = 10000, 
                         n.chains = 2, 
                         debug = FALSE)

results_second_FP_sim_14 <- second_FP_sim_14$summary
write.csv(results_second_FP_sim_14, here("02_data", "PFS_results_FP_sim14_108.csv"))

# export array for model averaging (only a sample to reduce computation time)
array_second_FP_sim_14 <- second_FP_sim_14$sims.array
time_point <- c("S[1,36]", "S[2,36]", "S[3,36]",
                "S[4,36]", "S[5,36]", "S[6,36]",
                "S[7,36]", "S[8,36]") # vector to extract timepoints of interest
array_second_FP_sim_14 <- array_second_FP_sim_14[, 1, time_point]
array_second_FP_sim_14 <- apply(array_second_FP_sim_14, 2, sample, size = 500) 
write.csv(array_second_FP_sim_14, here("02_data", "PFS_array_FP_sim14.csv"))

# Second-order FP, FE model, P1 = -1, P2 = 0.5 ----------------
second_FP$P1 <- -1
second_FP$P2 <- 0.5
second_FP$maxt <- 108 # increasing maximum time to extrapolate survival
second_FP_sim_11 <- bugs(data = second_FP, inits = second_inits, model.file = second_FP_model_FE, 
                         n.iter = 20000, parameters.to.save = "S", n.burnin = 10000, 
                         n.chains = 2, 
                         debug = FALSE)

results_second_FP_sim_11 <- second_FP_sim_11$summary
write.csv(results_second_FP_sim_11, here("02_data", "PFS_results_FP_sim11_108.csv"))

array_second_FP_sim_11 <- second_FP_sim_11$sims.array
array_second_FP_sim_11 <- array_second_FP_sim_11[, 1, time_point]
array_second_FP_sim_11 <- apply(array_second_FP_sim_11, 2, sample, size = 500)
write.csv(array_second_FP_sim_11, here("02_data", "PFS_array_FP_sim11.csv"))

# Second-order FP, FE model, P1 = 0, P2 = 0.5 ----------------
second_FP$P1 <- 0
second_FP$P2 <- 0.5
second_FP$maxt <- 108 # increasing maximum time to extrapolate survival
second_FP_sim_19 <- bugs(data = second_FP, inits = second_inits, model.file = second_FP_model_FE, 
                         n.iter = 20000, parameters.to.save = "S", n.burnin = 10000, 
                         n.chains = 2, 
                         debug = FALSE)

results_second_FP_sim_19 <- second_FP_sim_19$summary
write.csv(results_second_FP_sim_19, here("02_data", "PFS_results_FP_sim19_108.csv"))

array_second_FP_sim_19 <- second_FP_sim_19$sims.array
array_second_FP_sim_19 <- array_second_FP_sim_19[, 1, time_point]
array_second_FP_sim_19 <- apply(array_second_FP_sim_19, 2, sample, size = 500)
write.csv(array_second_FP_sim_19, here("02_data", "PFS_array_FP_sim19.csv"))

# Second-order FP, FE model, P1 = -1, P2 = 0.5 ----------------
second_FP$P1 <- -1
second_FP$P2 <- 0
second_FP$maxt <- 108 # increasing maximum time to extrapolate survival
second_FP_sim_10 <- bugs(data = second_FP, inits = second_inits, model.file = second_FP_model_FE, 
                         n.iter = 20000, parameters.to.save = "S", n.burnin = 10000, 
                         n.chains = 2, 
                         debug = FALSE)

results_second_FP_sim_10 <- second_FP_sim_10$summary
write.csv(results_second_FP_sim_10, here("02_data", "PFS_results_FP_sim10_108.csv"))

array_second_FP_sim_10 <- second_FP_sim_10$sims.array
array_second_FP_sim_10 <- array_second_FP_sim_10[, 1, time_point]
array_second_FP_sim_10 <- apply(array_second_FP_sim_10, 2, sample, size = 500)
write.csv(array_second_FP_sim_10, here("02_data", "PFS_array_FP_sim10.csv"))

# Second-order FP, FE model, P1 = -0.5, P2 = 0 ----------------
second_FP$P1 <- -0.5
second_FP$P2 <- 0
second_FP$maxt <- 108 # increasing maximum time to extrapolate survival
second_FP_sim_15 <- bugs(data = second_FP, inits = second_inits, model.file = second_FP_model_FE, 
                         n.iter = 20000, parameters.to.save = "S", n.burnin = 10000, 
                         n.chains = 2, 
                         debug = FALSE)

results_second_FP_sim_15 <- second_FP_sim_15$summary
write.csv(results_second_FP_sim_15, here("02_data", "PFS_results_FP_sim15_108.csv"))

array_second_FP_sim_15 <- second_FP_sim_15$sims.array
array_second_FP_sim_15 <- array_second_FP_sim_15[, 1, time_point]
array_second_FP_sim_15 <- apply(array_second_FP_sim_15, 2, sample, size = 500)
write.csv(array_second_FP_sim_15, here("02_data", "PFS_array_FP_sim15.csv"))

#8.2 ENZAMET as the trial of interest----------------------------------
# Second-order FP, FE model, P1 = -0.5, P2 = -0.5 ----------------
second_FP_model_FE_ENZA <- here("03_models", "second_FP_model_FE_ENZAMET.txt")
second_FP$P1 <- -0.5
second_FP$P2 <- -0.5
second_FP$maxt <- 108 # increasing maximum time to extrapolate survival
second_FP_sim_14_ENZA <- bugs(data = second_FP, inits = second_inits, model.file = second_FP_model_FE_ENZA, 
                         n.iter = 20000, parameters.to.save = "S", n.burnin = 10000, 
                         n.chains = 2, 
                         debug = FALSE)

results_second_FP_sim_14_ENZA <- second_FP_sim_14_ENZA$summary
write.csv(results_second_FP_sim_14_ENZA, here("02_data", "PFS_results_FP_sim14_ENZA.csv"))

# Second-order FP, FE model, P1 = -1, P2 = 0.5 ----------------
second_FP$P1 <- -1
second_FP$P2 <- 0.5
second_FP$maxt <- 108 # increasing maximum time to extrapolate survival
second_FP_sim_11_ENZA <- bugs(data = second_FP, inits = second_inits, model.file = second_FP_model_FE_ENZA, 
                         n.iter = 20000, parameters.to.save = "S", n.burnin = 10000, 
                         n.chains = 2, 
                         debug = FALSE)

results_second_FP_sim_11_ENZA <- second_FP_sim_11_ENZA$summary
write.csv(results_second_FP_sim_11_ENZA, here("02_data", "PFS_results_FP_sim11_ENZA.csv"))

# Second-order FP, FE model, P1 = 0, P2 = 0.5 ----------------
second_FP$P1 <- 0
second_FP$P2 <- 0.5
second_FP$maxt <- 108 # increasing maximum time to extrapolate survival
second_FP_sim_19_ENZA <- bugs(data = second_FP, inits = second_inits, model.file = second_FP_model_FE_ENZA, 
                         n.iter = 20000, parameters.to.save = "S", n.burnin = 10000, 
                         n.chains = 2, 
                         debug = FALSE)

results_second_FP_sim_19_ENZA <- second_FP_sim_19_ENZA$summary
write.csv(results_second_FP_sim_19_ENZA, here("02_data", "PFS_results_FP_sim19_ENZA.csv"))

# Second-order FP, FE model, P1 = -1, P2 = 0.5 ----------------
second_FP$P1 <- -1
second_FP$P2 <- 0
second_FP$maxt <- 108 # increasing maximum time to extrapolate survival
second_FP_sim_10_ENZA <- bugs(data = second_FP, inits = second_inits, model.file = second_FP_model_FE_ENZA, 
                         n.iter = 20000, parameters.to.save = "S", n.burnin = 10000, 
                         n.chains = 2, 
                         debug = FALSE)

results_second_FP_sim_10_ENZA <- second_FP_sim_10_ENZA$summary
write.csv(results_second_FP_sim_10_ENZA, here("02_data", "PFS_results_FP_sim10_ENZA.csv"))

# Second-order FP, FE model, P1 = -0.5, P2 = 0 ----------------
second_FP$P1 <- -0.5
second_FP$P2 <- 0
second_FP$maxt <- 108 # increasing maximum time to extrapolate survival
second_FP_sim_15_ENZA <- bugs(data = second_FP, inits = second_inits, model.file = second_FP_model_FE_ENZA, 
                         n.iter = 20000, parameters.to.save = "S", n.burnin = 10000, 
                         n.chains = 2, 
                         debug = FALSE)

results_second_FP_sim_15_ENZA <- second_FP_sim_15_ENZA$summary
write.csv(results_second_FP_sim_15_ENZA, here("02_data", "PFS_results_FP_sim15_ENZA.csv"))


#8.2 CHAARTED as the trial of interest----------------------------------
# Second-order FP, FE model, P1 = -0.5, P2 = -0.5 ----------------
second_FP_model_FE_CHAART <- here("03_models", "second_FP_model_FE_CHAARTED.txt")
second_FP$P1 <- -0.5
second_FP$P2 <- -0.5
second_FP$maxt <- 108 # increasing maximum time to extrapolate survival
second_FP_sim_14_CHAART <- bugs(data = second_FP, inits = second_inits, model.file = second_FP_model_FE_CHAART, 
                              n.iter = 20000, parameters.to.save = "S", n.burnin = 10000, 
                              n.chains = 2, 
                              debug = FALSE)

results_second_FP_sim_14_CHAART <- second_FP_sim_14_CHAART$summary
write.csv(results_second_FP_sim_14_CHAART, here("02_data", "PFS_results_FP_sim14_CHAART.csv"))

# Second-order FP, FE model, P1 = -1, P2 = 0.5 ----------------
second_FP$P1 <- -1
second_FP$P2 <- 0.5
second_FP$maxt <- 108 # increasing maximum time to extrapolate survival
second_FP_sim_11_CHAART <- bugs(data = second_FP, inits = second_inits, model.file = second_FP_model_FE_CHAART, 
                              n.iter = 20000, parameters.to.save = "S", n.burnin = 10000, 
                              n.chains = 2, 
                              debug = FALSE)

results_second_FP_sim_11_CHAART <- second_FP_sim_11_CHAART$summary
write.csv(results_second_FP_sim_11_CHAART, here("02_data", "PFS_results_FP_sim11_CHAART.csv"))

# Second-order FP, FE model, P1 = 0, P2 = 0.5 ----------------
second_FP$P1 <- 0
second_FP$P2 <- 0.5
second_FP$maxt <- 108 # increasing maximum time to extrapolate survival
second_FP_sim_19_CHAART <- bugs(data = second_FP, inits = second_inits, model.file = second_FP_model_FE_CHAART, 
                              n.iter = 20000, parameters.to.save = "S", n.burnin = 10000, 
                              n.chains = 2, 
                              debug = FALSE)

results_second_FP_sim_19_CHAART <- second_FP_sim_19_CHAART$summary
write.csv(results_second_FP_sim_19_CHAART, here("02_data", "PFS_results_FP_sim19_CHAART.csv"))

# Second-order FP, FE model, P1 = -1, P2 = 0.5 ----------------
second_FP$P1 <- -1
second_FP$P2 <- 0
second_FP$maxt <- 108 # increasing maximum time to extrapolate survival
second_FP_sim_10_CHAART <- bugs(data = second_FP, inits = second_inits, model.file = second_FP_model_FE_CHAART, 
                              n.iter = 20000, parameters.to.save = "S", n.burnin = 10000, 
                              n.chains = 2, 
                              debug = FALSE)

results_second_FP_sim_10_CHAART <- second_FP_sim_10_CHAART$summary
write.csv(results_second_FP_sim_10_CHAART, here("02_data", "PFS_results_FP_sim10_CHAART.csv"))

# Second-order FP, FE model, P1 = -0.5, P2 = 0 ----------------
second_FP$P1 <- -0.5
second_FP$P2 <- 0
second_FP$maxt <- 108 # increasing maximum time to extrapolate survival
second_FP_sim_15_CHAART <- bugs(data = second_FP, inits = second_inits, model.file = second_FP_model_FE_CHAART, 
                              n.iter = 20000, parameters.to.save = "S", n.burnin = 10000, 
                              n.chains = 2, 
                              debug = FALSE)

results_second_FP_sim_15_CHAART <- second_FP_sim_15_CHAART$summary
write.csv(results_second_FP_sim_15_CHAART, here("02_data", "PFS_results_FP_sim15_CHAART.csv"))

#8.3 LATITUDE as the trial of interest----------------------------------
# Second-order FP, FE model, P1 = -0.5, P2 = -0.5 ----------------
second_FP_model_FE_LATITUDE <- here("03_models", "second_FP_model_FE_LATITUDE.txt")
second_FP$P1 <- -0.5
second_FP$P2 <- -0.5
second_FP$maxt <- 108 # increasing maximum time to extrapolate survival
second_FP_sim_14_LATITUDE <- bugs(data = second_FP, inits = second_inits, model.file = second_FP_model_FE_LATITUDE, 
                                n.iter = 20000, parameters.to.save = "S", n.burnin = 10000, 
                                n.chains = 2, 
                                debug = FALSE)

results_second_FP_sim_14_LATITUDE <- second_FP_sim_14_LATITUDE$summary
write.csv(results_second_FP_sim_14_LATITUDE, here("02_data", "PFS_results_FP_sim14_LATITUDE.csv"))

# Second-order FP, FE model, P1 = -1, P2 = 0.5 ----------------
second_FP$P1 <- -1
second_FP$P2 <- 0.5
second_FP$maxt <- 108 # increasing maximum time to extrapolate survival
second_FP_sim_11_LATITUDE <- bugs(data = second_FP, inits = second_inits, model.file = second_FP_model_FE_LATITUDE, 
                                n.iter = 20000, parameters.to.save = "S", n.burnin = 10000, 
                                n.chains = 2, 
                                debug = FALSE)

results_second_FP_sim_11_LATITUDE <- second_FP_sim_11_LATITUDE$summary
write.csv(results_second_FP_sim_11_LATITUDE, here("02_data", "PFS_results_FP_sim11_LATITUDE.csv"))

# Second-order FP, FE model, P1 = 0, P2 = 0.5 ----------------
second_FP$P1 <- 0
second_FP$P2 <- 0.5
second_FP$maxt <- 108 # increasing maximum time to extrapolate survival
second_FP_sim_19_LATITUDE <- bugs(data = second_FP, inits = second_inits, model.file = second_FP_model_FE_LATITUDE, 
                                n.iter = 20000, parameters.to.save = "S", n.burnin = 10000, 
                                n.chains = 2, 
                                debug = FALSE)

results_second_FP_sim_19_LATITUDE <- second_FP_sim_19_LATITUDE$summary
write.csv(results_second_FP_sim_19_LATITUDE, here("02_data", "PFS_results_FP_sim19_LATITUDE.csv"))

# Second-order FP, FE model, P1 = -1, P2 = 0.5 ----------------
second_FP$P1 <- -1
second_FP$P2 <- 0
second_FP$maxt <- 108 # increasing maximum time to extrapolate survival
second_FP_sim_10_LATITUDE <- bugs(data = second_FP, inits = second_inits, model.file = second_FP_model_FE_LATITUDE, 
                                n.iter = 20000, parameters.to.save = "S", n.burnin = 10000, 
                                n.chains = 2, 
                                debug = FALSE)

results_second_FP_sim_10_LATITUDE <- second_FP_sim_10_LATITUDE$summary
write.csv(results_second_FP_sim_10_LATITUDE, here("02_data", "PFS_results_FP_sim10_LATITUDE.csv"))

# Second-order FP, FE model, P1 = -0.5, P2 = 0 ----------------
second_FP$P1 <- -0.5
second_FP$P2 <- 0
second_FP$maxt <- 108 # increasing maximum time to extrapolate survival
second_FP_sim_15_LATITUDE <- bugs(data = second_FP, inits = second_inits, model.file = second_FP_model_FE_LATITUDE, 
                                n.iter = 20000, parameters.to.save = "S", n.burnin = 10000, 
                                n.chains = 2, 
                                debug = FALSE)

results_second_FP_sim_15_LATITUDE <- second_FP_sim_15_LATITUDE$summary
write.csv(results_second_FP_sim_15_LATITUDE, here("02_data", "PFS_results_FP_sim15_LATITUDE.csv"))


#9. Export HRs for timepoints of interest--------------------------------
HR_params <- c("HR_2year", "HR_3year", "HR_4year")
# Second-order FP, FE model, P1 = -0.5, P2 = -0.5 ----------------
second_FP$P1 <- -0.5
second_FP$P2 <- -0.5
#second_FP$maxt <- 48 # pending to confirm this choice. At the moment ARASENS follow-up
second_FP_sim_HR <- bugs(data = second_FP, inits = second_inits, model.file = second_FP_model_FE, 
                         n.iter = 20000, parameters.to.save = HR_params, n.burnin = 10000, 
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


HRs <- cbind(HR_24, HR_36, HR_48)
write_csv(as_tibble(HRs), here("05_tables", "time-varying_HRs_PFS_v5.csv"))

#10. Code to extract the RMST at 108 months-----------------------------
# Second-order FP, FE model, P1 = 0, P2 = 0.5 (model with the lowest DIC)
second_FP$P1 <- -0.5
second_FP$P2 <- -0.5
second_FP_sim_14 <- bugs(data = second_FP, inits = second_inits, model.file = second_FP_model_FE, 
                         n.iter = 20000, parameters.to.save = "RMST", n.burnin = 10000, 
                         n.chains = 2, 
                         debug = FALSE)

# Absolute RSMT
rmst_daro    <- summary.stats(second_FP_sim_14$sims.array[,1,"RMST[1,108]"])
rmst_doc     <- summary.stats(second_FP_sim_14$sims.array[,1,"RMST[2,108]"])
rmst_enza    <- summary.stats(second_FP_sim_14$sims.array[,1,"RMST[3,108]"])
rmst_adt     <- summary.stats(second_FP_sim_14$sims.array[,1,"RMST[4,108]"])
rmst_abi     <- summary.stats(second_FP_sim_14$sims.array[,1,"RMST[5,108]"])
rmst_apa     <- summary.stats(second_FP_sim_14$sims.array[,1,"RMST[6,108]"])
rmst_enzadoc <- summary.stats(second_FP_sim_14$sims.array[,1,"RMST[7,108]"])
rmst_abidoc  <- summary.stats(second_FP_sim_14$sims.array[,1,"RMST[8,108]"])

rmst <- c(rmst_daro, rmst_doc, rmst_enza, rmst_adt, rmst_abi, rmst_apa,
          rmst_enzadoc, rmst_abidoc)

# Relative RSMT to reference treatment
diff_doc <- summary.stats(second_FP_sim_14$sims.array[,1,"RMST[2,108]"] - 
                          second_FP_sim_14$sims.array[,1,"RMST[1,108]"])
diff_enza <- summary.stats(second_FP_sim_14$sims.array[,1,"RMST[3,108]"] - 
                           second_FP_sim_14$sims.array[,1,"RMST[1,108]"])
diff_adt <- summary.stats(second_FP_sim_14$sims.array[,1,"RMST[4,108]"] - 
                          second_FP_sim_14$sims.array[,1,"RMST[1,108]"])
diff_abi <- summary.stats(second_FP_sim_14$sims.array[,1,"RMST[5,108]"] - 
                          second_FP_sim_14$sims.array[,1,"RMST[1,108]"])
diff_apa <- summary.stats(second_FP_sim_14$sims.array[,1,"RMST[6,108]"] - 
                          second_FP_sim_14$sims.array[,1,"RMST[1,108]"])
diff_enzadoc <- summary.stats(second_FP_sim_14$sims.array[,1,"RMST[7,108]"] - 
                              second_FP_sim_14$sims.array[,1,"RMST[1,108]"])
diff_abidoc <- summary.stats(second_FP_sim_14$sims.array[,1,"RMST[8,108]"] - 
                             second_FP_sim_14$sims.array[,1,"RMST[1,108]"])

diff <- c(NA, diff_doc, diff_enza, diff_adt, diff_abi, diff_apa,
          diff_enzadoc, diff_abidoc)

# Export results
write.csv(bind_cols(rmst, diff), here("05_tables", "rmst_pfs.csv"))
