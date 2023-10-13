# Based on P Orishaba proportional hazards NMA script for FE OS
# Adapted to include correlation between STAMPEDE hazard ratios
# Howard Thom 30-Dec-2023

# Load the R2OpenBUGS package
library(R2OpenBUGS)
library(readxl)
library(here)

# Load the data 
mHSPC_PFS_age <- read_excel(here("Meta-regressions/data", "mHSPC PFS age.xlsx"))


# Random effects model
model_normal_identity_re <- function()
{
  for(i in 1:ns){ # LOOP THROUGH STUDIES
    var[i] <- pow(se[i,2],2) # calculate variances
    prec[i] <- 1/var[i] # set precisions
    y[i,2] ~ dnorm(theta[i],prec[i]) # normal likelihood
    theta[i] <- delta[i] + (beta[t[i,2]]-beta[t[i,1]]) * (x[i]-mx)# model for linear predictor

    dev[i] <- (y[i,2]-theta[i])*(y[i,2]-theta[i])*prec[i] #Deviance contribution
    delta[i] ~ dnorm(md[i], tau) # trial-specific LOR distributions
    md[i] <- d[t[i,2]] - d[t[i,1]] # mean of treat effects distributions (with multi-arm trial correction)
  }
  totresdev <- sum(dev[]) #Total Residual Deviance
  d[1]<-0 # treatment effect is zero for reference treatment
  beta[1] <- 0 # covariate effect is zero for reference treatment
  for (k in 2:nt){ d[k] ~ dnorm(0,.0001)  # vague priors for treatment effects
    beta[k] <- B # common covariate effect
  }
  B ~ dnorm(0,.0001) # vague prior for covariate effect
  
  exp_beta <- exp(B) # exponential of beta
  

  # Estimate effects relative to treatment of interest
  for(k in 1:nt) {d_interest[k] <- d[k] - d[t_of_interest]}

  # Informative prior on heterogeneity variance
  het.var.prec <- pow(1.41, -2) # code for 1/(1.41*1.41)
  het.var ~ dlnorm(-4.18, het.var.prec) #lognormal distribution
  tau <- pow(het.var, -1)
}

model_mams_adjust_re <- function()
{
  # Non-MAMS designs analysed as usual
  for(i in 1:ns){ # LOOP THROUGH STUDIES
    
    var[i] <- pow(se[i,2],2) # calculate variances
    prec[i] <- 1/var[i] # set precisions
    y[i,2] ~ dnorm(theta[i],prec[i]) # normal likelihood
    theta[i] <- delta[i] + (beta[t[i,2]]-beta[t[i,1]]) * (x[i]-mx)# model for linear predictor
    dev[i] <- (y[i,2]-theta[i])*(y[i,2]-theta[i])*prec[i] #Deviance contribution
    
    delta[i] ~ dnorm(md[i], tau) # trial-specific LOR distributions
    md[i] <- d[t[i,2]] - d[t[i,1]] # mean of treat effects distributions (with multi-arm trial correction)
  }
  
  
  # MAMS designs
  y_mams[1:na_mams] ~ dmnorm(theta_mams[], prec_mams[ , ])
  for(k in 1:na_mams) {
    theta_mams[k] <- delta_mams[k] + (beta[t_mams[k,2]]-beta[t_mams[k,1]]) * (x_mams[k]-mx)
    delta_mams[k] ~ dnorm(md_mams[k], tau)
      md_mams[k] <- d[t_mams[k,2]]-d[t_mams[k,1]]
    dev_mams[k] <-
      (y_mams[k]-theta_mams[k]) * (y_mams[k]-theta_mams[k]) * prec_mams[k, k] #Deviance contribution
  }
  
  # Total residual deviance is sum of non-MAMS and MAMS contributions
  totresdev <- sum(dev[]) + sum(dev_mams[])
  
  # Treatment model same as in unadjusted analysis
  d[1]<-0 # treatment effect is zero for reference treatment
  beta[1] <- 0 # covariate effect is zero for reference treatment
  for (k in 2:nt){ d[k] ~ dnorm(0,.0001)  # vague priors for treatment effects
    beta[k] <- B # common covariate effect
  }
  B ~ dnorm(0,.0001) # vague prior for covariate effect
  
  exp_beta <- exp(B) # exponential of beta
  
  
  # Estimate effects relative to treatment of interest
  for(k in 1:nt) {d_interest[k] <- d[k] - d[t_of_interest]}

  # Informative prior on heterogeneity variance
  het.var.prec <- pow(1.41, -2) # code for 1/(1.41*1.41)
  het.var ~ dlnorm(-4.18, het.var.prec) #lognormal distribution
  tau <- pow(het.var, -1)
}

# Function to generate summary statistics using coda from MCMC simulations
summary_stats <- function(x, n_digits = 2, med = FALSE)
{
  if(med){
    return(paste0(format(median(x), digits = n_digits, nsmall = n_digits), " (", 
                  format(quantile(x, probs = c(0.025)), digits = n_digits, nsmall = n_digits), ", ", 
                  format(quantile(x, probs = c(0.975)), digits = n_digits, nsmall = n_digits), ")"))
  }else{
    return(paste0(format(mean(x), digits = n_digits, nsmall = n_digits)," (",
                  format(quantile(x, probs = c(0.025)), digits = n_digits, nsmall = n_digits),", ",
                  format(quantile(x, probs = c(0.975)), digits = n_digits, nsmall = n_digits),")"))
  }
}

# Function to generate cross relative treatment effects comparison table
# using the d (e.g. log odds, mean differences) from R2OpenBUGS NMA models
cross_effect <- function(bugs_object, t_names, med = FALSE, exp = TRUE)
{
  effect_table <- matrix(NA, nrow = length(t_names), ncol = length(t_names))
  rownames(effect_table) <- colnames(effect_table) <- t_names
  for(i in 2:length(t_names))
  {
    # If the results need to be exponentiated (e.g. for odds ratios or hazard ratios)
    if(exp) {
      # Comparisons with reference
      effect_table[i, 1] <- summary_stats(exp(bugs_object$sims.array[, , 
                                                                     grep("d", rownames(bugs_object$summary))[i - 1]]), med = med)
      effect_table[1, i] <- summary_stats(exp(-bugs_object$sims.array[, , 
                                                                      grep("d", rownames(bugs_object$summary))[i - 1]]), med = med)
      for(j in 2:length(t_names))
      {
        effect_table[i, j] <- summary_stats(exp(
          bugs_object$sims.array[, , grep("d", rownames(bugs_object$summary))[i - 1]]-
            bugs_object$sims.array[, , grep("d", rownames(bugs_object$summary))[j - 1]]
        ), med = med)			
      }  
    } else {
      # If results do not need to be exponentiated (e.g. for mean differences)
      effect_table[i, 1] <- summary_stats(bugs_object$sims.array[, , 
                                                                 grep("d", rownames(bugs_object$summary))[i - 1]], med = med)
      effect_table[1, i] <- summary_stats(-bugs_object$sims.array[, , 
                                                                  grep("d", rownames(bugs_object$summary))[i - 1]], med = med)
      for(j in 2:length(t_names))
      {
        effect_table[i, j] <- summary_stats(
          bugs_object$sims.array[, , grep("d", rownames(bugs_object$summary))[i - 1]]-
            bugs_object$sims.array[, , grep("d", rownames(bugs_object$summary))[j - 1]]
          , med = med)			
      }
    }
    
  }
  for(i in 1:length(t_names))effect_table[i, i] <- t_names[i]
  
  return(effect_table)
}


# Number of MCMC chains and samples
n_chains <- 3
num_sims <- 30000 * n_chains 
burn_in <- 30000 * n_chains	

# Define the bugs data 
# also, to get the correct number of dimensions is good to use a "comparator" arm with 0 for the lhr and the se
ns <- nrow(mHSPC_PFS_age)
t  <- array(c(mHSPC_PFS_age$t1, mHSPC_PFS_age$t2), dim = c(ns, 2)) 
nt <- max(t) 
y  <- array(c(rep(0, ns), mHSPC_PFS_age$y), dim = c(ns, 2))
se <- array(c(rep(0, ns), mHSPC_PFS_age$se), dim = c(ns, 2))
x <- mHSPC_PFS_age$x

# Set reference treatment to ADT, which is most common treatment
# This is treatment 4 in Excel sheet
# Temporary code for ADT
t[t == 4] <- 9999
# Temporary code for darolutamide
t[t == 1] <- 8888
# Set reference to ADT
t[t == 9999] <- 1
t[t == 8888] <- 4

study_names <- gsub("#", "", mHSPC_PFS_age$`#ID`)
rownames(t) <- rownames(y) <- rownames(se) <- study_names


# Separate out the STAMPEDE/MAMS studies
mams_indices <- grep("STAMPEDE", rownames(y))
# Data for the MAMS trials
y_mams <- y[mams_indices, 2]
se_mams <- se[mams_indices, ] # Not used as SE may not be correct and need covariance
x_mams <- x[mams_indices]
t_mams <- t[mams_indices, ]
na_mams <- length(y_mams)
# Data for the non-MAMS trials (analysed as usual)
y_non_mams <- y[-mams_indices, ]
se_non_mams <- se[-mams_indices, ]
x_non_mams <- x[-mams_indices]
t_non_mams <- t[-mams_indices, ] 
ns_non_mams <- dim(y_non_mams)[1]

# Covariance matrix for the MAMS trial
# From Section 4.3 of the SAP (FFS/PFS can be found there as well)
var_mams <- matrix(c(0.00810098, 0.00296, 0.00342,
                     0.00296, 0.00538584, 0.00131,
                     0.00342, 0.00131, 0.01905604),
                   nrow = 3)
# Inverse of covariance matrix is precision
prec_mams <- solve(var_mams)


# Bugs data for unadjusted model
bugs_data <- list(
  y = y,
  se = se,
  t = t,
  ns = ns, 
  nt = nt,
  x = x,
  mx = mean(x), # Average x
  t_of_interest = 4 # Darolutamide
  )

# Bugs data for adjusted model
bugs_data_mams <- list(
  y = y_non_mams,
  se = se_non_mams,
  x = x_non_mams,
  t = t_non_mams,
  ns = ns_non_mams, 
  y_mams = y_mams,
  prec_mams = prec_mams,
  t_mams = t_mams,
  x_mams = x_mams,
  na_mams = na_mams,
  nt = nt,
  mx = mean(x), # Average x (x is both x_non_mams and x_mams)
  t_of_interest = 4 # Darolutamide
  )

# Create initial values for MCMC simulation 
# initial values according to the number of parameters
# These are the same for both adjusted and unadjusted models
inits1 <- list(d=c( NA, rep(0, nt - 1)), het.var= 1, B = -1)
inits2 <- list(d=c( NA, rep(-1, nt - 1)), het.var = 2, B = 2)
inits3 <- list(d=c( NA, rep(2, nt - 1)), het.var = 0.5, B = 0.5)
bugs_inits <- list(inits1, inits2, inits3)

# Call OpenBUGS

bugs_object_re <- bugs(data = bugs_data, inits = bugs_inits,
                       parameters.to.save = c("d", "B", "exp_beta", "totresdev", "d_interest"),
                       model = model_normal_identity_re, clearWD = TRUE, 
                       summary.only = FALSE,
                       n.iter = (num_sims + burn_in), n.burnin = burn_in,
                       n.chains = n_chains, bugs.seed = 1 ,debug = TRUE)

bugs_object_re_mams <- bugs(data = bugs_data_mams, inits = bugs_inits,
                            parameters.to.save = c("d", "B", "exp_beta", "totresdev", "d_interest"),
                            model = model_mams_adjust_re, clearWD = TRUE, 
                            summary.only = FALSE,
                            n.iter = (num_sims + burn_in), n.burnin = burn_in,
                            n.chains = n_chains, bugs.seed = 1 ,debug = TRUE)


# Look at the raw bugs_object
bugs_object_re$DIC
# There is a difference but it is not huge
# Not that totresdev and DIC can't be compared as data are different
bugs_object_re$summary[grep("d_interest", rownames(bugs_object_re$summary)), c("mean", "2.5%", "97.5%")]
bugs_object_re_mams$summary[grep("d_interest", rownames(bugs_object_re$summary)), c("mean", "2.5%", "97.5%")]

# Format the odds ratios
cross_meandiff_re_nmr_pfs_age <- cross_effect(bugs_object = bugs_object_re_mams, t_names = c(1:8), med = TRUE, exp = TRUE)
write.csv(x = cross_meandiff_re_nmr_pfs_age, file = "cross_meandiff_re_nmr_pfs_age.csv")






