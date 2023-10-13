# Running a random effects, Normal likelihood, identity link NMA in R2OpenBUGS

# Load the R2OpenBUGS package
library(R2OpenBUGS)
library(readxl)
library(here)

# Load the data 
mHSPC_OS_data_low_volume <- read_excel(here("Subgroups/data", "mHSPC OS data low_volume.xlsx"))

# Normal likelihood, identity link, random effects
model_normal_identity_re <- function()
{
  for(i in 1:ns){ # LOOP THROUGH STUDIES
    #w[i,1] <- 0 # adjustment for multi-arm trials is zero for control arm
    delta[i,1] <- 0 # treatment effect is zero for control arm
    # HT: Main change was to remove the mu[] since theta[i] is a LHR and not log hazard
    # HT: The loop over k is redundant since no k index is used in lines 16-20
    #for (k in 1:2) { # LOOP THROUGH ARMS
      var[i] <- pow(se[i,2],2) # calculate variances
      prec[i] <- 1/var[i] # set precisions
      y[i,2] ~ dnorm(theta[i],prec[i]) # normal likelihood
      theta[i] <-  delta[i,2] # model for linear predictor
      dev[i] <- (y[i,2]-theta[i])*(y[i,2]-theta[i])*prec[i] #Deviance contribution
    #}
    # HT: As only one arm contributes, you don't need to sum over na[i]
    #resdev[i] <- sum(dev[i,1:na[i]]) # summed residual deviance contribution for this trial
    # HT: Again loop over k is redundant since code below should only be used for the
    # second arm. Note that the multi-arm correction is only needed if you are including
    # any trials with >=3 arms
    #for (k in 2:na[i]) { # LOOP THROUGH ARMS
    # HT: Your taud on line below was missing a second index
      delta[i,2] ~ dnorm(md[i],taud[i, 2]) # trial-specific LOR distributions
      md[i] <- d[t[i,2]] - d[t[i,1]] # mean of treat effects distributions (with multi-arm trial correction)
      # HT: Note that below simplifies to taud[i,2] = tau, as expected
      taud[i,2] <- tau *2*(2-1)/2 # precision of treat effects distributions (with multi-arm trial correction)
      #w[i,k] <- (delta[i,k] - d[t[i,k]] + d[t[i,1]]) # adjustment for multi-arm RCTs
      #sw[i,k] <- sum(w[i,1:k-1])/(k-1) # cumulative adjustment for multi-arm trials
    #}
  }
  # HT: Simplified as dev[i] is total residual deviance for each arm and resdev is redundant
  totresdev <- sum(dev[]) #Total Residual Deviance
  d[1]<-0 # treatment effect is zero for reference treatment
  for (k in 2:nt){ d[k] ~ dnorm(0,.0001) } # vague priors for treatment effects
  #sd ~ dunif(0,5) # vague prior for between-trial SD.
  #tau <- pow(sd,-2) # between-trial precision = (1/between-trial variance)
  # HT: I had to rename var to het.var as var[] is already used for trial
  # specific variance
  het.var.prec <- pow(1.41, -2) # code for 1/(1.41*1.41)
    het.var ~ dlnorm(-4.18, het.var.prec) #lognormal distribution
    tau <- pow(het.var, -1)
    
    for (k in 1:nt) {
      rk[k] <- nt+1-rank(d[],k) # assumes events are "good"
      #rk[k] <- rank(d[],k) # assumes events are "bad"
      best[k] <- equals(rk[k],1) #calculate probability that treat k is best
      for (h in 1:nt){ prob[h,k] <- equals(rk[k],h) } # calculates probability that treat k is h-th best
    }
}


model_normal_identity_re_ht <- function()
{
  for(i in 1:ns){ # LOOP THROUGH STUDIES
   
    var[i] <- pow(se[i,2],2) # calculate variances
    prec[i] <- 1/var[i] # set precisions
    y[i,2] ~ dnorm(theta[i],prec[i]) # normal likelihood
    theta[i] <- delta[i] # model for linear predictor
    dev[i] <- (y[i,2]-theta[i])*(y[i,2]-theta[i])*prec[i] #Deviance contribution
    
    delta[i] ~ dnorm(md[i], tau) # trial-specific LOR distributions
    md[i] <- d[t[i,2]] - d[t[i,1]] # mean of treat effects distributions (with multi-arm trial correction)
  }
  # HT: Simplified as dev[i] is total residual deviance for each arm and resdev is redundant
  totresdev <- sum(dev[]) #Total Residual Deviance
  d[1]<-0 # treatment effect is zero for reference treatment
  for (k in 2:nt){ d[k] ~ dnorm(0,.0001) } # vague priors for treatment effects
  #sd ~ dunif(0,5) # vague prior for between-trial SD.
  #tau <- pow(sd,-2) # between-trial precision = (1/between-trial variance)
  # HT: I had to rename var to het.var as var[] is already used for trial
  # specific variance
  het.var.prec <- pow(1.41, -2) # code for 1/(1.41*1.41)
  het.var ~ dlnorm(-4.18, het.var.prec) #lognormal distribution
  tau <- pow(het.var, -1)
  
  # ranking on relative scale                  #EK added
  for (k in 1:nt) {
    rk[k] <- nt+1-rank(d[],k) # assumes events are "good"
    #rk[k] <- rank(d[],k) # assumes events are "bad"
    best[k] <- equals(rk[k],1) #calculate probability that treat k is best
    for (h in 1:nt){ prob[h,k] <- equals(rk[k],h) } # calculates probability that treat k is h-th best
  }
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

calculate_sucra <- function(bugs_object, bugs_data, t_names)
{
  # Ranking probability indices
  prob.indices <- grep("prob", rownames(bugs_object$summary))
  
  # Probability that each treatment takes each rank plus cumulative
  rank_probs <- rank_cumprobs <- matrix(NA, nrow = bugs_data$nt, ncol = bugs_data$nt)
  # SUCRA	
  sucra <- rep(NA, bugs_data$nt)
  names(sucra) <- rownames(rank_cumprobs) <- rownames(rank_probs) <- t_names
  colnames(rank_cumprobs) <- colnames(rank_probs) <- paste("Rank", c(1:bugs_data$nt))
  for(i in 1:bugs_data$nt)
  {
    # prob[k, i] is probability treatment i is kth best
    rank_probs[i, ] <- bugs_object$summary[prob.indices, "mean"][c(0:(bugs_data$nt - 1))*bugs_data$nt+i]
    rank_cumprobs[i, ] <- cumsum(rank_probs[i, ])
    sucra[i] <- (1/(bugs_data$nt - 1))*sum(rank_cumprobs[i, 1:bugs_data$nt - 1])
  }
  
  return("sucra" = sucra)
}


# Number of MCMC chains and samples
n_chains <- 3
num_sims <- 50000 * n_chains 
burn_in <- 50000 * n_chains	

# HT: I needed to add below to actually load the data.
# Please ensure you load the Excel data each time, so that you aren't using
# an older version of the data!
# HT: I set working directory to source file location and put data in same directory
#mHSPC_OS_nice_network <- read_xlsx("mHSPC OS_nice network.xlsx")

# also, to get the correct number of dimensions is good to use a "comparator" arm with 0 for the lhr and the se
ns <- nrow(mHSPC_OS_data_low_volume)
t  <- array(c(mHSPC_OS_data_low_volume$t1, mHSPC_OS_data_low_volume$t2), dim = c(ns, 2)) 
# HT: Used max(t) here to be more general
nt <- max(t)
y  <- array(c(rep(0, ns), mHSPC_OS_data_low_volume$y), dim = c(ns, 2))
se <- array(c(rep(0, ns), mHSPC_OS_data_low_volume$se), dim = c(ns, 2))


# Load the data
#trial_data <- read.csv("Practical1_exercise1_data.csv")
# Store BUGS data in matrices
# Use as.matrix() to ensure data in correct format for bugs() to recognise
#y <- as.matrix(trial_data[, c("y1", "y2", "y3")])
#se <- as.matrix(trial_data[, c("se1", "se2", "se3")])
# Call the treatment matrix tr to avoid overwriting the transpose function t()
#tr <- as.matrix(trial_data[, c("t1", "t2", "t3")])
#na <- as.vector(trial_data[, "na"])

bugs_data <- list(
  y = y,
  se = se,
  t = t,
  #na = na, not needed in this model
  ns = ns, 
  nt = nt)


# Define the bugs data
#bugs_data <- list("y" = y,
                  #"se" = se,
                  #"t" = tr,
                  #"na" = na,
                  #"ns" = 7, "nt" = 5)

# Create initial values for MCMC simulation
# HT: I changed the sd to het.var as sd is not used in the model
# Note that I used het.var rather than tau, since the het.var are the random
# parameters in the BUGS model (i.e. defined by dlnorm)
# HT: I removed mu as no longer in model
# HT: I also changed d to automatically depend on number of treatments
inits1 <- list(d=c( NA, rep(0, nt - 1)), het.var= 1)
inits2 <- list(d=c( NA, rep(-1, nt - 1)), het.var = 2)
inits3 <- list(d=c( NA, rep(2, nt - 1)), het.var = 0.5)
bugs_inits <- list(inits1, inits2, inits3)

# Call OpenBUGS
bugs_object_re<-bugs(data = bugs_data, inits = bugs_inits,
                     parameters.to.save = c("d", "totresdev", "rk", "best", "prob"),
                     model = model_normal_identity_re, clearWD = TRUE, 
                     summary.only = FALSE,
                     n.iter = (num_sims + burn_in), n.burnin = burn_in,
                     n.chains = n_chains, bugs.seed = 1 ,debug = TRUE)

# Look at the raw bugs_object
bugs_object_re$DIC
bugs_object_re$summary

# Format the odds ratios
cross_meandiff_re_lowvol <- cross_effect(bugs_object = bugs_object_re, t_names = c("daro+doce+adt", "doce+adt", "enza+adt", "adt", "apa+adt", "abi+doce+adt", "abi+adt", "enza+doce+adt"), med = TRUE, exp = TRUE)
write.csv(x = cross_meandiff_re_lowvol, file = "cross_meandiff_re_lowvol.csv")

calculate_sucra(bugs_object = bugs_object_re, bugs_data = bugs_data, t_names = c("daro+doce+adt", "doce+adt", "enza+adt", "adt", "apa+adt", "abi+doce+adt", "abi+adt", "enza+doce+adt"))





