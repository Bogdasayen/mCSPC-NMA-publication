# Philip proportional hazards NMA script
# R2OpenBUGS 
# Fixed effects, Normal likelihood, identity link NMA in R2OpenBUGS

# Load the R2OpenBUGS package
library(R2OpenBUGS)
library(readxl)
library(here)

# Load the data 
mHSPC_OS_data_high_risk <- read_excel(here("Subgroups/data", "mHSPC OS data high risk.xlsx"))

# Normal likelihood, identity link, fixed effects
model_normal_identity_fe <- function()
{
  for(i in 1:ns){ # LOOP THROUGH STUDIES
    mu[i] ~ dnorm(0,.0001) # vague priors for all trial baselines
      var[i] <- pow(se[i, 2],2) # calculate variances 
      prec[i] <- 1/var[i] # set precisions
      y[i, 2] ~ dnorm(theta[i],prec[i]) # normal likelihood 
      theta[i] <- d[t[i,2]]-d[t[i,1]] # model for linear predictor 
      
      dev[i] <- (y[i,2]-theta[i])*(y[i,2]-theta[i])*prec[i] #Deviance contribution
    }
  
  totresdev <- sum(dev[]) 
  
  d[1]<-0 # treatment effect is zero for reference treatment
  for (k in 2:nt){ d[k] ~ dnorm(0,.0001) } # vague priors for treatment effects
  
  # ranking on relative scale                  
  for (k in 1:nt) {
    rk[k] <- nt+1-rank(d[],k) # assumes events are "good"
    #rk[k] <- rank(d[],k) # assumes events are "bad"
    best[k] <- equals(rk[k],1) #calculate probability that treat k is best
    for (h in 1:nt){ prob[h,k] <- equals(rk[k],h) } # calculates probability that treat k is h-th best
  }
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
num_sims <- 10000 * n_chains 
burn_in <- 10000 * n_chains	

# Define the bugs data 
# also, to get the correct number of dimensions is good to use a "comparator" arm with 0 for the lhr and the se
ns <- nrow(mHSPC_OS_data_high_risk)
nt <- max(mHSPC_OS_data_high_risk$t1)
t  <- array(c(mHSPC_OS_data_high_risk$t1, mHSPC_OS_data_high_risk$t2), dim = c(ns, 2)) 
y  <- array(c(rep(0, ns), mHSPC_OS_data_high_risk$y), dim = c(ns, 2))
se <- array(c(rep(0, ns), mHSPC_OS_data_high_risk$se), dim = c(ns, 2))

bugs_data <- list(
  y = y,
  se = se,
  t = t,
  ns = ns, 
  nt = nt)

# Create initial values for MCMC simulation DA: always check you have the right number of 
# initial values according to the number of parameters
# Also I've seen wrapping the lists in a function help BUGS to converge

bugs_inits <- function(){
  #chain 1
  list(d = c(NA, rep(0, nt - 1)))
  #chain 2
  list(d = c(NA, rep(-1, nt - 1)))
  #chain 3
  list(d = c(NA, rep(2, nt - 1)))
}

# Call OpenBUGS
bugs_object_fe<-bugs(data = bugs_data, inits = bugs_inits,
                     parameters.to.save = c("d", "totresdev", "rk", "best", "prob"),
                     model = model_normal_identity_fe, clearWD = TRUE, 
                     summary.only = FALSE,
                     n.iter = (num_sims + burn_in), n.burnin = burn_in,
                     n.chains = n_chains, bugs.seed = 1 ,debug = TRUE)

# Look at the raw bugs_object
bugs_object_fe$DIC
bugs_object_fe$summary

# Format the odds ratios
cross_meandiff_fe_OS_highrisk <- cross_effect(bugs_object = bugs_object_fe, t_names = c("daro+doce+adt", "doce+adt", "adt", "abi+adt", "apa+adt"), med = TRUE, exp = TRUE)
write.csv(x = cross_meandiff_fe_OS_highrisk, file = "cross_meandiff_fe_OS_highrisk.csv")

calculate_sucra(bugs_object = bugs_object_fe, bugs_data = bugs_data, t_names = c("daro+doce+adt", "doce+adt", "adt", "abi+adt", "apa+adt"))






