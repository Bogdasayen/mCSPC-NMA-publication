# Utils functions

# Function to extract fit statistics from a BUGS object
extract_fit <- function(bugs_fit){
  result_fit <- c(bugs_fit$mean$deviance, bugs_fit$pD, bugs_fit$DIC, bugs_fit$mean$totresdev)
  return(result_fit)
}

# Function to format results in CI style
format_results <- function(x, n_digits = 3) {
  return(paste0(format(x[1], digits = n_digits, nsmall = n_digits),
                " (", format(x[2], digits = n_digits, nsmall = n_digits),
                ", ", format(x[3], digits = n_digits, nsmall = n_digits), ")"))
}

# Function to extract summary statistics using CI style

summary.stats<-function(x,n.digits=2,med=FALSE)
{
  if(med){
    return(paste0(round(median(x),digits=n.digits)," (",round(quantile(x,probs=c(0.025)),digits=n.digits),", ",round(quantile(x,probs=c(0.975)),digits=n.digits),")"))
  }else{
    return(paste0(round(mean(x),digits=n.digits)," (",round(quantile(x,probs=c(0.025)),digits=n.digits),", ",round(quantile(x,probs=c(0.975)),digits=n.digits),")"))
  }
}