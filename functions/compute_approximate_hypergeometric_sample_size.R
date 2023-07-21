# Compute approximate sample size under hypergeometric sampling for a given 
  # sensitivity, design prevalence and consignment size. 
  # Uses the approximation formula from (Lane et al., 2018) to avoid the 
  # solution to bounce up and down due to integer rounding.
  
  # S: Sensitivity of the inspection (also called confidence level)
  # p: design prevalence (also called level of detection)
  # N: consignment size
  compute_approximate_hypergeometric_sample_size <- function(S, p, N){
    hypergeometric_n <- (1 - (1 - S)^(1 / (N * p))) * (N - (N * p - 1) / 2)
    return(hypergeometric_n)
  }