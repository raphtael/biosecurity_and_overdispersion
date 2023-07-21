# Analytical formula (Eq. 3 in manuscript) for mean propagule pressure per 
  # consignment under simple random sampling when the distribution of 
  # infestation rate among consignment is zero-inflated beta-binomial
  
  # a: alpha parameter of the Beta distribution for the pathway from which the consigenemnt comes from.
  # b: beta parameter of the Beta distribution for the pathway from which the consigenemnt comes from.
  # z: zero-inflated parameter. The proportion of consignments with infestation rate zero in the zero-inflated distribution. Default value of zero (i.e, default to the beta-binomial distribution)
  # N: vector of the number of units per consignement. One value for each consignment
  # n: inpsection sample size
  compute_propagule_pressure_per_cons <- function(a, b, n, N, z = 0){
    propagule_pressure <- 
      exp(lbeta(a, b + n) - lbeta(a, b)) * a / (a + b + n) * (N - n) * (1 - z)
    # Ensure propagule pressure cannot be negative (ensure cases where n > N 
    # will not be chosen in the greedy optimisation., since incrementing n past 
    # N can otherwise result in negative propagule_pressure).
    propagule_pressure[propagule_pressure < 0] <- 0 
    return(propagule_pressure)
  }