# Create function to compute propagule pressure (or leakage) per consignment when
# assuming no overdispersion (i.e., all consignment in the pathway are assumed 
# to have the same infestation rate)
  # p: vector of mean infestation rate. One value for each consignment
  # N: vector of the number of units per consignement. One value for each consignment
  # n: inpsection sample size
  compute_propagule_pressure_per_cons_no_overdispersion <- function(p, n, N){
    propagule_pressure <- (1 - p)^n * p * (N - n)
    # Ensure propagule pressure cannot be negative (ensure cases where n > N 
    # will not be chosen in the greedy optimisation., since incrementing n past 
    # N can otherwise result in negative propagule_pressure). 
    propagule_pressure[propagule_pressure < 0] <- 0 
    return(propagule_pressure)
  }
  
  
  