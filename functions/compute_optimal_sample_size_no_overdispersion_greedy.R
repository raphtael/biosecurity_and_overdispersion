################################################################################
#
# Compute optimal sample size assuming no overdispersion using greedy algorithm
#                                                                
################################################################################
                                                            
# Compute the optimal sample size per consignement (assuming no overdispersion) 
# that minimizes total leakage using a greedy algorithm
# Each row represents a consignment with its associated a, b, N values. 
# Requires loading compute_propagule_pressure_per_cons_no_overdispersion
  
  # p: vector of mean infestation rate. One value for each consignment
  # N: vector of the number of units per consignement. One value for each consignment
  # B: Total budget effort in terms of cumulative sample size across all consignments. Typically 600 times the number of consignements in the pathway.
  # n_increment_resolution: Increment resolution of the greedy algorithm at each step in terms of sample size

# Compute optimal sample size assuming no overdispsersion
compute_optimal_sample_size_no_overdispersion_greedy <- 
  function(p, N, B, n_increment_resolution = 10){
    # Define the pathway and consignement problem
      data_greedy <- as.data.table(cbind.data.frame(p = p, N = N, n = 0))
    # Compute initial leakage, and leakage if we added n + 1 to each consignement 
      leak_t1 <- data_greedy[, .(leakage = 
        compute_propagule_pressure_per_cons_no_overdispersion(
          p = p, 
          N = N, 
          n = n
        ))]
      leak_t2 <- data_greedy[, .(leakage = 
        compute_propagule_pressure_per_cons_no_overdispersion(
          p = p, 
          N = N, 
          n = n + n_increment_resolution
        ))]                                         
    # Distribute sampling effort among consignements, incrementing n by 
      # increments n by n_increment_resolution in every step of the loop
      n_loop_steps <- floor(B / n_increment_resolution)
      progress_bar <- txtProgressBar(min = 0, max = n_loop_steps, style = 3) 
      for (i in 1:n_loop_steps){ 
        # Progress bar
          setTxtProgressBar(progress_bar, i)
        # Find the consignment which would reduce leakage the most
          delta_leak <- leak_t1 - leak_t2
          cons_to_inspect <- which.max(delta_leak$leakage) 
        # Update n, leak_t1, and leak_t2
          data_greedy[cons_to_inspect, n := n + n_increment_resolution]
          leak_t1[cons_to_inspect, leakage := leak_t2[cons_to_inspect, leakage]]
          leak_t2[cons_to_inspect, leakage := data_greedy[cons_to_inspect, 
            compute_propagule_pressure_per_cons_no_overdispersion(
              p = p, 
              N = N, 
              n = n + n_increment_resolution
            )]]
      }
      close(progress_bar)
    # Return the optimal inspection size for each consignement
      return(data_greedy$n)
  }



