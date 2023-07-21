################################################################################
#
# This script contains custom R functions using numerical greedy algorithms
# to compute the optimal sample size that minimize propagule pressure or risk
#                                                                 
################################################################################
                                                            
# Compute the optimal sample size per consignement that minimizes total leakage using a greedy algorithm
  # Each consignement has its row, with its associated, a, b, N values. 
  # Requires loading the data.table package and custom compute_propagule_pressure_per_cons function
  
  # a: vector of alpha parameters of the Beta distribution. One value for each consignment
  # b: vector of beta parameters of the Beta distribution. One value for each consignment
  # N: vector of the number of units per consignement. One value for each consignment
  # B: Total budget effort in terms of cumulative sample size across all consignments. Typically 600 times the number of consignements in the pathway.
  # n_increment_resolution: Increment resolution of the greedy algorithm at each step in terms of sample size
  
# Compute optimal sample size  
compute_optimal_sample_size_greedy <- 
  function(a, b, z, N, B, n_increment_resolution = 10){
  # Define the pathway and consignement problem
    data_greedy <- as.data.table(
      cbind.data.frame(a = a, b = b, z = z, N = N, n = 0))
  # Compute initial leakage, and leakage if we added n+1 to each consignement 
    leak_t1 <- data_greedy[, 
      list(leakage = 
        compute_propagule_pressure_per_cons(a = a, 
                                            b = b, 
                                            z = z, 
                                            N = N, 
                                            n = n)
      )]
    leak_t2 <- data_greedy[, 
      list(leakage = 
        compute_propagule_pressure_per_cons(a = a, 
                                            b = b, 
                                            z = z, 
                                            N = N, 
                                            n = n + n_increment_resolution)
      )]                                          
  # Distribute sampling effort among consignements, 
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
        leak_t2[cons_to_inspect, leakage := 
          data_greedy[cons_to_inspect, 
            compute_propagule_pressure_per_cons(
              a = a, 
              b = b, 
              z = z, 
              N = N, 
              n = n + n_increment_resolution)]]
    }
    close(progress_bar)
  # Return the optimal inspection size for each consignement
    return(data_greedy$n)
}
 

################################################################################  
# Compute the optimal sample size per consignement that minimizes total risk using a greedy algorithm
  # Each consignement has its row, with its associated, a, b, N values. 
  # a: vector of alpha parameters of the Beta distribution. One value for each consignment
  # b: vector of beta parameters of the Beta distribution. One value for each consignment
  # N: vector of the number of units per consignement. One value for each consignment
  # B: Total budget effort in terms of cumulative sample size across all consignments. Typically 600 times the number of consignements in the pathway.
  # cj: cost per unit inspected
  # d: detectability of the inspection
  # R: relative risk associated with different pathways (linked to the likelihood of establishment, spread, and damage associated with the pest)

# Compute optimal sample size 
  compute_optimal_sample_size_greedy_with_varying_abcdNRz <- 
    function(a, b, N, B, cj = 1, d = 1, R = 1, z = 0, n_increment_resolution = 10){  
    # Define the pathway and consignement problem
      data_greedy <- as.data.table(
        cbind.data.frame(a = a, b = b, cj, d, N = N, R = R, z = z, n = 0))
    # Compute initial risk, and risk if we added cj + 1 to each consignement
      risk_t1 <- data_greedy[, 
        list(
          risk = compute_propagule_pressure_per_cons(
            a = a, 
            b = b, 
            z = z, 
            N = N, 
            n = n * d
          ) * R
        )] 
      
      risk_t2 <- data_greedy[, 
        list(
          risk = compute_propagule_pressure_per_cons(
            a = a, 
            b = b, 
            z = z, 
            N = N, 
            n = n * d + n_increment_resolution * d / cj
          ) * R
        )]                                    
    # Distribute sampling effort among consignements, 
      # increments n by n_increment_resolution in every step of the loop
      n_loop_steps <- floor(B / n_increment_resolution)
      progress_bar <- txtProgressBar(min = 0, max = n_loop_steps, style = 3)
      for (i in 1:n_loop_steps){ 
        # Progress bar
          setTxtProgressBar(progress_bar, i)
        # Find the consignment which would reduce risk the most
          delta_risk <- risk_t1 - risk_t2
          cons_to_inspect <- which.max(delta_risk$risk) 
        # Update n, risk_t1, and risk_t2
          data_greedy[cons_to_inspect, n := n + n_increment_resolution / cj]
          risk_t1[cons_to_inspect, risk := risk_t2[cons_to_inspect, risk] ]
          risk_t2[cons_to_inspect, risk := 
            data_greedy[cons_to_inspect, 
              compute_propagule_pressure_per_cons(
                a = a, 
                b = b, 
                z = z, 
                N = N, 
                n = n * d + n_increment_resolution * d / cj) * R ]]
      }
      close(progress_bar)
    # Return the optimal inspection size for each consignement
      return(data_greedy$n)
  }
    
  
  