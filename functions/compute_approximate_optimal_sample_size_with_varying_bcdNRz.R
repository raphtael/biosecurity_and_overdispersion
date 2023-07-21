################################################################################
# Compute approximate optimal sample size to minimize risk under varying parameters
# This corresponds to equation 8 in the manuscript
# Note: this solution assumes the same alpha for all pathways. If this is not 
# the case, then use the hybrid of a numerical and closed-form solution 
# (Eq.7 in the manuscript) or the numerical greedy algorithm solution.
  
  # Each row represents a consignment. 
  # a: vector of alpha parameters of the Beta distribution. Alpha values should 
    # be the same among consignments.
  # b: vector of beta parameters of the Beta distribution. Can vary among 
    # consignments.
  # z: vector of zero inflated parameters of the zero-inflated beta binomial. 
    # distribution (default value of zero, i.e., beta-binomial distribution).
    # Can vary among consignments.
  # N: vector of the number of units per consignement. One value for each consignment.
    # Can vary among consignments.
  # cj: vector of cost per sampled unit. Default value of 1 (i.e., total budget 
    # is defined in terms of sampling units). Can vary among consignments.
  # d: vector of detectability parameter (i.e., likelihood of detecting a 
    # biosecurity risk material (BRM) when present when inspecting the unit. 
    # Default value of 1 (i.e., perfect detectability). Can vary among consignments.
  # R: vector of relative risk parameter (damage caused by introducing a single 
    # unit of BRM into the country. Can vary among consignments.
  # B: Total budget effort in terms of cumulative sample size across all 
    # consignments. Typically 600 times the number of consignements in the pathway.

# Optimal sample size with varying bcdRNz values 
  compute_approximate_optimal_sample_size_with_varying_bcdNRz <- 
    function(a, b, N, B, cj = 1, d = 1, R = 1, z = 0){
      if (length(unique(a)) > 1) {
        stop("In compute_approximate_optimal_sample_size_varying_bcdNRz, 
          alpha is assumed to have the same value for all consignments. Please use 
          the hybrid of numerical and closed form solution (lambda) or the greedy
          algorithm solution if alpha varies among consignments.")
      }
      data_input <- cbind.data.frame(a = a, 
                                     b = b, 
                                     cj = cj, 
                                     d = d, 
                                     N = N, 
                                     R = R, 
                                     z = z, 
                                     B = B)
      data_input <- as.data.table(data_input)
      data_input[, row_id := 1:.N]
      data_input[, optimal_n := 
        - b / d + (sum(cj / d * b) + B) / 
        (sum(cj / d * (d / cj * b^a * N * R * (1 - z))^(1/(a + 2))) ) * 
        (d / cj * b^a * N * R * (1 - z))^(1/(a + 2)) * 1/d] # Compute optimal sample size per consignment
      # Deal with eventual sample size < 0 values
      while( sum(data_input$optimal_n < 0) > 0){
        # Set to zero               
        data_input[optimal_n < 0, optimal_n := 0]
        # Select consignments with non-zero sample size                      
        data_input_filtered <- data_input[optimal_n > 0, ]   
        # Update optimal sample size
        data_input_filtered[, updated_optimal_n := 
          - b / d + (sum(cj / d * b) + B) / 
          (sum(cj / d * (d / cj * b^a * N * R * (1 - z))^(1/(a + 2))) ) * 
          (d / cj * b^a * N * R * (1 - z))^(1/(a + 2)) * 1/d]
        # Merge zero and non zero sample size
        data_input <- merge(
          data_input[, .(row_id, a, b, cj, d, N, R, z, B, optimal_n)], 
          data_input_filtered[, .(row_id, updated_optimal_n)], 
          by = 'row_id', 
          all = TRUE) 
        # Update non-zero sample size values with their new optimal value  
        data_input[optimal_n > 0, optimal_n := updated_optimal_n] 
        data_input[, updated_optimal_n := NULL]
      }
      setkey(data_input, row_id) # Sort by row
      return(data_input$optimal_n)
    }
    