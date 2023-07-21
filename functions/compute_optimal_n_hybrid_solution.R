# Compute optimal sample size using the hybrid of numerical and closed-form 
  # solution (Eq.7 in the manuscript).
  # The Lagrangian multiplier (lambda) value is found numerically as the value
  # that obey the constraint outlined in Eq.S33 in section S4.2.5 of the 
  # supplementary material. 
  # risk and return optimal sample size. This solves Eq.7 in the manuscript
  
  # lambda: lagrange multiplier. Input function
  # a: vector of alpha parameters of the Beta distribution. Can vary among 
    # consignments.
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
    # consignments. Typically 600 times the number of consignments in the pathway.
  compute_optimal_n_hybrid_solution <- 
    function(lambda, a, b, N, B, cj = 1, d = 1, R = 1, z = 0){ 
      # Set dataframe of inputs
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
      # Compute initial conditions
        data_input[, lambda := 
          optimize(
            f = function(x){
              compute_constraint_for_hybrid_solution(
                lambda = x, 
                a = a, 
                b = b, 
                N = N, 
                B = B[1], 
                cj = cj, 
                d = d, 
                R = R, 
                z = z)^2
              }, # We square the constraint to optimise a convex function
            lower = -100, 
            upper = 0, 
            tol = 10^-10)$minimum]
        data_input[, optimal_n := 
          -b / d + (-lambda)^(-1/(a+2)) * 
          (d / cj * a * (a + 1) * b^a * N * R * (1 - z))^(1/(a+2)) / d]                                                                                       
      # Deal with eventual sample size < 0 values
        while( sum(data_input$optimal_n < 0) > 0){
          # Set to zero               
          data_input[optimal_n < 0, optimal_n := 0]
          # Select consignments with non-zero sample size                      
          data_input_filtered <- data_input[optimal_n > 0, ]   
          # Update optimal sample size
          data_input_filtered[, updated_lambda := 
            optimize(
              f = function(x) compute_constraint_for_hybrid_solution(
                lambda = x, 
                a = a, 
                b = b, 
                N = N, 
                B = B[1], 
                cj = cj, 
                d = d, 
                R = R, 
                z = z)^2, # We square the constraint to optimise a convex function
              lower = -100, 
              upper = 0, 
              tol = 10^-10)$minimum]      
          data_input_filtered[, updated_optimal_n := 
            -b / d + (-updated_lambda)^(-1/(a+2)) * 
            (d / cj * a * (a + 1) * b^a * N * R * (1 - z))^(1/(a+2)) / d]
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
      # Return solution
        setkey(data_input, row_id) # Sort by row
        return(data_input$optimal_n)
     }
            