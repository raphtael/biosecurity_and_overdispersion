################################################################################
# Compute approximate optimal sample size to minimize propagule pressure 
# This corresponds to equation 4 in the manuscript

# Note: this solution assumes the same alpha and beta for all consignments. If 
# this is not the case, then use either the numerical greedy algorithm solution, 
# the hybrid of a numerical and closed-form solution (Eq.7 in the manuscript), 
# or the approximate closed-form solution (Eq.8 in the manuscript, if alpha is 
# fixed) to the generalised risk equation.
  
  # Each row represents a consignment. 
  # a: alpha parameter of the Beta distribution.
  # b: beta parameter of the Beta distribution.
  # N: vector of the number of units per consignement. One value for each consignment.
    # Can vary among consignments.
  # B: Total budget effort in terms of cumulative sample size across all 
    # consignments. E.g., 600 times the number of consignements in the pathway.
  
# Compute optimal sample size but deal with the n > = 0 constraint
  compute_approximate_optimal_sample_size <- 
    function(a, b, N, B){
      if (length(unique(a)) > 1 | length(unique(b)) > 1) {
        stop("Alpha and beta should be fixed within a pathway")
      }
      data_input <- cbind.data.frame(a = a, 
                                     b = b,
                                     N = N, 
                                     B = B)
      data_input <- as.data.table(data_input)
      # Compute the number of consignments in the pathway
      data_input[, C := .N]
      data_input[, row_id := 1:.N] 
      # Compute optimal sample size per consignment
      data_input[, optimal_n := 
        - b + (b * C + B) / sum(N^(1/(a + 2))) * N^(1/(a + 2))] 
      # Deal with eventual sample size < 0 values
      while( sum(data_input$optimal_n < 0) > 0){
        # Set to zero
        data_input[optimal_n < 0, optimal_n := 0]            
        # Select consignments with non-zero sample size
        data_input_filtered <- data_input[optimal_n > 0, ]
        # Update optimal sample size
        data_input_filtered[, updated_optimal_n := 
          - b + (b * C + B) / sum(N^(1/(a + 2))) * N^(1/(a + 2))]
        # Merge zero and non zero sample size
        data_input <- merge(
          data_input[, .(row_id, a, b, N, B, optimal_n)], 
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