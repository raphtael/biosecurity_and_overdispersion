# Function to compute sample size with hypergeometric sampling under budget 
# constraint.

# Start with an initial sensitivity and design prevalence (e.g., 95% and 0.5%, 
# to compare with a fixed sample size of 600). Compared to a fixed sample size,
# the reduced sample size for small consignments under hypergeometric sampling 
# frees up a few samples that can be reallocated to larger consignments. This 
# reallocation was implemented by iteratively increasing the sensitivity of the 
# inspection until reaching our total budget constraint. 

# B: Total budget effort in terms of cumulative sample size across all 
  # consignments. Typically 600 times the number of consignements in the pathway.
# N: consignment size
# S0: Initial sensitivity of the inspection (also called confidence level). 
  # Default value of 0.95
  # p: design prevalence (also called level of detection). Default value of 0.005
# S_increment_resolution: resolution of the increment in terms of sensitivity

# Compute sample size with hypergeometric sampling under budget constraint.
  compute_hypergeometric_sample_size_under_budget_constraint <- 
    function(B, N, p = 0.005, S0 = 0.95, S_increment_resolution = 0.001){
      hypergeometric_n <- 
        compute_approximate_hypergeometric_sample_size(S = S0, p = p, N = N)
      Si <- S0 # Initialise sensitivity at e.g., 0.95
      while(sum(hypergeometric_n) < B | Si == 1){
        # Increase the sensitivity of the test
          Si <- Si + S_increment_resolution
          hypergeometric_n <- 
            compute_approximate_hypergeometric_sample_size(S = Si, p = p, N = N)
      }
      # Recalculate sample size using the previous step's sensitivity where we 
      # were below the total budget
      S_to_use <- Si - S_increment_resolution 
      hypergeometric_n <- 
        compute_approximate_hypergeometric_sample_size(S = S_to_use, p = p, N = N)
      return(cbind.data.frame(hypergeometric_n, S_to_use)) 
    }  
