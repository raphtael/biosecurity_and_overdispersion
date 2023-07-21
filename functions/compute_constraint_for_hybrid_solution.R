# Compute the constraint Eq.S33 outlined in section S4.2.5 of the supplementary 
  # material. Computing this constraint is one of the necessary step to find the 
  # optimal sample size using the hybrid of numerical and closed-form solution 
  # (Eq.7 in the manuscript). This  fuction is one of the component of the
  # compute_optimal_n_from_lambda_with_varying_abcdNRz function. 
  
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
  compute_constraint_for_hybrid_solution <- function(lambda, a, b, N, B, cj = 1, d = 1, R = 1, z = 0){
    constraint <- - sum(cj / d * b) - B + sum((-lambda)^(-1/(a+2)) * (a * (a + 1) * d / cj * b^a * N * R * (1 - z))^(1 / (a+2)) * cj / d)
    return(constraint)
  }  
