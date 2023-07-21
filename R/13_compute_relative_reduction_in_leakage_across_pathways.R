################################################################################
#
#    Compute relative reduction in propagule pressure across all pathways
#
################################################################################

# Set working directory
  setwd("C:/Users/rtrouve/Dropbox/overdispersion_paper")
  
# Load libraries
  library(data.table)
  source('./functions/compute_propagule_pressure_per_cons.R')
  source('./functions/extract_bb_parameters_per_pathway_from_brms_fit.R')  
  source('./functions/extract_zibb_parameters_per_pathway_from_brms_fit.R')  
  source('./functions/inv_logit.R')  
  source('./functions/compute_approximate_optimal_sample_size.R') 
  source('./functions/compute_optimal_sample_size_greedy.R')
  source('./functions/compute_approximate_optimal_sample_size_with_varying_bcdNRz.R')
  source('./functions/compute_optimal_n_hybrid_solution.R')
  source('./functions/compute_constraint_for_hybrid_solution.R')                 
                      
################################################################################ 
#              Optimisation of sampling efforts across pathways
################################################################################ 

# Read the summary table and models
  data_all_pathways <- fread('./data/data_all_pathways.csv') 
  load('./outputs/fit_beta_binomial_hierarchical.Rdata')
  load('./outputs/fit_zibb_hierarchical_fixed_a.Rdata')

# Extract coefficients of the beta-binomial and ZIBB and merge with data
  data_coef_summary_bb <- 
    extract_bb_parameters_per_pathway_from_brms_fit(
      model = fit_beta_binomial_hierarchical
    )
  data_coef_summary_zibb_fixed_a <- 
    extract_zibb_parameters_per_pathway_from_brms_fit(
      model = fit_zibb_hierarchical_fixed_a
    )    
  data_all_with_coef <- 
    merge(
      data_coef_summary_bb[, 
        list(
          group_id, 
          a = alpha_mean, 
          b = beta_mean
        )
      ], 
      data_all_pathways, 
      by = 'group_id'
    )
  data_all_with_coef <- 
    merge(
      data_coef_summary_zibb_fixed_a[, 
        list(
          group_id, 
          a_zibb_fixed_a = alpha_mean, 
          b_zibb_fixed_a = beta_mean, 
          z_zibb_fixed_a = zi_mean
        )
      ], 
      data_all_with_coef, 
      by = 'group_id'
    )
    
    
# Select fresh produce pathway data
  data_fresh_with_coef <- data_all_with_coef[data_type %in% c('fresh_produce'),]
# Compute total budget 
  # B = number of consignments * 600 samples per consignment
  data_fresh_with_coef[, B := .N * 600] 
  
# Compute sample size 
  # Approximate optimal sample size that optimise based solely on variation in N
    # Ignore variation in alpha and beta among pathways. Eq. 4 in manuscript.
    data_fresh_with_coef[, optimal_n_based_on_N_variation_only := 
      compute_approximate_optimal_sample_size(
        a = median(a), 
        b = median(b), 
        N = N, 
        B = B
      )]
  # Approximate optimal sample size proportional to sqrt(N) 
    # (assumes alpha = 0, neglect beta offset)   
    data_fresh_with_coef[, optimal_n_prop_sqrt_N := 
      B / sum(sqrt(N)) * sqrt(N)] 
  # Hybrid of numerical and closed-form solution 
    # Allows alpha and beta to vary. Eq. 7 in manuscript
    data_fresh_with_coef[, optimal_n_hybrid_solution := 
      compute_optimal_n_hybrid_solution(
        a = a, 
        b = b, 
        N = N, 
        B = B
      )]
  # Compute optimal sample size using the closed-form solution for a simplified 
    # problem. Uses the same alpha for all pathways, but let the z parameter 
    # absorb the variation in non-zero infestation rate among pathways.
    # Eq. 8 in manuscript. 
    data_fresh_with_coef[, optimal_n_zibb_fixed_a := 
      compute_approximate_optimal_sample_size_with_varying_bcdNRz(
        a = a_zibb_fixed_a, 
        b = b_zibb_fixed_a, 
        z = z_zibb_fixed_a, 
        N = N, 
        B = B
      )]
  # Greedy algorithm solution
    ptm <- proc.time()
      data_fresh_with_coef[, optimal_n_greedy := 
        compute_optimal_sample_size_greedy(
          a = a, 
          b = b, 
          z = 0, 
          N = N, 
          B = B[1], 
          n_increment_resolution = 10 # Use 1 for improved (but slower) results
        )]
    proc.time() - ptm
  
# Compute propagule pressure
  data_fresh_with_coef[, propagule_pressure_fixed_n := 
    compute_propagule_pressure_per_cons(a = a, b = b, n = 600, N = N)]
  data_fresh_with_coef[, propagule_pressure_optimal_n_based_on_N_variation_only := 
    compute_propagule_pressure_per_cons(a = a, 
                                        b = b, 
                                        n = optimal_n_based_on_N_variation_only, 
                                        N = N)]
  data_fresh_with_coef[, propagule_pressure_optimal_n_prop_sqrt_N := 
    compute_propagule_pressure_per_cons(a = a, 
                                        b = b, 
                                        n = optimal_n_prop_sqrt_N, 
                                        N = N)]
  data_fresh_with_coef[, propagule_pressure_optimal_n_hybrid_solution := 
    compute_propagule_pressure_per_cons(a = a, 
                                        b = b, 
                                        n = optimal_n_hybrid_solution, 
                                        N = N)]
  data_fresh_with_coef[, propagule_pressure_optimal_n_zibb_fixed_a := 
    compute_propagule_pressure_per_cons(a = a, 
                                        b = b, 
                                        n = optimal_n_zibb_fixed_a, 
                                        N = N)]
  data_fresh_with_coef[, propagule_pressure_optimal_n_greedy := 
    compute_propagule_pressure_per_cons(a = a, 
                                        b = b, 
                                        n = optimal_n_greedy, 
                                        N = N)]
# Write dataframe
  fwrite(
    data_fresh_with_coef, 
    './outputs/data_fresh_with_coef_with_added_leakage_across_pathways.csv'
  )  

                                       
# Compute total propagule pressure across pathways for different sampling regimes                                        
  data_fresh_with_coef <- 
    fread('./outputs/data_fresh_with_coef_with_added_leakage_across_pathways.csv')
  
  data_propagule_pressure_summary <- data_fresh_with_coef[, 
    list(
      propagule_pressure_fixed_n = sum(propagule_pressure_fixed_n),
      propagule_pressure_optimal_n_based_on_N_variation_only = 
        sum(propagule_pressure_optimal_n_based_on_N_variation_only), 
      propagule_pressure_optimal_n_prop_sqrt_N = 
        sum(propagule_pressure_optimal_n_prop_sqrt_N),
      propagule_pressure_optimal_n_hybrid_solution = 
        sum(propagule_pressure_optimal_n_hybrid_solution),
      propagule_pressure_optimal_n_zibb_fixed_a = 
        sum(propagule_pressure_optimal_n_zibb_fixed_a),
      propagule_pressure_optimal_n_greedy = 
        sum(propagule_pressure_optimal_n_greedy)
    )]                                        
                                          
# Relative propagule pressure reduction compared to a fixed sample size  
  data_propagule_pressure_summary[, 
    list(
      round(propagule_pressure_optimal_n_based_on_N_variation_only / propagule_pressure_fixed_n, 3),
      round(propagule_pressure_optimal_n_prop_sqrt_N / propagule_pressure_fixed_n, 3),
      round(propagule_pressure_optimal_n_hybrid_solution / propagule_pressure_fixed_n, 3),
      round(propagule_pressure_optimal_n_zibb_fixed_a / propagule_pressure_fixed_n, 3),
      round(propagule_pressure_optimal_n_greedy / propagule_pressure_fixed_n, 3)
    )]
    
 