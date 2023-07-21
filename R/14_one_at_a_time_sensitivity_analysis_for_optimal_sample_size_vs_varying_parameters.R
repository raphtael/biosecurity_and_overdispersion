################################################################################
#
#  One at a time sensitivity analysis of optimal sample size vs. varying parameter
#
################################################################################

# Set working directory
  setwd("C:/Users/rtrouve/Dropbox/overdispersion_paper")
  
# Load libraries
  library(ggplot2)
  library(data.table)
  library(latex2exp)
  library(cowplot)
  source('./functions/compute_propagule_pressure_per_cons.R')
  source('./functions/extract_bb_parameters_per_pathway_from_brms_fit.R')  
  source('./functions/extract_zibb_parameters_per_pathway_from_brms_fit.R')  
  source('./functions/inv_logit.R')
  source('./functions/compute_optimal_sample_size_greedy.R')
  source('./functions/compute_optimal_n_hybrid_solution.R')
  source('./functions/compute_constraint_for_hybrid_solution.R')                 
  
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
   
# Compute median and range per parameter of the sensitivity analysis
  a_median <- median(data_fresh_with_coef$a)
  a_min <- min(data_fresh_with_coef$a)
  a_max <- max(data_fresh_with_coef$a)
  sd_a <- sd(log(data_fresh_with_coef$a))
  b_median <- median(data_fresh_with_coef$b) 
  b_min <- min(data_fresh_with_coef$b)
  b_max <- max(data_fresh_with_coef$b)
  N_median <- median(data_fresh_with_coef$N) 
  N_min <- min(data_fresh_with_coef$N)
  N_max <- max(data_fresh_with_coef$N)
  z_median <- median(data_fresh_with_coef$z_zibb_fixed_a) # Proportion of excess zero
  z_min <- min(data_fresh_with_coef$z_zibb_fixed_a)
  z_max <- max(data_fresh_with_coef$z_zibb_fixed_a)

################################################################################   
# Sensitivity analysis on the alpha parameter
  # Define sensitivity analysis sampling design
    set.seed(1001)
    a_distr <- exp(runif(n = 100, min = log(a_min), max = log(a_max + 0.1))) 
    data_sim_varying_a <- expand.grid(a = a_distr, 
                                      b = b_median, 
                                      z = 0, 
                                      N = N_median)
    data_sim_varying_a <- as.data.table(data_sim_varying_a)
  # Compute sample size
    ptm <- proc.time()
      data_sim_varying_a[, optimal_n_greedy := 
        compute_optimal_sample_size_greedy(a = a, 
                                           b = b, 
                                           z = z, 
                                           N = N, 
                                           B = .N * 600, 
                                           n_increment_resolution = 10)]
    proc.time() - ptm
    data_sim_varying_a[, optimal_n_hybrid_solution := 
      compute_optimal_n_hybrid_solution(a = a, b = b, N = N, B = .N * 600)]
  # Plot optimal sample size vs. varying alpha
    p_a <- ggplot(data = data_sim_varying_a) + 
      geom_point(aes(a, optimal_n_greedy)) + 
      geom_line(aes(a, optimal_n_hybrid_solution), col = 'dodgerblue') + 
      ylab(TeX('$n_i^*$')) +
      xlab(TeX('$\\alpha_i$')) +
      theme_bw() +
      theme(panel.grid.minor = element_blank())
    p_a
  
# Sensitivity analysis on beta parameter
  # Define sensitivity analysis sampling design
    set.seed(1002)
    b_distr <- exp(runif(n = 100, min = log(b_min), max = log(b_max))) # 
    data_sim_varying_b <- expand.grid(a = a_median, 
                                      b = b_distr, 
                                      z = 0, 
                                      N = N_median)
    data_sim_varying_b <- as.data.table(data_sim_varying_b)
  # Compute sample size  
    ptm <- proc.time()
      data_sim_varying_b[, optimal_n_greedy := 
        compute_optimal_sample_size_greedy(a = a, 
                                           b = b, 
                                           z = z, 
                                           N = N, 
                                           B = .N * 600, 
                                           n_increment_resolution = 4)]
    proc.time() - ptm
    data_sim_varying_b[, optimal_n_hybrid_solution := 
      compute_optimal_n_hybrid_solution(a = a, b = b, N = N, B = .N * 600)]
  # Plot optimal sample size vs. varying beta
    p_b <- ggplot(data = data_sim_varying_b, aes(b, optimal_n_greedy)) + 
      geom_point() + 
      geom_line(aes(b, optimal_n_hybrid_solution), col = 'dodgerblue') +
      ylab(TeX('$n_i^*$')) +
      xlab(TeX('$\\beta_i$')) +
      theme_bw() +
      theme(panel.grid.minor = element_blank())
    p_b
  
# Sensitivity analysis on consignment size N
  # Define sensitivity analysis sampling design
    N_distr <- exp(runif(n = 100, min = log(N_min), max = log(N_max))) # 
    data_sim_varying_N <- expand.grid(a = a_median, b = b_median, z = 0, N = N_distr)
    data_sim_varying_N <- as.data.table(data_sim_varying_N)
  # Compute sample size
    ptm <- proc.time()
      data_sim_varying_N[, optimal_n_greedy := 
        compute_optimal_sample_size_greedy(a = a, 
                                           b = b, 
                                           z = z, 
                                           N = N, 
                                           B = .N * 600, 
                                           n_increment_resolution = 10)]
    proc.time() - ptm
    data_sim_varying_N[, optimal_n_hybrid_solution := 
      compute_optimal_n_hybrid_solution(a = a, b = b, N = N, B = .N * 600)]
  # Plot optimal sample size vs. varying N    
    p_N <- ggplot(data = data_sim_varying_N, aes(N, optimal_n_greedy)) + 
      geom_point() + 
      geom_line(aes(N, optimal_n_hybrid_solution), col = 'dodgerblue') +
      ylab(TeX('$n_i^*$')) +
      xlab(TeX('$N_i$')) +
      theme_bw() +
      theme(panel.grid.minor = element_blank())
    p_N  
 
# Sensitivity analysis on the zero-inflated parameter z
  # Define sensitivity analysis sampling design
    set.seed(1001)
    z_distr <- exp(runif(n = 100, min = log(z_min), max = log(z_max))) 
    data_sim_varying_z <- expand.grid(a = a_median, b = b_median, z = z_distr, N = N_median)
    data_sim_varying_z <- as.data.table(data_sim_varying_z)
  # Compute sample size
    ptm <- proc.time()
      data_sim_varying_z[, optimal_n_greedy := 
        compute_optimal_sample_size_greedy(a = a, 
                                           b = b, 
                                           z = z, 
                                           N = N, 
                                           B = .N * 600, 
                                           n_increment_resolution = 10)]
    proc.time() - ptm
    data_sim_varying_z[, optimal_n_hybrid_solution := 
      compute_optimal_n_hybrid_solution(a = a, b = b, N = N, z = z, B = .N * 600)]
  # Plot optimal sample size vs. varying z 
    p_z <- ggplot(data = data_sim_varying_z, aes(1 - z, optimal_n_greedy)) + 
      geom_point() + 
      geom_line(aes(1 - z, optimal_n_hybrid_solution), col = 'dodgerblue') +
      ylab(TeX('$n_i^*$')) +
      xlab(TeX('$q_i = 1 - \\z_i$')) +
      theme_bw() +
      theme(panel.grid.minor = element_blank())
    p_z  

# Sensitivity analysis on relative risk R
  # Define sensitivity analysis sampling design
    set.seed(1001)
    R_distr <- exp(runif(n = 100, min = log(0.1), max = log(10)))
    data_sim_varying_R <- expand.grid(a = a_median, b = b_median, z = 0, N = N_median, R = R_distr)
    data_sim_varying_R <- as.data.table(data_sim_varying_R)
  # Compute sample size
    ptm <- proc.time()
      data_sim_varying_R[, optimal_n_greedy := 
        compute_optimal_sample_size_greedy_with_varying_abcdNRz(
          a = a, 
          b = b, 
          z = z, 
          N = N, 
          R = R,
          B = .N * 600, 
          n_increment_resolution = 10
        )]
    proc.time() - ptm
    data_sim_varying_R[, optimal_n_hybrid_solution := 
      compute_optimal_n_hybrid_solution(a = a, b = b, N = N, R = R, B = .N * 600)]
  # Plot optimal sample size vs. varying R
    p_R <- ggplot(data = data_sim_varying_R, aes(R, optimal_n_greedy)) + 
      geom_point() + 
      geom_line(aes(R, optimal_n_hybrid_solution), col = 'dodgerblue') +
      ylab(TeX('$n_i^*$')) +
      xlab(TeX('$R_i$')) +
      theme_bw() +
      theme(panel.grid.minor = element_blank())
    p_R 
      
# Sensitivity analysis on cost per sample c
  # Define sensitivity analysis sampling design
    set.seed(1001)
    c_distr <- exp(runif(n = 100, min = log(0.1), max = log(10)))
    data_sim_varying_c <- expand.grid(a = a_median, b = b_median, z = 0, cj = c_distr, N = N_median)
    data_sim_varying_c <- as.data.table(data_sim_varying_c)
  # Compute sample size
    ptm <- proc.time()
      data_sim_varying_c[, optimal_n_greedy := 
        compute_optimal_sample_size_greedy_with_varying_abcdNRz(
          a = a, 
          b = b, 
          z = z, 
          cj = cj,
          N = N, 
          B = .N * 600, 
          n_increment_resolution = 10
        )]
    proc.time() - ptm
    data_sim_varying_c[, optimal_n_hybrid_solution := 
      compute_optimal_n_hybrid_solution(a = a, b = b, N = N, cj = cj, B = .N * 600)]
  # Plot optimal sample size vs. varying c 
    p_c <- ggplot(data = data_sim_varying_c, aes(1 / cj, optimal_n_greedy)) + 
      geom_point() + 
      geom_line(aes(1 / cj, optimal_n_hybrid_solution), col = 'dodgerblue') +
      ylab(TeX('$n_i^*$')) +
      xlab(TeX('$\\frac{1}{c_i}$')) +
      theme_bw() +
      theme(panel.grid.minor = element_blank())
    p_c 
  
# Sensitivity analysis on detectability d
  # Define sensitivity analysis sampling design
    set.seed(1001)
    d_distr <- exp(runif(n = 100, min = log(0.1), max = log(1)))
    data_sim_varying_d <- expand.grid(a = a_median, b = b_median, z = 0, d = d_distr, N = N_median)
    data_sim_varying_d <- as.data.table(data_sim_varying_d)
  # Compute sample size
    ptm <- proc.time()
      data_sim_varying_d[, optimal_n_greedy := 
        compute_optimal_sample_size_greedy_with_varying_abcdNRz(
          a = a, 
          b = b, 
          z = z, 
          d = d,
          N = N, 
          B = .N * 600, 
          n_increment_resolution = 10
        )]
    proc.time() - ptm
    data_sim_varying_d[, optimal_n_hybrid_solution := 
      compute_optimal_n_hybrid_solution(a = a, b = b, N = N, d = d, B = .N * 600)]
  # Plot optimal sample size vs. varying d
    p_d <- ggplot(data = data_sim_varying_d, aes(d, optimal_n_greedy)) + 
      geom_point() + 
      geom_line(aes(d, optimal_n_hybrid_solution), col = 'dodgerblue') +
      ylab(TeX('$n_i^*$')) +
      xlab(TeX('$d_i$')) +
      theme_bw() +
      theme(panel.grid.minor = element_blank())
    p_d  

  # Plot effective optimal sample size vs. varying d
    p_d_eff <- ggplot(data = data_sim_varying_d, aes(d, optimal_n_greedy * d)) + 
      geom_point() + 
      geom_line(aes(d, optimal_n_hybrid_solution * d), col = 'dodgerblue') +
      ylab(TeX('$n_{eff_i}^* = n_i^* d_i$')) +
      xlab(TeX('$d_i$')) +
      theme_bw() +
      theme(panel.grid.minor = element_blank())
    p_d_eff 
  
  
# Combine plots into a single figure 
  pdf('./figs/optimal_n_vs_varying_parameters_one_at_a_time.pdf', width = 14, height = 7)
    plot_grid(p_a, p_b, p_z, p_N, p_R, p_c, p_d, p_d_eff, ncol = 4, align = 'hv')
  dev.off()
 
  