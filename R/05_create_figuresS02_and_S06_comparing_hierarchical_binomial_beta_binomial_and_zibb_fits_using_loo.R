################################################################################
#
#Compare hierarchical Beta-Binomial vs. binomial and ZIBB models, group by group 
#                        
################################################################################

# Set working directory
  setwd("C:/Users/rtrouve/Dropbox/overdispersion_paper")

# Load packages
  library(ggplot2)
  library(data.table)
  library(brms)
  library(latex2exp)
  source('./functions/brms_beta_binomial_familly.R')
  source('./functions/brms_zero_inflated_beta_binomial_familly.R')
  source('./functions/compute_loo_per_group.R')
  tick_break_log_scale <- c(10, 30, 100, 300, 1000, 3000, 10000)
   
################################################################################
# Read data
  data_all_pathways <- fread('./data/data_all_pathways.csv')

################################################################################  
# Compare binomial and beta-binomial models
  # Load models
    load('./outputs/fit_binomial_hierarchical.Rdata')
    load('./outputs/fit_beta_binomial_hierarchical.Rdata')
    expose_functions(fit_beta_binomial_hierarchical, vectorize = TRUE)
  
  # Compute loo for whole dataset
    loo_binom <- brms::loo(fit_binomial_hierarchical)
    loo_beta_binom <- brms::loo(fit_beta_binomial_hierarchical)
    loo_binom
    loo_beta_binom
    
  # Compute and store loo per group for binomial and beta binomial model 
    data_loo_results <- data_all_pathways[, 
      list(
        loo_binom = compute_loo_per_group(
          model = fit_binomial_hierarchical, 
          k = k, 
          n = n, 
          group_id = group_id
        ),
        loo_beta_binom = compute_loo_per_group(
          model = fit_beta_binomial_hierarchical, 
          k = k, 
          n = n, 
          group_id = group_id
        )
      ), 
      by = .(data_type, group_id)
    ]    
    summary(data_loo_results[, .(loo_beta_binom - loo_binom)])
    
  # Plot results              
    p1 <- ggplot(data = data_loo_results, aes(loo_binom, loo_beta_binom)) +
      geom_point() + 
      scale_x_log10(breaks = tick_break_log_scale) + 
      scale_y_log10(breaks = tick_break_log_scale) + 
      geom_abline() + 
      ylab('LOO Beta-Binomial model') +
      xlab('LOO Binomial model') + 
      coord_fixed(ratio = 1) +
      theme_bw() + 
      theme(panel.grid.minor=element_blank())
    ggsave(
      plot = p1, 
      filename = './figs/loo_beta_binom_vs_binom_model_comparison.pdf', 
      height = 8, 
      width = 8, 
      units = 'cm'
    )  
      
################################################################################
# Compare BB and ZIBB
  # Load model
    load('./outputs/fit_zibb_hierarchical.Rdata')
    expose_functions(fit_zibb_hierarchical, vectorize = TRUE)
    
  # Compute loo for whole dataset
    loo_zibb <- brms::loo(fit_zibb_hierarchical)
    loo_beta_binom <- brms::loo(fit_beta_binomial_hierarchical)
    loo_zibb
    loo_beta_binom
    
  # Compute and store loo per group for binomial and beta binomial model 
    data_loo_results_zibb <- data_all_pathways[, 
      list(
        loo_zibb = compute_loo_per_group(
          model = fit_zibb_hierarchical, 
          k = k, 
          n = n, 
          group_id = group_id
        ),
        loo_beta_binom = compute_loo_per_group(
          model = fit_beta_binomial_hierarchical, 
          k = k, 
          n = n, 
          group_id = group_id
        )
      ), 
      by = .(data_type, group_id)
    ]
    summary(data_loo_results_zibb[, .(loo_beta_binom - loo_zibb)])
    
  # Plot results              
    p2 <- ggplot(data = data_loo_results_zibb, aes(loo_zibb, loo_beta_binom)) +
      geom_point() + 
      scale_x_log10(breaks = tick_break_log_scale) + 
      scale_y_log10(breaks = tick_break_log_scale) + 
      geom_abline() + 
      ylab('LOO Beta-Binomial model') +
      xlab('LOO ZIBB model') + 
      coord_fixed(ratio = 1) + 
      ggtitle('a.') +
      theme_bw() + 
      theme(panel.grid.minor=element_blank())
    ggsave(
      plot = p2, 
      filename = './figs/loo_beta_binom_vs_zibb_model_comparison.pdf', 
      height = 8, 
      width = 8, 
      units = 'cm'
    )

################################################################################
# Compare BB and ZIBB with fixed alpha
  # Load model  
    load('./outputs/fit_zibb_hierarchical_fixed_a.Rdata')
    expose_functions(fit_zibb_hierarchical_fixed_a, vectorize = TRUE)
  
  # Compute loo for whole dataset
    loo_zibb_fixed_a <- brms::loo(fit_zibb_hierarchical_fixed_a)
    loo_beta_binom <- brms::loo(fit_beta_binomial_hierarchical)
    loo_zibb_fixed_a
    loo_beta_binom
    
  # Compute and store loo per group for binomial and beta binomial model 
    data_loo_results_zibb_fixed_a <- data_all_pathways[, 
      list(
        loo_zibb_fixed_a = compute_loo_per_group(
          model = fit_zibb_hierarchical_fixed_a, 
          k = k, 
          n = n, 
          group_id = group_id
        ),
        loo_beta_binom = compute_loo_per_group(
          model = fit_beta_binomial_hierarchical, 
          k = k, 
          n = n, 
          group_id = group_id
        )
      ), 
      by = .(data_type, group_id)
    ] 
    summary(
      data_loo_results_zibb_fixed_a[, .(loo_beta_binom - loo_zibb_fixed_a)]
    )
    
  # Plot results              
    p3 <- ggplot(data = data_loo_results_zibb_fixed_a, 
             aes(loo_zibb_fixed_a, loo_beta_binom)) +
      geom_point() + 
      scale_x_log10(breaks = tick_break_log_scale) + 
      scale_y_log10(breaks = tick_break_log_scale) + 
      geom_abline() + 
      ylab('LOO Beta-Binomial model') +
      xlab(TeX('LOO ZIBB with fixed $\\alpha$ model')) + 
      coord_fixed(ratio = 1) +
      ggtitle('b.') +
      theme_bw() + 
      theme(panel.grid.minor=element_blank())
    ggsave(
      plot = p3, 
      filename = './figs/loo_beta_binom_vs_zibb_fixed_a_model_comparison.pdf', 
      height = 8, 
      width = 8, 
      units = 'cm'
    )
    
    
    