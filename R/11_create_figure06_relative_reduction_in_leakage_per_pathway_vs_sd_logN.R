################################################################################
#
#  Create Fig.06: relative reduction in propagule pressure vs. sigma_log(N)
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
  source('./functions/compute_approximate_optimal_sample_size.R') 
  source('./functions/compute_optimal_sample_size_no_overdispersion_greedy.R')
  source('./functions/compute_propagule_pressure_per_cons_no_overdispersion.R')
  source('./functions/compute_hypergeometric_sample_size_under_budget_constraint.R') 
  source('./functions/compute_approximate_hypergeometric_sample_size.R')
                   
################################################################################ 
#              Optimisation of sampling efforts within pathways
################################################################################ 
# Read data
  data_all_pathways <- fread('./data/data_all_pathways.csv')

# Load models  
  load('./outputs/fit_beta_binomial_hierarchical.Rdata')

# Extract coefficients of the beta-binomial and ZIBB and merge with data
  data_coef_summary_bb <- 
    extract_bb_parameters_per_pathway_from_brms_fit(
      model = fit_beta_binomial_hierarchical
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
  
# Compute sigma log(N) and mean N per pathway  
  data_fresh_with_coef <- data_all_with_coef[data_type %in% c('fresh_produce'),]
  data_fresh_with_coef[, sd_logN := sd(log(N)), by = group_id]
  data_fresh_with_coef[, mean_N := mean(N), by = group_id]

# Compute sample size per consignment per pathway for different sampling regimes
  data_fresh_with_coef[, optimal_n := 
    compute_approximate_optimal_sample_size(
      a = a, 
      b = b, 
      N = N, 
      B = .N * 600
    ), 
    by = .(group_id)]
  data_fresh_with_coef[, hypergeometric_n := 
    compute_hypergeometric_sample_size_under_budget_constraint(
      B = .N * 600, 
      p = 0.005, 
      N = N, 
      S0 = 0.95,                     
      S_increment_resolution = 0.00001
    )[, 'hypergeometric_n'], 
    by = group_id]
  data_fresh_with_coef[, prop_n := 600 * .N * N / sum(N), by = .(group_id)]
  data_fresh_with_coef[, optimal_n_no_overdisp :=         # Take longest
    compute_optimal_sample_size_no_overdispersion_greedy(
      p = a / (a + b), 
      N = N, 
      B = .N * 600, 
      n_increment_resolution = 1
    ), 
    by = .(group_id)]
     
# Compute propagule pressure per pathway
  data_fresh_with_coef[, propagule_pressure_optimal_n := 
    compute_propagule_pressure_per_cons(a = a, b = b, n = optimal_n, N = N)]
  data_fresh_with_coef[, propagule_pressure_fixed_n := 
    compute_propagule_pressure_per_cons(a = a, b = b, n = 600, N = N)]
  data_fresh_with_coef[, propagule_pressure_prop_n := 
    compute_propagule_pressure_per_cons(a = a, b = b, n = prop_n, N = N)]
  data_fresh_with_coef[, propagule_pressure_optimal_n_no_overdisp := 
    compute_propagule_pressure_per_cons(a = a, 
                                        b = b, 
                                        n = optimal_n_no_overdisp, 
                                        N = N)]
  data_fresh_with_coef[, propagule_pressure_hypergeometric_n := 
    compute_propagule_pressure_per_cons(a = a, 
                                        b = b, 
                                        n = hypergeometric_n, 
                                        N = N)]
  fwrite(data_fresh_with_coef, 
         './outputs/data_fresh_with_coef_with_added_leakage.csv')

# Compute propagule pressure summary per pathway
  data_fresh_with_coef <- 
    fread('./outputs/data_fresh_with_coef_with_added_leakage.csv')
  
  data_propagule_pressure_summary <- data_fresh_with_coef[, 
    list(
      propagule_pressure_optimal_n = sum(propagule_pressure_optimal_n), 
      propagule_pressure_fixed_n = sum(propagule_pressure_fixed_n),
      propagule_pressure_hypergeometric_n = 
        sum(propagule_pressure_hypergeometric_n),
      propagule_pressure_prop_n = sum(propagule_pressure_prop_n),
      propagule_pressure_optimal_n_no_overdisp = 
        sum(propagule_pressure_optimal_n_no_overdisp),
      a = mean(a),
      b = mean(b),
      sd_logN = mean(sd_logN),
      C = .N
    ), 
    by = .(data_type, group_id)]
   
# Create figure 06 comparing the ratio of propagule pressure under optimal vs. 
  # alternative sample size and how this ratio varies with pathway sd_logN
    data_sim <- expand.grid(sd_logN = seq(0, 1.25, length = 100), 
                            a = median(data_propagule_pressure_summary$a))
  # Expected relative propagule pressure reduction if consignment size is 
    # lognormally distributed (Eq. 5 in manuscript)
    data_sim$pred <- with(data_sim, exp( (-1/2 + 1 / (2 * (2 + a))) * sd_logN^2))
  
# Optimal sampling vs. fixed sample size
  p1 <- 
    ggplot(
      data = data_propagule_pressure_summary, 
      aes(sd_logN, propagule_pressure_optimal_n / propagule_pressure_fixed_n)
    ) + 
    geom_point(size = .8) + 
    geom_line(data = data_sim, aes(sd_logN, pred), col = 'dodgerblue') + 
    scale_x_continuous(breaks = seq(0, 2, by = 0.25)) +
    scale_y_continuous(breaks = seq(0, 1, by = 0.2), limits = c(0, 1)) + 
    xlab(TeX('$\\sigma_{log(N)}$')) +
    ylab('Relative propagule \n pressure reduction') +
    facet_wrap(~'Optimal n vs. fixed n') +
    annotate(
      'text', 
      x = 0.8, 
      y = 0.60, 
      label = TeX('$exp(-0.26 x^2)$', output = 'character'), 
      parse = T, 
      col = 'dodgerblue'
    ) +
    theme_bw() +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = c(.7, .2), 
      plot.margin = unit(c(.1, .1, .1, .1), "cm")
    )
  p1
# Optimal sampling vs. hypergeometric sampling
  p2 <- 
    ggplot(
      data = data_propagule_pressure_summary, 
      aes(
        sd_logN, 
        propagule_pressure_optimal_n / propagule_pressure_hypergeometric_n
      )
    ) + 
    geom_point(size = .8) + 
    geom_line(data = data_sim, aes(sd_logN, pred), col = 'dodgerblue') + 
    scale_x_continuous(breaks = seq(0, 2, by = 0.25)) +
    scale_y_continuous(breaks = seq(0, 1, by = 0.2), limits = c(0, 1)) + 
    xlab(TeX('$\\sigma_{log(N)}$')) +
    ylab('') +
    facet_wrap(~'Optimal vs. hypergeometric n') +
    theme_bw() +
    theme(
      panel.grid.minor = element_blank(),
      axis.title.y = element_blank(), 
      axis.ticks.y = element_blank(), 
      axis.text.y = element_blank(),
      plot.margin = unit(c(.1, .1, .1, .1), "cm")
    )
  p2
# Optimal sampling vs. sampling proportioal to consignment size
  p3 <- 
    ggplot(
      data = data_propagule_pressure_summary, 
      aes(sd_logN, propagule_pressure_optimal_n / propagule_pressure_prop_n)
    ) + 
    geom_point(size = .8) + 
    geom_line(data = data_sim, aes(sd_logN, pred), col = 'dodgerblue') + 
    scale_x_continuous(breaks = seq(0, 2, by = 0.25)) +
    scale_y_continuous(breaks = seq(0, 1, by = 0.2), limits = c(0, 1)) + 
    xlab(TeX('$\\sigma_{log(N)}$')) +
    ylab('') +
    facet_wrap(~'vs. n proportional to N') +
    theme_bw() +
    theme(
      panel.grid.minor = element_blank(),
      axis.title.y = element_blank(), 
      axis.ticks.y = element_blank(), 
      axis.text.y = element_blank(),
      plot.margin = unit(c(.1, .1, .1, .1), "cm")
    )
  p3
# Optimal sampling vs. optimal sampling asuming no overdispsersion
  p4 <- 
    ggplot(
      data = data_propagule_pressure_summary, 
      aes(
        sd_logN, 
        propagule_pressure_optimal_n / propagule_pressure_optimal_n_no_overdisp
      )
    ) + 
    geom_point(size = .8) + 
    geom_line(data = data_sim, aes(sd_logN, pred), col = 'dodgerblue') + 
    scale_x_continuous(breaks = seq(0, 2, by = 0.25)) +
    scale_y_continuous(breaks = seq(0, 1, by = 0.2), limits = c(0, 1)) + 
    xlab(TeX('$\\sigma_{log(N)}$')) +
    ylab('') +
    facet_wrap(~'vs. optimal n assuming no-overdisp') +
    theme_bw() +
    theme(
      panel.grid.minor = element_blank(),
      axis.title.y = element_blank(), 
      axis.ticks.y = element_blank(), 
      axis.text.y = element_blank(),
      plot.margin = unit(c(.1, .1, .1, .1), "cm")
    )
  p4
  
# Combine plots and save figure 06 
  pdf('./figs/fig06_relative_propagule_pressure_reduction_vs_sd_logN.pdf', 
      height = 6/2.54, 
      width = 24/2.54)
    plot_grid(p1, p2, p3, p4, ncol = 4, rel_widths = c(.374, .3, .3, .3))
  dev.off() 