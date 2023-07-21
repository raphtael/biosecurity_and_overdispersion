################################################################################
#
#         Posterior predictive check for the hierarchical models
#
################################################################################

# set working directory
  setwd("C:/Users/rtrouve/Dropbox//overdispersion_paper")

# load packages
  library(ggplot2)
  library(data.table)
  library(brms)
  library(RColorBrewer)
  library(latex2exp)
  library(cowplot)
  source('./functions/brms_beta_binomial_familly.R')
  source('./functions/brms_zero_inflated_beta_binomial_familly.R')
  source('./functions/extract_bb_parameters_per_pathway_from_brms_fit.R')  
  source('./functions/inv_logit.R')
  source('./functions/compute_prop_zero.R')    
  tick_break_log_scale <- c(
    0.00001, 0.00003, 0.0001, 0.0003, 0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 1, 3, 
    10, 30, 100, 300, 1000, 3000, 10000, 30000, 10^5
  )
################################################################################
# Read data
  data_all_pathways <- fread('./data/data_all_pathways.csv')
  data_all_pathways[, row_id := 1:.N]   

################################################################################
# pp_check for the beta-binomial model
# Load model
  load('./outputs/fit_beta_binomial_hierarchical.Rdata')
  expose_functions(fit_beta_binomial_hierarchical, vectorize = TRUE)
  
# pp_check on the empirical cumulative distribution function (ecdf)
  # Posterior sampling and merge with dataset
    n_mcmc_reps <- 50   
    data_sim <- posterior_predict(
      fit_beta_binomial_hierarchical, 
      ndraws = n_mcmc_reps)
    data_sim <- t(data_sim)
    data_sim <- as.data.table(data_sim)
    colnames(data_sim) <- as.character(1:ncol(data_sim))
    data_sim[, row_id := 1:.N]
    data_sim_long <- melt(data_sim, 
                          id.vars = 'row_id', 
                          variable.name = 'mcmc_reps_id')
    data_sim_long <- merge(
      data_sim_long, 
      data_all_pathways[, .(row_id, data_type, group_id, n)], 
      by = 'row_id') 
  # Create figure  
    p1 <- 
      ggplot(
        data = data_sim_long, 
        aes(value / n, group = mcmc_reps_id, 
        col = data_type)
      ) + 
      stat_ecdf(size=.1) + 
      stat_ecdf(
        data = data_all_pathways, 
        aes(k / n, group = NULL), 
        size=.2, 
        col = 'black'
      ) +
      facet_wrap(~group_id, scales = 'free_x', ncol = 7) +
      scale_colour_brewer(palette = 'Pastel2', name = '') +
      theme_bw() +
      theme(
        panel.grid.minor = element_blank(),
        axis.title = element_blank(),
        legend.position = 'bottom', 
        strip.text.x = element_text(margin = margin(.08, 0, .08, 0, "cm"))
      ) +
      xlab('Infestation rate pj among different consignments of the pathway') +
      ylab('Empirical cumulative distribution function of pj')
  
  # Save figure
    ggsave(
      plot = p1, 
      filename = './figs/pp_check_ecdf_bb.pdf', 
      width = 22, 
      height = 24, 
      units = 'cm'
    )
  
# Posterior predictive check on the proportion of zeroes per pathway
  # Posterior sampling
    n_mcmc_reps <- 100   
    data_sim <- posterior_predict(
      fit_beta_binomial_hierarchical, 
      ndraws = n_mcmc_reps
    )  
    data_sim <- t(data_sim)
    data_sim <- as.data.table(data_sim)
    colnames(data_sim) <- as.character(1:ncol(data_sim))
    data_sim[, row_id := 1:.N]
  
  # Get long format and merge it with pathway details
    data_sim_long <- melt(
      data_sim, 
      id.vars = 'row_id', 
      variable.name = 'mcmc_reps_id'
    )
    data_sim_long <- 
      merge(
        data_sim_long, 
        data_all_pathways[, .(row_id, data_type, group_id, k, n)], 
        by = 'row_id'
      ) 
  # Compute summary
    data_sim_summary <- data_sim_long[, 
      list(
        prop_zero_pred = compute_prop_zero(value),
        prop_zero_obs = compute_prop_zero(k)
      ), 
      by = list(data_type, group_id, mcmc_reps_id)
    ]
  
  # Plot prop zero figures 
    p2 <- ggplot(data = data_sim_summary, aes(prop_zero_pred, fill = data_type)) + 
      geom_histogram() + 
      geom_vline(aes(xintercept = prop_zero_obs)) + 
      facet_wrap(~group_id, ncol = 7) + 
      scale_fill_brewer(palette = 'Pastel2', name = '') +
      theme_bw() +
      theme(
        panel.grid.minor = element_blank(),
        legend.position = 'bottom', 
        strip.text.x = element_text(margin = margin(.08, 0, .08, 0, "cm")),
        axis.title = element_blank()
      ) +
      xlim(c(0, 1)) +
      xlab('Proportion of compliant consignments')
    ggsave(plot = p2, filename = './figs/pp_check_prop_compliant_bb.pdf', 
           width = 22, height = 24, units = 'cm')

# Check prop. of zeroes and mean infestation rate per pathway obs. vs. predicted
  # Extract parameter values per pathways and merge with data
    data_coef_summary <- 
      extract_bb_parameters_per_pathway_from_brms_fit(
        model = fit_beta_binomial_hierarchical
      )
    data_all_pathways_summary <- data_all_pathways[, 
      list(
        prop_zero_obs = compute_prop_zero(k), 
        p_mean_obs = mean(k/n), 
        n_mean = mean(n)
      ), 
      by = .(data_type, group_id)]
    data_coef_summary <- merge(data_coef_summary, 
                               data_all_pathways_summary, 
                               by = 'group_id')
  # Compute predicted mean infestation rate and proportion of zero per pathways
    data_coef_summary[, p_mean_pred := 
      alpha_mean / (alpha_mean + beta_mean)] # Mean infestation rate formula
    data_coef_summary[, prop_zero_pred :=    # Pr(X = 0 | alpha, beta, n)
      exp(lbeta(alpha_mean, beta_mean + n_mean) - lbeta(alpha_mean, beta_mean)), 
      by = group_id]
 
  # Compute p_min and p_max
    p_vector <- c(data_coef_summary$p_mean_pred, data_coef_summary$p_mean_obs)
    p_min <- floor(min(p_vector) * 10^6) / 10^4
    p_max <- ceiling(max(p_vector) * 10^6) / 10^4 
  
  # Plot figures 
    p3 <- 
      ggplot(
        data = data_coef_summary, 
        aes(p_mean_pred * 100, p_mean_obs * 100, col = data_type)
      ) + 
      geom_point() + 
      geom_abline() +
      scale_colour_brewer(palette = 'Set2', name = '') +
      ylab(TeX('$p_{obs}$ (\\%)')) +
      xlab(TeX('$p_{pred}$ (\\%)')) +
      scale_x_log10(breaks = tick_break_log_scale, limits = c(p_min, p_max)) + 
      scale_y_log10(breaks = tick_break_log_scale, limits = c(p_min, p_max)) +
      theme_bw() +
      theme(
        legend.position = 'none',
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)
      )
    p4 <- 
      ggplot(
        data = data_coef_summary, 
        aes(prop_zero_pred, prop_zero_obs, col = data_type)
      ) + 
      geom_point() + 
      geom_abline() +
      scale_colour_brewer(palette = 'Set2', name = '') +
      ylab('Proportion of compliant consignments obs.') +
      xlab('Proportion of compliant consignments pred.') +
      xlim(0, 1) +
      ylim(0, 1) +
      theme_bw() +
      theme(
        legend.position = c(.75, .25), 
        legend.background = element_blank(),
        panel.grid.minor = element_blank(),
      )
  # Save figure  
    pdf(file = './figs/p_mean_and_prop_zero_obs_vs_pred_bb.pdf', 
        width = 20/2.54, 
        height = 10/2.54)
      plot_grid(p3, p4, ncol = 2, align = 'hv')
    dev.off()  
     
################################################################################
# Posterior predictive check for the zero-inflated BB model 
# Load model  
  load('./outputs/fit_zibb_hierarchical.Rdata') 
  expose_functions(fit_zibb_hierarchical, vectorize = TRUE)
  
# pp_check on the empirical cumulative distribution function (ecdf)
  # Posterior sampling and merge with dataset
    n_mcmc_reps <- 50   
    data_sim <- posterior_predict(fit_zibb_hierarchical, ndraws = n_mcmc_reps)
    data_sim <- t(data_sim)
    data_sim <- as.data.table(data_sim)
    colnames(data_sim) <- as.character(1:ncol(data_sim))
    data_sim[, row_id := 1:.N]
    data_sim_long <- melt(data_sim, 
                          id.vars = 'row_id', 
                          variable.name = 'mcmc_reps_id')
    data_sim_long <- merge(
      data_sim_long, 
      data_all_pathways[, .(row_id, data_type, group_id, n)], 
      by = 'row_id') 
  # Create figure  
    p5 <- 
      ggplot(
        data = data_sim_long, 
        aes(value / n, group = mcmc_reps_id, 
        col = data_type)
      ) + 
      stat_ecdf(size=.1) + 
      stat_ecdf(
        data = data_all_pathways, 
        aes(k / n, group = NULL), 
        size=.2, 
        col = 'black'
      ) +
      facet_wrap(~group_id, scales = 'free_x', ncol = 7) +
      scale_colour_brewer(palette = 'Pastel2', name = '') +
      theme_bw() +
      theme(
        panel.grid.minor = element_blank(),
        axis.title = element_blank(),
        legend.position = 'bottom', 
        strip.text.x = element_text(margin = margin(.08, 0, .08, 0, "cm"))
      ) +
      xlab('Infestation rate pj among different consignments of the pathway') +
      ylab('Empirical cumulative distribution function of pj')
  
  # Save figure
    ggsave(
      plot = p5, 
      filename = './figs/pp_check_ecdf_zibb.pdf', 
      width = 22, 
      height = 24, 
      units = 'cm'
    ) 
 
################################################################################
# Posterior predictive check for the zero-inflated BB with fixed alpha model 
# Load model  
  load('./outputs/fit_zibb_hierarchical_fixed_a.Rdata') 
  expose_functions(fit_zibb_hierarchical_fixed_a, vectorize = TRUE)
  
# pp_check on the empirical cumulative distribution function (ecdf)
  # Posterior sampling and merge with dataset
    n_mcmc_reps <- 50   
    data_sim <- posterior_predict(
      fit_zibb_hierarchical_fixed_a, 
      ndraws = n_mcmc_reps
    )
    data_sim <- t(data_sim)
    data_sim <- as.data.table(data_sim)
    colnames(data_sim) <- as.character(1:ncol(data_sim))
    data_sim[, row_id := 1:.N]
    data_sim_long <- melt(data_sim, 
                          id.vars = 'row_id', 
                          variable.name = 'mcmc_reps_id')
    data_sim_long <- merge(
      data_sim_long, 
      data_all_pathways[, .(row_id, data_type, group_id, n)], 
      by = 'row_id') 
  # Create figure  
    p6 <- 
      ggplot(
        data = data_sim_long, 
        aes(value / n, group = mcmc_reps_id, 
        col = data_type)
      ) + 
      stat_ecdf(size = .1) + 
      stat_ecdf(
        data = data_all_pathways, 
        aes(k / n, group = NULL), 
        size=.2, 
        col = 'black'
      ) +
      facet_wrap(~group_id, scales = 'free_x', ncol = 7) +
      scale_colour_brewer(palette = 'Pastel2', name = '') +
      theme_bw() +
      theme(
        panel.grid.minor = element_blank(),
        axis.title = element_blank(),
        legend.position = 'bottom', 
        strip.text.x = element_text(margin = margin(.08, 0, .08, 0, "cm"))
      ) +
      xlab('Infestation rate pj among different consignments of the pathway') +
      ylab('Empirical cumulative distribution function of pj')
  
  # Save figure
    ggsave(
      plot = p6, 
      filename = './figs/pp_check_ecdf_zibb_fixed_a.pdf', 
      width = 22, 
      height = 24, 
      units = 'cm'
    )
                       