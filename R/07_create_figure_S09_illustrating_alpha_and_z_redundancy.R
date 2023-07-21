################################################################################
#
# Illustrate the redundancy between the alpha and z parameters (Fig. S9)
#                    
################################################################################

# Set working directory
  setwd("C:/Users/rtrouve/Dropbox/overdispersion_paper")

# Load packages
  library(ggplot2)
  library(data.table)
  library(brms)
  library(latex2exp)
  source('./functions/extract_bb_parameters_per_pathway_from_brms_fit.R')  
  source('./functions/extract_zibb_parameters_per_pathway_from_brms_fit.R')  
  source('./functions/inv_logit.R')  
  
################################################################################
# Load models 
   load('./outputs/fit_beta_binomial_hierarchical.Rdata')  
   load('./outputs/fit_zibb_hierarchical.Rdata')  
   load('./outputs/fit_zibb_hierarchical_fixed_a.Rdata')
   
# Extract coefficients of the beta-binomial and zero-inflated BB distributions
  data_coef_summary_bb <- 
    extract_bb_parameters_per_pathway_from_brms_fit(
      model = fit_beta_binomial_hierarchical
    )
  data_coef_summary_zibb <- 
    extract_zibb_parameters_per_pathway_from_brms_fit(
      model = fit_zibb_hierarchical
    )
  data_coef_summary_zibb_fixed_a <- 
    extract_zibb_parameters_per_pathway_from_brms_fit(
      model = fit_zibb_hierarchical_fixed_a
    )

# Merge datasets
  data_coef_summary <- merge(
    data_coef_summary_bb[, list(
      group_id,
      a_bb = alpha_mean,
      b_bb = beta_mean
    )],
    data_coef_summary_zibb[, list(
      group_id, 
      a_zibb = alpha_mean, 
      b_zibb = beta_mean, 
      z_zibb = zi_mean
    )],
    by = 'group_id'
  )

  data_coef_summary <- merge(
    data_coef_summary_zibb_fixed_a[, list(
      group_id, 
      a_zibb_fixed_a = alpha_mean, 
      b_zibb_fixed_a = beta_mean, 
      z_zibb_fixed_a = zi_mean
    )], 
    data_coef_summary, 
    by = 'group_id'
  )
  
# Plot figure
  p1 <- ggplot(data = data_coef_summary, aes(a_bb, z_zibb_fixed_a)) + 
    geom_point() +
    xlab(TeX('$\\alpha_j$ per pathway in the Beta-Binomial model')) +
    ylab(TeX('$z_j$ per pathway in the ZIBB model with fixed $\\alpha$')) + 
    theme_bw() + 
    theme(panel.grid.minor=element_blank())
  ggsave(
    plot = p1, 
    filename = './figs/a_vs_z_in_BB_vs_ZIBB_with_fixed_a_models.pdf', 
    height = 10, 
    width = 10, 
    units = 'cm'
  )
    
