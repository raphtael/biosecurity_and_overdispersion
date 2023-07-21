################################################################################
#
#     Write a table with parameter estimates associated with each pathway
#                                                    
################################################################################

# Set working directory
  setwd("C:/Users/rtrouve/Dropbox/overdispersion_paper")

# Load packages
  library(data.table)
  library(brms)
  source('./functions/extract_bb_parameters_per_pathway_from_brms_fit.R')
  source('./functions/extract_zibb_parameters_per_pathway_from_brms_fit.R')
  source('./functions/inv_logit.R')
  
################################################################################
# Read data
  data_all_pathways <- fread('./data/data_all_pathways.csv')
  data_summary_per_pathway <- data_all_pathways[, 
    list(
      C = .N, 
      median_n = median(n)
    ), 
    by = list(data_type, group_id)]
    
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

# Merge coefficient summary tables for different models  
  data_coef_summary <- 
    merge(
      data_summary_per_pathway, 
      data_coef_summary_bb, 
      by = 'group_id'
    )
   
  data_coef_summary <- 
    merge(
      data_coef_summary, 
      data_coef_summary_zibb[, 
        list(
          group_id, 
          a_zibb = alpha_mean, 
          b_zibb = beta_mean, 
          z_zibb = zi_mean
        )], 
      by = 'group_id'
    )
  
  data_coef_summary <- 
    merge(
      data_coef_summary, 
      data_coef_summary_zibb_fixed_a[, 
        list(
          group_id, 
          a_zibb_fixed_a = alpha_mean, 
          b_zibb_fixed_a = beta_mean, 
          z_zibb_fixed_a = zi_mean
        )], 
      by = 'group_id'
    )
  
  data_coef_summary[, p_mean_pred_perc :=  # Compute mean infestation rate
    alpha_mean / (alpha_mean + beta_mean) * 100] 
  data_coef_summary
 
# Format data
  names_character <- 
    colnames(data_coef_summary)[-which(colnames(data_coef_summary) %in% 
      c("group_id", "C", "median_n", "beta_mean", "b_zibb", "b_zibb_fixed_a"))]
    for(col in names_character) 
    set(data_coef_summary, 
        j = col, 
        value = 
          formatC(
            data_coef_summary[[col]], 
            dig = 2, 
            format = "f", 
            width = 5
          )
    )
  data_coef_summary[, C := 
    formatC(C, dig = 0, format = "f", width = 3) ]
  data_coef_summary[, median_n := 
    formatC(median_n, dig = 0, format = "f", width = 5) ]
  data_coef_summary[, beta_mean := 
    formatC(beta_mean, dig = 1, format = "f", width = 4) ]
  data_coef_summary[, b_zibb := 
    formatC(b_zibb, dig = 1, format = "f", width = 4) ]
  data_coef_summary[, b_zibb_fixed_a := 
    formatC(b_zibb_fixed_a, dig = 1, format = "f", width = 4) ]
  
# Write table for Appendix
  # csv format   
    fwrite(data_coef_summary, './outputs/table_pathway_parameters.csv')               

               