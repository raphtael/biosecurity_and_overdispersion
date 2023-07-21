## Function extracting alpha and beta parameters per pathway from 
 # a hierarchical beta binomial distribution based on brms fit
 # Inputs are model fit. Outputs are alpha and beta parameters 
 # mean and 95% credible intervals per pathway.
  extract_bb_parameters_per_pathway_from_brms_fit <- 
    function(model){
      # Alpha
        data_log_alpha <- coef(model)$`group_id`[, , 'Intercept']
        group_id_name <- rownames(data_log_alpha)
        data_log_alpha <- as.data.table(data_log_alpha)
        data_log_alpha[, group_id := group_id_name]
        data_log_alpha[, alpha_mean := exp(Estimate) * exp(Est.Error^2/2) ]
        data_log_alpha[, alpha_min := exp(Q2.5) ]
        data_log_alpha[, alpha_max := exp(Q97.5) ]
      # Beta
        data_log_beta <- coef(model)$`group_id`[, , 'beta_Intercept'] 
        data_log_beta <- as.data.table(data_log_beta)
        data_log_beta[, group_id := group_id_name]
        data_log_beta[, beta_mean := exp(Estimate) * exp(Est.Error^2/2) ]
        data_log_beta[, beta_min := exp(Q2.5) ]
        data_log_beta[, beta_max := exp(Q97.5) ]
      # Merge datasets together
        data_brms_summary <- 
          cbind.data.frame(
            data_log_alpha[, .(group_id, alpha_mean, alpha_min, alpha_max)], 
            data_log_beta[, .(beta_mean, beta_min, beta_max)]
          )
      # Return summary statistics
        data_brms_summary <- as.data.table(data_brms_summary)
        data_brms_summary
    } 

   