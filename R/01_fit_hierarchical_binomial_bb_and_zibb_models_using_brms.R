################################################################################
#
# Use a hierarchical models to estimate the alpha and beta parameters of the 
# Beta-Binomial distribution associated with each pathway.
# Also fit alternative Binomial, and zero-inflated Beta-Binomial models      
#                                                     
################################################################################

# Set working directory
  setwd("C:/Users/rtrouve/Dropbox/overdispersion_paper")

# Load packages
  library(data.table)
  library(brms)
  source('./functions/brms_beta_binomial_familly.R')
  source('./functions/brms_zero_inflated_beta_binomial_familly.R')

################################################################################
# Read data
  data_all_pathways <- fread('./data/data_all_pathways.csv')

# Set priors for hierarchical models
  my_priors_binomial_hierarchical <- 
    c(
      prior(prior = constant(7), class = 'df', group = 'group_id'),
      prior(prior = normal(0, 2), class = 'sd')
    )
  my_priors_beta_binomial_hierarchical <- 
    c(
      prior(prior = constant(7), class = 'df', group = 'group_id'),
      prior(prior = normal(log(0.3) - 1^2/2, 1), class = 'Intercept'),
      prior(prior = student_t(3, log(50) - 2^2/2, 2), class = 'Intercept', dpar = 'beta'),
      prior(prior = normal(0, 2), class = 'sd'), 
      prior(prior = normal(0, 2), class = 'sd', dpar ='beta')
    )
  my_priors_zibb_hierarchical <- 
    c(
      prior(prior = constant(7), class = 'df', group = 'group_id'),
      prior(prior = normal(log(0.3) - 1^2/2, 1), class = 'Intercept'),
      prior(prior = student_t(3, log(50) - 2^2/2, 2), class = 'Intercept', dpar = 'beta'),
      prior(prior = normal(0, 2), class = 'sd'), 
      prior(prior = normal(0, 2), class = 'sd', dpar ='beta'), 
      prior(prior = normal(0, 2), class = 'sd', dpar ='zi')
    )
  my_priors_zibb_hierarchical_fixed_a <- 
    c(
      prior(prior = constant(7), class = 'df', group = 'group_id'),
      prior(prior = normal(log(0.3) - 1^2/2, 1), class = 'Intercept'),
      prior(prior = student_t(3, log(50) - 2^2/2, 2), class = 'Intercept', dpar = 'beta'),
      prior(prior = normal(0, 2), class = 'sd', dpar ='beta'), 
      prior(prior = normal(0, 2), class = 'sd', dpar ='zi')
    )                                          
   
# Fit hierarchical models
  fit_binomial_hierarchical <- 
    brm(
      bf(k | trials(n) ~ 1 + (1|gr(group_id, dist = "student"))),
      data = data_all_pathways,
      family = binomial, 
      prior = my_priors_binomial_hierarchical,
      chain = 1,
      iter = 200,
      init = 0,
      control = list(adapt_delta = 0.99)
    )
  save(fit_binomial_hierarchical, 
       file = './outputs/fit_binomial_hierarchical.Rdata')
  
  fit_beta_binomial_hierarchical <- 
    brm(
      bf(k | vint(n) ~ 1 + (1|group_id),
      beta ~ 1 + (1|gr(group_id, dist = "student"))),
      data = data_all_pathways,
      family = beta_binomial2, 
      stanvars = stanvars,
      prior = my_priors_beta_binomial_hierarchical,
      chain = 1,
      iter = 2000,
      init = 0,
      control = list(adapt_delta = 0.99)
    )
  save(fit_beta_binomial_hierarchical, 
       file = './outputs/fit_beta_binomial_hierarchical.Rdata')
  
  fit_zibb_hierarchical <- 
    brm(
      bf(k | vint(n) ~ 1 + (1|group_id),
      beta ~ 1 + (1|gr(group_id, dist = "student")), 
      zi ~ 1 + (1|gr(group_id, dist = "student"))),
      data = data_all_pathways,
      family = zero_inflated_beta_binomial, 
      stanvars = stanvars_zibb,
      prior = my_priors_zibb_hierarchical,
      chain = 2,
      iter = 2000,
      init = 0,
      control = list(adapt_delta = 0.99)
    )
   save(fit_zibb_hierarchical, 
        file = './outputs/fit_zibb_hierarchical.Rdata')
   
   fit_zibb_hierarchical_fixed_a <- 
     brm(
       bf(k | vint(n) ~ 1,
       beta ~ 1 + (1|gr(group_id, dist = "student")), 
       zi ~ 1 + (1|gr(group_id, dist = "student"))),
       data = data_all_pathways,
       family = zero_inflated_beta_binomial, 
       stanvars = stanvars_zibb,
       prior = my_priors_zibb_hierarchical_fixed_a,
       chain = 2,
       iter = 2000,
       init = 0,
       control = list(adapt_delta = 0.99)
     )
   save(fit_zibb_hierarchical_fixed_a, 
        file = './outputs/fit_zibb_hierarchical_fixed_a.Rdata')
    