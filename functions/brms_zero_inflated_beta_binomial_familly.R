################################################################################
#
#   Define a custom zero-inflated beta-binomial distribution for use in brms
#
################################################################################
# The procedure to define custom response distributions in brms is explained here:
# https://cran.r-project.org/web/packages/brms/vignettes/brms_customfamilies.html

# We also defined functions for posterior predictive check and model comparisons

# mu represents the alpha parameter
# beta is the beta parameter
# zi is the zero-inflated probability parameter
# vint1 is the sample size of the inspection
zero_inflated_beta_binomial <- custom_family(
  "zero_inflated_beta_binomial", dpars = c("mu", "beta", "zi"),
  links = c("log", "log", "logit"), lb = c(0, 0, 0), ub = c(NA, NA, 1),
  type = "int", vars = "vint1[n]"   
)

stan_funs_zibb <- "
 real zero_inflated_beta_binomial_lpmf(int y, real mu, real beta, real zi, int T) { 
    if (y == 0) { 
      return log_sum_exp(bernoulli_lpmf(1 | zi), 
                         bernoulli_lpmf(0 | zi) + 
                         beta_binomial_lpmf(y | T, mu, beta)); 
    } else { 
      return bernoulli_lpmf(0 | zi) +  
             beta_binomial_lpmf(y | T, mu, beta); 
    } 
  } 
"

stanvars_zibb <- stanvar(scode = stan_funs_zibb, block = "functions")

log_lik_zero_inflated_beta_binomial <- function(i, prep) {
  zi <- get_dpar(prep, "zi", i = i)
  mu <- brms::get_dpar(prep, "mu", i = i)
  beta <- brms::get_dpar(prep, "beta", i = i)
  trials <- prep$data$vint1[i]
  y <- prep$data$Y[i]
  zero_inflated_beta_binomial_lpmf(y, mu, beta, zi, trials)
}
 
posterior_predict_zero_inflated_beta_binomial <- function(i, prep, ...) {
  zi <- get_dpar(prep, "zi", i = i)
  mu <- get_dpar(prep, "mu", i = i)
  beta <- get_dpar(prep, "beta", i = i)
  trials <- prep$data$vint1[i]
  zi_sample <- runif(prep$ndraws, 0, 1)
  ifelse(
    zi_sample < zi, 0, 
    emdbook::rbetabinom(n = prep$ndraws, shape1 = mu, shape2 = beta, size = trials)
  )
}