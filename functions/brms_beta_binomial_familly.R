################################################################################
#
#        Define a custom beta-binomial distribution for use in brms
#
################################################################################
# The procedure to define custom response distributions in brms is explained here:
# https://cran.r-project.org/web/packages/brms/vignettes/brms_customfamilies.html

# The distribution that we defined is parameterized directly in terms of alpha 
# and beta instead of mean and dispersion as is usual in brms
# We also defined functions for posterior predictive check and model comparisons


# mu represents the alpha parameter
# beta is the beta parameter
# vint1 is the sample size of the inspection
beta_binomial2 <- custom_family(
  "beta_binomial2", dpars = c("mu", "beta"), 
  links = c("log", "log"), lb = c(0, 0), ub = c(NA, NA),
  type = "int", vars = "vint1[n]"   
)

stan_funs <- "
  real beta_binomial2_lpmf(int y, real mu, real beta, int T) {
    return beta_binomial_lpmf(y | T, mu, beta);
  }
  int beta_binomial2_rng(real mu, real beta, int T) {
    return beta_binomial_rng(T, mu, beta);
  }
"

stanvars <- stanvar(scode = stan_funs, block = "functions")

posterior_predict_beta_binomial2 <- function(i, prep, ...) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  beta <- brms::get_dpar(prep, "beta", i = i)
  trials <- prep$data$vint1[i]
  beta_binomial2_rng(mu, beta, trials)
}
 
log_lik_beta_binomial2 <- function(i, prep) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  beta <- brms::get_dpar(prep, "beta", i = i)
  trials <- prep$data$vint1[i]
  y <- prep$data$Y[i]
  beta_binomial2_lpmf(y, mu, beta, trials)
}
 