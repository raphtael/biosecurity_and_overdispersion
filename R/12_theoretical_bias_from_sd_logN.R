################################################################################
#
# Simulated and analytical bias (reduction) in average propagule pressure per 
# consignment when we increase SD log(N). Also test the influence of consignment 
# size distribution on this bias.
#                                               
################################################################################

# The bias should provide an analytical expression for the reduction in leakage 
# when using the optimal sample size formula 
# compared to using the same sample size for all consignments (this is similar to the bias in log transformed variables).

# Set working directory
  setwd('C:/Users/rtrouve/Dropbox/overdispersion_paper')

# Load libraries
  library(data.table)
  library(ggplot2)
  library(cowplot)
  library(latex2exp)
  source('./functions/compute_optimal_sample_size_no_overdispersion_greedy.R')
  source('./functions/compute_propagule_pressure_per_cons_no_overdispersion.R')
  source('./functions/compute_approximate_optimal_sample_size.R')
  source('./functions/compute_propagule_pressure_per_cons.R')
      
################################################################################
# Theoretical prop. pressure reduction for optimal sample size vs. fixed sample 
# size under log-normally distributed N

# Experimental design                                                  
  data_sim <- expand.grid(SD = seq(0, 4, length = 20), 
                          a = 10^(seq(-5, 2, by = 1)), 
                          b = 8, 
                          C = 10^5, 
                          B = 10^5 * 600, 
                          mean_N = 10^5)
  data_sim <- as.data.table(data_sim)
 
# Compute propagule pressure reduction for optimal vs. fixed sample size
  compute_rel_propagule_pressure_reduction <- function(a, b, C, B, mean_N, SD){
    N1 <- exp(rnorm(C, log(mean_N) - SD^2/2, SD))
    n1 <- compute_approximate_optimal_sample_size(a = a, b = b, N = N1, B = B)
    n2 <- B /C    # Fixed sample size
    propagule_pressure1 <- 
      sum(compute_propagule_pressure_per_cons(a = a, b = b, n = n1, N = N1))
    propagule_pressure2 <- 
      sum(compute_propagule_pressure_per_cons(a = a, b = b, n = n2, N = N1))
    return(propagule_pressure1 / propagule_pressure2) # Ratio
  }
  data_sim[, leak_ratio := 
    compute_rel_propagule_pressure_reduction(a = a, 
                                             b = b, 
                                             C = C, 
                                             B = B, 
                                             mean_N = mean_N, 
                                             SD = SD), 
    by = .(a, SD)]
  data_sim[, leak_ratio_pred := exp((- 1/2 + 1 / (2 * (a + 2))) * SD^2) ] 
  
# Plot
  p1 <- ggplot(data = data_sim, aes(SD, leak_ratio, group = a)) + 
    geom_point(col = 'red3') + 
    geom_line(aes(SD, leak_ratio_pred), col = 'black') + 
    facet_wrap(~a, label = 'label_both', ncol = 4) +
    theme_bw() + 
    theme(panel.grid.minor = element_blank()) +
    xlab(TeX('$\\sigma_{log(N)}$')) +
    ylab('Relative reduction in leakage')
  p1
  ggsave(
    plot = p1, 
    filename = './figs/theoretical_bias_vs_sd_logN.pdf', 
    width = 20, 
    height = 10, 
    units = 'cm'
  )    

################################################################################
# Theoretical prop. pressure reduction for optimal sample size vs. sample size
# proportional to consignment size under log-normally distributed N

# Experimental design                                                  
  data_sim <- expand.grid(SD = seq(0, 4, length = 20), 
                          a = 10^(seq(-5, 2, by = 1)), 
                          b = c(1, 10, 100), 
                          C = 10^5, 
                          B = 10^5 * 600, 
                          mean_N = 10^5)
  data_sim <- as.data.table(data_sim)

# Compute propagule pressure reduction for optimal vs. sample size prop. to N 
  compute_rel_propagule_pressure_reduction <- function(a, b, C, B, mean_N, SD){
    N1 <- exp(rnorm(C, log(mean_N) - SD^2/2, SD))
    N1 <- N1 * mean_N / mean(N1) # Get exact mean
    n1 <- compute_approximate_optimal_sample_size(a = a, b = b, N = N1, B = B)
    n2 <- N1 / sum(N1) * B
    propagule_pressure1 <- 
      sum(compute_propagule_pressure_per_cons(a = a, b = b, n = n1, N = N1))
    propagule_pressure2 <- 
      sum(compute_propagule_pressure_per_cons(a = a, b = b, n = n2, N = N1))
    return(propagule_pressure1 / propagule_pressure2)
  }
  data_sim[, leak_ratio := 
    compute_rel_propagule_pressure_reduction(a = a, 
                                             b = b, 
                                             C = C, 
                                             B = B, 
                                             mean_N = mean_N, 
                                             SD = SD), 
    by = .(a, b, SD)]
  data_sim[, leak_ratio_pred := exp((- 1/2 + 1 / (2 * (a + 2))) * SD^2)]
           
# Plot
  p2 <- ggplot(data = data_sim, aes(SD, leak_ratio, group = paste(a, b))) + 
    geom_point(aes(col = as.factor(b))) + 
    geom_line(aes(SD, leak_ratio_pred), col = 'black') + 
    facet_wrap(~a, label = 'label_both', ncol = 4) +
    theme_bw() + 
    theme(panel.grid.minor = element_blank()) + 
    xlab(TeX('$\\sigma_{log(N)}$')) +
    ylab('Relative reduction in leakage')+
    ylim(-0.01, 1.01) +
    guides(col = guide_legend(title = TeX("$\\beta =$")))
  p2
  ggsave(
    plot = p2, 
    filename = './figs/theoretical_bias_vs_sd_logN_sample_size_prop_to_N.pdf', 
    width = 20, 
    height = 10, 
    units = 'cm'
  )   

################################################################################
# Theoretical prop. pressure reduction for optimal n vs. optimal n assuming
# no-overdispersion under log-normally distributed N

# Experimental design.
  # Note, very slow. C = 10^5 can take weeks on a single computer. 
  # Split in 10^4 batches and merge at the end                                       
  data_sim <- expand.grid(SD = seq(0, 4, length = 20), 
                          a = 10^(seq(-3, -1, by = 1)), 
                          b = c(1, 10, 100), 
                          C = 10^2,  # Increase this number for reliable result
                          B = 10^2 * 600, 
                          mean_N = 10^5)
  data_sim <- as.data.table(data_sim)
  
# Compute propagule pressure reduction for optimal vs. optimal sample size 
# assuming no overdispsersion
  compute_leakage_reduction_non_overdisp <- function(a, b, C, B, mean_N, SD){
    N1 <- exp(rnorm(C, log(mean_N) - SD^2/2, SD))
    n1 <- compute_approximate_optimal_sample_size(a = a, b = b, N = N1, B = B)
    n2 <- 
      compute_optimal_sample_size_no_overdispersion_greedy(
        p = a / (a + b), 
        N = N1, 
        B = B, 
        n_increment_resolution = 10  # Use 1 for improved (but slower) results
      )
    propagule_pressure1 <- 
      sum(compute_propagule_pressure_per_cons(a = a, b = b, n = n1, N = N1))
    propagule_pressure2 <- 
      sum(compute_propagule_pressure_per_cons(a = a, b = b, n = n2, N = N1))
    return(propagule_pressure1 / propagule_pressure2)
  }
  data_sim[, leak_ratio := # Extremely slow (days) due to greedy algorithm
    compute_leakage_reduction_non_overdisp(a = a, 
                                           b = b, 
                                           C = C, 
                                           B = B, 
                                           mean_N = mean_N, 
                                           SD = SD), 
    by = .(a, b, SD)] 
  data_sim[, leak_ratio_pred := exp((- 1/2 + 1 / (2 * (a + 2))) * SD^2) ] 
     
# Plot
  p3 <- ggplot(data = data_sim, aes(SD, leak_ratio, group = paste(a, b))) + 
    geom_point(aes(col = as.factor(b))) + 
    geom_line(aes(SD, leak_ratio_pred), col = 'black') + 
    facet_wrap(~a, label = 'label_both', ncol = 4) +
    theme_bw() + 
    theme(panel.grid.minor = element_blank()) + # + write down expression
    xlab(TeX('$\\sigma_{log(N)}$')) +
    ylab('Relative reduction in leakage')+
    ylim(-0.01, 1.01) +
    guides(col = guide_legend(title = TeX("$\\beta =$")))
  p3
  ggsave(
    plot = p3, 
    filename = './figs/theoretical_bias_vs_sd_logN_optimal_sample_size_assuming_no_overdispersion.pdf', 
    width = 20, 
    height = 10, 
    units = 'cm'
  )    


################################################################################
# Testing alternatives to lognormal distribution for N
# lognormal distribution  
  data_sim <- expand.grid(SD = seq(0, 4, length = 20), 
                          a = 0.1, 
                          b = 8, 
                          C = 10^6, 
                          B = 10^6 * 600, 
                          mean_N = 10^5)
  data_sim <- as.data.table(data_sim) 
# Compute relative propagule pressure reduction  
  compute_rel_propagule_pressure_reduction <- function(a, b, C, B, mean_N, SD){
    N1 <- exp(rnorm(n = C, log(mean_N) - SD^2/2, SD))   
    N1 <- N1 * mean_N / mean(N1)
    n1 <- compute_approximate_optimal_sample_size(a = a, b = b, N = N1, B = B)
    n2 <- B / C
    propagule_pressure1 <- 
      sum(compute_propagule_pressure_per_cons(a = a, b = b, n = n1, N = N1))
    propagule_pressure2 <- 
      sum(compute_propagule_pressure_per_cons(a = a, b = b, n = n2, N = N1))
    return(propagule_pressure1 / propagule_pressure2)
  }
  data_sim[, leak_ratio := 
    compute_rel_propagule_pressure_reduction(a = a, 
                                             b = b, 
                                             C = C, 
                                             B = B, 
                                             mean_N = mean_N, 
                                             SD = SD), 
    by = .(a, SD)]
  data_sim[, leak_ratio_pred := exp((- 1/2 + 1 / (2 * (a + 2))) * SD^2)]
# Plot
  p4 <- ggplot(data = data_sim, aes(SD, leak_ratio, group = a)) + 
    geom_point(col = 'red3') + 
    geom_line(aes(SD, leak_ratio_pred), col = 'black') + 
    facet_wrap(~'log(N_i) ~ Normal') +
    theme_bw() + 
    theme(panel.grid.minor = element_blank()) + # + write down expression
    xlab(TeX('$\\sigma_{log(N)}$')) +
    ylab('Relative reduction in leakage') +
   theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), 'cm'))
  p4
  
# Bimodal distribution.
  data_sim <- expand.grid(SD = seq(0, 4, length = 20), 
                          a = 0.1, 
                          b = 8, 
                          C = 10^5, 
                          B = 10^5 * 600, 
                          mean_N = 10^5, 
                          prop_cons_in_large_size_class = c(0.5)) 
  data_sim <- as.data.table(data_sim)  
  compute_rel_propagule_pressure_reduction <- function(a, b, C, B, prop_cons_in_large_size_class, mean_N, SD){
    # Bimodal on log scale. We can use p = 0.2; sqrt(p * (1-p)) * 1/sqrt(p * (1-p))* SD to rescale for various p
    N1 <- rbinom(n = C, size = 1, prob = prop_cons_in_large_size_class) * 
      1/sqrt(prop_cons_in_large_size_class * (1-prop_cons_in_large_size_class)) * SD #
    N1 <- exp(N1)   
    N1 <- N1 * mean_N / mean(N1) # Rescale
    n1 <- compute_approximate_optimal_sample_size(a = a, b = b, N = N1, B = B)
    n2 <- B / C
    propagule_pressure1 <- 
      sum(compute_propagule_pressure_per_cons(a = a, b = b, n = n1, N = N1))
    propagule_pressure2 <- 
      sum(compute_propagule_pressure_per_cons(a = a, b = b, n = n2, N = N1))
    return(propagule_pressure1 / propagule_pressure2)
  }
  data_sim[, leak_ratio := 
    compute_rel_propagule_pressure_reduction(
      a = a, 
      b = b, 
      C = C, 
      B = B, 
      mean_N = mean_N, 
      SD = SD, 
      prop_cons_in_large_size_class = prop_cons_in_large_size_class
    ), 
    by = .(prop_cons_in_large_size_class, a, b, SD)]
  data_sim[, leak_ratio_pred := exp((- 1/2 + 1 / (2 * (a + 2))) * SD^2) ] 
  
  data_sim[, tail(leak_ratio, 1), by = prop_cons_in_large_size_class]  #! The effect of prop is stricking and so regular  
  plot(data_sim[, tail(leak_ratio, 1), by = prop_cons_in_large_size_class]); abline(0, 1) # But this depends on alpha and also beta  
# Plot
  p5 <- 
    ggplot(
      data = data_sim, 
      aes(
        x = SD, 
        y = leak_ratio, 
        group = prop_cons_in_large_size_class
      )
    ) + 
    geom_point(col = 'red3') + 
    geom_line(aes(SD, leak_ratio_pred), col = 'black') + 
    facet_wrap(~'log(N_i) ~ Bimodal') +
    theme_bw() + 
    theme(
      panel.grid.minor = element_blank(), 
      plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), 'cm'), 
      axis.text.y = element_blank(), 
      axis.ticks.y = element_blank()
    ) +
    xlab(TeX('$\\sigma_{log(N)}$')) + 
    ylab('') 
  p5
  
# Student t distribution  
  set.seed(1001)
  data_sim <- expand.grid(SD = seq(0, 4, length = 20), 
                          a = 0.1, 
                          b = 8, 
                          C = 10^6, 
                          B = 10^6 * 600, 
                          mean_N = 10^5)
  data_sim <- as.data.table(data_sim) 
  compute_rel_propagule_pressure_reduction <- function(a, b, C, B, mean_N, SD){
    DF <- 10
    rt_N <- rt(n = C, df = DF)
    rt_N_scaled <- rt_N / sqrt(DF / (DF - 2)) * SD # sd = sqrt(df / (df - 2))
    N1 <- exp(rt_N_scaled)   
    N1 <- N1 * mean_N / mean(N1)
    n1 <- compute_approximate_optimal_sample_size(a = a, b = b, N = N1, B = B)
    n2 <- B / C
    propagule_pressure1 <- 
      sum(compute_propagule_pressure_per_cons(a = a, b = b, n = n1, N = N1))
    propagule_pressure2 <- 
      sum(compute_propagule_pressure_per_cons(a = a, b = b, n = n2, N = N1))
    return(propagule_pressure1 / propagule_pressure2)
  }
  data_sim[, leak_ratio := 
    compute_rel_propagule_pressure_reduction(a = a, 
                                             b = b, 
                                             C = C, 
                                             B = B, 
                                             mean_N = mean_N, 
                                             SD = SD
    ), 
    by = .(a, SD)]
  data_sim[, leak_ratio_pred := exp((- 1/2 + 1 / (2 * (a + 2))) * SD^2) ]  
# Plot
  p6 <- ggplot(data = data_sim, aes(SD, leak_ratio, group = a)) + 
    geom_point(col = 'red3') + 
    geom_line(aes(SD, leak_ratio_pred), col = 'black') + 
    facet_wrap(~'log(N_i) ~ Student-t') +
    theme_bw() + 
    theme(
      panel.grid.minor = element_blank(),
      plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), 'cm'), 
      axis.text.y = element_blank(), 
      axis.ticks.y = element_blank()
    ) +
    xlab(TeX('$\\sigma_{log(N)}$')) + 
    ylab('')
  p6
  
# Right skewed distribution  
  set.seed(1002)
  data_sim <- expand.grid(SD = seq(0, 4, length = 20), a = 0.1, b = 8, C = 10^6, B = 10^6 * 600, mean_N = 10^4)
  data_sim <- as.data.table(data_sim) 
  compute_rel_propagule_pressure_reduction <- function(a, b, C, B, mean_N, SD){                      
    rgamma_N <- rgamma(n = C, shape = 2, scale = SD / sqrt(2)) # With gamma, SD = scale * sqrt(shape)
    N1 <- exp(rgamma_N)   
    N1 <- N1 * mean_N / mean(N1)
    n1 <- compute_approximate_optimal_sample_size(a = a, b = b, N = N1, B = B)
    n2 <- B / C
    propagule_pressure1 <- 
      sum(compute_propagule_pressure_per_cons(a = a, b = b, n = n1, N = N1))
    propagule_pressure2 <- 
      sum(compute_propagule_pressure_per_cons(a = a, b = b, n = n2, N = N1))
    return(propagule_pressure1 / propagule_pressure2)
  }
  data_sim[, leak_ratio := 
    compute_rel_propagule_pressure_reduction(a = a, 
                                             b = b, 
                                             C = C, 
                                             B = B, 
                                             mean_N = mean_N, 
                                             SD = SD                                        
    ), 
    by = .(a, SD)]
  data_sim[, leak_ratio_pred := exp((- 1/2 + 1 / (2 * (a + 2))) * SD^2) ]  
# Plot
  p7 <- ggplot(data = data_sim, aes(SD, leak_ratio, group = a)) + 
    geom_point(col = 'red3') + 
    geom_line(aes(SD, leak_ratio_pred), col = 'black') + 
    facet_wrap(~'log(N_i) ~ Gamma') +
    theme_bw() + 
    theme(panel.grid.minor = element_blank()) + 
    theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), 'cm'), axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
    xlab(TeX('$\\sigma_{log(N)}$')) + 
    ylab('')
  p7

 pdf(file = './figs/theoretical_bias_vs_sd_logN_alternative_distr.pdf', 
     width = 20/2.54, 
     height = 5/2.54)
   plot_grid(p4, p5, p6, p7, ncol = 4, rel_widths = c(.9, .9, .9, .9))
 dev.off()
