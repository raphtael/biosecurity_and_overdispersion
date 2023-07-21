################################################################################
#
#   Create figure 4 showing a toy example of propagule pressure optimisation 
#
################################################################################

# set working directory
  setwd("C:/Users/rtrouve/Dropbox/overdispersion_paper")

# load packages
  library(ggplot2)
  library(data.table)
  library(cowplot)
  library(latex2exp)
  source('./functions/compute_propagule_pressure_per_cons.R')
  source('./functions/compute_propagule_pressure_per_cons_no_overdispersion.R')
  tick_break_log_scale <- c(1, 3, 10, 30, 100)
 
################################################################################
# Varying consignment size. Pathway alpha and beta parameters (FP02).
  data_sim <- expand.grid(a1 = 0.18,   
                          a2 = 0.18, 
                          b1 = 87.4, 
                          b2 = 87.4, 
                          n1 = 0:1200, 
                          N1 = 10^4, 
                          N2 = 10^5) 
  data_sim <- as.data.table(data_sim)
  data_sim[, n2 := max(n1) - n1]
  
# Compute mean leakage using closed form solutions    
  data_sim[, leak_pathway1 := 
    compute_propagule_pressure_per_cons(a = a1, b = b1, n = n1, N = N1)] 
  data_sim[, leak_pathway2 := 
    compute_propagule_pressure_per_cons(a = a2, b = b2, n = n2, N = N2)] 
  data_sim[, leak_no_overdisp_pathway1 := 
    compute_propagule_pressure_per_cons_no_overdispersion(
      p = a1 / (a1 + b1), 
      n = n1, 
      N = N1
    )] 
  data_sim[, leak_no_overdisp_pathway2 := 
    compute_propagule_pressure_per_cons_no_overdispersion(
      p = a2 / (a2 + b2), 
      n = n2, 
      N = N2
    )] 
  data_sim[, leak_total := leak_pathway1 + leak_pathway2]
  data_sim[, leak_no_overdisp_total := 
    leak_no_overdisp_pathway1 + leak_no_overdisp_pathway2]
 
# Compute sample size under different sampling regimes 
  optimal_n1 <- data_sim[leak_total == min(leak_total), n1]
  fixed_n1 <- 600
  # In hypergeometric sampling, consignments have the same design prevalence 
  # (e.g., 0.5%) and same sensitivity and rather than the same sample size
  data_sim[, S_hypergeometric1 := 
    1 - dhyper(x = 0, m = 0.005 * N1, n = N1 - 0.005 * N1, k = n1)]
  data_sim[, S_hypergeometric2 := 
    1 - dhyper(x = 0, m = 0.005 * N2, n = N2 - 0.005 * N2, k = n2)]
  hypergeometric_n1 <-  
    data_sim[which.min((S_hypergeometric1 - S_hypergeometric2)^2), n1]
  prop_n1 <- data_sim[n1 == round(N1 / (N1 + N2) * max(n1)), n1]
  optimal_no_overdisp_n1 <- 
    data_sim[leak_no_overdisp_total == min(leak_no_overdisp_total), n1]
  
# Compute propagule pressure under different sampling regimes  
  leak_optimal_n1 <- data_sim[n1 == optimal_n1, leak_total]
  leak_fixed_n1 <- data_sim[n1 == fixed_n1, leak_total]
  leak_hypergeometric_n1 <- data_sim[n1 == hypergeometric_n1, leak_total] 
  leak_prop_n1 <- data_sim[n1 == prop_n1, leak_total]
  leak_min_zero_overdisp_n1 <- data_sim[n1 == optimal_no_overdisp_n1, leak_total]
  leak_optimal_n1
  leak_fixed_n1
  leak_hypergeometric_n1
  leak_prop_n1
  leak_min_zero_overdisp_n1
   
# Plot figure                                
  p1 <- ggplot(data = data_sim, aes(n1, leak_pathway1)) + 
    geom_line(col = "#0368C4", lwd = 1.2) + 
    geom_line(aes(n1, leak_pathway2), col = "#00CDCA", lwd = 1.2) +
    geom_line(aes(n1, leak_pathway1 + leak_pathway2), col = 'black', lwd = 1.2)+ 
    annotate(geom = 'point', x = optimal_n1, y = leak_optimal_n1) + 
    annotate(
      geom = 'text', 
      x = optimal_n1, 
      y = leak_optimal_n1 + 3, 
      label = round(leak_optimal_n1, 1), 
      size = 3
    ) +
    annotate(geom = 'point', x = fixed_n1, y = leak_fixed_n1) +
    annotate(
      geom = 'text', 
      x = fixed_n1, 
      y = leak_fixed_n1 + 8, 
      label = round(leak_fixed_n1, 1), 
      size = 3
    ) +
    annotate(geom = 'point', x = hypergeometric_n1, y = leak_hypergeometric_n1)+ 
    annotate(
      geom = 'text', 
      x = hypergeometric_n1 - 30, 
      y = leak_hypergeometric_n1 + 3, 
      label = round(leak_hypergeometric_n1, 1), 
      size = 3
    ) + 
    annotate(geom = 'point', x = prop_n1, y = leak_prop_n1) + 
    annotate(
      geom = 'text', 
      x = prop_n1 + 30, 
      y = leak_prop_n1 + 5, 
      label = round(leak_prop_n1, 1), 
      size = 3
    ) + 
    annotate(
      geom = 'point', 
      x = optimal_no_overdisp_n1, 
      y = leak_min_zero_overdisp_n1
    ) + 
    annotate(
      geom = 'text', 
      x = optimal_no_overdisp_n1, 
      y = leak_min_zero_overdisp_n1 + 8, 
      label = round(leak_min_zero_overdisp_n1, 1), 
      size = 3
    ) +
    annotate(
      geom = 'point', 
      x = fixed_n1, 
      y = data_sim[n1 == fixed_n1, leak_pathway1],
      col = "#0368C4", 
    ) + 
    annotate(
      geom = 'text', 
      x = fixed_n1 + 50, 
      y = data_sim[n1 == fixed_n1, leak_pathway1] + 0.5, 
      label = round(data_sim[n1 == fixed_n1, leak_pathway1], 1), 
      size = 3,
      col = "#0368C4"
    ) + 
    annotate(
      geom = 'point', 
      x = fixed_n1, 
      y = data_sim[n1 == fixed_n1, leak_pathway2],
      col = "#00CDCA", 
    ) + 
    annotate(
      geom = 'text', 
      x = fixed_n1 + 50, 
      y = data_sim[n1 == fixed_n1, leak_pathway2] - 3, 
      label = round(data_sim[n1 == fixed_n1, leak_pathway2], 1), 
      size = 3,
      col = "#00CDCA"
    ) + 
    theme_bw() +
    theme(panel.grid.minor = element_blank()) +
    scale_y_log10(breaks = tick_break_log_scale) +
    scale_x_continuous(
      breaks = seq(0, 1200, by = 300), 
      sec.axis = sec_axis(
        ~ 1200 - ., 
        name = TeX("$n_2$"), 
        breaks = seq(0, 1200, by = 300)
      )
    ) +
    xlab(TeX("$n_1$")) +
    ylab('Propagule pressure')
  p1
  
# Save figure
  ggsave(
    plot = p1, 
    filename = './figs/fig04_toy_example_propagule_pressure_optimisation_FP02.pdf', 
    height = 6 * 1.5, 
    width = 8 * 1.5, 
    units = 'cm'
  )
  
    