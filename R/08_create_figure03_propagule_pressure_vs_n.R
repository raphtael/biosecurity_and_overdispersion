################################################################################
#
#        Create figure 3 showing propagule pressure vs. sample size 
#                                                    
################################################################################

# Set working directory
  setwd("C:/Users/rtrouve/Dropbox/overdispersion_paper")
  
# Load libraries
  library(ggplot2)
  library(data.table)
  tick_break_log_scale <- c(0.1, 0.3, 1, 3, 10, 30, 100)
  source('./functions/compute_propagule_pressure_per_cons.R')
  source('./functions/compute_propagule_pressure_per_cons_no_overdispersion.R')
  
################################################################################
# Simulate data. Same mean infestation rate p = 0.5% but varying overdispersion
  data_sim1 <- expand.grid(n = 0:1200, a = 0.05, b = 9.95, N = 10^4)
  data_sim2 <- expand.grid(n = 0:1200, a = 0.5, b = 99.5, N = 10^4)

# Compute propagule pressure  
  data_sim1 <- as.data.table(data_sim1)
  data_sim1[, propagule_pressure1 := 
    compute_propagule_pressure_per_cons(
      a = a, 
      b = b, 
      n = n, 
      N = N
    )]
  data_sim1[, propagule_pressure2 := 
    compute_propagule_pressure_per_cons_no_overdispersion(
      p = a / (a + b),
      n = n,
      N = N
    )]
  data_sim2 <- as.data.table(data_sim2)
  data_sim2[, propagule_pressure1 := 
    compute_propagule_pressure_per_cons(
      a = a, 
      b = b, 
      n = n, 
      N = N
    )]
  data_sim2[, propagule_pressure2 := 
    compute_propagule_pressure_per_cons_no_overdispersion(
      p = a / (a + b), 
      n = n, 
      N = N
    )]
  
# Plot figure
  p <- ggplot() +
    geom_line(
      data = data_sim1, 
      aes(n, propagule_pressure1, col = 'High overdispersion; p_j ~ Beta(0.05, 9.95)')
    ) +
    geom_line(
      data = data_sim2,
      aes(n, propagule_pressure1, col = 'Low overdispersion; p_j ~ Beta(0.5, 99.5)')
    ) +
    geom_line(
      data = data_sim1,
      aes(n, propagule_pressure2, col = 'No overdispersion; p = 0.005')
    ) + 
    scale_y_log10(
      breaks = tick_break_log_scale,
      labels = scales::comma
    ) +
    scale_x_continuous(
      breaks = seq(0, 1200, by = 200),
      labels = scales::comma
    ) +
    coord_cartesian(ylim = c(0.1, 100)) +
    xlab('Inspection sample size n') +
    ylab('Propagule pressure per consignement') + 
    scale_colour_brewer(palette = 'Blues', direction = -1) +
    theme_bw() + 
    theme(
      panel.grid.minor = element_blank(),
      legend.position = c(.57, .9), 
      legend.title = element_blank(), 
      legend.direction = 'vertical', 
      legend.background = element_blank(),
      legend.text = element_text(size = 8.5)
    )
      
# Save figure
  ggsave(
    plot = p, 
    filename = './figs/fig03_propagule_pressure_vs_n.pdf', 
    width = 10, 
    height = 10, 
    units = 'cm'
  )
  
   
   