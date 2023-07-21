################################################################################
#
# Create figure 1 showing a scatterplot of the alpha and beta parameter values 
# of the beta distribution associated with each pathways
#                                                                
################################################################################

# Set working directory
  setwd("C:/Users/rtrouve/Dropbox/overdispersion_paper")

# Load packages
  library(ggplot2)
  library(data.table)
  library(brms)
  library(latex2exp) # For LateX expresions in ggplot 
  library(RColorBrewer)
  library(cowplot)
  source('./functions/extract_bb_parameters_per_pathway_from_brms_fit.R')  
  tick_break_log_scale <- 
    c(0.01, 0.03, 0.1, 0.3, 1, 3, 10, 30, 100, 300, 1000, 3000, 10000)
  
################################################################################
# Read data
  data_all_pathways <- fread('./data/data_all_pathways.csv')
     
# Posterior predictive check
  load('./outputs/fit_beta_binomial_hierarchical.Rdata')
  
# Plot coefficients 
  data_coef_summary <- 
    extract_bb_parameters_per_pathway_from_brms_fit(
      model = fit_beta_binomial_hierarchical
    )
  data_coef_summary <- 
    merge(
      data_coef_summary, 
      data_all_pathways[, .(data_type = data_type[1]), by = group_id], 
      by = 'group_id'
    )
  
# Scatterplot alpha vs. beta 
  plot_main <- 
    ggplot(
      data = data_coef_summary, 
      aes(alpha_mean, beta_mean, col = data_type)
    ) +
    geom_point() +
    scale_x_log10(breaks = tick_break_log_scale) +
    scale_y_log10(breaks = tick_break_log_scale) +
    scale_colour_brewer(palette = 'Dark2', name = '') +
    xlab(TeX('$\\alpha_j$')) +
    ylab(TeX('$\\beta_j$')) +
    theme_bw() +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = c(0.2, 0.83),
      legend.title = element_blank(),
      legend.background = element_blank()
    ) +
    annotate(
      geom = 'segment',
      x = 0.175,
      xend = 0.175,
      y = 580,
      yend = 100,
      arrow = arrow(length = unit(0.15, "cm")),
      col = 'black',
      size = 0.4
    ) +
    annotate(geom = 'point', x = c(.35, 0.009), y = c(94.1 + .35, 3.3062), col = 'black')
  plot_main
    
# Create p_j histogram to be used as inset. 
  # Use the 'FP02' pathway which is at the top and center of the figure  
  select_group_id <- c('FP02')
  data_FP02 <- data_all_pathways[group_id %in% select_group_id, ]
  data_coef_summary_FP02 <- data_coef_summary[group_id %in% select_group_id, ] 
  plot_inset <- ggplot(data = data_FP02, aes(k / n)) +
  stat_bin(aes(y = ..density..), fill = 'grey70', bins = 20) +
  stat_function(
    fun = dbeta,
    args = list(
      shape1 = data_coef_summary_FP02$alpha_mean,
      shape2 = data_coef_summary_FP02$beta_mean
    ),
    col = 'black'
  ) +
  scale_y_sqrt() +
  xlab(TeX('$p_j$')) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 6),
    panel.grid = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.x = element_text(size = 8)
  )
  plot_inset

# Add both figures together and save
  plot_with_inset <- ggdraw() +
    draw_plot(plot_main) +
    draw_plot(plot_inset, x = 0.62, y = .62, width = .33, height = .35)
  ggsave(
    plot_with_inset, 
    filename = './figs/fig01_scatterplot_of_beta_vs_alpha_values_with_literature_values.pdf', 
    width = 10, 
    height = 10, 
    units = 'cm'
  )
  
    
             