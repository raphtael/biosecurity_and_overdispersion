################################################################################
#
#  Create figure 5 showing numerical and theoretical optimal sample size 
#  vs. consignment size for case study pathway FP02
#                                                   
################################################################################

# Set working directory
  setwd("C:/Users/rtrouve/Dropbox/overdispersion_paper")
  
# Load libraries
  library(ggplot2)
  library(data.table)
  library(cowplot)
  source('./functions/extract_bb_parameters_per_pathway_from_brms_fit.R')
  source('./functions/compute_optimal_sample_size_greedy.R')
  source('./functions/compute_propagule_pressure_per_cons.R')
  
################################################################################  
# Read data
  data_all_pathways <- fread('./data/data_all_pathways.csv') 

# Load model  
  load('./outputs/fit_beta_binomial_hierarchical.Rdata')

# Extract coefficients of the beta-binomial and zero-inflated BB distributions
  data_coef_summary_bb <- 
    extract_bb_parameters_per_pathway_from_brms_fit(
      model = fit_beta_binomial_hierarchical
    )
  
# Select pathway FP02 (~median pathway, also used as inset in figure 1) 
  data_pathway_FP02 <- data_all_pathways[group_id == 'FP02', ] # n_units_per_consignement
  data_coef_summary_bb_FP02 <- data_coef_summary_bb[group_id == 'FP02', ]

# Extract alpha and beta parameters
  data_pathway_FP02[, alpha_mean := data_coef_summary_bb_FP02$alpha_mean] 
  data_pathway_FP02[, beta_mean := data_coef_summary_bb_FP02$beta_mean] 

# Compute optimal sample size using greedy algorithm   
  ptm <- proc.time()  # Takes ~3 mn
    data_pathway_FP02[, optimal_n_greedy := 
      compute_optimal_sample_size_greedy(
        a = alpha_mean, 
        b = beta_mean, 
        z = 0, 
        N = N, 
        B = .N * 600, 
        n_increment_resolution = 1
      )
    ] 
  proc.time() - ptm                          
    
# Compute optimal sample size using theory 
  data_pred_FP02 <- expand.grid(
    N = seq(0, max(data_pathway_FP02$N), length = 10000), 
    a = data_pathway_FP02$alpha_mean[1], 
    b = data_pathway_FP02$beta_mean[1], 
    C = nrow(data_pathway_FP02), 
    B = nrow(data_pathway_FP02) * 600)
  data_pred_FP02 <- as.data.table(data_pred_FP02)
  data_pred_FP02[, constant :=  # See Eq.4 in manuscript
    (C * b + B) / sum(data_pathway_FP02$N^(1 / (mean(a) + 2))) ] 
  data_pred_FP02[, optimal_n_theory := - b + constant * N^(1/(a + 2))]
  data_pred_FP02[optimal_n_theory < 0, optimal_n_theory := 0]
  
# Create main plot   
  plot_main <- ggplot(data = data_pathway_FP02, aes(N, optimal_n_greedy)) + 
    geom_point(size = .5) + 
    geom_line(
      data = data_pred_FP02, 
      aes(N, optimal_n_theory), 
      col = 'dodgerblue'
    ) + 
    scale_x_continuous(
      breaks = seq(0, 10^6, by = 10^5), labels = scales::comma
    ) +
    scale_y_continuous(
      breaks = seq(0, 10^4, by = 200), labels = scales::comma
    ) +
    coord_cartesian(xlim = c(0, 450000), ylim = c(0, 1600)) + 
    xlab('Number of units per consignment') + 
    ylab('Optimal inspection sample size') + 
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.minor=element_blank()
    )
  plot_main
   
# Create zoomed-in inset 
  data_pred_FP02_inset <- expand.grid(
    N = seq(0, 12000, length = 1000), 
    a = data_pathway_FP02$alpha_mean[1], 
    b = data_pathway_FP02$beta_mean[1], 
    C = nrow(data_pathway_FP02), 
    B = nrow(data_pathway_FP02) * 600)
  data_pred_FP02_inset <- as.data.table(data_pred_FP02_inset)
  data_pred_FP02_inset[, constant := 
    (C * b + B) / sum(data_pathway_FP02$N^(1/(mean(a) + 2))) ]
  data_pred_FP02_inset[, optimal_n_theory := - b + constant * N^(1/(a + 2))]
  data_pred_FP02_inset[optimal_n_theory < 0, optimal_n_theory := 0]
  plot_inset <- ggplot() + 
    geom_line(
      data = data_pred_FP02_inset, 
      aes(N, optimal_n_theory), 
      col = 'dodgerblue'
    ) + 
    xlab('') + 
    ylab('') +
    geom_point(
      data = data_pathway_FP02[N < 12000, ], 
      aes(N, optimal_n_greedy), 
      size = .5
    ) +
    scale_x_continuous(
      breaks = seq(0, 12000, by = 3000), labels = scales::comma
    ) + 
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text = element_text(size = 6),
      axis.title = element_blank(),
      panel.grid.minor=element_blank()
    )
  plot_inset
  
# Combine both plots  
  plot_main_plus_inset <- plot_main + 
    annotation_custom(
      ggplotGrob(plot_inset), 
      xmin = 2 * 10^5, 
      xmax = 4.5 * 10^5, 
      ymin = -20, 
      ymax = 800
    )
  plot_main_plus_inset

# Save figure   
  ggsave(
    plot = plot_main_plus_inset, 
    filename = './figs/fig05_optimal_sample_size_vs_N_for_pathway_FP02.pdf', 
    width = 10, 
    height = 10, 
    units = 'cm'
  )   


