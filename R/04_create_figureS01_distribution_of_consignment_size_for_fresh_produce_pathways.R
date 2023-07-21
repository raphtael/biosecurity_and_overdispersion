################################################################################
#
#     Plot the distribution of consignment size for fresh produce pathways
#                                                   
################################################################################

# Set working directory
  setwd("C:/Users/rtrouve/Dropbox/overdispersion_paper")

# Load packages
  library(ggplot2)
  library(data.table)
  library(latex2exp)
  tick_break_log_scale <- 
    c(1000, 3000, 10000, 30000, 10^5, 3 * 10^5, 10^6, 3 * 10^6, 10^7)
  
################################################################################
# Read data
  data_all_pathways <- fread('./data/data_all_pathways.csv')

# Filter data for fresh produce (have reliable consignment size)
  data_fresh_produce <- data_all_pathways[data_type == 'fresh_produce', ]
  
# Compute mean N and SD log(N) per pathway  
  data_fresh_produce[, sd_N := sd(log(N)), by = group_id]
  data_fresh_produce[, mean_N := mean(N), by = group_id]
  data_fresh_produce[, sd_N_label := 
    paste0('$\\sigma_{log(N)} = ', format(x = round(sd_N, 2), nsmall = 2), '$')]
  data_fresh_produce[, mean_N_label := 
    paste0('$\\bar{N} = ', format(x = round(mean_N), big.mark = ","), '$')]         
  data_fresh_produce[, group_label := 
    paste0(group_id, ': ', mean_N_label, '; ', sd_N_label)]
  
# Plot distribution of consignment size per pathway
  p <- ggplot(data = data_fresh_produce, aes(N)) + 
    geom_histogram() + 
    facet_wrap(~TeX(group_label, output = "character"), 
      scale = 'free_y', 
      ncol = 5, 
      labeller = label_parsed
    ) + 
    scale_x_log10(breaks = tick_break_log_scale, label = scales::comma) + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
  ggsave(
    plot = p, 
    filename = './figs/hist_logN_fresh_produces.pdf', 
    width = 12, 
    height = 12
  )
  
  