# Compute leave-one-out (LOO) cross-validation metrics for each group present in 
# the hierarchcial model
  # model: brms model fit
  # k: Vector of number of infested units found in each inspection
  # n: Vector of sample size associated with each inspection
  # group_id: Grouping variable (pathway ID) associated with each inspection
compute_loo_per_group <- function(model, k, n, group_id){
  selected_data <- cbind.data.frame(k = k, n = n, group_id = group_id)
  loo_for_selected_data <- brms::loo(x = model, newdata = selected_data)
  return(loo_for_selected_data$estimates['looic', 'Estimate'])
}