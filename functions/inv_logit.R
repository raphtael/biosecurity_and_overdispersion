# Inverse logit function 
  inv_logit <- function(x) {
    1 / (1 + exp(-x))
  }