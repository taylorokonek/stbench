#' Compute acceptance probability in M-H algorithm
#' 
#' @param old_thetas a vector of old fitted values
#' @param new_thetas a vector of proposed fitted values
#' @param old_betas a vector of old fixed effects
#' @param new_betas a vector of proposed fixed effects
#' @param z a vector of national level estimates, in order arrange(time)
#' @param intercept_means a vector of fixed effect means
#' @param var_z a vector of national level variances, in order arrange(time)
#' @param var_plus a vector of fixed effect variances
#' @return the acceptance probability
#' 
#' @author Taylor Okonek
#' @noRd
#' @keywords internal
A_full <- function(old_thetas, new_thetas, old_betas, new_betas, 
                   z, intercept_means, var_z, var_plus) {
  num <- 1
  denom <- 1
  for (i in 1:length(z)) {
    num <- num * dnorm(z[i], new_thetas[i], sqrt(var_z[i]))
    denom <- denom * dnorm(z[i], old_thetas[i], sqrt(var_z[i]))
  }
  for (i in 1:length(old_betas)) {
    num <- num * dnorm(old_betas[i], intercept_means[i], 
                       sqrt(var_plus[i]))
    denom <- denom * dnorm(new_betas[i], intercept_means[i], 
                           sqrt(var_plus[i]))
  }
  return(min(1, num/denom))
}