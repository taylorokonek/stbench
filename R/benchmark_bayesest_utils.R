#' Project a sample from a posterior into the constrained benchmarked space
#' using exact benchmarking and MSE loss
#' 
#' @param weight_vec a vector of weights, typically standardized population counts
#' @param bayes_est a vector of posterior means computed from the posterior samples
#' @param t the national level estimate for a given year
#' @return the projected sample
#' 
#' @author Taylor Okonek
#' @noRd
#' @keywords internal

benchmark_project <- function(weight_vec, bayes_est, t) {
  bayes_est + (1/sum(weight_vec^2)) * (t - sum(weight_vec * bayes_est)) * weight_vec
}

#' Vectorized version of benchmark_project 
#' 
#' @param weight_vec a vector of weights, typically standardized population counts
#' @param bayes_est a vector of posterior means computed from the posterior samples
#' @param samp_j sample
#' @param weight_j weight
#' @param t the national level estimate for a given year
#' @return the projected sample
#' 
#' @author Taylor Okonek
#' @noRd
#' @keywords internal
# 
benchmark_project_samp <- function(weight_vec, bayes_est, t, samp_j, weight_j) {
  samp_j + (1/sum(weight_vec^2)) * (t - sum(weight_vec * bayes_est)) * weight_j
}
Vectorize(benchmark_project_samp, vectorize.args = "samp_j")