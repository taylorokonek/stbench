#' Benchmark posterior draws via a constrained Bayes estimator approach 
#' 
#' Benchmarks draws from a distribution using the method described in
#' the paper "Bayesian benchmarking with applications to small area estimation" by
#' G.S. Datta, M. Ghosh, R. Steorts, and J. Maples. Takes in a matrix of posterior draws,
#' and national-level benchmark data.
#' 
#' @param posterior_samps a matrix containing posterior draws at the area-level. Rows should be
#' arranged in order arrange(region, time), with each column containing an independent draw from
#' the posterior.
#' @param nregion the number of regions you have estimates for. Defaults to the number of rows
#' in \code{posterior_samps}.
#' @param ntime the number of time points you have estimates for. Defaults to 1.
#' @param natl a vector of national level estimates, arranged in order of time if you have multiple
#' time points in your data.
#' @param pop_weights a vector of population weights for use in the benchmarking constraint. Must
#' sum to one at each time point, and be in order arrange(region, time)
#' @return A matrix of benchmarked samples in order arrange(region, time)
#' 
#' @author Taylor Okonek
#' @export benchmark_bayesest

benchmark_bayesest <- function(posterior_samps, 
                              nregion = NA, 
                              ntime = 1, 
                              natl,
                              pop_weights = NA) {
  
  # error handling
  
  # make sure dimensions of posterior_samps align with nregion and ntime
  if (is.na(nregion)) {
    nregion <- nrow(posterior_samps)
  }
  
  # if nregion > 1 and pop_weights is NA, throw error. otherwise, set pop_weights = rep(1, ntime)
  if (nregion > 1) {
    if (length(pop_weights) == 1) {
      if (is.na(pop_weights)) {
        stop("must specify pop_weights if nregion > 1")
      }
    }
  } else {
    pop_weights <- rep(1, ntime)
  }
  
  if ((ntime * nregion) != nrow(posterior_samps)) {
    stop("ntime * nregion must equal the number of rows in posterior_samps")
  }
  
  # construct region_time_id_df 
  region_time_id_df <- expand.grid(region = 1:nregion, time = 1:ntime) %>%
    arrange(region, time)
  
  # create matrix for benchmarked samples
  bench_mat <- matrix(NA, nrow = nregion * ntime, ncol = ncol(posterior_samps))
  
  for (i in 1:ntime) {
    # get smaller samp_mat only containing the year of interest
    temp_mat <- posterior_samps[which(region_time_id_df$time == i),]
    
    # compute posterior means for each area
    bayes_est <- apply(temp_mat, 1, mean)
    
    # get national-level estimate to benchmark to for this year
    t <- natl[i]
    
    # get weights for this time point
    temp_pop <- pop_weights[which(region_time_id_df$time == i)]
    
    # use benchmark_project to get projected means for all regions in this year
    bench_ests <- benchmark_project(temp_pop, bayes_est, t)
    
    # use benchmark_project_samp to get projected samples for all regions in this year
    for (j in 1:nrow(temp_mat)) {
      bench_mat[which(region_time_id_df$region == j & region_time_id_df$time == i),] <- unlist(benchmark_project_samp(temp_pop, bayes_est, t, samp_j = temp_mat[j,], weight_j = temp_pop[j]))
    }
    
  }
  
  # return matrix of benchmarked samples
  return(bench_mat)
  
}