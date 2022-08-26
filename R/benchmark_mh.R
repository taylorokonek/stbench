#' Benchmark posterior draws via a Metropolis-Hastings algorithm
#' 
#' Benchmarks draws from a posterior distribution using a Metropolis-Hastings approach. Takes in 
#' a matrix of posterior draws, a matrix of fixed effect draws, and national-level benchmark data.
#' 
#' @param posterior_samps a matrix containing posterior draws at the area-level. Rows should be
#' arranged in order arrange(region, time), with each column containing an independent draw from
#' the posterior.
#' @param posterior_samps_fe a matrix containing posterior draws of fixed effects (typically intercepts)
#' that have had prior means and variances shifted in the unbenchmarked model, with each column
#' containing an independent draw from the posterior.
#' @param nregion the number of regions you have estimates for. Defaults to the number of rows
#' in \code{posterior_samps}.
#' @param ntime the number of time points you have estimates for. Defaults to 1.
#' @param natl a vector of national level estimates, arranged in order of time if you have multiple
#' time points in your data.
#' @param natl_sd a vector of standard deviations for national level estimates, arranged in order of
#' time if you have multiple time points in your data.
#' @param pop_weights a vector of population weights for use in the benchmarking constraint. Must
#' sum to one at each time point, and be in order arrange(region, time)
#' @param fe_prior_means a vector of prior means for the fixed effects specified in \code{posterior_samps_fe}. Must
#' be in row order of posterior_samps_fe.
#' @param fe_prior_variances a vector of prior variances for the fixed effects specified in \code{posterior_samps_fe} 
#' Must be in row order of posterior_samps_fe.
#' @return A list containing: 
#' \itemize{
#' \item fitted_list: a list of matrices of benchmarked posterior samples of fitted values in 
#' order arrange(time). Each matrix will have rows arranged in order arrange(region).
#' \item natl_list: a list of vectors containing aggregated national-level samples that were 
#' accepted, in order arrange(time)
#' \item prop_accepted: the proportion of samples accepted during sampling. 
#' } 
#' 
#' @author Taylor Okonek
#' @export benchmark_mh

benchmark_mh <- function(posterior_samps, 
                         posterior_samps_fe,
                         nregion = NA, 
                         ntime = 1, 
                         natl,
                         natl_sd, 
                         pop_weights = NA,
                         fe_prior_means,
                         fe_prior_variances) {
  
  # error handling
  
  # make sure dimensions of posterior_samps align with nregion and ntime
  if (is.na(nregion)) {
    nregion <- nrow(posterior_samps)
  }
  
  # make sure dimensions of fe_prior_means and fe_prior_variances align with posterior_samps_fe
  if (length(fe_prior_means) != length(fe_prior_variances)) {
    stop("must specify the same number of fixed effect prior means and variances in fe_prior_means and fe_prior_variances")
  }
  if (length(fe_prior_means) != nrow(posterior_samps_fe)) {
    stop("must specify prior means and variances for all fixed effects in posterior_samps_fe")
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
  
  # prepare matrix for q at each time point
  q_mat <- matrix(0, nrow = ntime, ncol = ncol(posterior_samps))
  
  # fill in q_mat
  if (nregion == 1) {
    q_mat <- posterior_samps
  } else {
    for (i in 1:ntime) {
      temp_samps <- posterior_samps[region_time_id_df$time == i,]
      temp_pop <- pop_weights[region_time_id_df$time == i]
      
      for (j in 1:ncol(posterior_samps)) {
        q_mat[i,j] <- sum(temp_samps[,j] * temp_pop)
      }
    }
  }
  
  ## MH algorithm
  
  # set-up
  old_thetas <- q_mat[,1]
  old_betas <- posterior_samps_fe[,1]
  accepted_ids <- c()
  num_accepted <- 0
  
  for (i in 2:ncol(q_mat)) {
    new_thetas <- q_mat[, i]
    new_betas <- posterior_samps_fe[, i]
    accept_prob <- A_full(old_thetas = old_thetas, old_betas = old_betas, 
                          new_thetas = new_thetas, new_betas = new_betas, 
                          z = natl, intercept_means = fe_prior_means, 
                          var_z = natl_sd^2, var_plus = fe_prior_variances)
    accept_yn <- rbinom(n = 1, size = 1, prob = accept_prob)
    if (accept_yn) {
      accepted_ids <- c(accepted_ids, i)
      old_betas <- new_betas
      old_thetas <- new_thetas
      num_accepted <- num_accepted + 1
    }
    else {
      if (i > 2) {
        accepted_ids <- c(accepted_ids, tail(accepted_ids, 1))
      }
    }
  }
  
  # set up fitted_list and prop_accepted for return
  fitted_list_combined <- list()
  q_list_combined <- list()
  
  for (i in 1:ntime) {
    fitted_list_combined[[i]] <- posterior_samps[region_time_id_df$time == i, accepted_ids]
    q_list_combined[[i]] <- q_mat[i, accepted_ids]
  }
  
  # compute proportion of samples accepted
  nsamps <- ncol(posterior_samps)
  prop_accepted <- num_accepted/nsamps
  
  return(list(fitted_list = fitted_list_combined,
    natl_list = q_list_combined,
    prop_accepted = prop_accepted))
  
}