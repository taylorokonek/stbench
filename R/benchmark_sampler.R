#' Benchmark posterior draws via a rejection sampler
#' 
#' Benchmarks draws from a posterior distribution using a rejection sampling approach. Takes in 
#' a matrix of posterior draws, and national-level benchmark data.
#' 
#' @param posterior_samps a matrix containing posterior draws at the area-level. Rows should be
#' arranged in order arrange(region, time), with each column containing an independent draw from
#' the posterior.
#' @param nregion the number of regions you have estimates for. Defaults to the number of rows
#' in \code{posterior_samps}.
#' @param ntime the number of time points you have estimates for. Defaults to 1.
#' @param natl a vector of national level estimates, arranged in order of time if you have multiple
#' time points in your data.
#' @param natl_sd a vector of standard deviations for national level estimates, arranged in order of
#' time if you have multiple time points in your data.
#' @param pop_weights a vector of population weights for use in the benchmarking constraint. Must
#' sum to one at each time point, and be in order arrange(region, time)
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
#' @export benchmark_sampler

benchmark_sampler <- function(posterior_samps, 
                              nregion = NA, 
                              ntime = 1, 
                              natl,
                              natl_sd, 
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
  
  
  ## rejection sampler
  
  # generate Us
  U <- runif(ncol(posterior_samps), 0, 1)
  
  # set up fitted_list and prop_accepted for return
  fitted_list <- list()
  prop_accepted <- 0
  
  # COULD BE MORE EFFICIENT
  accept_ratio <- matrix(0, nrow = nrow(q_mat), ncol = ncol(q_mat))
  
  # this does the same thing as below but with matrices
  # accept_vec <- c()
  # sigma_inv <- diag(ntime) * (1/(natl_sd)^2)
  # for (j in 1:ncol(q_mat)) {
  #   accept_vec[j] <- exp(- ((q_mat[,j] - natl) %*% sigma_inv %*% (q_mat[,j] - natl))/2)
  # }
  
  for (i in 1:ntime) {
    for (j in 1:ncol(q_mat)) {
      accept_ratio[i, j] <- exp((-1 / (2 * natl_sd[i] ^ 2)) * (q_mat[i, j] - natl[i]) ^ 2)
    }
  }
  
  # get accepted samples if we take all time points at once
  multi_accepted_samps <- rep(FALSE, ncol(accept_ratio))
  for (i in 1:ncol(accept_ratio)) {
    multi_accepted_samps[i] <- ( U[i] < prod(accept_ratio[,i]) )
  }
  
  # get accepted samples for each time point separately
  U_uni <- runif(ncol(posterior_samps) * ntime, 0, 1) %>%
    matrix(nrow = ntime)
  
  uni_accepted_samps <- matrix(FALSE, ncol = ncol(accept_ratio), nrow = nrow(accept_ratio))
  for (i in 1:nrow(accept_ratio)) {
    for (j in 1:ncol(accept_ratio)) {
      uni_accepted_samps[i, j] <- U_uni[i, j] < accept_ratio[i, j]
    }
  }
  
  fitted_list_combined <- list()
  fitted_list_separate <- list()
  q_list_combined <- list()
  q_list_separate <- list()
  
  for (i in 1:ntime) {
    fitted_list_combined[[i]] <- posterior_samps[region_time_id_df$time == i, multi_accepted_samps]
    fitted_list_separate[[i]] <- posterior_samps[region_time_id_df$time == i, uni_accepted_samps[i,]]
    q_list_combined[[i]] <- q_mat[i, multi_accepted_samps]
    q_list_separate[[i]] <- q_mat[i, uni_accepted_samps[i,]]
  }
  
  # compute proportion of samples accepted
  nsamps <- ncol(posterior_samps)
  prop_accepted <- length(q_list_combined[[1]])/nsamps
  
  return(list(#fitted_list_combined = fitted_list_combined,
              #fitted_list_separate = fitted_list_separate,
              fitted_list = fitted_list_combined,
              #q_list_combined = q_list_combined,
              #q_list_separate = q_list_separate,
              natl_list = q_list_combined,
              prop_accepted = prop_accepted))
              #multi_accepted_samps = multi_accepted_samps,
              #uni_accepted_samps = uni_accepted_samps))
  
}